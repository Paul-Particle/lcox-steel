"""Shared logging setup for the lcox-steel pipeline.

One-stop configuration so every Snakemake script and standalone run produces the
same structured output. Adapted from PyPSA-Eur's `_helpers.configure_logging`.

Usage in a script
-----------------
    import logging
    from common._logging import configure_logging, progress

    if "snakemake" not in globals():
        from common._stubs import snakemake

    configure_logging(snakemake)
    log = logging.getLogger(__name__)

    log.info("starting %s", area)
    for item in progress(items, desc="cells"):
        ...

Behaviour
---------
* Under Snakemake (rule has `log: "logs/.../foo.log"`): file handler at that
  path, plus a stderr stream handler. Engine messages still go to the top-level
  `--log logs/snakemake.log` if you pass that flag.
* Standalone (`python workflow/scripts/.../foo.py`): falls back to
  `logs/<script_stem>.log` (best effort) + stderr.
* `sys.excepthook` is overridden so uncaught exceptions land in the log file
  rather than bypassing it.

Verbosity
---------
Driven by `config.yaml`:

    logging:
      level: INFO                       # DEBUG for verbose runs
      third_party_level: WARNING        # atlite/cdsapi/pypsa/etc.
      format: "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"

Override on the command line: ``snakemake --config 'logging={level: DEBUG}'``.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Iterable, Iterator, TypeVar

# Loggers that flood at INFO level — kept at WARNING by default.
_NOISY_LOGGERS = (
    "atlite",
    "cdsapi",
    "urllib3",
    "entsoe",
    "pypsa",
    "linopy",
    "fiona",
    "rasterio",
    "matplotlib",
    "numexpr",
    "asyncio",
)

_DEFAULT_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"
_DEFAULT_DATEFMT = "%Y-%m-%d %H:%M:%S"

_CONFIGURED = False


def _snakemake_config(snakemake: Any) -> dict:
    """Return `snakemake.config['logging']` (or empty dict) safely."""
    cfg = getattr(snakemake, "config", None)
    if not isinstance(cfg, dict):
        return {}
    section = cfg.get("logging", {})
    return section if isinstance(section, dict) else {}


def _snakemake_log_path(snakemake: Any) -> Path | None:
    """Resolve the Snakemake-injected log file path, if any.

    Snakemake exposes `snakemake.log` as a NamedList. We accept either a named
    'python' entry (PyPSA-Eur convention) or the first positional entry.
    """
    log_attr = getattr(snakemake, "log", None)
    if not log_attr:
        return None
    # NamedList supports .get() for named outputs; fall back to indexing.
    try:
        named = log_attr.get("python")  # type: ignore[union-attr]
    except (AttributeError, TypeError):
        named = None
    if named:
        return Path(named)
    try:
        return Path(log_attr[0])
    except (IndexError, TypeError):
        return None


def _fallback_log_path() -> Path | None:
    """Pick a per-script log file when running standalone (no snakemake.log)."""
    main = sys.modules.get("__main__")
    stem = None
    if main is not None and hasattr(main, "__file__") and main.__file__:
        stem = Path(main.__file__).stem
    if not stem:
        return None
    # Repo root is two levels above workflow/common/_logging.py.
    repo_root = Path(__file__).resolve().parent.parent.parent
    return repo_root / "logs" / f"{stem}.log"


def configure_logging(
    snakemake: Any = None,
    *,
    level: str | int | None = None,
    log_path: str | Path | None = None,
) -> None:
    """Wire root logger to a file handler and stderr, honouring config overrides.

    Calling twice is safe — handlers are replaced rather than appended.
    """
    global _CONFIGURED

    cfg = _snakemake_config(snakemake) if snakemake is not None else {}
    resolved_level = level if level is not None else cfg.get("level", "INFO")
    fmt = cfg.get("format", _DEFAULT_FORMAT)
    third_party_level = cfg.get("third_party_level", "WARNING")

    # File path priority: explicit arg > snakemake.log > script-stem fallback.
    file_path: Path | None
    if log_path is not None:
        file_path = Path(log_path)
    elif snakemake is not None:
        file_path = _snakemake_log_path(snakemake) or _fallback_log_path()
    else:
        file_path = _fallback_log_path()

    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if file_path is not None:
        try:
            file_path.parent.mkdir(parents=True, exist_ok=True)
            handlers.insert(0, logging.FileHandler(file_path))
        except OSError:
            # Read-only filesystem or similar — stay on stderr only.
            pass

    formatter = logging.Formatter(fmt=fmt, datefmt=_DEFAULT_DATEFMT)
    for h in handlers:
        h.setFormatter(formatter)

    root = logging.getLogger()
    for old in list(root.handlers):
        root.removeHandler(old)
        old.close()
    for h in handlers:
        root.addHandler(h)
    root.setLevel(resolved_level)

    tame_third_party(third_party_level)

    if not _CONFIGURED:
        sys.excepthook = _excepthook
        _CONFIGURED = True


def tame_third_party(level: str | int = "WARNING") -> None:
    """Quiet known-noisy upstream loggers."""
    for name in _NOISY_LOGGERS:
        logging.getLogger(name).setLevel(level)


def _excepthook(exc_type, exc_value, exc_traceback) -> None:
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.getLogger().error(
        "uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
    )


# ── Progress bars ─────────────────────────────────────────────────────────────

T = TypeVar("T")


def progress(
    iterable: Iterable[T],
    *,
    desc: str,
    total: int | None = None,
    unit: str = "it",
    leave: bool = False,
) -> Iterator[T]:
    """tqdm wrapper that plays nicely with Snakemake log files.

    Writes to stderr with `disable=None` — tqdm auto-disables when stderr isn't
    a TTY (i.e. when Snakemake is capturing it into a log file). On disable a
    single start/finish line is logged instead, so non-interactive runs still
    leave a trace.
    """
    from tqdm import tqdm

    log = logging.getLogger("progress")

    try:
        total_n = total if total is not None else len(iterable)  # type: ignore[arg-type]
    except TypeError:
        total_n = None

    is_tty = sys.stderr.isatty()
    if not is_tty:
        if total_n is not None:
            log.info("%s — starting (n=%d)", desc, total_n)
        else:
            log.info("%s — starting", desc)

    bar = tqdm(
        iterable,
        desc=desc,
        total=total_n,
        unit=unit,
        leave=leave,
        dynamic_ncols=True,
        mininterval=0.2,
        smoothing=0.1,
        disable=None,  # auto-off when stderr isn't a TTY
        file=sys.stderr,
    )
    try:
        for item in bar:
            yield item
    finally:
        bar.close()
        if not is_tty:
            log.info("%s — done", desc)
