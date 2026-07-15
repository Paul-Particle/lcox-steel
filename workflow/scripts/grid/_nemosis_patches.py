"""Monkey-patches for NEMOSIS to work around AEMO API changes.

Applied on import via `import _nemosis_patches` in download_nem.py.
Delete this file and that import once the upstream package ships fixes.

Patch 1: AEMO now blocks NEMOSIS's stale Chrome 80 User-Agent with HTTP 403.

Patch 2: MMSDM filenames from Aug-2024 contain literal '#' characters
  (e.g. PUBLIC_ARCHIVE#DISPATCHPRICE#FILE01#202408010000.zip). requests cannot
  fetch these — it strips '#' as a URL fragment — so we send a pre-encoded path
  via http.client instead.

  How '#' must be encoded depends on AEMO's current front-end, and it has changed:
  the direct nemweb.com.au archive wants a single '%23', but an earlier
  Azure-fronted endpoint rejected '%23' (HTTP 400) and only accepted a double
  '%2523' (which Azure decoded one layer back to '%23'). Rather than bet on one,
  we try each encoding in turn and use the first that returns a valid zip, so the
  patch survives AEMO flipping the front-end again.
"""

import http.client
import io
import ssl
import zipfile

import nemosis.downloader as _dl

_dl.USR_AGENT_HEADER.clear()
_dl.USR_AGENT_HEADER["User-Agent"] = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
)


def _fetch_zip(url_encoded: str) -> bytes:
    """GET a pre-encoded URL via http.client (bypassing requests' fragment stripping)."""
    from urllib.parse import urlsplit
    u = urlsplit(url_encoded)
    path = u.path + (f"?{u.query}" if u.query else "")
    conn = http.client.HTTPSConnection(
        u.hostname, u.port or 443, context=ssl.create_default_context(), timeout=120  # pyright: ignore
    )
    try:
        conn.request("GET", path,
                     headers={**_dl.USR_AGENT_HEADER, "Host": u.hostname, "Accept": "*/*"})  # pyright: ignore
        resp = conn.getresponse()
        status, reason = resp.status, resp.reason
        body = resp.read() if status == 200 else b""
    finally:
        conn.close()
    if status != 200:
        raise IOError(f"GET {url_encoded} -> {status} {reason}")
    return body


def _download_unzip_csv_patched(url: str, down_load_to: str) -> None:
    """Download and extract a NEMOSIS MMSDM zip, trying each '#' encoding (Patch 2)."""
    literal = url.replace("%2523", "#").replace("%23", "#")  # normalise to a literal '#'
    # Single '%23' (current nemweb) first, then double '%2523' (historical Azure front-end).
    candidates = [literal.replace("#", "%23"), literal.replace("#", "%2523")]
    last_err: Exception | None = None
    for candidate in candidates:
        try:
            body = _fetch_zip(candidate)
        except IOError as e:
            last_err = e
            continue
        try:
            zipfile.ZipFile(io.BytesIO(body)).extractall(down_load_to)
            return
        except zipfile.BadZipFile as e:  # 200 but not a real zip (e.g. an error page)
            last_err = IOError(f"GET {candidate} -> 200 but not a zip: {e}")
    raise last_err  # both encodings failed; NEMOSIS's run() logs "not downloaded"


_dl.download_unzip_csv = _download_unzip_csv_patched
