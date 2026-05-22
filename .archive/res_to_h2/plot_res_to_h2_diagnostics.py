import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import matplotlib.pyplot as plt

CONFIG_PATH = Path("config_hannah.yaml")
DATA_DIR = Path("data/res_cf")
RESULTS_DIR = Path("results")
PLOTS_DIR = RESULTS_DIR / "plots"


# ----------------------------
# Helpers
# ----------------------------
def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merges override into base (without mutating inputs)."""
    out = dict(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def is_bad_fstring_name(name_from_cfg) -> bool:
    """YAML can't evaluate f-strings. Catch common 'f\"{...}\"' / braces patterns."""
    if not name_from_cfg:
        return True
    s = str(name_from_cfg)
    return ("{" in s) or s.startswith("f\"") or s.startswith("f'")


def scenario_to_name(scenario: dict) -> str:
    country = scenario["country"]
    tech = scenario["tech"]
    variant = scenario["variant"]
    year = int(scenario["year"])
    name_from_cfg = scenario.get("name")
    if is_bad_fstring_name(name_from_cfg):
        return f"{country}_{tech}_{variant}_{year}"
    return str(name_from_cfg)


def load_atlite_cf_with_time(country: str, tech: str, variant: str, year: int) -> pd.Series:
    """
    Returns a pandas Series indexed by timestamp with values in [0,1] (CF).
    """
    cc = country.lower()
    if variant == "avg":
        filename = f"{cc}_cf_{year}.csv"
    elif variant == "bestsite_p95":
        filename = f"{cc}_cf_{year}_bestsite_p95.csv"
    else:
        raise ValueError(f"Unknown variant: {variant}")

    path = DATA_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Could not find Atlite CF file: {path}")

    df = pd.read_csv(path, parse_dates=["time"])
    col = f"{tech}_cf"
    if col not in df.columns:
        raise KeyError(f"Column {col} not found in {filename}")

    s = pd.Series(df[col].astype(float).to_numpy(), index=pd.to_datetime(df["time"]))
    if not np.isfinite(s.to_numpy()).all():
        raise ValueError("CF profile contains NaN or inf values.")
    if s.max() <= 0:
        raise ValueError("CF profile has non-positive maximum.")
    return s


def read_summary_csv(scenario_name: str) -> pd.Series:
    """
    Reads the 1-row summary CSV produced by run_res_to_h2_optimisation.py.
    Returns a Series (row).
    """
    p = RESULTS_DIR / f"{scenario_name}_res_to_h2_summary.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"Missing summary CSV for scenario '{scenario_name}'. Expected: {p}\n"
            "Run run_res_to_h2_optimisation.py first (with summary export enabled)."
        )
    df = pd.read_csv(p)
    if len(df) != 1:
        raise ValueError(f"Expected exactly 1 row in {p}, found {len(df)}")
    return df.iloc[0]


def get_summary_field(row: pd.Series, candidates: list[str], required: bool = True):
    for c in candidates:
        if c in row.index and pd.notna(row[c]):
            return row[c]
    if required:
        raise KeyError(f"None of these columns were found (or all NaN): {candidates}")
    return np.nan


def make_power_series(cf: pd.Series, res_capacity_mw: float) -> pd.Series:
    return cf * float(res_capacity_mw)


def plot_window(
    *,
    scenario_name: str,
    window_label: str,
    cf_window: pd.Series,
    electrolyser_mw: float,
    res_capacity_mw: float,
    ratio: float,
    design_label: str,
    out_path: Path,
):
    """
    Plots P_RES(t), P_used(t), and P_EL line for a given window,
    plus deficit duration + cumulative deficit energy plots.
    """

    # ---- compute power series ----
    p_res = cf_window * float(res_capacity_mw)          # MW, indexed by time
    p_el = float(electrolyser_mw)                       # MW constant
    p_used = np.minimum(p_res.to_numpy(), p_el)         # MW array

    # --------------------------------------------------
    # Deficit diagnostics (hours vs energy) for THIS WINDOW
    # --------------------------------------------------
    timestep_hours = 1.0  # your CF is hourly

    deficit_mw = np.maximum(p_el - p_res.to_numpy(), 0.0)     # MW array
    deficit_mwh = deficit_mw * timestep_hours                 # MWh array

    duration_deficit_mw = np.sort(deficit_mw)[::-1]           # MW sorted
    cum_deficit_mwh = np.cumsum(np.sort(deficit_mwh)[::-1])   # MWh cumulative

    # NOTE: for a *window* plot, annual baseload is misleading.
    # Use window baseload so the dashed line is comparable within the window.
    window_baseload_mwh = p_el * len(p_res) * timestep_hours
    allowed_deficit_mwh = 0.05 * window_baseload_mwh

    # -----------------------------
    # 1) Time-series plot
    # -----------------------------
    plt.figure(figsize=(12, 4))

    plt.plot(p_res.index, p_res.to_numpy(), label="P_RES(t)")
    plt.plot(p_res.index, p_used, label="P_used(t)")
    plt.axhline(p_el, linestyle="--", label="P_EL (electrolyser MW)")

    if "day" in window_label.lower():
        hours = p_res.index.hour
        plt.xticks(
            ticks=p_res.index[::3],
            labels=[f"{h:02d}" for h in hours[::3]],
        )
        plt.xlabel("Hour of day")
    else:
        plt.xticks(rotation=45)
        plt.xlabel("Time")

    plt.title(
        f"{scenario_name} | {design_label} | {window_label}\n"
        f"ratio={ratio:.3f} | RES={res_capacity_mw:,.0f} MW | EL={p_el:,.0f} MW"
    )
    plt.ylabel("Power (MW)")
    plt.legend(loc="best")
    plt.grid(True)
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()

    # -----------------------------
    # 2) Deficit duration curve
    # -----------------------------
    dur_path = out_path.with_name(out_path.stem + "__deficit_duration.png")
    plt.figure(figsize=(6, 4))
    plt.plot(duration_deficit_mw)
    plt.xlabel("Hours (sorted)")
    plt.ylabel("Power deficit (MW)")
    plt.title(f"{scenario_name} | {design_label} | deficit duration\n{window_label}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(dur_path, dpi=200)
    plt.close()

    # -----------------------------
    # 3) Cumulative deficit energy
    # -----------------------------
    cum_path = out_path.with_name(out_path.stem + "__cum_deficit.png")
    plt.figure(figsize=(6, 4))
    plt.plot(cum_deficit_mwh, label="Cumulative deficit (MWh)")
    plt.axhline(allowed_deficit_mwh, linestyle="--", label="Allowed deficit (5% of window baseload)")
    plt.xlabel("Hours (sorted)")
    plt.ylabel("Energy deficit (MWh)")
    plt.title(f"{scenario_name} | {design_label} | cumulative deficit\n{window_label}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(cum_path, dpi=200)
    plt.close()
# ----------------------------
# Main
# ----------------------------
def main():
    print("Loading config from:", CONFIG_PATH.resolve())
    cfg = load_config(CONFIG_PATH)

    defaults = cfg.get("defaults", {})
    scenarios = cfg.get("scenarios", [])
    if not isinstance(scenarios, list) or len(scenarios) == 0:
        raise ValueError("config_hannah.yaml must contain a non-empty list under `scenarios:`")

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Fixed diagnostic windows (label by dates to avoid hemisphere confusion)
    jan_week_start = pd.Timestamp("2023-01-15 00:00:00")
    jan_week_end = pd.Timestamp("2023-01-21 23:00:00")
    jan_day_start = pd.Timestamp("2023-01-18 00:00:00")
    jan_day_end = pd.Timestamp("2023-01-18 23:00:00")

    jul_week_start = pd.Timestamp("2023-07-01 00:00:00")
    jul_week_end = pd.Timestamp("2023-07-07 23:00:00")
    jul_day_start = pd.Timestamp("2023-07-04 00:00:00")
    jul_day_end = pd.Timestamp("2023-07-04 23:00:00")

    for i, scenario_raw in enumerate(scenarios, start=1):
        scenario = deep_merge(defaults, scenario_raw)
        scenario_name = scenario_to_name(scenario)
        country = scenario["country"]
        tech = scenario["tech"]
        variant = scenario["variant"]
        year = int(scenario["year"])

        print("\n" + "=" * 80)
        print(f"Plotting scenario {i}/{len(scenarios)}: {scenario_name}")
        print("=" * 80)

        # Load CF with timestamps
        cf = load_atlite_cf_with_time(country=country, tech=tech, variant=variant, year=year)

        # Load summary (design points)
        summary = read_summary_csv(scenario_name)

        # Handle both old and new column naming (robust)
        electrolyser_mw = float(
            get_summary_field(summary, ["electrolyser_mw", "electrolyser_mw_mw"], required=True)
        )

        best_ratio = float(
            get_summary_field(summary, ["best_ratio_res_per_el", "best_ratio"], required=True)
        )

        firm_achievable = bool(
            get_summary_field(summary, ["firm_achievable"], required=False)
        ) if "firm_achievable" in summary.index else False

        firm_ratio_feasible = get_summary_field(
            summary,
            ["firm_ratio_feasible_res_per_el", "firm_ratio_feasible"],
            required=False,
        )

        # Decide which design points to plot
        design_points: list[tuple[str, float]] = [("cost_optimal", best_ratio)]
        if firm_achievable and np.isfinite(firm_ratio_feasible):
            design_points.append(("firm95", float(firm_ratio_feasible)))
        else:
            print("Note: firm95 not achievable (or missing); plotting BEST only.")

        # Build windows
        windows = [
            ("Jan week (2023-01-15..21)", cf.loc[jan_week_start:jan_week_end]),
            ("Jan day (2023-01-18)", cf.loc[jan_day_start:jan_day_end]),
            ("Jul week (2023-07-01..07)", cf.loc[jul_week_start:jul_week_end]),
            ("Jul day (2023-07-04)", cf.loc[jul_day_start:jul_day_end]),
        ]

        # Plot each design point × window
        for design_label, ratio in design_points:
            res_capacity_mw = ratio * electrolyser_mw

            for window_label, cf_window in windows:
                if cf_window.empty:
                    print(f"Warning: empty CF window for {scenario_name} | {window_label}. Skipping.")
                    continue

                safe_window = (
                    "jan_week" if "Jan week" in window_label else
                    "jan_day" if "Jan day" in window_label else
                    "jul_week" if "Jul week" in window_label else
                    "jul_day"
                )

                out_name = f"{scenario_name}__{design_label}__{safe_window}.png"
                out_path = PLOTS_DIR / out_name

                plot_window(
                    scenario_name=scenario_name,
                    window_label=window_label,
                    cf_window=cf_window,
                    electrolyser_mw=electrolyser_mw,
                    res_capacity_mw=res_capacity_mw,
                    ratio=ratio,
                    design_label=design_label,
                    out_path=out_path,
                )

                print("Saved:", out_path)

    print("\nDone. Plots in:", PLOTS_DIR.resolve())


if __name__ == "__main__":
    main()