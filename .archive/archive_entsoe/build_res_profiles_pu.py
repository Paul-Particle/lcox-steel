from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

# Model scope
AREAS = ["DE_LU", "FR", "ES"]
RES_TECHS = ["solar", "wind_onshore", "wind_offshore"]

GEN_PATH = ROOT / "data" / "processed_data.feather"
CAPS_PATH = ROOT / "data" / "installed_capacities_latest.csv"
OUT_PATH = ROOT / "data" / "res_profiles_pu.feather"


def main():
    if not GEN_PATH.exists():
        raise FileNotFoundError(f"Missing generation file: {GEN_PATH}")
    if not CAPS_PATH.exists():
        raise FileNotFoundError(f"Missing capacities file: {CAPS_PATH}")

    gen = pd.read_feather(GEN_PATH)
    if not isinstance(gen.columns, pd.MultiIndex) or gen.columns.nlevels != 2:
        raise ValueError("Expected generation columns as MultiIndex (area, tech).")

    caps = pd.read_csv(CAPS_PATH).set_index("area")

    # Check required techs exist in caps file
    missing_caps_cols = [t for t in RES_TECHS if t not in caps.columns]
    if missing_caps_cols:
        raise ValueError(f"Capacities CSV missing columns: {missing_caps_cols}")

    out = {}
    for area in AREAS:
        if area not in gen.columns.get_level_values(0):
            raise ValueError(f"Area {area} not found in generation data.")
        if area not in caps.index:
            raise ValueError(f"Area {area} not found in capacities CSV.")

        for tech in RES_TECHS:
            col = (area, tech)
            if col not in gen.columns:
                raise ValueError(f"Missing generation column {col} in processed_data.feather")

            installed_mw = float(caps.loc[area, tech])
            series_mw = gen[col].astype(float)

            # Guard for technologies with zero installed capacity (e.g. offshore ES)
            if installed_mw <= 0:
                out[col] = 0.0 * series_mw
            else:
                out[col] = series_mw / installed_mw

    df_pu = pd.DataFrame(out, index=gen.index)

    # Optional: keep it tidy
    df_pu.columns = pd.MultiIndex.from_tuples(df_pu.columns, names=["area", "tech"])
    df_pu = df_pu.sort_index(axis=1)

    # -----------------------------------------------------------------------------
# Cap extreme per-unit RES outliers (soft guardrail)
#
# We normalise generation (MW) by installed capacity (MW) to obtain per-unit
# profiles. Because installed capacity is taken as a single annual snapshot
# (here: 2024) while generation spans real-world hourly data, occasional
# inconsistencies can arise (e.g. commissioning within the year, reporting
# artefacts, or category mismatches in ENTSO-E data).
#
# To prevent rare, non-representative spikes from dominating optimisation
# results (e.g. electrolyser sizing), we softly cap per-unit values at a high
# percentile (99.9%) per (area, technology). This preserves the full profile
# shape while removing extreme tails.
# -----------------------------------------------------------------------------

    CAP_PCT = 0.999  # 99.9th percentile cap per (area, tech)

    for col in df_pu.columns:
        cap = df_pu[col].quantile(CAP_PCT)
        if pd.notna(cap) and cap > 0:
            df_pu[col] = df_pu[col].clip(lower=0, upper=cap)


    df_pu.to_feather(OUT_PATH)
    print(f"Saved per-unit RES profiles to: {OUT_PATH}")

    # Quick sanity print
    summary = (
        df_pu
        .stack("area", future_stack=True)
        .describe(percentiles=[0.5])
        .loc[["min", "50%", "max"]]
    )
    print("\nPer-unit profile summary (min/median/max):")
    print(summary.round(3))


if __name__ == "__main__":
    main()