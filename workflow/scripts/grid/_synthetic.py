"""Synthetic source implementation, called by retrieve_grid_data.py.

Makes up an hourly, Spain-like market series for the demo, so an export furnace
can be priced after a fresh clone with no API key and no third-party data in the
repository. The shapes are plausible, not fitted: a solar day scaled by season,
wind as a persistent random walk, fixed nuclear with a refuelling dip, gas
filling what a load profile leaves over, and a price that follows the gas share.
Nothing about it describes a real market.

Variants
--------
dayahead   "price" only
emissions  "price" and one column per carrier, named as the ENTSO-E downloader
           names them so the emission factor table reads it unchanged
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from _helpers_grid import iso

# Module-level logger only — retrieve_grid_data.py installs the handlers.
log = logging.getLogger(__name__)

# Fixed, so a regenerated file solves to the same demo numbers.
RANDOM_SEED = 2025

# Carriers the synthetic zone has none of, kept so its columns match a real
# `emissions` series one for one.
ABSENT_CARRIERS = ["brown_coal", "coal_gas", "oil_shale", "peat", "geothermal",
                   "marine", "wind_offshore", "energy_storage"]


def retrieve(snakemake, market_area: str) -> None:
    """Make the requested (variant, date range) window up and write it out."""
    area       = snakemake.wildcards.area
    variant    = snakemake.wildcards.variant
    start_date = snakemake.wildcards.start_date
    end_date   = snakemake.wildcards.end_date
    if variant not in ("dayahead", "emissions"):
        raise ValueError(
            f"the synthetic market makes `dayahead` and `emissions` series only, "
            f"not {variant!r}: it has no load or cross-border flows to put in `full`."
        )

    rng = np.random.default_rng(RANDOM_SEED)
    time = pd.date_range(iso(start_date), f"{iso(end_date)} 23:00", freq="h", name="time")
    local_hour = (time.hour + 1) % 24
    season = np.cos(2 * np.pi * (time.dayofyear - 172) / 365)  # +1 midsummer, -1 midwinter

    day_length = 12 + 2.5 * season
    sunrise = 13.5 - day_length / 2
    daylight = np.clip(np.sin(np.pi * (local_hour - sunrise) / day_length), 0, None)
    daily_cloud = np.repeat(rng.uniform(0.55, 1.0, len(time) // 24 + 1), 24)[:len(time)]
    solar = 24000 * daylight * (0.75 + 0.25 * season) * daily_cloud

    wind_walk = np.zeros(len(time))
    for i in range(1, len(time)):
        wind_walk[i] = 0.97 * wind_walk[i - 1] + rng.normal(0, 0.25)
    wind_onshore = 19000 / (1 + np.exp(-(wind_walk - 0.3 * season - 0.2)))

    load = (28000 + 3500 * np.sin(np.pi * np.clip(local_hour - 7, 0, 16) / 16)
            + 2500 * np.abs(season))
    refuelling = np.isin(time.month, [4, 10])
    nuclear = np.where(refuelling, 5000.0, 7000.0)
    must_run = pd.DataFrame({
        "biomass": rng.normal(420, 40, len(time)),
        "hard_coal": np.clip(rng.normal(150, 120, len(time)), 0, None),
        "oil": rng.normal(28, 8, len(time)),
        "hydro_river": 950 - 350 * season,
        "hydro_reservoir": 2500 + 1500 * np.sin(np.pi * np.clip(local_hour - 8, 0, 14) / 14),
        "pumped_storage": 1800 * np.clip(np.sin(np.pi * (local_hour - 18) / 5), 0, None),
        "nuclear": nuclear,
        "other": rng.normal(7, 3, len(time)),
        "other_re": rng.normal(69, 9, len(time)),
        "solar": solar,
        "waste": rng.normal(190, 30, len(time)),
        "wind_onshore": wind_onshore,
    }, index=time).clip(lower=0.0)
    gas = np.clip(load - must_run.sum(axis=1), 1500, 16500)
    price = np.clip(0.011 * (gas - 1500) + 10 + rng.normal(0, 8, len(time)), -5, 240)

    out_df = pd.DataFrame({"price": price}, index=time)
    if variant == "emissions":
        out_df = pd.concat(
            [out_df, must_run.assign(gas=gas, **{carrier: 0.0 for carrier in ABSENT_CARRIERS})],
            axis=1,
        )

    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out_path, index=True)
    log.info(f"wrote synthetic {area} {out_path} ({len(out_df)} rows × {out_df.shape[1]} cols)")
