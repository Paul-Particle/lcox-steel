"""NEM download primitives — imported by retrieve_nem.py.

Provides four download functions (price, generation, load, crossborder) that
pull data from NEMOSIS and return DataFrames in a shape compatible with their
ENTSO-E counterparts in download_entsoe.py.

The AEMO compatibility patches in _nemosis_patches.py must be applied before
any NEMOSIS call; they are activated by the `import _nemosis_patches` line
below, which must remain the first import in this file.
"""

import _nemosis_patches  # noqa: F401  — applies AEMO compatibility patches on import; must be first
import logging
from pathlib import Path

import pandas as pd
from nemosis import dynamic_data_compiler

# Module-level logger only — the rule script (retrieve_nem.py) installs handlers.
log = logging.getLogger(__name__)


def download_price(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads price data or gets it from the cached feather files or csv files if rebuild=True."""
    log.info("fetching prices")
    prices = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCHPRICE",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "REGIONID", "RRP"], # regional reference price in Australian Dollar
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    prices["SETTLEMENTDATE"] = pd.to_datetime(prices["SETTLEMENTDATE"])
    prices = prices.pivot_table(
        index="SETTLEMENTDATE", columns="REGIONID", values="RRP"
    )
    prices.columns = pd.MultiIndex.from_tuples([(col, 'price') for col in prices.columns])
    prices.index.name = None
    return prices


def _resolve_generator_excel(cache_dir: Path) -> Path:
    """Return path to the NEM Registration and Exemption List, fetching via NEMOSIS if absent.

    A snapshot is committed to the repo as .xlsx; NEMOSIS saves it as .xls (its
    historical filename) even though AEMO now serves XLSX content. Either extension
    is accepted and read with engine='openpyxl' downstream.
    """
    stem = "NEM Registration and Exemption List"
    for ext in (".xlsx", ".xls"):
        candidate = cache_dir / f"{stem}{ext}"
        if candidate.exists():
            return candidate

    from nemosis import defaults as _nemosis_defaults
    from nemosis.data_fetch_methods import static_downloader_map

    table_name = "Generators and Scheduled Loads"
    target = cache_dir / _nemosis_defaults.names[table_name]
    url = _nemosis_defaults.static_table_url[table_name]
    log.info(f"generator excel not found; downloading via NEMOSIS to {target}")
    cache_dir.mkdir(parents=True, exist_ok=True)
    try:
        static_downloader_map[table_name](url, str(target))
    except Exception as e:
        raise FileNotFoundError(
            f"NEMOSIS download of {stem!r} failed: {e}. "
            f"Fetch manually from {url} (served as XLSX) and place it in "
            f"{cache_dir} as .xls or .xlsx — see README §External data files."
        ) from e
    return target


def download_generation(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads and processes generation data."""
    generator_file_path = _resolve_generator_excel(cache_dir)

    generator_info = pd.read_excel(
        generator_file_path, sheet_name="PU and Scheduled Loads", engine="openpyxl"
    )
    generator_info.columns = generator_info.columns.str.strip().str.replace("\n", " ")

    def determine_gen_type(df):
        """Assign a generator type to each row using ordered regex patterns.

        First match wins. Dict order is intentional — see inline comments.
        The Excel contains typos in fuel/technology strings, so some patterns
        intentionally cover misspellings (e.g. 'natrual gas').
        """
        # Dict order is load-bearing: first match wins.
        # Ordering rules (each intentional for emissions classification):
        #   - gas and biomass before oil: natural-gas/diesel → gas; bagasse+diesel → biomass
        #   - oil before coal_gas and waste: coal-seam-methane and landfill-methane classify
        #     as oil (the 'ethane' substring of 'methane' matches the oil regex on purpose —
        #     these generators are treated as dirty/leaky for emissions estimation)
        #   - other_re before pumped_storage: sewerage burning sludge wins over pump storage
        #   - wind_onshore before energy_storage: co-located wind+battery treated as wind
        type_regex = {
            'hard_coal'      : 'black coal',
            'brown_coal'     : 'brown coal',
            'gas'            : 'natural gas|natrual gas',
            'biomass'        : 'bagasse|biomass|biogas',
            'oil'            : 'diesel|ethane|kerosene',
            'coal_gas'       : 'coal seam methane|coal mine gas',
            'waste'          : 'methane',
            'solar'          : 'solar',
            'other_re'       : 'sewerage',
            'wind_onshore'   : 'wind - onshore|wind',
            'energy_storage' : 'battery',
            'pumped_storage' : 'pump storage|^- -$',
            'hydro_river'    : 'run of river',
            'hydro'          : 'hydro',
        }
        gen_info = df['gen_info']
        gen_type = pd.Series('', index=df.index, dtype=str)
        matched = pd.Series(False, index=df.index)
        for t, pattern in type_regex.items():
            m = (~matched) & gen_info.str.contains(pattern, na=False)
            gen_type[m] = t
            matched[m] = True
        return pd.DataFrame({'gen_info': gen_info, 'gen_type': gen_type})

    nem_gen_names = dict(
        generator_info.copy()
        .loc[
            :, ["DUID", "Fuel Source - Descriptor", "Technology Type - Descriptor"]
        ]  # DUID = dispatchable unit ID
        .dropna()
        .pivot_table(
            index=["Fuel Source - Descriptor"],
            columns=["Technology Type - Descriptor"],
            values="DUID",
            aggfunc="count",
            fill_value=0,
        )
        .stack()  # getting all combinations of Fuel Source and Technology Type
        .loc[lambda s: s > 0]  # filter for existing combinations
        .reset_index()  # put the columns back into the df
        .assign(
            gen_info=lambda df: (
                df["Fuel Source - Descriptor"].str.lower()
                + " "
                + df["Technology Type - Descriptor"].str.lower()
            )
        )
        .drop(
            columns=["Technology Type - Descriptor", "Fuel Source - Descriptor", 0]
        )  # 0 column is number of plants for each combination, i.e. not relevant
        .pipe(determine_gen_type)  # apply regex-based categorization
        .values
    )
    # Ensure that there are no nan values in the dict, to catch changes to future versions of the excel.
    # This can happen if a new 'Fuel Source - Descriptor' and 'Technology Type - Descriptor' combination
    # is not covered by the regex in `determine_gen_type`.
    assert not any(pd.isna(v) for v in nem_gen_names.values()), "NaN values found in nem_gen_names mapping. Please update the regex in 'determine_gen_type'."

    # final mapping of generator ID to generator type
    map_df = (generator_info.copy()
        .loc[:, ["DUID", "Region", "Fuel Source - Descriptor", "Technology Type - Descriptor"]]
        .assign(gen_info = lambda df: df["Fuel Source - Descriptor"].str.lower() + ' ' + df["Technology Type - Descriptor"].str.lower())
        .drop(columns=["Technology Type - Descriptor", "Fuel Source - Descriptor"])
        .assign(gen_type = lambda df: df["gen_info"].map(nem_gen_names))
    )

    log.info("fetching generation")
    unit_generation = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCH_UNIT_SCADA",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "DUID", "SCADAVALUE"], # SCADAVALUE = generation in MW
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    unit_generation["SETTLEMENTDATE"] = pd.to_datetime(unit_generation["SETTLEMENTDATE"])

    generation = (
        unit_generation.merge(map_df.dropna(subset=['gen_type']), on="DUID", how="left")
        .groupby(["SETTLEMENTDATE", "Region", "gen_type"])["SCADAVALUE"]
        .sum()
        .reset_index()
        .pivot_table(
            index="SETTLEMENTDATE",
            columns=["Region", "gen_type"],
            values="SCADAVALUE",
        )
        .fillna(0)
    )
    generation = generation.clip(lower=0)
    generation.index.name = None
    return generation


def download_load(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads and processes load data."""
    log.info("fetching load")
    load = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCHREGIONSUM",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "REGIONID", "TOTALDEMAND"],
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    load["SETTLEMENTDATE"] = pd.to_datetime(load["SETTLEMENTDATE"])

    # Pivot load
    load = load.pivot_table(
        index="SETTLEMENTDATE", columns="REGIONID", values="TOTALDEMAND"
    )
    load.columns = pd.MultiIndex.from_tuples([(col, 'load') for col in load.columns])
    load.index.name = None
    return load


def download_crossborder(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads and processes cross-border flow data.

    NEM records one signed METEREDMWFLOW value per interconnector:
      positive → flow from `from_region` to `to_region`
      negative → reverse

    ENTSO-E records two unsigned values per border from each area's perspective:
      to_{neighbour}   (exports)
      from_{neighbour} (imports)
    and in rare cases both can be non-zero simultaneously (bilateral trading).

    To produce an output with the same shape as the ENTSO-E crossborder table —
    so that VIC1 and DE_LU can be compared side-by-side — each NEM interconnector
    row is expanded into four rows, one per (region, direction) perspective:

      (a) from_region exporting: to_{to_region}      = |flow|  (main direction)
      (b) to_region importing:   from_{from_region}  = |flow|  (same physical flow)
      (c) from_region importing: from_{to_region}    = 0       (reverse simultaneous,
      (d) to_region exporting:   to_{from_region}    = 0        always 0 in NEM)

    The zero-valued rows (c, d) exist so that any aggregation or comparison code
    that expects both directions per border (as in ENTSO-E) finds the expected
    columns rather than NaN.
    """
    log.info("fetching cross-border flows")
    interconnector = dynamic_data_compiler(
        start_time,
        end_time,
        "DISPATCHINTERCONNECTORRES",
        cache_dir,
        select_columns=["SETTLEMENTDATE", "INTERCONNECTORID", "METEREDMWFLOW", "MWFLOW"],
        fformat="feather",
        keep_csv=True,
        rebuild=rebuild,
    )
    interconnector["SETTLEMENTDATE"] = pd.to_datetime(interconnector["SETTLEMENTDATE"])

    # Normalise INTERCONNECTORID names to match ENTSO-E border naming convention
    ic_names = {
        k: v
        for k, v in zip(
            ["N-Q-MNSP1", "NSW1-QLD1", "T-V-MNSP1", "V-S-MNSP1", "V-SA", "VIC1-NSW1"],
            ["NSW1-QLD1", "NSW1-QLD1", "TAS1-VIC1", "VIC1-SA1", "VIC1-SA1", "VIC1-NSW1"],
        )
    }

    regions = interconnector['INTERCONNECTORID'].map(ic_names).str.split('-', expand=True)
    interconnector['from_region'] = regions[0]
    interconnector['to_region'] = regions[1]

    # Build the four flow columns (see docstring for the (a)-(d) mapping)
    flow = interconnector.METEREDMWFLOW
    crossborder = interconnector.assign(
        # (a) from_region perspective: exporting to to_region
        flow_type_for_region    = flow.apply(lambda v: 'to_'   if v >= 0 else 'from_') + interconnector.to_region,
        # (b) to_region perspective: importing from from_region
        flow_type_for_other     = flow.apply(lambda v: 'from_' if v >= 0 else 'to_'  ) + interconnector.from_region,
        # (c) from_region perspective: reverse simultaneous import = 0
        flow_type_r_for_region  = flow.apply(lambda v: 'from_' if v >= 0 else 'to_'  ) + interconnector.to_region,
        # (d) to_region perspective: reverse simultaneous export = 0
        flow_type_r_for_other   = flow.apply(lambda v: 'to_'   if v >= 0 else 'from_') + interconnector.from_region,
        flow_value_for_region   = flow.abs(),
        flow_value_for_other    = flow.abs(),
        flow_value_r_for_region = 0.0,   # NEM has no concurrent reverse flows
        flow_value_r_for_other  = 0.0,
    )

    # Rename columns to a common schema so all four perspectives can be stacked
    common_names = {
        'from_region':           'region',
        'to_region':             'region',
        'flow_type_for_region':  'flow_type',
        'flow_type_for_other':   'flow_type',
        'flow_type_r_for_region':'flow_type',
        'flow_type_r_for_other': 'flow_type',
        'flow_value_for_region': 'flow_value',
        'flow_value_for_other':  'flow_value',
        'flow_value_r_for_region':'flow_value',
        'flow_value_r_for_other': 'flow_value',
    }

    base_cols = ['SETTLEMENTDATE', 'INTERCONNECTORID', 'METEREDMWFLOW']
    perspectives = [
        ('from_region', 'flow_type_for_region',   'flow_value_for_region'),    # (a)
        ('to_region',   'flow_type_for_other',    'flow_value_for_other'),     # (b)
        ('from_region', 'flow_type_r_for_region', 'flow_value_r_for_region'),  # (c) zeros
        ('to_region',   'flow_type_r_for_other',  'flow_value_r_for_other'),   # (d) zeros
    ]
    parts = [
        crossborder[base_cols + [region_col, flow_type_col, flow_value_col]]
        .rename(columns=common_names)
        for region_col, flow_type_col, flow_value_col in perspectives
    ]

    crossborder = pd.concat(parts, ignore_index=True)
    crossborder = crossborder.pivot_table(index='SETTLEMENTDATE', values='flow_value', columns=['region', 'flow_type'], aggfunc='sum')
    return crossborder


DOWNLOADERS = {
    "price": download_price,
    "generation": download_generation,
    "load": download_load,
    "crossborder": download_crossborder,
}
