import http.client
import io
import ssl
import zipfile
from datetime import datetime
from pathlib import Path

import nemosis.downloader as _dl
import pandas as pd
from nemosis import dynamic_data_compiler

# Patch 1: AEMO now blocks NEMOSIS's stale Chrome 80 User-Agent with HTTP 403.
_dl.USR_AGENT_HEADER.clear()
_dl.USR_AGENT_HEADER["User-Agent"] = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
)

# Patch 2: MMSDM filenames from Aug-2024 contain literal '#' characters.
# requests percent-encodes '#' → '%23', which AEMO's Azure endpoint rejects (HTTP 400).
# Sending a literal '#' also fails (proxies strip it as a fragment). The only working
# encoding is double-encoding '%2523'; Azure decodes one layer to '%23' as expected.
# requests would re-normalise '%2523' → '%23', so we bypass it entirely with http.client.
def download_unzip_csv_patched(url: str, down_load_to: str) -> None:
    from urllib.parse import urlsplit
    url_fixed = url.replace("#", "%2523").replace("%23", "%2523")
    u = urlsplit(url_fixed)
    path = u.path + (f"?{u.query}" if u.query else "")
    conn = http.client.HTTPSConnection(
        u.hostname, u.port or 443, context=ssl.create_default_context(), timeout=120 # pyright:ignore
    )
    try:
        conn.request("GET", path,
                     headers={**_dl.USR_AGENT_HEADER, "Host": u.hostname, "Accept": "*/*"}) # pyright: ignore
        resp = conn.getresponse()
        if resp.status != 200:
            raise IOError(f"GET {url} -> {resp.status} {resp.reason}")
        body = resp.read()
    finally:
        conn.close()
    zipfile.ZipFile(io.BytesIO(body)).extractall(down_load_to)

_dl.download_unzip_csv = download_unzip_csv_patched

if "snakemake" not in globals():
    from common._stubs import snakemake


def download_price_data(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads price data or gets it from the cached feather files or csv files if rebuild=True."""
    print("Fetching prices...")
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
    print(f"Generator excel not found; downloading via NEMOSIS to {target}...")
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


def download_generation_data(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads and processes generation data."""
    generator_file_path = _resolve_generator_excel(cache_dir)

    generator_info = pd.read_excel(
        generator_file_path, sheet_name="PU and Scheduled Loads", engine="openpyxl"
    )
    generator_info.columns = generator_info.columns.str.strip().str.replace("\n", " ")

    
    def determine_gen_type(df):
        dfi = df['gen_info']
        type_regex = {
            'hard_coal' : 'black coal',
            'brown_coal' : 'brown coal',
            'coal_gas' : 'coal seam methane|coal mine gas',
            'waste' : 'methane',
            'oil' : 'diesel|ethane|kerosene',
            'biomass' : 'bagasse|biomass|biogas',
            'gas' : 'natural gas|natrual gas',
            'solar' : 'solar',
            'energy_storage' : 'battery',
            'pumped_storage' : 'pump storage|^- -$',
            'other_re' : 'sewerage',
            'wind_onshore' : 'wind - onshore|wind',
            'hydro_river' : 'run of river',
            'hydro' : 'hydro',
        }
        for t, r in type_regex.items():
            m = dfi.str.contains(r)
            df.loc[m, 'gen_type'] = t
            df.loc[m, 'gen_info_temp'] = df.loc[m, 'gen_info']
            df.loc[m, 'gen_info'] = ''
        df = df.drop(columns=['gen_info'])
        df.columns = ['gen_type', 'gen_info']
        return df.loc[:,['gen_info', 'gen_type']]


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


    print("Fetching generation...")
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


def download_load_data(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads and processes load data."""
    print("Fetching load...")
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

def download_crossborder_data(start_time: str, end_time: str, cache_dir: Path, rebuild: bool) -> pd.DataFrame:
    """Downloads and processes cross-border flow data."""
    print("Fetching cross-border flows...")
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

    # convert fromatting to match ENTSO-E
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

    crossborder = (interconnector.copy()
        .assign(flow_type_for_region    = lambda df: df.METEREDMWFLOW.apply(lambda v: 'to_'   if v >= 0 else 'from_') + df.to_region  )
        .assign(flow_type_for_other     = lambda df: df.METEREDMWFLOW.apply(lambda v: 'from_' if v >= 0 else 'to_'  ) + df.from_region)
        .assign(flow_type_r_for_region  = lambda df: df.METEREDMWFLOW.apply(lambda v: 'from_' if v >= 0 else 'to_'  ) + df.to_region  )
        .assign(flow_type_r_for_other   = lambda df: df.METEREDMWFLOW.apply(lambda v: 'to_'   if v >= 0 else 'from_') + df.from_region)
        .assign(flow_value_for_region   = lambda df: df.METEREDMWFLOW.apply(lambda v: v       if v >= 0 else -v))
        .assign(flow_value_for_other    = lambda df: df.METEREDMWFLOW.apply(lambda v: v       if v >= 0 else -v))
        .assign(flow_value_r_for_region = lambda df: df.METEREDMWFLOW.apply(lambda v: 0.0))
        .assign(flow_value_r_for_other  = lambda df: df.METEREDMWFLOW.apply(lambda v: 0.0))

    )
    common_names = {
        'from_region':'region',
        'to_region':'region',
        'flow_type_for_region': 'flow_type',
        'flow_type_for_other': 'flow_type',
        'flow_type_r_for_region': 'flow_type',
        'flow_type_r_for_other': 'flow_type',
        'flow_value_for_region': 'flow_value',
        'flow_value_for_other': 'flow_value',
        'flow_value_r_for_region': 'flow_value',
        'flow_value_r_for_other': 'flow_value',
    }

    cols_a = ['SETTLEMENTDATE', 'from_region','flow_type_for_region', 'flow_value_for_region', 'INTERCONNECTORID', 'METEREDMWFLOW']
    df_a = crossborder[cols_a].rename(columns=common_names)

    cols_b = ['SETTLEMENTDATE', 'to_region', 'flow_type_for_other', 'flow_value_for_other',  'INTERCONNECTORID',  'METEREDMWFLOW']
    df_b = crossborder[cols_b].rename(columns=common_names)

    cols_c = ['SETTLEMENTDATE', 'from_region', 'flow_type_r_for_region', 'flow_value_r_for_region',  'INTERCONNECTORID',  'METEREDMWFLOW']
    df_c = crossborder[cols_c].rename(columns=common_names)

    cols_d = ['SETTLEMENTDATE', 'to_region', 'flow_type_r_for_other', 'flow_value_r_for_other',  'INTERCONNECTORID',  'METEREDMWFLOW']
    df_d = crossborder[cols_d].rename(columns=common_names)

    crossborder = pd.concat([df_a, df_b, df_c, df_d], ignore_index=True)
    crossborder = crossborder.pivot_table(index='SETTLEMENTDATE', values='flow_value', columns=['region', 'flow_type'], aggfunc='sum')
    return crossborder


def process_area(snakemake) -> None:
    """Build one NEM region's price (hourly) + full (native 5-min) parquets.

    NEMOSIS caches the national MMSDM tables under data/nem_cache; refresh by
    deleting those files. Dates/area come from output-path wildcards.
    """
    cache_dir = Path("data/nem_cache")
    area = snakemake.wildcards.area
    start_time = (
        datetime.strptime(snakemake.wildcards.start_date, "%Y%m%d").strftime("%Y/%m/%d")
        + " 00:00:00"
    )
    end_time = (
        datetime.strptime(snakemake.wildcards.end_date, "%Y%m%d").strftime("%Y/%m/%d")
        + " 23:59:59"
    )

    prices      = download_price_data(      start_time, end_time, cache_dir, rebuild=False)
    generation  = download_generation_data( start_time, end_time, cache_dir, rebuild=False)
    load        = download_load_data(       start_time, end_time, cache_dir, rebuild=False)
    crossborder = download_crossborder_data(start_time, end_time, cache_dir, rebuild=False)

    df_combined = pd.concat([prices, load, generation, crossborder], axis=1)
    if df_combined.index.tz is not None:
        df_combined.index = df_combined.index.tz_localize(None)

    # Single-region rich frame (native 5-min): this region's metrics + derived residual.
    load_s = df_combined[(area, "load")]
    wind = df_combined.get((area, "wind_onshore"), pd.Series(0, index=df_combined.index))
    solar = df_combined.get((area, "solar"), pd.Series(0, index=df_combined.index))
    derived = pd.DataFrame(index=df_combined.index)
    derived["wind"] = wind  # NEM has no offshore, so wind == onshore
    derived["residual"] = load_s - (wind + solar)
    full = pd.concat([df_combined[area], derived], axis=1)
    full.index.name = "time"

    # Price (model input): hourly mean of the 5-min RRP, converted AUD→EUR so it
    # matches the EUR-denominated ENTSO-E prices the optimisation consumes.
    # (The _full frame stays raw AUD.)
    eur_per_aud = snakemake.params.eur_per_aud
    price = (df_combined[(area, "price")].resample("1h").mean() * eur_per_aud).rename("price")
    price.index.name = "time"

    full_out = Path(snakemake.output.full)
    full_out.parent.mkdir(parents=True, exist_ok=True)
    full.to_parquet(full_out, index=True)
    price.to_frame().to_parquet(snakemake.output.prices, index=True)
    print(f"Wrote {snakemake.output.prices} ({len(price)} hourly prices, "
          f"{price.index.min()} .. {price.index.max()})")
    print(f"Wrote {full_out} ({full.shape[0]} rows × {full.shape[1]} cols)")


if __name__ == "__main__":
    process_area(snakemake)
