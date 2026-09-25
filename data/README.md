# Data in this repository

The repository LICENSE covers the code. It does not cover the third-party files
listed here: they came from their publishers on those publishers' terms, and every
one of them asks to be credited. Each entry says where a file came from and what
the credit is.

No market price or generation series is shipped. The demo's export twins melt their
iron in `ES_SYN`, a Spain-like market that `workflow/scripts/grid/_synthetic.py`
makes up on the first run. Its numbers are invented: they describe no real market
and must not be quoted as Spain's.

## `cutouts/VIC1_20250101_20251231_backup.nc`

A Victoria (AUS) slice of ERA5 hourly reanalysis for 2025, prepared by atlite, so
the demo needs no CDS download.

- **Source:** Hersbach, H. et al. (2023), ERA5 hourly data on single levels from 1940
  to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS).
  DOI: [10.24381/cds.adbb2d47](https://doi.org/10.24381/cds.adbb2d47).
- **Licence:** [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). The whole
  Climate Data Store has been under it since 2 July 2025.
- **Credit:** Contains modified Copernicus Climate Change Service information (2025).
  Neither the European Commission nor ECMWF is responsible for any use made of the
  Copernicus information or data it contains.
- **Changes:** cut to Victoria's bounding box, with variables derived by atlite.

## `resources/shapes/*_geo.parquet` and `*_offshore_geo.parquet`

Area geometries, committed because building them needs the Natural Earth and EEZ
zips that the repository does not ship.

- **Onshore** (`{AREA}_geo.parquet`): countries come from Natural Earth 1:110m
  Admin 0. NEM regions, Canadian provinces and Brazilian submarkets are dissolved
  from Natural Earth 1:10m Admin 1 states and provinces
  (`scratch/cut_canada_provinces.py`, `scratch/cut_brazil_submarkets.py`). Natural
  Earth is [public domain](https://www.naturalearthdata.com/about/terms-of-use/),
  so no credit is required, but it is given here anyway.
- **Offshore** (`{AREA}_offshore_geo.parquet`): the area's EEZ, minus the onshore
  geometry above and clipped to a coastal buffer. The EEZ comes from Flanders
  Marine Institute (2023), Maritime Boundaries Geodatabase:
  Maritime Boundaries and Exclusive Economic Zones (200NM), version 12,
  <https://www.marineregions.org/>, DOI:
  [10.14284/632](https://doi.org/10.14284/632). Licence:
  [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: clipped to
  one area's zone and coastal buffer.

## `data/nem_cache/NEM Registration and Exemption List.xlsx`

An unmodified snapshot of the Australian Energy Market Operator's (AEMO) register of
NEM participants and generating units, committed because AEMO's own hosting is
unreliable.

- **Source and author:** AEMO,
  <https://www.aemo.com.au/-/media/Files/Electricity/NEM/Participant_Information/NEM-Registration-and-Exemption-List.xls>.
- **Terms:** AEMO's [copyright permissions](https://www.aemo.com.au/privacy-and-legal-notices/copyright-permissions)
  give "general permission for anyone to use AEMO Material for any purpose, but only
  with accurate and appropriate attribution of the relevant AEMO Material and AEMO as
  its author."

## `data/entsoe_cache/entsoe_bidding_zones.csv`

A hand-maintained list of bidding-zone codes, written for this repository and
covered by its LICENSE. It holds no data from the ENTSO-E Transparency Platform.

## Everything else under `data/`

The market caches (`entsoe_cache/`, `nem_cache/`, `ons_cache/`, `canada_cache/`) and the
downloaded shapes are ignored by version control. They are local working copies of
third-party data, are not redistributed here, and are rebuilt from their sources on
demand. Whoever downloads them is bound by those sources' terms. For ENTSO-E, only
the items on its CC-BY 4.0 "List of Data Available for Free Re-Use" may be
redistributed freely, and day-ahead prices and generation per type are not among
them.
