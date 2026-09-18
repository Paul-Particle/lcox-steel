# Data in this repository

The repository LICENSE covers the code. It does not cover anything described here:
these files came from third parties on their terms, and those terms are not ours to
pass on.

## `resources/timeseries/ESP_grid_emissions_20250101_20251231.parquet`

Hourly generation by carrier for the Spanish bidding zone (ES), calendar year 2025,
8760 rows. Committed because the demo's `-export` routes melt their iron in Spain and
price it against Spain's own market, so without this file the repo's default target
cannot be solved after a fresh clone.

- **Source:** ENTSO-E Transparency Platform, <https://transparency.entsoe.eu>,
  retrieved through the platform's RESTful API.
- **Data item:** Aggregated Generation per Type, Article 16.1(b) & (c) of Commission
  Regulation (EU) No 543/2013.
- **Form:** a reshaped and aggregated derivative, not a verbatim API response.

**Licence status: re-used under the Terms of Use, not open data.** ENTSO-E publishes a
"List of Data Available for Free Re-Use" under CC-BY 4.0, and that list covers only an
enumerated subset of the platform — load forecasts, unavailability, cross-zonal
capacity, physical flows, balancing. It contains no Article 16 data, so this series is
**not** CC-BY and must not be labelled as such. It is redistributed under clause 3.1 of
the Terms of Use, which permits re-use subject to naming the platform as the source.
Any copyright or related right in the underlying data may be held by the Primary Owner
of Data, Red Eléctrica de España.

ENTSO-E does not endorse this work and exercises no control over the accuracy of what
is published here.

Two things that mislead people about this, both worth knowing before repeating them:
the widespread claim that ENTSO-E generation data is CC-BY comes from secondary
summaries rather than the list itself; and the restrictive "non-commercial, personal
use" wording often quoted is from the separate entsoe.eu corporate disclaimer, which
does not govern the Transparency Platform. There is no commercial-use restriction in
the platform's terms.

If this repository ever ships a broad multi-country, multi-year cache rather than this
one series, the "small slice" argument weakens and the question is worth putting to
ENTSO-E in writing.

## `data/assumptions/Assumptions_yaml_inputs.xlsx`

The techno-economic inputs `config/assumptions.yaml` is generated from, by
`scratch/generate_assumptions.py`. One row per value: the dotted YAML key it
writes to, the value, its unit, the key it replaced, and the row ID it was
confirmed against.

- **Source:** compiled by Future Cleantech Architects from approximately 150
  third-party sources. This file is a single extract tab from a larger internal
  workbook where each figure is attributed to its own source; the
  `Confirmed row ID` column is the reference back into it.
- **Form:** a compilation. The individual figures are facts drawn from their
  respective publishers; the selection and arrangement are ours.

**Licence status: not established per-source here.** The underlying publications
carry their own terms, which are recorded in the master workbook rather than in
this repository. Before quoting a figure externally, check its source entry
there. Nothing in this file should be treated as open data on the strength of
being committed.

## Everything else under `data/`

The market caches (`entsoe_cache/`, `nem_cache/`, `ons_cache/`, `canada_cache/`) and the
downloaded shapes are gitignored. They are local working copies of third-party data,
are not redistributed here, and are rebuilt from their sources on demand.
