# Two gaps from the September test run, filled

The 13–14 September run left an Ontario-shaped hole in the grid cases and a MOE
turndown sensitivity that existed only islanded. Both are now run. 48 networks:
22 Canadian grid on 2024, 26 MOE turndown on grid over the thirteen market areas.

Provenance is marked inline, as in `testrun-2026-09-13.md`: *(measured)* off the
solved networks, *(derived)* arithmetic on them, *(sourced)* from a document,
*(unestablished)* outside what this run can settle.

---

## 1. Canada on 2024

Ontario's HOEP series ends 2025-04-30, where the Market Renewal Program replaced
it with locational marginal prices that nothing here reads. 2024 is the last full
year, so the Canadian grid case is a 2024 scenario. Alberta rides along on the
same window so the two provinces are read on one year; Alberta's 2025 case stays
in `standard-grid`.

**The price is HOEP alone, and that is the thing to know before reading any
number below.** HOEP is commodity-only and excludes the Global Adjustment, which
is most of what an Ontario consumer actually pays. 2024 averages 20.6 EUR/MWh
*(measured, off the IESO file)* against 89.3 for Germany *(measured)*. That gap is
the missing charge, not cheap power. So the Ontario LCOS here is a wholesale
energy component, **not a delivered cost**, and does not belong in a ranking
against the all-in prices every other area carries. Adding a GA term is what
would make it comparable; whether Class A or Class B is the right basis is a real
question, because Class A is billed on contribution to system peak and would
therefore interact with the flexibility the model is optimising.

LCOS EUR/t, domestic routes *(measured)*:

| area | ew-eaf | h2-dri-eaf | mix/ng-dri-eaf | moe-eaf |
|---|---|---|---|---|
| ONT | **635** | 785 | 702 | 670 |
| AB  | **673** | 818 | 712 | 724 |

Export twins run 677–733 (ONT) and 733–913 (AB).

Ontario coming out below Alberta on every route is what the missing Global
Adjustment buys it, so the ordering between the two provinces carries no
information *(unestablished)*. Alberta's own 2024 number against its 2025 one
(673 vs 638 for `ew-eaf`) is a like-for-like year comparison and does mean
something.

Emission intensity is blank for both provinces, as designed: AESO publishes no
per-carrier generation and IESO's needs a separate XML report, so the retrieve
rule serves the price alone and says so.

---

## 2. MOE turndown, on grid

`moe-turndown-70` carried no `grid` tech row, so it resolved to the seven
islanded areas and the grid half was never run. Nothing failed — that is the same
silent shape as the missing overlay file, and worth noting as a second instance:
a scenario's *tech rows* decide whether it is a grid case at all, and getting
that wrong produces a smaller correct-looking run rather than an error.

All 26 pairs have `inputs_hash` values distinct from their `standard-grid` twins
*(measured)*, so the overlay applied.

**On grid, dropping the MOE cell from 0.95 to 0.70 is worth almost nothing:
0.0 to 6.6 EUR/t, 0.0–0.9%** *(measured)*, against **55–131 EUR/t, 6–12%,
islanded** *(measured, 13 September)*.

| area | base 0.95 | turndown 0.70 | Δ | | area | base 0.95 | turndown 0.70 | Δ |
|---|---|---|---|---|---|---|---|---|
| TAS1 | 769.4 | 762.7 | −6.6 | | VIC1 | 738.0 | 736.7 | −1.3 |
| AB   | 685.3 | 681.6 | −3.7 | | QLD1 | 726.8 | 726.2 | −0.6 |
| SA1  | 715.9 | 712.5 | −3.4 | | FRA  | 766.6 | 766.3 | −0.3 |
| NSW1 | 743.0 | 740.8 | −2.3 | | BR_S, BR_SE | | | −0.1 |
| DEU  | 827.2 | 825.1 | −2.1 | | BR_N, BR_NE | | | −0.0 |
| ESP  | 797.1 | 795.2 | −1.9 | | | | | |

### Why, established rather than assumed

The cost decomposition separates the two cases cleanly *(measured, MEUR/yr,
turndown minus base)*:

| | renewables | battery | total |
|---|---|---|---|
| islanded, mean over 7 areas | **−48.5** | −29.6 to −88.5 | −55 to −131 |
| grid, mean over 13 areas | **+1.7** | −0.0 to −4.9 | −0.03 to −6.7 |

Islanded, the inflexible cell is what forces the renewable fleet and the battery
to be oversized: holding 95% through every calm hour has to be paid for in
capacity that sits idle the rest of the year. Letting the cell fall to 70% lets
the optimiser shrink both, and that is the whole of the saving — renewables fall
by 48.5 MEUR/yr on average and the battery by 30–88.

On grid the connection already provides that firmness, so there is nothing to
shrink. Renewables barely move and change sign by area; what little saving there
is comes mostly from a smaller battery, partly offset by *more* renewable and
grid spend. The flexibility has no scarcity left to be valuable against.

**A first guess that did not survive.** The obvious reading — that on grid the
cell would earn its turndown by following the hourly price — is not supported.
Correlation between the saving and each area's price volatility is +0.43 on
price standard deviation and +0.33 on the p95−p5 spread across 13 areas
*(measured)*, which with that n is not a relationship. South Australia has the
widest price spread of any area (σ = 228 EUR/MWh) and saves 3.4 EUR/t; Tasmania
has half that spread and saves the most. Whatever orders these thirteen numbers,
it is not price volatility alone *(unestablished)*.

### The steel store, again

`cost_steel_store_meur` changes by +0.00 to +0.05 MEUR/yr in every one of the 26
grid runs, as it did islanded. The standing note that this sensitivity is
meaningless without the finished-steel store is wrong in both configurations: the
store is never built, and islanded the saving arrives anyway by a different
route. The grid result does not rescue that note either — it makes the
sensitivity nearly worthless on grid for an unrelated reason.

---

## What this run does not settle

- No Global Adjustment, so no Ontario delivered price and no defensible
  comparison of Ontario against any other area.
- Alberta's capacity factors still come from `scratch/cf_for_landlocked_area.py`,
  not the pipeline: `bestsite_p95` computes all three technologies on every call
  and `pick_p95_cell` fails on an all-NaN distance array for an area with no
  offshore geometry, taking the well-defined onshore and solar outputs with it.
- The four Brazilian submarkets are four hardlinks to one whole-Brazil cutout,
  placed by hand and outside the keyed cache, so a fresh worktree that rebuilds
  them goes to CDS. Hit once while setting this run up.
- `make_area_geometry` reads Natural Earth Admin 0, so it silently produces the
  whole country for a top-level sub-national area — Alberta came out as all of
  Canada in a fresh worktree. Sub-national shapes come from `scratch/cut_*.py`
  and the committed parquets; only AB and ONT lack a committed copy. Issue #64.
