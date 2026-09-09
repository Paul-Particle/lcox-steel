# Grid costs: modelling practice and the policy/physics separation

Desk research for the WP4 grid-cost component. Two questions:

1. How do PyPSA-based models treat grid costs?
2. Is there a defensible way to separate policy effects from technological/physical
   effects in what a plant actually pays?

Status: literature/desk survey. Nothing here is wired into `assumptions.yaml` yet.

---

## 0. Where this model stands today

`build_network._add_grid_import()` adds an unconstrained generator at
`marginal_cost = price_series`, where `price_series` is the raw ENTSO-E day-ahead
or NEM spot price. That is the **wholesale energy component only**. A real
grid-connected plant pays wholesale + network charges + taxes/levies, minus
whatever exemptions it qualifies for.

So the current grid pathway is not "the grid pathway" — it is "the grid pathway
for a consumer facing zero network charges and zero levies", which is one specific
(and in Germany, briefly real — see §118(6) EnWG below) policy corner.

---

## 1. How PyPSA models handle grid costs

There are two entirely different things called "grid cost", and the literature
splits cleanly along that line.

### 1.1 System models: grid cost is endogenous capacity expansion

PyPSA-Eur / PyPSA-DE / PyPSA-Earth build transmission as an investment decision.
Line and link extensions carry an annualised capex in EUR/MW/km from the
`PyPSA/technology-data` repository (mostly Danish Energy Agency figures), scaled
by length and a route factor, and the optimiser trades grid build-out against
generation and storage siting. There is no tariff in the model: the cost of the
grid is a term in the objective function, and its dual (the nodal price spread)
is the economic signal.

Consequences worth knowing:

- **No distribution grid.** PyPSA-DE states explicitly that "a detailed
  representation of the electricity distribution grid is outside the scope of the
  model", while noting distribution investment needs are of a magnitude similar to
  transmission. Anything connecting below the transmission level is therefore
  under-costed by construction.
- **No taxes, no levies, no cross-subsidy.** The objective is resource cost. There
  is no VAT, no renewables surcharge, no exemption regime, no regulated return
  above WACC.
- **Tariffs, where they appear at all, are a post-processing division.** PyPSA-DE
  reports "average grid tariffs" falling by 7.5 EUR2020/MWh under integrated
  planning versus the national grid development plan — derived from total
  annualised transmission cost (191 bn EUR2020 in their run) over demand. It is a
  reported quantity, not a modelled charge, and no consumer in the model responds
  to it.

### 1.2 Project / LCOX models: grid cost is an exogenous adder

Project-scale tools (model.energy and its descendants, the green-H2 and green-iron
cost literature) either model an **islanded** plant with no grid at all, or a
**grid-connected** plant whose electricity price is wholesale plus a flat
EUR/MWh network charge. That adder is almost always a single number taken from
statistics rather than derived.

Representative figure: Hydrogen Europe's 2024 *Clean Hydrogen Production Pathways*
report uses 80 EUR/MWh wholesale + **29.3 EUR/MWh average network fees** for
grid-connected alkaline electrolysis, and the same 29.3 EUR/MWh on top of a
60 EUR/MWh PPA when the PPA is wheeled over the grid. Direct coupling avoids the
adder entirely, which is quoted as worth 10-30% of production cost.

This is the pattern our model would follow. The honest framing is that the adder
is a **scenario lever**, not a datum.

### 1.3 The transmission-cost number itself is contested

`technology-data` issue #125 tabulates AC overhead transmission cost assumptions,
and the spread is the single most useful thing in this memo:

| Source | AC overhead, EUR/MW/km |
|---|---|
| CIGRE (min) | 130 |
| CIGRE (max) | 250 |
| ACER 2023 average | 371 |
| **PyPSA-Eur** | **442** |
| ACER 2023 Q3 | 481 |
| NEP 2021 (DE) | 736 |
| NEP 2023 (DE) | 1325-1384 |

Same technology, same continent, a factor of ~10 end to end and a factor of 3
between PyPSA-Eur and the German network development plan. DC underground is
worse: 870-1159 EUR/MW/km in the academic estimates versus 3300-3800 in NEP 2023.

Note this against our own `transmission.cost_per_mw_per_km_eur: 400`, commented
"HVDC overhead, order-of-magnitude (PyPSA-Eur-ish)" — that sits on the low,
engineering-cost end of the range. Defensible, but it should be labelled as a
*physical* cost estimate rather than a *German-project* cost estimate.

---

## 2. Separating policy from physics

Short answer: **yes, but only approximately, and the clean separations are
institutional rather than statistical.** Four workable approaches, in descending
order of rigour-per-effort for our purposes.

### 2.1 The component decomposition (standard, cheap, leaky)

Eurostat (Regulation (EU) 2016/1952, series `nrg_pc_205_c`) splits every reported
price into three components:

- `DP_ES` — energy and supply
- `DP_NC` — network costs
- `DP_TL` — taxes and levies

with a further breakdown of taxes and levies by purpose (renewables/CHP support,
nuclear, social, security of supply, energy efficiency, environmental and excise,
VAT). The Ecofys/Fraunhofer *Prices and costs of EU energy* study builds its whole
analysis on this taxonomy, and adds the distinction that matters most here:

- **Direct impact** — taxes, levies and certificate obligations that change the
  retail price by construction.
- **Indirect impact** — policies that change the generation mix, supply routes or
  market design, and therefore move the wholesale price. Same policy, invisible in
  the tax column.

The naive reading — "energy + network = physics, taxes + levies = policy" — is
wrong in both directions, and the report says so:

- Network charges are *not* pure physics. Renewables connection costs land in
  distribution networks but "are often covered by special levies, or they are
  included in renewable energy surcharges", so the same physical cost appears in
  different columns in different countries. And the allocation *between* consumer
  classes is a pure policy choice (below).
- The energy component is not policy-free either — that is the entire
  merit-order-effect literature.

So the decomposition is the right *reporting* frame and a bad *causal* frame. Use
it for transparency, not for attribution.

### 2.2 The regulatory separation: LRMC vs residual (best available)

This is the strongest line, because the split is written into the tariff-setting
rules themselves, in both our jurisdictions.

Network revenue requirements are based on **embedded (sunk) costs**, which are
systematically higher than **long-run marginal cost**. Efficient pricing says
charge LRMC; cost recovery says charge more. The gap is the **residual**, and — in
the words of the tariff-design literature — residual costs are sunk, no change in
user behaviour affects them, and **how the residual is allocated across load,
generation and storage is mostly a policy choice**.

That gives a defensible two-way split of the network charge:

```
network charge  =  LRMC component      (forward-looking, physical, behaviour-responsive)
                +  residual component  (sunk recovery, allocation is policy)
```

Australia institutionalises exactly this. The National Electricity Rules require
network businesses to set tariffs based on the long-run marginal cost of the
service, then adjust to recover the residual subject to stand-alone and avoidable
cost bounds on each tariff class — a bounded Ramsey adjustment. AER-published LRMC
studies (e.g. the Energeia LRMC report) give the forward-looking number directly,
so for a NEM connection point the LRMC and the residual are separately observable,
at least at tariff-class granularity.

Germany does not publish the split, but the same decomposition is conceptually
available and the residual share is large.

**This is the recommendation for the model**: parameterise the grid charge as
`lrmc_component + residual_component` and sweep the residual, rather than sweeping
one opaque EUR/MWh number. It makes the sensitivity say something.

### 2.3 The counterfactual / natural-experiment approach (strong where it applies)

Where a cost driver is a discrete, dated political decision, the counterfactual is
credible and the effect is directly attributable. The cleanest example available
is the one already visible in the table above.

In December 2015 the German cabinet mandated **underground cabling priority**
(*Erdkabelvorrang*) for new HVDC corridors, on political rather than technical
grounds — a concession to Bavarian opposition to overhead lines. Underground
transmission costs several times overhead depending on ground conditions; the
official additional-cost estimate at the time was 3-8 bn EUR and proved
substantially optimistic. SuedLink alone is now around 10 bn EUR.

That decision is a large part of why NEP 2023 says 1325-1384 EUR/MW/km where CIGRE
says 130-250. **The gap between an engineering cost database and a national
network development plan is a first-order estimate of the policy premium on grid
infrastructure.** It is crude, it conflates permitting, compensation, undergrounding
mandates and route factors, and it says nothing about which of those a modeller
should treat as avoidable. But it is a real number with a real mechanism behind it,
and it bounds the question.

The equivalent method on the energy component is the merit-order-effect
literature's synthetic-supply counterfactual: rebuild the supply curve without the
policy-supported capacity and re-clear. Doable, but a project in its own right.

### 2.4 The exemption stack (decisive for *this* plant, and pure policy)

For an H2-DRI plant the biggest single driver is not the tariff level — it is
whether the plant pays it at all. Every item below is a policy instrument with an
expiry date, and none of them is a physical cost.

**Germany:**

- **§118(6) EnWG** — electrolysers commissioned before **4 August 2029** are exempt
  from network charges for **20 years**. The Bundesnetzagentur has stated that a
  full exemption is not tenable under EU law and not expedient in energy-policy
  terms, and is consulting on replacing it with conditional reduced charges. So
  the exemption is real today, contested, and dated.
- **§19(2) StromNEV** — two industrial privileges:
  - *Atypical grid use*: an individual charge no lower than 20% of the published
    one; reductions up to ~90%.
  - *Bandlastprivileg* (7000-hour rule): >10 GWh/yr consumption and >7000 full-load
    hours qualifies for a heavily reduced individual charge.

  Combined, these were worth >1 bn EUR/yr across ~400 band-load and ~4200 atypical
  customers, and the regime is slated for replacement.
- **Levy and tax relief**: electricity-tax reductions for industry (~1.7 bn EUR in
  2023), offshore/CHP levy reductions (~1 bn EUR), indirect ETS cost compensation
  (~3 bn EUR). Energy-intensive industry is 18% of EU energy consumption but 2% of
  energy-tax revenue.
- **2026 transmission-charge subsidy**: a 6.5 bn EUR federal subsidy cuts the
  headline TSO charge sharply for 2026, contingent on legislation. A one-year,
  budget-financed price. *(Reported headline figures for this vary a lot between
  voltage levels and sources — the industrial HV energy charge is on the order of
  0.4-0.7 ct/kWh plus a capacity component, an order of magnitude below the
  low-voltage pass-through figure. Do not mix the two. Check against the TSOs' own
  2026 price sheets before use.)*

Directly relevant to us: `plant.availability_target: 0.95` means >8300 full-load
hours, so a modelled H2-DRI plant **comfortably clears the 7000-hour Bandlast
threshold**. The model's own operating assumption places it in the privileged
class.

**RFNBO rules** are the other policy lever on the grid pathway, and they act on
*eligibility* rather than price: additionality, and monthly temporal correlation
tightening to **hourly from 1 January 2030**. Modelled cost impacts of strict
temporal correlation run up to ~6.6 EUR/kg H2, with hourly matching forcing
oversized wind+electrolysis configurations. For a grid-connected green-steel
plant, this is the difference between "buy cheap grid power" and "buy hour-matched
certified power", and it is entirely regulatory.

**Australia:** no equivalent electrolyser network-charge exemption. Large loads
near 40 GWh/yr get individually-set locational TUOS at their connection point.
The DE/AUS contrast is therefore not just a resource contrast — it is a large
policy-regime contrast sitting inside the same EUR/MWh number.

---

## 3. What I'd propose for the model

1. **Split the grid-import price into named terms** rather than one adder:
   `wholesale (from ENTSO-E/NEM) + network_lrmc + network_residual + levies_and_taxes`,
   each defaulting to a documented source, each sweepable. The sum is what enters
   `marginal_cost`; the split is what makes the sensitivity interpretable.
2. **Make the exemption regime an explicit scenario dimension**, not a value baked
   into the network term. At minimum three points: full published tariff /
   Bandlast-privileged / §118(6)-exempt. These are discrete legal states, so a
   discrete scenario axis is the honest representation; a continuous sweep over
   EUR/MWh hides the fact that the intermediate values do not exist.
3. **Keep `transmission.cost_per_mw_per_km_eur` explicitly labelled as an
   engineering cost**, and add a policy-premium multiplier (1.0 for the CIGRE/
   PyPSA-Eur view, ~3 for the NEP view) rather than silently picking one. The
   issue #125 table is the citation.
4. **Report the decomposition in the outputs.** If the grid pathway's LCOS is
   quoted, quote it with the component stack visible, because the headline number
   is dominated by which policy corner was assumed.

Open question for Pauline: whether WP4 wants the *observed* charge (what a plant
pays today, exemptions included) or the *cost-reflective* charge (LRMC, no
exemptions). These differ by more than the technology uncertainty does, and the
answer determines whether the grid pathway is being compared to the RES pathway on
a like-for-like basis.

---

## Sources

- PyPSA/technology-data issue #125, transmission cost comparison table —
  https://github.com/PyPSA/technology-data/issues/125
- PyPSA-DE: Open-source German energy system model reveals savings from integrated
  planning — https://arxiv.org/abs/2510.09414
- PyPSA-Eur documentation — https://pypsa-eur.readthedocs.io/
- Ariadne PyPSA model documentation — https://ariadneprojekt.de/en/model-documentation-pypsa/
- Hydrogen Europe, Clean Hydrogen Production Pathways Report 2024 —
  https://hydrogeneurope.eu/wp-content/uploads/2024/06/2024_H2E_CleanH2ProductionPathwaysReport.pdf
- Eurostat, Electricity price statistics (components methodology) —
  https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Electricity_price_statistics
- Ecofys/Fraunhofer ISI, Prices and costs of EU energy, Final Report (2016) —
  https://www.isi.fraunhofer.de/content/dam/isi/dokumente/ccx/2016/report_ecofys2016.pdf
- Bruegel, Europe's under-the-radar industrial policy: intervention in electricity
  pricing — https://www.bruegel.org/policy-brief/europes-under-radar-industrial-policy-intervention-electricity-pricing
- smartEn / FTI Consulting, A Roadmap for Cost-Reflective Electricity Network
  Tariffs in the EU (2025) —
  https://smarten.eu/wp-content/uploads/2025/03/FTI-Consulting-Report_smartEn_03-2025_DIGITAL_V2.pdf
- Simshauser, Efficient tariff structures for distribution network services —
  https://www.sciencedirect.com/science/article/abs/pii/S0313592615300552
- AER / Energeia, Long-Run Marginal Cost Final Report —
  https://www.aer.gov.au/system/files/PWC%20-%2011.05%20-%20Energeia%20-%20LRMC%20Report%20-%2031%20Jan%202023%20-%20Public.PDF
- AEMC, Economic Concepts for Pricing Electricity Network Services (NERA) —
  https://www.aemc.gov.au/sites/default/files/content/f2475394-d9f6-497d-b5f0-8d59dabf5e1c/NERA-Economic-Consulting-%E2%80%93-Network-pricing-report.PDF
- AEMO, Pricing methodology for prescribed shared transmission services —
  https://www.aemo.com.au/-/media/files/electricity/nem/participant_information/fees/2023/revised-pricing-methodology-for-1-july-2022-to-30-june-2027.pdf
- GvW, Charges, levies, taxes and fees for electricity for hydrogen production in
  electrolysers (§118(6) EnWG) —
  https://www.gvw.com/en/news/blog/detail/umlagen-abgaben-steuern-und-entgelte-bei-dem-strombezug-fuer-die-wasserstofferzeugung-in-elektrolyseuren
- Gleiss Lutz, Bundesnetzagentur discussion paper on industrial network charges —
  https://www.gleisslutz.com/en/know-how/shape-industrial-network-charges-come-bundesnetzagentur-publishes-discussion-paper-general-electricity-framework
- Amprion, Surcharge for special grid utilisation (§19 StromNEV levy) —
  https://www.amprion.net/Market/Levies/Surcharge-for-special-grid-utilisation/
- Effects of the delegated act on RFNBO on production costs and distribution of
  hydrogen production capacities across Europe —
  https://www.sciencedirect.com/science/article/pii/S0301421526001126
- Hydrogen in the European power sector — impacts of regulatory frameworks for
  green hydrogen — https://www.sciencedirect.com/science/article/pii/S0301421526001345
- Estimating the merit-order effect using coarsened exact matching (counterfactual
  synthetic supply) — https://www.sciencedirect.com/science/article/pii/S0301421523005165
