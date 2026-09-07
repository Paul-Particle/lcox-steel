# EU ETS treatment of an iron and steel installation

Research date: 7 September 2026. All figures and quotes verified against primary EU
texts (EUR-Lex, Commission guidance, DEHSt, ecologie.gouv.fr) unless flagged otherwise.

**Headline caveat on currency of sources.** Three things changed after the 2023 revision
that most secondary write-ups still get wrong, and all three matter for your model:

1. Annex I activity 5 no longer says *pig* iron — a hydrogen-DRI or electrolytic iron
   plant **is** an ETS installation, even with zero direct emissions.
2. From 1 January 2026 the "exchangeability of fuel and electricity" correction was
   **deleted** (Art. 22 of Reg. 2019/331). EU free allocation for EAF/DRI steel now
   covers part of the plant's *electricity* emissions.
3. DRI (sponge iron) was folded into the **hot metal benchmark** at 1.248 tCO2e per
   **tonne of DRI**, while the CBAM free-allocation adjustment for imported DRI is only
   0.295–0.397 tCO2e/t. That asymmetry is the whole answer to question 8.

---

## 1. SCOPE — which emissions are covered

### 1.1 Direct emissions only, as a matter of law

The ETS covers emissions "from sources in an installation". Electricity consumed by an
installation is not an emission of that installation; it is an emission of the generator,
which is itself an ETS installation. There is no scope-2 obligation anywhere in
Directive 2003/87/EC.

Directive (EU) 2023/959 did change Article 3(b): the words "into the atmosphere" were
deleted from the definition of 'emissions'. The Commission's own guidance says this
does **not** widen installation boundaries:

> "The deletion of 'into the atmosphere' in the definition does not impact whether an
> installation should be included in the EU ETS or the installation boundaries. However,
> it paves the way to a consistent treatment of CO2 transfer and carbon capture and
> utilisation (CCU) activities (Article 12(3b) of the EU ETS Directive)."
> — *Guidance on Interpretation of Annex I of the EU ETS Directive*, §2.3.5
> https://climate.ec.europa.eu/document/download/edc93136-82a0-482c-bf47-39ecaf13b318_en?filename=policy_ets_gd0_annex_i_euets_directive_en.pdf

The mirror-image confirmation on the free-allocation side, for the EAF carbon steel
benchmark (pre-2026 wording, Guidance Document 9):

> "Emissions related to the production of the consumed electricity are excluded from the
> system boundaries." (hot metal)
> "Allocation should however be based on direct emissions only."  (EAF carbon steel)
> — *Sector-specific guidance* GD9, §4 and §5,
> https://climate.ec.europa.eu/system/files/2019-07/p4_gd9_sector_specific_guidance_en.pdf

### 1.2 Which sources count for an EAF

Benchmark system boundary (Delegated Reg. (EU) 2024/873 amending 2019/331, Annex I),
**EAF carbon steel**:

> "All processes directly or indirectly linked to the process units electric arc furnace,
> secondary metallurgy, casting and cutting, post-combustion unit, dedusting unit,
> vessels heating stands, casting ingots preheating stands, scrap drying and scrap
> preheating are included. Processes downstream of casting are not included."
> https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=OJ:L_202400873

Reportable source streams inside that boundary (MRR, Reg. (EU) 2018/2066, plus GD9):

| Source | Type | Typical magnitude |
|---|---|---|
| Graphite **electrode consumption** (oxidised C) | process | 1–2 kg electrode/t steel → ~4–7 kg CO2/t |
| **Carbon/coal injection** and charge carbon | process | ~10–15 kg C/t → ~40–55 kg CO2/t |
| **Carbon-bearing charge materials** (pig iron, HBI/DRI carbon, cast iron, carburiser) | process | route-dependent |
| **Carbonates / fluxes** — limestone, dolomite, calcined at temperature | process | ~10–30 kg CO2/t |
| **Natural gas / oxy-fuel burners**, ladle and tundish preheating, scrap drying/preheating | combustion | ~10–40 kg CO2/t |
| Post-combustion of CO from the furnace off-gas | combustion | counted in the above |

GD9 makes the accounting basis explicit: the EAF benchmark accounts for "carbon from
electrodes and scrap that is oxidised in the electric arc furnace".

Scrap itself is **not** a fuel and carries no embedded charge — note that CN 7204
(ferrous waste and scrap) is expressly listed in CBAM Annex II with the ferro-alloy
exclusions, and scrap-based steel gets the lowest benchmarks in the system.

### 1.3 Which sources count for a DRI shaft

Hot metal benchmark system boundary as amended by 2024/873 (emphasis on the new items):

> "All processes directly or indirectly linked to the process units blast furnace, hot
> metal treatment units, blast furnace blowers, blast furnace hot stoves, **direct reduced
> iron reactor, electric arc furnace and electric smelting furnace for sponge iron**,
> basic oxygen furnace, secondary metallurgy units, vacuum ladles, casting units
> (including cutting), slag treatment unit, burden preparation, BF and other gas treatment
> units, dedusting units, scrap pre-heating, coal drying for PCI, vessels preheating
> stands, casting ingots preheating stands, compressed air production, dust treatment unit
> (briquetting), sludge treatment unit (briquetting), steam injection in BF unit, steam
> generation plant, converter BOF gas cooling and miscellaneous are included."

Sources that count for a shaft furnace:
- **Reducing gas carbon** — for a natural-gas MIDREX/Energiron shaft, the CH4/CO in the
  reformed gas that ends up as CO2 in the top gas. This is the dominant term (~0.5–0.65
  tCO2/t DRI for gas-based).
- **Fuels** — reformer burners, top-gas fuel, shaft preheating, HBI briquetting machine
  drives that burn fuel.
- **Carbonates** — limestone/dolomite in the pellet or added as flux, and any calcination.
- **Carbon retained in the DRI** is *not* an emission at the shaft; it is emitted later in
  the EAF and is captured there as a carbon-bearing charge material.

For a **100 % hydrogen** shaft, the reducing gas carries no carbon. Residual direct
emissions are confined to any residual natural gas used for preheating/trim, flux
calcination, and the briquetting line. In the limit this is close to zero — which is
exactly why the "is it in the ETS at all" question in section 2 matters.

### 1.4 Direct emission intensities — sourced numbers

**Regulatory numbers (top-decile, direct only)** from Commission Implementing Regulation
(EU) 2025/2620, Annex point 5 (CBAM benchmarks). Column A is the *process step only*;
Column B is *cumulative including precursors*. These are derived from "the average direct
emissions of the 10 % best installations under these ETS benchmarks in the new baseline
years 2021 and 2022" (recital 17).
https://eur-lex.europa.eu/eli/reg_impl/2025/2620/oj

| CN code / good | Route | Column A (process) tCO2e/t | Column B (cumulative) tCO2e/t |
|---|---|---|---|
| 7203 10 00 — DRI/sponge iron | gas-based DR | **0.295** | **0.397** |
| 7201 — pig iron | BF | 1.089 | 1.210 |
| 7206 — ingots, iron/non-alloy steel | (C) BF/BOF | 0.150 | 1.288 |
| 7206 — ingots | (D) DRI/EAF | 0.027 | 0.424 |
| 7206 — ingots | (E) Scrap/EAF | **0.027** | **0.027** |
| 7207 — semi-finished | (C) BF/BOF | 0.188 | 1.364 |
| 7207 — semi-finished | (D) DRI/EAF | 0.065 | 0.475 |
| 7207 — semi-finished | (E) Scrap/EAF | 0.065 | 0.066 |

Route legend is given verbatim in the same Annex: "(C) Carbon Steel based on BF/BOF;
(D) Carbon Steel based on DRI/EAF; (E) Carbon Steel based on Scrap/EAF; (J) High alloy
Steel (based on EAF)".

So **the top-decile EU scrap-EAF has direct emissions of ~0.027 tCO2/t of crude steel at
the ingot stage, ~0.065–0.072 tCO2/t at the hot-rolled semis stage.**

CBAM default values (IR (EU) 2025/2621, Dec 2025) for finished steel goods are quoted in
trade press as BF-BOF **1.370**, DRI-EAF **0.481**, scrap-EAF **0.072** tCO2e/t — i.e. the
Column-B ladder plus one further downstream step. (Secondary source: EUROMETAL,
https://eurometal.net/eu-commission-finalizes-cbam-benchmarks-default-values-ahead-of-january-2026-launch/ ;
cbamguide https://cbamguide.com/sectors/steel/benchmarks/ — I could not retrieve the OJ
text of 2025/2621 itself, so treat the exact figures as second-hand.)

**Industry / literature numbers (typical, not top-decile)**

| Route | Direct (scope 1) tCO2/t crude steel | Indirect (scope 2) |
|---|---|---|
| Scrap-EAF | **0.06–0.1** | ~0.4 (grid-dependent) |
| Integrated BF-BOF | 1.8–3.0 | small |
| Gas-DRI-EAF | 0.7–1.2 (total) | — |
| Coal-DRI-EAF | 2–3 | — |

Source: SteelOnTheNet, *Steelmaking CO2 emissions by process step*,
http://www.steelonthenet.com/resources/kb/co2-emissions.html (attributes the EU plant
figures to European Commission data).

**H2-DRI-EAF.** There is no single authoritative regulatory number. Reported residual
direct emissions are **0.05–0.2 tCO2/t crude steel**, arising from electrode consumption,
carburising carbon added at the EAF, flux/lime calcination, and any residual natural gas
for process heat. For EU modelling, the defensible anchors are:
- lower bound ≈ the CBAM DRI/EAF *process* benchmark 0.027 tCO2/t (EAF step only, top
  decile) plus a near-zero H2 shaft;
- a working central value of **0.05–0.10 tCO2/t** if you add carburiser and lime;
- upper bound ~0.16 tCO2/t where the literature includes lime calcination plus residual
  NG (SteelWatch/EFI-type accounting).
Cross-check: Sandbag cites JRC data giving actual EU DRI at **0.39 tCO2/t DRI** — but that
is the *gas-based* EU fleet, not H2.

For a **scrap-EAF in the EU** the safest single number for a scope-1 ETS obligation is
**0.06–0.10 tCO2/t crude steel**, with 0.027–0.072 tCO2/t being the top-decile regulatory
benchmark values.

---

## 2. IS A HYDROGEN-DRI OR ELECTROLYTIC IRON PLANT AN ETS INSTALLATION?

**Yes. Unambiguously, since 1 January 2024.** Two amendments in Directive (EU) 2023/959
did it.

### 2.1 Activity 5 no longer says "pig iron"

Directive (EU) 2023/959, Annex, point (1)(c)(iii):

> "the fifth row is replaced by the following:
> 'Production of **iron** or steel (primary or secondary fusion) including continuous
> casting, with a capacity exceeding 2,5 tonnes per hour — Carbon dioxide'"
> https://eur-lex.europa.eu/eli/dir/2023/959/oj

The Commission's Annex I guidance spells out the intent:

> "In the iron and steel sector, 'iron production' is not limited anymore to 'pig iron'.
> Therefore, production routes leading to other forms of iron, in particular sponge iron
> (also called DRI (Direct Reduced Iron) or HBI (Hot briquetted Iron)) are covered. **This
> enables inclusion (with free allocation) of steel production routes using hydrogen or
> even electrolytic iron reduction processes.**"
> — *Guidance on Interpretation of Annex I*, §3.3

That sentence covers **molten oxide electrolysis and aqueous electrowinning** explicitly.

### 2.2 A zero-emission installation is still an installation

Article 2(1) was rewritten from "This Directive shall apply to *emissions from* the
activities listed in Annexes I and III" to "This Directive shall apply to **the
activities** listed in Annexes I and III". Guidance §2.3.4, headed "New: zero-emission
installations possible":

> "This clarifies that the activity carried out in an installation (and meeting its
> threshold, if applicable) is the relevant criterion for inclusion in the EU ETS, not
> whether the installation actually emits greenhouse gases covered by the Directive."

The guidance footnote also records that this overturns the reading in ECJ case C-577/16
(*Trinseo Deutschland*), which had held that an installation with no direct emissions
fell outside the ETS. **Any secondary source relying on Trinseo is out of date.**

### 2.3 The threshold

The 2.5 t/h capacity threshold applies to the *iron or steel production* activity. There
is no 20 MW combustion condition attached to activity 5 (unlike activity 6, "production or
processing of ferrous metals", which needs combustion units >20 MW). So a hydrogen shaft
or an electrolytic iron cell house above 2.5 t/h (≈21 kt/y at full load, i.e. essentially
any commercial plant) is in.

Note the guidance's boundary between activities 5 and 6:
> "'Production of iron and steel including continuous casting' ends at iron or steel in
> primary forms (slabs, ingots etc.), and it relates to the product benchmarks 'hot
> metal', 'EAF carbon steel' or 'EAF high alloy steel'." (§3.4.2)

### 2.4 The on-site electrolyser is a *second* Annex I activity

Directive 2023/959 also replaced the 24th row:
> "'Production of hydrogen (H2) and synthesis gas with a production capacity exceeding
> **5 tonnes per day** — Carbon dioxide'"

Two changes, per Guidance §3.3:
> "The limitation to the production processes 'reforming or partial oxidation' has been
> removed. Together with the removal of the need for GHG emissions in the installation
> itself (section 2.3.4) this means that **all kinds of electrolysis processes will be
> included.** The threshold has been reduced from 25 to 5 tonnes production per day."

A 1 Mt/y H2-DRI plant consumes roughly 50–60 kg H2/t DRI ≈ 150–180 t H2/day — thirty times
the threshold. So the electrolyser is itself an Annex I activity.

**⚠ Open question worth resolving for your model.** The hydrogen product benchmark is
defined on hydrogen "exported from the sub-installation concerned … as net saleable
product", and the ammonia benchmark expressly absorbs its captive hydrogen ("All processes
directly or indirectly linked to the production of the ammonia **and the intermediate
product hydrogen** are included"). No equivalent carve-out was written for hot metal / DRI.
Whether captive H2 fed to an on-site shaft earns the hydrogen benchmark (**7.98 tCO2e/t H2**
for 2026-2030 — worth ~0.44 tCO2e per tonne of steel on top of the hot metal allocation) or
is absorbed into the hot metal sub-installation is not settled by the texts I could reach.
This is a large swing and should be checked with DEHSt/the relevant competent authority.

---

## 3. FREE ALLOCATION AND BENCHMARKS

### 3.1 Benchmark values

2021–2025: **Commission Implementing Regulation (EU) 2021/447**,
https://eur-lex.europa.eu/eli/reg_impl/2021/447/oj

2026–2030: **Commission Implementing Regulation (EU) 2026/1412** of 26 June 2026, OJ
29 June 2026, in force 30 June 2026.

| Product benchmark | Unit | 2021–2025 | 2026–2030 |
|---|---|---|---|
| Coke | tCO2e/t | 0.217 | **0.143** |
| Sintered ore → renamed **agglomerated iron ore** | tCO2e/t | 0.157 | **0.086** |
| **Hot metal** (now incl. DRI) | tCO2e/t | 1.288 | **1.248** |
| **EAF carbon steel** | tCO2e/t | 0.215 | **0.142** |
| **EAF high alloy steel** | tCO2e/t | 0.268 | **0.176** |
| Iron casting | tCO2e/t | 0.282 | **0.163** |
| **Hydrogen** | tCO2e/t | 8.85 (as amended by 2024/873) | **7.98** |
| Heat benchmark | tCO2e/TJ | 47.3 | **31.2** |
| Fuel benchmark | tCO2e/TJ | 42.6 | **28.1** |

Note the 2026–2030 table splits benchmarks into those "without electricity data
collection" (coke, agglomerated iron ore, hot metal) and those "with electricity data
collection" (EAF carbon steel, EAF high alloy steel, iron casting, hydrogen) — the tell
that the electrified benchmarks now carry an indirect component (see §3.3).

Values cross-checked against the OJ text and against cbamguide's IR 2026/1412 reference
table (https://cbamguide.com/carbon/ets-benchmarks/) and EUROMETAL
(https://eurometal.net/european-commission-publishes-delayed-lower-ferrous-ets-benchmarks/).

### 3.2 Is there a DRI benchmark? — No. DRI sits inside hot metal.

Commission Delegated Regulation (EU) **2024/873** rewrote the hot metal definition:

> "Iron produced from iron ores for primary steelmaking including (a) liquid iron
> saturated with carbon for further processing, considered as product of blast furnaces,
> and expressed in tonnes of liquid iron at the exit point of the blast furnace, excluding
> liquid iron produced from sponge iron under (b); **(b) sponge iron at the exit point of a
> direct reduced iron reactor**"

and the two EAF benchmarks now both end with:

> "**Steel produced from iron sponge already covered under the hot metal benchmark is not
> covered by this benchmark.**"

Consequences:
- **Activity level for a DRI plant is tonnes of DRI at the shaft exit**, not tonnes of
  steel. Sandbag confirms this reading and quantifies it:
  > "the 2024 reform of the Free Allocation Regulation further extended this benchmark to
  > the production of direct reduced iron (DRI), which has much lower emissions (0.39
  > tCO2/t in the EU, according to the JRC, compared to 1.248 EUA for the benchmark)
  > creating windfall profits of **0.858 EUA per tonne of DRI produced**."
  > https://sandbag.be/2026/07/07/eu-steelmaking-the-ets-money-is-coming/
- The Commission excluded DRI-route data from the 2026–2030 recalculation, so 1.248 is
  still a blast-furnace-derived number applied to DRI.
- Also new in 2024/873: pellets folded into the (renamed) agglomerated iron ore benchmark,
  and electrolysis folded into the hydrogen benchmark.

### 3.3 The 2026 change that most write-ups miss: exchangeability deleted

Under the 2021–2025 rules the EAF benchmarks were built on *total* (direct + indirect)
emissions, but allocation was scaled back to the direct share (GD9 §5):

> F(p,k) = [Em_direct + Em_NetHeatImport] / [Em_direct + Em_NetHeatImport + Em_indirect]
> × BM_p × HAL_p × CLEF(p,k),  with Em_indirect = Elec.use × **0.376** tCO2/MWh

Delegated Regulation (EU) 2024/873, Art. 1, point (11): **"Article 22 is deleted"**.
Recital 9:

> "In order to incentivise the electrification of industrial processes to significantly
> reduce emissions from such processes, it is necessary to remove the rules for the
> exchangeability of fuel and electricity. Consequently, highly or entirely electrified
> processes covered by the EU ETS should benefit from free allocation in the same way as
> processes with high direct emissions. Therefore, the amount of free allocation should be
> determined regardless of the share of direct and indirect emissions for installations
> falling under the same benchmark. **Even though free allocation for those processes will
> cover also indirect emissions**, it does not necessarily imply that carbon leakage risks
> … have been fully addressed … **In turn, financial measures to compensate indirect costs
> passed on in electricity prices should not compensate the same indirect costs covered by
> free allocation.**"

Applicable date, recital 40: for incumbent installations these provisions "should apply to
allocations relating to the period **from 1 January 2026**".

Confirmed independently in IR (EU) 2025/2620, recital 17:
> "The rules for the exchangeability of fuel and electricity have been removed for the
> determination of free allocation under the EU ETS starting in 2026. This means that free
> allocation granted under some ETS product benchmarks in the steel sector **will cover
> indirect emissions to a certain extent.** As the CBAM scope currently only covers direct
> emissions in the steel sector, only the direct emission share of the respective ETS
> benchmarks should be considered when determining the corresponding CBAM benchmarks."

Arithmetic implication: EAF carbon steel benchmark 2026-2030 = **0.142** tCO2e/t, of which
the direct share (the CBAM benchmark for scrap/EAF ingots) is **0.027** tCO2e/t. So roughly
**80 % of an EU EAF's free allocation is now compensation for electricity emissions it does
not itself owe allowances for.**

### 3.4 CBAM factor / free-allocation phase-out — exact schedule

Verbatim, ETS Directive Art. 10a(1a) as inserted by Directive (EU) 2023/959:

> "Subject to the application of Regulation (EU) 2023/956, no free allocation shall be
> given in relation to the production of goods listed in Annex I to that Regulation.
> By way of derogation … A factor reducing the free allocation for the production of those
> goods shall be applied (CBAM factor). The CBAM factor shall be equal to 100 % for the
> period between the entry into force of that Regulation and the end of 2025 and … shall be
> equal to **97,5 % in 2026, 95 % in 2027, 90 % in 2028, 77,5 % in 2029, 51,5 % in 2030,
> 39 % in 2031, 26,5 % in 2032 and 14 % in 2033. From 2034, no CBAM factor shall apply.**"

| Year | CBAM factor (= free allocation retained) | Share phased out |
|---|---|---|
| ≤2025 | 100 % | 0 % |
| 2026 | 97.5 % | 2.5 % |
| 2027 | 95 % | 5 % |
| 2028 | 90 % | 10 % |
| 2029 | 77.5 % | 22.5 % |
| 2030 | 51.5 % | 48.5 % |
| 2031 | 39 % | 61 % |
| 2032 | 26.5 % | 73.5 % |
| 2033 | 14 % | 86 % |
| 2034 | 0 % | 100 % |

Note it is a cliff, not a ramp, in the last step: 14 % → 0 %.

Sectors covered: iron and steel, cement, aluminium, fertilisers, hydrogen, electricity
(CBAM Reg. Annex I).

Allowances freed by the reduction go to the Innovation Fund (Art. 10a(1a) 4th subpara and
Art. 10a(8)).

### 3.5 ⚠ The schedule is under active revision (as of Sept 2026)

**COM(2026) 616**, published 17 July 2026 (the ETS "Phase 5" review), proposes to
**reintroduce 15 % of the phased-out free allocation from 2028** and to **extend the
phase-out endpoint from 2034 to 2038**, with free allocation from 2031 conditional on a
verified decarbonisation plan (80 % on approval, 20 % on demonstrated reductions), and
clawback of allowances if activity relocates outside the EU. It is a co-decision proposal;
Parliament/Council agreement is targeted for early 2027. Not yet law.
- ICAP: https://icapcarbonaction.com/en/news/eu-commission-publishes-eu-ets-review-proposal
- Carbon Brief Q&A: https://www.carbonbrief.org/qa-what-the-eus-carbon-market-review-means-for-climate-action
- EUROFER reaction: https://eurometal.net/ec-proposes-extending-free-ets-allowances-for-cbam-sectors-until-2038-eurofer-raises-concerns/

Separately, **COM(2026) 619** (17 July 2026) is a stand-alone proposal to raise free
allocation determined by the **heat and fuel fallback benchmarks** for 2026–2030, to be
paid out with the 2027 allocation for the 2026 year.
https://climate.ec.europa.eu/document/download/672f0f66-a1c2-41d4-a6d6-ff82fea1eb62_en

**Recommendation for your model: run the statutory schedule as base case and the
COM(2026) 616 schedule as a sensitivity.**

---

## 4. INDIRECT COST COMPENSATION

### 4.1 Legal basis

ETS Directive Art. 10a(6), as replaced by Directive (EU) 2023/959:

> "Member States **should adopt financial measures** … in favour of sectors or subsectors
> which are exposed to a genuine risk of carbon leakage due to significant indirect costs
> that are actually incurred from greenhouse gas emission costs passed on in electricity
> prices, provided that such financial measures are in accordance with State aid rules …
> **The financial measures adopted should not compensate indirect costs covered by free
> allocation in accordance with the benchmarks established pursuant to paragraph 1** …
> Where a Member State spends an amount higher than the equivalent of 25 % of the auction
> revenues … for the year in which the indirect costs were incurred, it shall set out the
> reasons for exceeding that amount."

Note this is *permissive* ("should adopt"), not mandatory — which is why coverage varies
by Member State.

### 4.2 The State aid instrument

*Guidelines on certain State aid measures in the context of the system for greenhouse gas
emission allowance trading post-2021*, **OJ C 317, 25.9.2020, p. 5** (2020/C 317/04), as
supplemented in **2022** and again in **December 2025**.
Consolidated text (EFTA/ESA mirror, includes both supplements):
https://www.eftasurv.int/cms/sites/default/files/documents/2020_ETS_Guidelines_and_2022_and_2025_supplements_consolidated_0.pdf

**Formula (point 28):**

- Where a product-specific electricity consumption efficiency benchmark exists (Annex II):
  `A_max,t = A_i × C_t × P_(t-1) × E × AO_t`
- Otherwise (fall-back):
  `A_max,t = A_i × C_t × P_(t-1) × EF × AEC_t`

where
- `A_i` = aid intensity (fraction),
- `C_t` = applicable regional CO2 emission factor, or a market-based CO2 emission factor,
  in tCO2/MWh,
- `P_(t-1)` = **EUA forward price at year t−1** in EUR/tCO2 (so the compensation lags the
  carbon price by a year),
- `E` = product-specific electricity consumption efficiency benchmark (MWh/t),
- `EF` = fall-back benchmark = **"80 per cent of actual electricity consumption"**, reduced
  by **1.09 % per year** from t = 2022,
- `AO_t` = actual output; `AEC_t` = actual electricity consumption.

**Aid intensity.** Originally 75 %. The **December 2025 supplement raised it to 80 %** for
the core sectors (point 27):
> "The aid is proportionate and has a sufficiently limited negative effect on competition
> and trade if it does not exceed **80 % of the indirect emission costs incurred for the
> sectors listed in Table 1 of Annex I and 75 % for the sectors listed in Table 2** of
> Annex I or any further sectors considered eligible pursuant to the procedure set out in
> point (21)."

The same amendment added ~20 sectors and 2 subsectors (organic chemicals, parts of
ceramics/glass/battery value chains), allows Member States to notify unlisted sectors, and
introduces a green-investment condition for large beneficiaries. Eligible from costs
incurred **1 January 2025**.
- Commission/press coverage: https://europeansting.com/2025/12/24/commission-amends-ets-state-aid-guidelines-to-tackle-carbon-leakage-for-more-energy-intensive-industries/

Optional Member State cap (point 31): indirect costs payable per undertaking may be limited
to **1.5 % of gross value added**.

### 4.3 Is steel eligible? — Yes

Annex I, entry 7: **NACE 24.10 "Manufacture of basic iron and steel and of ferro-alloys"**.
Also relevant: 24.51 (casting of iron, all product categories), 20.11.11.50 (hydrogen),
24.42 aluminium, 24.43/24.44/24.45 non-ferrous.

**Important for EAFs:** Annex II contains an electricity efficiency benchmark for *basic
oxygen steel* (0.03385 MWh/t crude cast steel) and for ferro-alloys, but **there is no
Annex II benchmark for EAF steel**. EAF steelmakers therefore fall under the **fall-back**
approach: 80 % of actual electricity consumption, degressive at 1.09 %/yr. DEHSt confirms
this empirically — the single largest fall-back element in Germany is "sector '2410 —
Manufacture of basic iron and steel and of ferro-alloys'" (13 % of 2023 aid, 12 % of 2024).

### 4.4 Regional CO2 emission factors (Annex III, tCO2/MWh)

| Zone | Factor |
|---|---|
| **Germany, Luxembourg** | **0.73** |
| **France** | **0.43** |
| **Spain, Portugal** | **0.47** |
| Belgium | 0.37 |
| Netherlands | 0.44 |
| Italy | 0.44 |
| Austria | 0.33 |
| Poland | 0.78 |
| Sweden | 0.60 |
| Czechia | 0.89 |
| Bulgaria, Romania | 0.91 |

A plant on a corporate PPA can instead use a **market-based CO2 emission factor**, which for
a genuinely renewable PPA drives `C_t` — and hence the compensation — toward zero.

### 4.5 Which Member States operate a scheme, and how much

Schemes approved/operating include **Germany, France, Netherlands, Belgium (incl. Wallonia),
Spain, Italy, Greece, Finland, Slovakia, Lithuania, Poland, Luxembourg, Austria** (€233 m,
approved Sept 2023). Multi-year Commission approvals: Germany **€27.5 bn** for 2021-2030
(SA.100559, Dec 2022); France **€13.5 bn** for 2021-2030 (Dec 2022).
- https://ec.europa.eu/commission/presscorner/detail/es/ip_22_4925
- https://france.representation.ec.europa.eu/informations/la-commission-autorise-un-regime-francais-de-135-milliards-deuros-visant-compenser-les-couts-des-2022-12-01-0_fr
- https://energy.ec.europa.eu/news/state-aid-commission-approves-eu233-million-austrian-scheme-compensate-energy-intensive-companies-2023-09-21_en

**Germany — verified against DEHSt's own reports** (Strompreiskompensation / Electricity
Price Compensation, EPC). Payments are made in the year after the cost year.
https://www.dehst.de/EN/Topics/SPK/spk_node.html

| Cost/accounting year | Total aid | Companies | Installations | EUA reference price used | Iron & steel share |
|---|---|---|---|---|---|
| 2022 | **€1.64 bn** | — | — | €54.06 | — |
| 2023 | **€2.395 bn** | 351 | 707 | €83.59 | **€644.53 m** (26.9 %, 106 applications) |
| 2024 | **€2.784 bn** | 366 | 729 | €89.29 | **€716.89 m** (25.8 %, 109 applications) |

Sources: DEHSt *EPC Report 2023* (as of 18/12/2024) and *EPC Report 2024*.
- https://www.dehst.de/SharedDocs/downloads/EN/spk/Auswertungsbericht_2023_Englische_Version.pdf
- https://www.dehst.de/SharedDocs/downloads/EN/spk/Auswertungsbericht_2024_Englische_Version.pdf

Total budget available in Germany's Climate and Transformation Fund (KTF): €3.896 bn for
2023, €2.85 bn for 2024 — no pro-rata reduction was needed in either year.

**➜ Your €1.6 bn figure for Germany is the *2022* cost year, not 2023.** The 2023 figure is
**€2.40 bn** and 2024 is **€2.78 bn**. Update the model input.

**France — verified against the official Art. 10a(6) reporting to the Commission**
(*Rapportage France coûts indirects au titre de 2022*, published 2024):
https://www.ecologie.gouv.fr/sites/default/files/documents/Rapportage%20France%20couts%20indirects%20au%20titre%20de%202022%20(2024).pdf

For costs incurred in **2022**:

| Sector | Compensation (1a+1b) | +1.5 % GVA top-up | Total 2022 |
|---|---|---|---|
| **24.10 Sidérurgie** | €130,794,148 | €16,096,123 | **€146,890,271** |
| 24.42 Aluminium | €125,012,374 | €33,841,654 | €158,854,028 |
| 20.13 Inorganic chemicals | €111,850,268 | €12,514,302 | €124,364,570 |
| 17.12 Paper/board | €71,292,561 | €4,298,718 | €75,591,279 |
| … | | | |
| **TOTAL** | €531,384,540 | €72,766,199 | **€604,150,739** |

Advance paid during 2023 on account of 2023 costs: €203,318,759.

France also notes the 2022 compensation was **33.2 %** of its non-aviation auction revenues
(€1,820,455,395), i.e. above the 25 % soft ceiling, and explains why (France gets few
auction rights relative to its industrial base because allocation keys were set on early-ETS
verified emissions, which are low for France's decarbonised power system).

**➜ Your "~€600 m France" figure is right in magnitude but is the *2022* cost year total
(€604 m), of which steel was €147 m.** I could not retrieve the 2023-cost-year report.

### 4.6 Is this the stated reason CBAM excludes indirect emissions? — **Yes, explicitly.**

CBAM Regulation (EU) 2023/956, **recital 19**:

> "The CBAM should also apply to indirect emissions. Those indirect emissions are the
> emissions arising from the generation of electricity used to produce the goods to which
> this Regulation applies. … **Indirect emissions should, however, not be taken into
> account initially for the goods in respect of which financial measures apply in the Union
> that compensate for indirect emissions costs incurred from greenhouse gas emission costs
> passed on in electricity prices. Those goods are identified in Annex II to this
> Regulation.** Future revisions of the EU ETS … and, in particular, revisions of the
> compensation measures of the indirect costs should be appropriately reflected as regards
> the scope of application of the CBAM."

CBAM **Annex II** is headed *"List of goods for which only direct emissions are to be taken
into account, pursuant to Article 7(1)"* and covers:
- **Chapter 72 — Iron and steel** in full, *except* ferro-silicon, ferro-silico-manganese,
  ferro-silico-chromium, ferro-molybdenum, ferro-tungsten, ferro-titanium, ferro-vanadium,
  ferro-niobium, ferro-phosphorus, ferro-silico-magnesium, other ferro-alloys, and **7204
  ferrous waste and scrap**;
- headings 7301–7311, 7318, 7326 (downstream steel articles);
- **all aluminium** headings;
- **2804 10 00 — Hydrogen**.

So the causal chain is stated in law: *EU steelmakers get indirect-cost compensation →
therefore CBAM does not charge importers for indirect emissions on those goods.*
https://eur-lex.europa.eu/eli/reg/2023/956/oj

Directive (EU) 2023/959 recital 11 completes the loop from the other side: after the 2024
benchmark redefinitions, "it is necessary to ensure that producers do not receive double
compensation for the same emissions with both free allocation and indirect costs
compensation, and thus to adjust accordingly the financial measures to compensate indirect
costs passed on in electricity prices."

---

## 5. CARBON LEAKAGE AND EXPORTS

### 5.1 What the law says today

There is **no export rebate and no export exemption**. EU steel exported out of the EU has
already had its direct emissions charged under the ETS (net of free allocation) and gets
nothing back at the border.

Directive 2003/87/EC Art. 10a(1a), 5th subparagraph, contains only a review clause:

> "By 31 December 2024 … the Commission shall assess the carbon leakage risk for goods
> subject to CBAM and produced in the Union **for export to third countries which do not
> apply the EU ETS or a similar carbon pricing mechanism.** … Where the report concludes
> that there is a carbon leakage risk … the Commission shall, where appropriate, submit a
> legislative proposal to address that carbon leakage risk **in a manner that is compliant
> with the rules of the World Trade Organization, including Article XX of the General
> Agreement on Tariffs and Trade 1994**, and takes into account the decarbonisation of
> installations in the Union."

And Art. 30(2) as amended: "Before 1 January 2028, and every two years thereafter … the
Commission shall assess the impact of CBAM on the risk of carbon leakage, **including in
relation to exports**."

The WTO qualifier is the crux: a straight export rebate looks like a prohibited export
subsidy / border tax adjustment on a non-product tax, which is why the Commission has
repeatedly declined to propose one.

### 5.2 State of the debate (as of Sept 2026)

- **3 July 2025** — Commission announced it would propose a measure so that CBAM sectors get
  "equal treatment for all goods, whether produced and sold in the EU, imported into the EU
  or exported", with a proposal by end-2025.
  https://taxation-customs.ec.europa.eu/news/cbam-commission-announces-plan-mitigate-carbon-leakage-risk-exporters-2025-07-03_en
- **16 December 2025** — the Commission's CBAM package delivered **not** an export rebate but
  a **Temporary Decarbonisation Fund**: **25 % of CBAM revenues in 2028–2029**, reimbursing
  exporters for part of the CO2 cost incurred as free allocation is withdrawn in **2026 and
  2027**, conditional on decarbonisation investment. Same package proposed extending CBAM to
  ~**180 downstream steel and aluminium products from 1 January 2028**.
  https://icapcarbonaction.com/en/news/eu-cbam-enters-compliance-phase-and-outlines-path-ahead
- **June 2026** — ECOFIN adopted the Council position on the CBAM review; EUROFER says
  "loopholes on circumvention, downstream and exports remain".
  https://www.eurofer.eu/press-releases/cbam-eu-ministers-make-progress-but-loopholes-on-circumvention-downstream-and-exports-remain-warns-steel-industry
- **17 July 2026** — COM(2026) 616 addresses residual leakage by *slowing the phase-out*
  (15 % reinstated from 2028, endpoint 2038) rather than by an export mechanism.

**Bottom line: as of today, an export rebate has been proposed by industry, repeatedly
studied, and repeatedly not proposed by the Commission.** The de facto substitutes are (i)
the slower free-allocation phase-out and (ii) the two-year Temporary Decarbonisation Fund.

---

## 6. INTERACTION WITH CBAM — how the importer's certificates are reduced

### 6.1 Enabling provision

CBAM Regulation Art. 31:

> "1. The CBAM certificates to be surrendered in accordance with Article 22 of this
> Regulation **shall be adjusted to reflect the extent to which EU ETS allowances are
> allocated free of charge** in accordance with Article 10a of Directive 2003/87/EC to
> installations producing, within the Union, the goods listed in Annex I to this Regulation.
> 2. The Commission is empowered to adopt implementing acts laying down detailed rules for
> the calculation of the adjustment … taking account of the different benchmarks used in the
> EU ETS for free allocation with a view to combining those benchmarks into corresponding
> values for the goods concerned, and **taking into account relevant input materials
> (precursors)**."

### 6.2 The mechanism — "SEFA", Implementing Regulation (EU) 2025/2620

*Commission Implementing Regulation (EU) 2025/2620 of 18 December 2025 … as regards the
calculation of the free allocation adjustment to the number of CBAM certificates*,
OJ L, 22.12.2025. https://eur-lex.europa.eu/eli/reg_impl/2025/2620/oj

**Equation 1 — the adjustment**

    FAA_g = SEFA_(g,y) × M_g

where `FAA_g` = free allocation adjustment for good g, `SEFA_(g,y)` = specific embedded free
allocation of good g in year y [tCO2e/t], `M_g` = mass imported.

**Equation 2 — process-level specific free allocation (actual data)**

    SFAProc_(g,y) = CBAM_y × CSCF_y × BM*_g

- `CBAM_y` = the CBAM factor from ETS Art. 10a(1a) (97.5 % in 2026 …)
- `CSCF_y` = the cross-sectoral correction factor under Art. 14(6) of Reg. 2019/331
- `BM*_g` = the process-related CBAM benchmark, Annex point 5 **Column A**

**Equation 3 — simple good:** `SEFA_(g,y) = SFAProc_(g,y)`

**Equation 4 — complex good (recursive over precursors):**

    SEFA_(g,y) = SFAProc_(g,y) + Σ_i [ m_(i,y) × SEFA_(i,y') ]

with `m_(i,y) = M_(i,y) / AL_(i,y)` (Equation 5) — precursor mass per tonne of good, using
the activity levels defined in IR (EU) 2025/2547. "Where precursors are themselves complex
goods, the calculation of SEFA_i shall be repeated recursively using Equations 2, 3 and 4,
as appropriate, until no more precursors are relevant."

**Equation 6 — default-value route:** `SEFA_(g,y) = CBAM_y × CSCF_y × BM_g`, with `BM_g`
from **Column B** (cumulative, precursors already included). Selection depends on country of
origin, CN code, alloy grade, and the default production route for that origin per IR (EU)
2025/2621. Where several alloy grades share a CN code, "the highest benchmark value given
for the relevant production year is used".

Electricity: "The free allocation adjustment for electrical energy (CN code 2716 00 00)
shall be zero" (Art. 2(2)), because Art. 10a(1) forbids free allocation for power
generation.

The final obligation is therefore, per good:
`certificates due = (embedded direct emissions × M_g) − FAA_g`, floored at zero and with no
refund of a negative balance.

### 6.3 The deliberate anti-DRI adjustment — read this one carefully

IR 2025/2620 **recital 16** is the single most important paragraph for your Australian-HBI
case:

> "Currently, **direct reduced iron (DRI) is covered by the hot metal benchmark of the EU
> ETS** and, without further differentiation, imports of steel based on natural-gas DRI
> would receive a free allocation adjustment that exceeds their embedded emissions for the
> first years in which a CBAM obligation is due, which means that no CBAM certificates would
> be due for DRI-based goods. Compared to this, secondary steel imports would face a CBAM
> obligation, despite having lower actual embedded emissions than DRI. In addition, the
> potential free allocation adjustment stemming from the hot metal benchmark would create a
> situation in which more carbon-intensive natural gas based DRI imports would receive more
> free allocation adjustment than secondary steel producers … **a dedicated CBAM benchmark
> for natural gas-based DRI should be created.** Taking into account the relative level of
> embedded emissions, the level of the DRI benchmark should be chosen to ensure that the
> CBAM obligation for primary natural gas-based DRI imports is lower than for primary blast
> furnace steel, but higher than for secondary steel."

Result: **imported DRI (CN 7203 10 00) is adjusted at 0.295 tCO2e/t (actual data) or 0.397
tCO2e/t (defaults), not at the 1.248 tCO2e/t hot metal benchmark that an EU DRI producer
receives.** The Commission consciously broke the mirror to protect EU scrap-EAF operators —
and in doing so opened a ~0.85–0.95 tCO2e/t gap in favour of DRI made inside the EU.

---

## 7. EUA PRICE

### 7.1 Spot / front-year, today

**€84.98 /tCO2 on 7 September 2026** (+0.95 % on the day, +3.29 % on the month, +10.12 % on
the year). All-time high €105.73 (February 2023).
Source: Trading Economics, EU Carbon Permits, https://tradingeconomics.com/commodity/carbon

Corroborating 2026 marks:
- Dec-2026 contract traded €74–77 in May 2026, monthly average €74.04.
- 29 June – 28 July 2026: high €86.91, low €76.92, average €80.96.
- ~€84.2/t at end-2025 (Statista series
  https://www.statista.com/statistics/1322214/carbon-prices-european-union-emission-trading-scheme/).

Regulatory reference prices actually used for indirect cost compensation (`P_(t-1)`, DEHSt):
2022 → **€54.06**; 2023 → **€83.59**; 2024 → **€89.29**.

### 7.2 Forward projections

| Source | Date | 2030 | 2035 | Note |
|---|---|---|---|---|
| **GMK Center consensus** (median of BNEF, ABN Amro, Refinitiv, ICIS, S&P Global, Aurora, PIK) | 2 Dec 2025 | **€126/t** (range €80–147) | — | also 2026 €85/t, 2027 €100/t |
| **ABN AMRO** baseline | 2026 | **€145/t** | **€200/t** | ESG Economist / Carbon Market Strategist |
| **Veyt** | 2025/26 | — | "+€100 by 2035" | forecasts €160/t as early as 2027 |
| Trading Economics model | 7 Sep 2026 | — | — | €84.71 end-Q3 2026; ~€92.08 in 12 months |

**⚠ Correct a common misquote.** BloombergNEF's headline **€149/t by 2030 is for ETS2**
(road transport, buildings, small industry, launching 2027) — *not* ETS1. Published
6 March 2025.
https://about.bnef.com/insights/commodities/europes-new-emissions-trading-system-expected-to-have-worlds-highest-carbon-price-in-2030-at-e149-bloombergnef-forecast-reveals/
Consensus: https://gmk.center/en/infographic/carbon-price-in-the-eu-ets-to-hit-e126-t-by-2030/

**Suggested modelling values (ETS1):** €85/t for 2026, €100/t for 2027, **€120–145/t for
2030**, **€170–200/t for 2035**, with a low case that holds €80–100 flat if COM(2026) 616's
looser free allocation and a weaker 2040 target both land.

---

## 8. SYNTHESIS — (a) H2-DRI-EAF in Germany vs (b) identical plant in Australia exporting HBI to a German EAF

Illustrative assumptions: 1 t crude steel from ~1.08 t HBI/DRI; year 2026 (CBAM factor
97.5 %); EUA €85/t; CSCF taken as 1.0 (the 2026–2030 value was not yet published in the
sources I reached — it has historically been ~0.87–1.0 and it scales *both* EU free
allocation and the importer's adjustment, so it partly cancels).

### (a) H2-DRI-EAF, Germany — one integrated ETS installation

| Item | Basis | tCO2e per t crude steel | € per t steel at €85 |
|---|---|---|---|
| ETS surrender obligation — EAF direct (electrodes, carburiser, lime, burners) | MRR source streams | +0.05 to +0.10 | +€4 to +€9 |
| ETS surrender obligation — H2 shaft direct | residual NG / flux only | +0.00 to +0.03 | +€0 to +€3 |
| ETS surrender obligation — electricity | **not in scope** | 0 | €0 |
| **Free allocation, hot metal benchmark** | 1.248 × 1.08 t DRI × 0.975 | **−1.314** | **−€112** |
| Free allocation, on-site H2 benchmark (**disputed**, see §2.4) | 7.98 × ~0.055 t H2 × 0.975 | −0.43 (if granted) | −€36 (if granted) |
| **Net ETS position** | | **≈ −1.20 (surplus)** | **≈ −€100/t (income)** |
| Electricity ETS cost, passed through in the power price | ~3.5–4.5 MWh/t × marginal generator CO2 | implicit, not an ETS liability | grid: potentially €50–150/t; PPA/own RES: ≈€0 |
| Indirect cost compensation available | NACE 24.10, fall-back = 80 % of AEC, `C_t` = 0.73 (DE), `A_i` = 0.80, `P_(t-1)` = prior-year EUA | credit | up to ~0.80 × 0.73 × 85 × 0.80 × MWh/t ≈ €40/MWh·t — **but** must not double-count what free allocation already covers (Art. 10a(6); 2024/873 rec. 9), and a renewable PPA lets `C_t` be market-based ≈ 0 |

**The German plant is a net receiver of carbon value in 2026, not a payer.** It surrenders
~0.05–0.13 EUAs and receives ~1.31 (possibly ~1.74) per tonne of steel. That surplus decays
with the CBAM factor: 1.314 → 1.28 (2027) → 1.21 (2028) → 1.04 (2029) → **0.69 (2030)** → 0
(2034), unless COM(2026) 616 passes and stretches it to 2038.

### (b) Identical H2-DRI plant in Australia, HBI exported to a merchant German EAF

| Item | Who pays | tCO2e per t crude steel | € per t steel |
|---|---|---|---|
| Australian plant, EU ETS | — | **0 — outside the EU ETS entirely** | €0 |
| Australian plant, EU free allocation | — | **0 — no allocation is possible outside the ETS** | €0 |
| Australian plant, Australian carbon | Safeguard Mechanism applies only above 100 kt CO2e/yr; a genuine H2-DRI plant will normally be below | ~0 | ~€0 |
| **CBAM on HBI import** (CN 7203 10 00) | importer | direct embedded only (Annex II); verified H2 route ≈ 0.00–0.05 | €0–4 |
| CBAM free-allocation adjustment on HBI | importer | −0.295 × 1.08 × 0.975 ≈ **−0.31** — but the balance is **floored at zero, not refunded** | €0 |
| German merchant EAF, direct emissions | EU EAF | +0.05 to +0.10 | +€4 to +€9 |
| German merchant EAF, free allocation | EU EAF | **EAF carbon steel 0.142 × 0.975 ≈ −0.138** (⚠ see interpretive risk below) | **−€12** |
| German merchant EAF, electricity | not in ETS scope; ~0.6–0.9 MWh/t | implicit | small; compensable |
| **Net EU carbon cost of route (b)** | | **≈ −0.04 to −0.09 (small surplus)** | **≈ −€3 to −€8/t** |

**⚠ Interpretive risk on the merchant EAF's benchmark.** Both EAF benchmarks now exclude
"steel produced from iron sponge already covered under the hot metal benchmark", while the
hot metal benchmark's activity level is *tonnes of sponge iron at the DRI reactor exit* —
which a merchant EAF does not produce. The natural reading is that imported HBI is **not**
"already covered under the hot metal benchmark", so the EAF keeps the EAF carbon steel
benchmark (0.142). A stricter categorical reading would leave a merchant HBI-melting EAF
with **no product benchmark at all**. This is worth confirming with DEHSt before you rely on
it — it is ±€12/t of steel and it decides whether route (b) is roughly neutral or mildly
negative.

### Where they differ

1. **The gap is entirely in free allocation, not in emissions pricing.** Both routes owe
   almost nothing in carbon *charges*; the German integrated plant collects ~1.31 EUAs/t of
   steel of *free allowances* it does not need, and the Australian producer collects nothing.
   At €85/t that is a **~€100–110 per tonne of steel advantage to producing the iron inside
   the EU**, in 2026, decaying to ~€60 by 2030 and to zero in 2034 (or 2038 under COM(2026)
   616). This dwarfs typical HBI ocean-freight differentials and is, on today's law, the
   dominant policy signal in your Australia-vs-Germany comparison.
2. **The Commission built the asymmetry on purpose.** IR 2025/2620 recital 16 says so in
   terms: imported DRI is adjusted at 0.295/0.397, not at 1.248, because a full mirror would
   have zeroed out CBAM for DRI-based goods and disadvantaged EU scrap-EAF.
3. **The adjustment is a floor, not a transfer.** A zero-carbon Australian HBI producer
   gains nothing from the free-allocation adjustment because it has no positive obligation to
   reduce. The adjustment only helps *dirty* importers. So CBAM gives a green exporter no
   credit at all, while the identical EU plant is paid.
4. **Nobody is priced for electricity carbon intensity — on either side.** The EU plant's
   electricity emissions are not its ETS liability (they belong to the generator), and 2026
   free allocation now *reimburses* it for part of them (2024/873 recital 9; IR 2025/2620
   recital 17). The Australian plant's electricity emissions are excluded from CBAM by Annex
   II. **So an Australian H2-DRI plant on a coal grid and one on 100 % solar face the same
   (zero) EU carbon cost.** Your model should not expect CBAM or the ETS to discriminate
   between them; only voluntary/contractual green-iron premiums do.
5. **The one place electricity carbon intensity does bite is the compensation formula, and it
   bites backwards.** `C_t` = 0.73 for Germany vs 0.43 for France vs 0.47 for Spain/Portugal
   means a German EAF is entitled to *more* indirect cost compensation than a French one, for
   the same MWh — because German power is dirtier and therefore its ETS pass-through is
   larger. That is a compensation for cost, not a penalty for carbon.
6. **Exports get nothing.** If either plant's steel leaves the EU, no rebate exists (§5).

### Sensitivities that would change the ranking

- **CSCF** for 2026–2030 (scales EU free allocation and the CBAM adjustment together).
- Whether the **on-site electrolyser earns the 7.98 tCO2e/t hydrogen benchmark** (§2.4) —
  ±€36/t of steel.
- Whether the **merchant EAF keeps the EAF carbon steel benchmark** (above) — ±€12/t.
- **COM(2026) 616**: +15 % free allocation from 2028 and endpoint 2038 keeps the EU-location
  advantage alive roughly four years longer.
- **EUA price**: every €10/t on the EUA moves the EU-location advantage by ~€13/t of steel in
  2026 and ~€7/t in 2030.

---

## Where sources conflict or are out of date

| Claim commonly seen | Status |
|---|---|
| "A plant with no direct emissions is outside the ETS" (ECJ C-577/16 *Trinseo*) | **Superseded** by the Art. 2(1) amendment in Dir. 2023/959; Commission guidance §2.3.4 says so explicitly. |
| "Annex I covers *pig* iron production" | **Out of date** since 1 Jan 2024 — the word "pig" was deleted. |
| "Hydrogen threshold is 25 t/day, reforming/partial oxidation only" | **Out of date** — now 5 t/day, all processes including electrolysis. |
| "EAF free allocation is scaled by the direct/total emissions ratio" (GD9 formula, 0.376 tCO2/MWh) | **Correct to 2025, wrong from 2026** — Art. 22 of Reg. 2019/331 deleted by 2024/873. |
| "Hot metal benchmark = 1.288" | 2021–2025 value. **1.248** for 2026–2030 (IR 2026/1412). |
| "EAF carbon steel benchmark = 0.209" | Wrong even for 2021-25 — the value in IR 2021/447 is **0.215**; 2026–2030 is **0.142**. |
| "There is a DRI benchmark" | **No.** DRI sits inside the hot metal benchmark (2024/873), measured in tonnes of sponge iron. A separate *CBAM* benchmark for DRI does exist (0.295/0.397) but that is not an ETS product benchmark. |
| "Aid intensity for indirect cost compensation is 75 %" | **80 %** for Annex I Table 1 sectors (incl. steel) from cost year 2025, per the Dec 2025 supplement; 75 % remains for Table 2. |
| "Germany paid €1.6 bn for 2023" | That is the **2022** cost year. 2023 = **€2.395 bn**; 2024 = **€2.784 bn** (DEHSt). |
| "France paid ~€600 m for 2023" | €604 m is the **2022** cost year total; steel's share was €146.9 m. |
| "BloombergNEF sees €149/t EUA in 2030" | That is **ETS2**, not ETS1. |
| "Free allocation ends in 2034" | True in current law; **COM(2026) 616 (17 July 2026) proposes 2038** with 15 % reinstated from 2028. Not yet adopted. |
| "CBAM will get an export rebate" | Repeatedly requested by EUROFER, **not proposed**. The Commission's answer (Dec 2025) is a two-year Temporary Decarbonisation Fund funded by 25 % of 2028–29 CBAM revenues, covering 2026–27 production. |

---

## Source list

**Primary EU law**
- Directive (EU) 2023/959 (ETS revision) — https://eur-lex.europa.eu/eli/dir/2023/959/oj
- Directive 2003/87/EC consolidated — https://eur-lex.europa.eu/legal-content/EN/ALL/?uri=celex:32003L0087
- Regulation (EU) 2023/956 (CBAM) — https://eur-lex.europa.eu/eli/reg/2023/956/oj
- Delegated Regulation (EU) 2019/331 (FAR)
- Delegated Regulation (EU) 2024/873 (FAR amendment; DRI into hot metal; Art. 22 deleted) — https://eur-lex.europa.eu/eli/reg_del/2024/873/oj
- Implementing Regulation (EU) 2021/447 (benchmarks 2021–2025) — https://eur-lex.europa.eu/eli/reg_impl/2021/447/oj
- Implementing Regulation (EU) 2026/1412 (benchmarks 2026–2030) — https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=OJ:L_202601412
- Implementing Regulation (EU) 2025/2620 (CBAM free allocation adjustment / SEFA) — https://eur-lex.europa.eu/eli/reg_impl/2025/2620/oj
- Implementing Regulation (EU) 2025/2547 (CBAM embedded emissions methods)
- Implementing Regulation (EU) 2025/2621 (CBAM default values)
- Regulation (EU) 2018/2066 (MRR)

**Commission guidance and State aid**
- *Guidance on Interpretation of Annex I of the EU ETS Directive* (Dec 2024) — https://climate.ec.europa.eu/document/download/edc93136-82a0-482c-bf47-39ecaf13b318_en?filename=policy_ets_gd0_annex_i_euets_directive_en.pdf
- *Guidance Document 9 — Sector-specific guidance* — https://climate.ec.europa.eu/system/files/2019-07/p4_gd9_sector_specific_guidance_en.pdf
- *ETS State aid Guidelines* 2020/C 317/04 + 2022 and 2025 supplements (consolidated) — https://www.eftasurv.int/cms/sites/default/files/documents/2020_ETS_Guidelines_and_2022_and_2025_supplements_consolidated_0.pdf
- COM(2026) 619 (fallback benchmark uplift) — https://climate.ec.europa.eu/document/download/672f0f66-a1c2-41d4-a6d6-ff82fea1eb62_en
- Free allocation to industrial installations — https://climate.ec.europa.eu/areas-action/carbon-markets/eu-emissions-trading-system-eu-ets/free-allocation/allocation-industrial-installations_en
- CBAM export leakage announcement (3 Jul 2025) — https://taxation-customs.ec.europa.eu/news/cbam-commission-announces-plan-mitigate-carbon-leakage-risk-exporters-2025-07-03_en
- German scheme approval (€27.5 bn) — https://ec.europa.eu/commission/presscorner/detail/es/ip_22_4925
- French scheme approval (€13.5 bn) — https://france.representation.ec.europa.eu/informations/la-commission-autorise-un-regime-francais-de-135-milliards-deuros-visant-compenser-les-couts-des-2022-12-01-0_fr

**Member State reporting**
- DEHSt EPC Report 2023 — https://www.dehst.de/SharedDocs/downloads/EN/spk/Auswertungsbericht_2023_Englische_Version.pdf
- DEHSt EPC Report 2024 — https://www.dehst.de/SharedDocs/downloads/EN/spk/Auswertungsbericht_2024_Englische_Version.pdf
- DEHSt SPK overview — https://www.dehst.de/EN/Topics/SPK/spk_node.html
- France, *Rapportage coûts indirects au titre de 2022* — https://www.ecologie.gouv.fr/sites/default/files/documents/Rapportage%20France%20couts%20indirects%20au%20titre%20de%202022%20(2024).pdf

**Analysis, prices and trade press (secondary)**
- Sandbag, *EU steelmaking: the ETS money is coming* (7 Jul 2026) — https://sandbag.be/2026/07/07/eu-steelmaking-the-ets-money-is-coming/
- ICAP, *EU Commission publishes EU ETS review proposal* — https://icapcarbonaction.com/en/news/eu-commission-publishes-eu-ets-review-proposal
- ICAP, *EU CBAM enters compliance phase* — https://icapcarbonaction.com/en/news/eu-cbam-enters-compliance-phase-and-outlines-path-ahead
- Carbon Brief ETS review Q&A — https://www.carbonbrief.org/qa-what-the-eus-carbon-market-review-means-for-climate-action
- EUROMETAL on 2026-2030 ferrous benchmarks — https://eurometal.net/european-commission-publishes-delayed-lower-ferrous-ets-benchmarks/
- EUROMETAL on CBAM benchmarks/default values — https://eurometal.net/eu-commission-finalizes-cbam-benchmarks-default-values-ahead-of-january-2026-launch/
- EUROFER on CBAM export loopholes — https://www.eurofer.eu/press-releases/cbam-eu-ministers-make-progress-but-loopholes-on-circumvention-downstream-and-exports-remain-warns-steel-industry
- JRC, *Greenhouse gas intensities of the EU steel industry and its trading partners* (JRC129297) — https://publications.jrc.ec.europa.eu/repository/bitstream/JRC129297/JRC129297_01.pdf
- SteelOnTheNet CO2 by process step — http://www.steelonthenet.com/resources/kb/co2-emissions.html
- Trading Economics EU Carbon Permits — https://tradingeconomics.com/commodity/carbon
- GMK Center EUA consensus forecast (2 Dec 2025) — https://gmk.center/en/infographic/carbon-price-in-the-eu-ets-to-hit-e126-t-by-2030/
- BloombergNEF ETS2 €149/t 2030 (6 Mar 2025) — https://about.bnef.com/insights/commodities/europes-new-emissions-trading-system-expected-to-have-worlds-highest-carbon-price-in-2030-at-e149-bloombergnef-forecast-reveals/
- ABN AMRO carbon market research — https://www.abnamro.com/research/en/our-research/carbon-market-strategist-carbon-prices-heat-up-in-2026
- cbamguide benchmark tables — https://cbamguide.com/carbon/ets-benchmarks/ , https://cbamguide.com/sectors/steel/benchmarks/ , https://cbamguide.com/carbon/free-allocation/
