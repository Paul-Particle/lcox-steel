# Green iron trade, "renewables pull", and carbon-policy treatment of exported low-carbon iron

Literature and policy review compiled 2026-09-07.

Scope: academic and think-tank work relevant to modelling H2-DRI / MOE / aqueous-electrowinning
green iron in Australia and Brazil, exported as HBI into an EU EAF, versus production in
Germany / France / Spain — with emphasis on the CBAM indirect-emissions exclusion.

---

## 1. Johnson, Åhman, Nilsson & Li (2025), Nature Communications — the closest paper to the question

**Citation.** Constantin Johnson, Max Åhman, Lars J. Nilsson, Zhenxi Li, "Emerging green steel
markets surrounding the EU emissions trading system and carbon border adjustment mechanism",
*Nature Communications* 16, published **13 October 2025**. DOI 10.1038/s41467-025-64440-9.
Division of Environmental and Energy Systems Studies, Lund University.

- Paywalled/redirecting at nature.com; **open full text at PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC12518805/**
- PubMed record: https://pubmed.ncbi.nlm.nih.gov/41083494/
- Nature landing page: https://www.nature.com/articles/s41467-025-64440-9
- ResearchGate copy: https://www.researchgate.net/publication/396457227

### 1.1 Scope 2 exclusion under CBAM — what the paper actually says

> "scope 2 emissions originating from electricity production are not included in the CBAM for
> iron, steel, and hydrogen due to World Trade Organization (WTO) compatibility and the conflict
> with the current European indirect cost compensation (ICC) scheme."

This confirms the framing in the modelling context: the exclusion is *deliberate*, and its stated
justification is (a) WTO non-discrimination — the EU cannot charge importers for indirect emissions
while it *compensates* its own producers for them — and (b) internal consistency with the ICC
(state-aid indirect cost compensation) scheme.

**The load-bearing nuance the paper adds:** ICC is not applied EU-wide. It covers only
**15 of 27 Member States**, and even there it only *partially* covers the marginal ETS-induced
electricity cost. The paper therefore argues:

> "scope 2 emissions can in theory be partially included in the CBAM whilst respecting WTO
> principles."

That is a genuinely new legal argument — the WTO objection is weaker than usually assumed
because the EU's own indirect-cost shield is partial and patchy.

### 1.2 The stated risk of emission-intensive HBI / H2-DRI-EAF imports

Explicit, and phrased almost exactly as the modelling context frames it:

> "there is a risk of emission-intensive imports, such as HBI and H2-DRI-EAF steel, where
> hydrogen is produced by fossil electricity."

### 1.3 Modelling approach

- Dynamic **cost-optimisation model built in MATLAB**; two time-sequentially coupled sub-models:
  an **investment model** (capacity additions with a 3-year construction lag) and a
  **production-system / market model** (linear programme, dual-simplex).
- Geography: **EU27 + China, India, Brazil, Russia, Japan, South Korea, Türkiye, USA, Australia**
  — >90% of global crude steel.
- Routes: BF-BOF (brownfield / greenfield / reline), scrap-EAF (100% scrap),
  **H2-DRI-EAF with 12.5% scrap charge**, NG-DRI-EAF (modelled, then dropped as non-competitive).
- Scenario axes: **CBAM scope 1 only vs CBAM scope 1 & 2**; EU ETS **€80–150/tCO2**;
  hydrogen at **€4.12 and €5.50 /kgH2**.

Cost assumptions (€2022), useful as external calibration points:

| Item | Value |
|---|---|
| BF-BOF brownfield CAPEX | €228 /tCS |
| BF-BOF greenfield CAPEX | €592 /tCS |
| H2-DRI-EAF CAPEX (100% DRI) | €716 /tCS |
| Reduction shaft | €308 /tCS |
| Electrolyser | €0.585 /W (doubled to €1.17/W to cover a 10-yr replacement over a 20-yr plant life) |
| WACC | 8% |
| O&M | 2% of investment |
| BF-BOF life | 19 yr greenfield/brownfield, 10.5 yr reline |
| H2-DRI-EAF life | 20 yr |

### 1.4 Quantification

**The single most useful number for the indirect-emissions question:**

> inclusion of scope 2 emissions "further mitigates imported embedded emissions
> (**12.9 Mt CO2**) due to the higher grid-emission intensity of third countries compared to
> the EU27 average."

i.e. moving CBAM from scope-1-only to scope-1-and-2 avoids ~**12.9 Mt CO2** of imported embedded
emissions in their EU steel-import scenario. That is the closest thing in the literature to a
direct quantification of the CBAM indirect-emissions gap for steel.

Other numbers:

- EU secondary (scrap-EAF) production 54 Mt (2025) → **68 Mt (2035)**.
- BF-BOF becomes non-competitive from **2030** at ETS >€150/tCO2, from **2033** at >€90/tCO2.
- Production cost per tonne crude steel: EU BF-BOF **€508**; China BF-BOF **€461**;
  India BF-BOF **€456**; EU H2-DRI-EAF **€636–762** (at €60–100/MWh electricity);
  NG-DRI-EAF **€604–685**.
- Emission intensities used: H2-DRI-EAF low-emission **64 kgCO2/tCS**; H2-DRI-EAF
  **high-emission 437 kgCO2/tCS** (this is the fossil-electricity-hydrogen case);
  NG-DRI-EAF 1.13 tCO2/tCS low, **1.67 tCO2/tCS** with upstream methane.
- Marginal abatement: hydrogen would need to fall by **€2.12–3.95 /kgH2** between 2025 and 2026
  to make H2-DRI-EAF beat BF-BOF at ETS €80–150.
- **Transport cost €4.44–39.27 /tonne** depending on distance; the paper notes shipping is
  cheaper than beneficiating iron ore when the Fe-content gap exceeds ~7 pp. This matters for
  a Brazil/Australia → EU HBI freight assumption.
- Competitive H2-DRI-EAF exporters absent carbon pricing: **China, India, Brazil, Australia**.
- Competitive secondary-steel exporters under CBAM: Brazil, USA, Japan, China.
- At risk: Türkiye (unless the grid decarbonises), Japan and South Korea (scrap-scarce).

### 1.5 Interaction with ETS free-allocation phase-out — an under-appreciated subsidy

The paper's sharpest structural point: the 2026 ETS benchmark revisions **extend free allocation
to the low-emission intermediates used in H2-DRI-EAF**:

- **Sintered ore benchmark** from 2026 now includes **pelletising** → free allocation for pellets.
- **Hot metal benchmark** from 2026 now includes **DRI reactors and sponge iron** → free
  allocation for DRI production.
- **Hydrogen benchmark** from 2026 now includes **electrolysis** → free allocation for green H2.

> "The revisions of the EU ETS will effectively subsidise the production of low-emission
> intermediate products used in H2-DRI-EAF steelmaking."

Phase-out schedule: free allocation declines with the CBAM factor 2026–2034, with additional
reduction years 2031–2033 for the sintered-ore, hydrogen, hot-metal and EAF-carbon-steel
benchmarks; **full phase-out in 2034** (coke and lime benchmarks continue). Annual benchmark
reduction rates rise from 0.2–1.6%/yr (2021–2025) to **0.3–2.5%/yr** after 2025.

Critical modelling finding: *"The competitiveness of H2-DRI-EAF relative to BF-BOF is invariant
to the benchmark value"* — a lower benchmark removes net-negative ETS revenue from the H2 route
by exactly as much as it adds carbon cost to the BF route.

**Implication for the model in the context:** an EU-located H2-DRI plant receives free allocation
on DRI, pellets *and* hydrogen from 2026 to 2033 — a real, quantifiable revenue stream that an
Australian or Brazilian exporter does not get. This partly offsets the CBAM scope-2 asymmetry
running the other way. Both need to be in the same scenario to get the sign right.

### 1.5b HBI in the paper — and a big stated limitation

**The model does not actually trade HBI.** The authors say so:

> "by not accounting for HBI-trade in the model, the results for several countries highlight
> barriers for domestic H2-DRI-EAF investments which could be overcome by importing HBI."

So the paper's EU-import results are a *lower bound* on green-iron trade. Its HBI statements are
discursive rather than modelled:

- "HBI-trade could enable the EU to retain its current steel production volume with the phase-out
  of free allocation" — and likewise for Japan and South Korea.
- Cost-competitive HBI export candidates from the model: **China, India, the USA, Brazil,
  Australia**; from the literature it also cites Chile, Canada, and MENA.
- No HBI price is given.

### 1.5c Free-riding / resource shuffling — the paper names it

> "countries can invest in emission-reducing technologies for **single installations targeted at
> exports**, and thereby **retain export volumes while most of the installed steel capacity remains
> emission-intensive**."

It also flags verification risk:

> "The challenge of verification and detection of inaccurate reporting of emissions or
> circumvention through, e.g., recirculation of carbon prices paid in the country of origin, may
> also affect the net competitiveness between producers."

Model caveat stated by the authors: "The model assumes perfect markets and foresight ... which
negates global overproduction with perfect pricing of steel and emissions" — so the reshuffling
risk is identified but **not quantified**.

### 1.5d Scope 1 vs scope 1&2 scenario differences

- "**EU secondary and H2-DRI-EAF steel production is therefore at competitive disadvantage with a
  CBAM scope 1 coverage only.**" — the exclusion actively penalises the *electricity-intensive*
  EU routes, which is the exact asymmetry in the modelling context.
- "No change in production or trade of BF-BOF was found as a result of electricity prices"
  between the two scenarios — BF-BOF barely uses electricity, so scope 2 coverage only bites on
  the EAF/DRI/H2 routes.
- Scope 2 coverage "favours countries with low grid-emission intensities" — which for the two
  modelled exporters cuts opposite ways: **Brazil's grid would clear it easily; an
  Australian coal-grid plant would not** (an islanded RE plant would).

### 1.5e Carbon prices in the model

- EU ETS: linear from **€80.32/tCO2 (2022) to €150/tCO2 (2035)**.
- CBAM cost = embedded emissions × (EU ETS price − effective carbon price paid in the third
  country), per **Article 9 of Regulation (EU) 2023/956**.
- CBAM factor phases in **2026–2034**; the paper flags "considerable stress post year 2029
  attributed to the non-linear decline of the CBAM factor."
- Modelled exporters that get H2-DRI-EAF investment **without any significant domestic carbon
  price**: "China, India, Brazil, and Australia ... due to their potential for affordable
  renewable electricity."
- Countervailing Brazil finding: "competitive BF-BOF steel and increasing global demand resulted
  in additional brownfield and greenfield BF-BOF investments in **Brazil and India**, further
  increasing existing carbon lock-ins." Brazil's secondary-steel potential is "limited by domestic
  scrap availability."

### 1.6 Policy recommendations (all of them)

1. **Three relocation pathways acknowledged as legitimate:** reconfigure within the EU; import
   HBI from established steelmaking countries; outsource H2-DRI to developing countries.
2. Extend the **EU ETS to upstream and midstream methane**, and extend **CBAM to natural gas and
   methane**; base the methane extension on **GWP20**, not GWP100.
3. **ICC reform:** "ICC should only compensate renewable electricity or be phased out and scope 2
   emissions fully included in the CBAM."
4. CBAM plus **climate trade diplomacy** to break carbon path-dependency.
5. **International revenue recycling** of CBAM certificate revenue into bilateral decarbonisation
   programmes.
6. Trade-dependent partner countries gain most from adopting their own carbon pricing.
7. **Common CBAMs (CCBAMs)** / green materials clubs among similar-ambition countries.
8. Technical: reduce the **Actual Activity Level update threshold**, or update Historical Activity
   Level more often (ex-ante vs ex-post distortion).
9. Warn on **NG-DRI lock-in**: ETS costs would push NG-DRI-EAF to NG-DRI-EAF-CCS soon after
   deployment, creating fossil lock-in with no zero-emission endpoint.
10. Frame resilience as **diversified supply chains and clean trade partnerships**, not domestic
    production mandates.
11. Avoid subsidy races: "a renewables pull can enable rapid deep decarbonisation of the European
    steel industry."

---

## 1.7 The legal hook, verified in the regulation itself

Regulation (EU) 2023/956, **Recital 19** (https://eur-lex.europa.eu/eli/reg/2023/956/oj/eng):

> "Indirect emissions should, however, **not be taken into account initially** for the goods in
> respect of which financial measures apply in the Union that compensate for indirect emissions
> costs incurred from greenhouse gas emission costs passed on in electricity prices. Those goods
> are identified in **Annex II** to this Regulation."

Two things worth carrying into a model:

1. The exclusion is **conditional on the existence of EU indirect cost compensation** — it is not
   a WTO-driven absolute. Remove or narrow ICC and the legal justification for the exclusion
   evaporates. This is exactly the lever Johnson et al. (2025) recommend pulling
   (recommendation 3, §1.6).
2. The word is "**initially**". The Regulation itself contemplates later inclusion, so a
   long-horizon model that treats scope-2 exclusion as permanent is making a policy assumption,
   not reading the law. A scenario in which indirect emissions are phased in — most plausibly
   alongside ICC phase-out around the 2034 full free-allocation phase-out — is defensible.

---

## 2. The "renewables pull" literature

### 2.1 Verpoort, Gast, Hofmann & Ueckerdt (2024), Nature Energy — the canonical quantification

**Citation.** P. C. Verpoort, L. Gast, A. Hofmann, F. Ueckerdt, "Impact of global heterogeneity
of renewable energy supply on heavy industrial production and green value chains",
*Nature Energy* **9**, 491–503, published **24 April 2024**. DOI 10.1038/s41560-024-01492-z.
PIK Potsdam, Research Department 3 — Transformation Pathways.

- Nature: https://www.nature.com/articles/s41560-024-01492-z
- **Open PDF: https://publications.pik-potsdam.de/pubman/item/item_29811_4/component/file_29872/29811oa.pdf**
- Press: https://www.eurekalert.org/news-releases/1042076

**Framework.** `Relocation savings = energy-cost savings (renewables pull) − transport penalty −
financing penalty`. Financing penalty = WACC 5% (RE-scarce) → 8% (RE-rich); found to be small.
Technology parameters represent **2040**. Commodities: hot rolled coil (steel), urea, ethylene.

**Four import cases** by depth of relocation:
- Base case: full domestic production.
- **Case 1A**: import H2 by ship (H2 transport €50/MWh).
- **Case 1B**: import H2 by pipeline (H2 transport €15/MWh).
- **Case 2**: import the **intermediate** — for steel this is **DRI**.
- **Case 3**: import the (semi)finished product — HRC.

**Headline table (Table 1) — relocation savings for full relocation (case 3):**

| Price case | Electrolysis price RE-rich / RE-scarce | Baseload RE-rich / RE-scarce | Δ (€/MWh) | Steel | Urea | Ethylene |
|---|---|---|---|---|---|---|
| Weak pull | 30 / 50 | 50 / 70 | 20 | **8.7%** | 14.1% | 20.6% |
| Medium pull | 30 / 70 | 50 / 90 | 40 | **18.3%** | 32.1% | 37.6% |
| Strong pull | 15 / 85 | 35 / 105 | 70 | **31.5%** | 55.0% | 60.0% |

**The result that matters most for an HBI-trade model — savings by depth of relocation
(medium pull, Δ€40/MWh):**

| Depth | Steel | Urea | Ethylene |
|---|---|---|---|
| Case 1A — ship hydrogen | **1%** | 2% | 2% |
| Case 1B — pipeline hydrogen | **9%** | 19% | 19% |
| **Case 2 — import the intermediate (DRI)** | **13%** | 25% | 37% |
| Case 3 — import finished product | **18%** | 32% | 38% |

Authors' own conclusion: importing the intermediate is the **"sweet spot" of relocation** —
it captures ~72% of the full steel relocation saving (13 of 18 points) while leaving the EAF,
casting and rolling value-added in the importing region. Shipping hydrogen instead is
"a potentially expensive and risky strategy" and "challenges the H2 import strategies of some
RE-scarce regions."

Steel savings are the *smallest* of the three commodities "where raw-material costs (iron ore and
so on) are high" — energy is a smaller share of the cost base than for urea or ethylene.

**German green-relocation-protection case study.** The subsidy Germany would have to pay to
prevent relocation of steel + urea + ethylene: **€6–18 bn/yr** (scenario 1, full retention) or
**€3–9 bn/yr** (scenario 2, partial retention, roughly the capacity Germany's IPCEIs and CCfDs
envisage transforming by 2030). Compared in the paper to German federal ministry budgets.

**They name "resource shuffling" explicitly** — but in a *different* sense than the CBAM usage:
as a risk to the **exporting** region, "using RE potentials only for exports instead of domestic
climate mitigation." Worth distinguishing from the CBAM/California sense.

**Honest caveat the authors state:** for DRI, "the emergence of a global market is unclear, yet
existent dependencies on iron-ore imports raise the question whether switching to DRI imports
would create much difference" — i.e. they think the security-of-supply objection to HBI imports
is weak because Europe already imports the ore.

### 2.2 Samadi, Fischer & Lechtenböhmer (2023) — the concept paper, and a direct resource-shuffling hit

**Citation.** Sascha Samadi (Wuppertal Institute), Andreas Fischer (German Economic Institute,
Cologne), Stefan Lechtenböhmer (Wuppertal Institute), "The renewables pull effect: How regional
differences in renewable energy costs could influence where industrial production is located in
the future", *Energy Research & Social Science* **104** (2023) 103257.

- Open PDF: https://sci4climate.nrw/wp-content/uploads/2023/12/Samadi_et-al_2023_The-renewables-pull-effect.pdf
- Also: https://d-nb.info/1302339036/34
- Institute page: https://wupperinst.org/en/a/wi/a/s/ad/8256/

Uses **DRI and ammonia** as the two worked examples. Key findings:

- With renewable energy, energy costs are expected to reach **~50% of total DRI production cost**
  and **~90% of ammonia** by 2035 (their Fig. 2, production costs in €/t, 2035).
- Transport costs are low for both DRI and ammonia relative to energy costs.
- DRI relocation is favoured because good RE conditions **coincide with iron ore** in
  **South Africa, Canada, Brazil and Northern Sweden** — which are already the four main ore
  suppliers to the German steel industry.
- Existing blast-furnace assets are of **limited value** for H2-DRI, so asset stranding is not a
  strong brake (unlike ammonia, where existing synthesis loops carry over).
- Survey evidence: **32% of metal production/processing companies** and **24% of chemical**
  companies expect relocation due to renewables pull, against a **20% all-industry average**.

**The CBAM/resource-shuffling passage (footnote 9) — a direct hit on priority item 5:**

> "In reality, it is possible that a country with abundant renewable energy sources is already
> today using low-cost renewable energy sources (such as hydropower) to a certain extent for the
> production of industrial goods. In the case of another country enacting more stringent climate
> policies in combination with a CBAM, producers in the renewable-rich country may have an
> incentive to use existing renewable energy sources to produce goods for this particular export
> market, while using fossil fuels for the production of goods for the domestic market or other
> export markets. Such **'resource shuffling'** could negate or minimise the policy's climate
> mitigation effects."

Note the sting for a Brazil case specifically: Brazil's grid is ~85–90% renewable already, largely
hydro. A Brazilian HBI plant contracting existing hydro for EU-bound exports is the textbook
resource-shuffling case — no new renewable capacity, EU-bound tonnes look clean, marginal
generation elsewhere on the grid turns to gas.

Main text (their Fig. 1b) also formalises the mechanism the modelling context relies on: under
unilateral tightening in country A **plus** a CBAM, a green producer in RE-rich country B
"would not be faced with a (significant) border adjustment price" and therefore exports into A.
Renewables pull is *amplified*, not suppressed, by CBAM.

### 2.3 Nykvist, Gong, Algers & Åhman (2025) — **the dissent; read this one**

**Citation.** Björn Nykvist, Jindan Gong, Jonas Algers, Max Åhman, "Renewables pull and strategic
push – What drives hydrogen-based steel relocation?", *Applied Energy* **395** (2025) 126189.
DOI 10.1016/j.apenergy.2025.126189. CC-BY (hybrid OA). SEI + Lund.

- https://www.sciencedirect.com/science/article/pii/S0306261925009195
- https://www.sei.org/publications/renewables-pull-and-strategic-push/
- https://lup.lub.lu.se/search/publication/20b1e6b9-a499-4e3f-acfb-23780fc1005a
- SSRN preprint: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5174276

This paper **contradicts** the Verpoort/Samadi consensus:

- Compares **five value-chain configurations** on harmonised assumptions (previous studies used
  incompatible ones, making cross-comparison meaningless).
- "**The renewables pull effect is sensitive to assumptions and weaker than previously found.**"
- "The cost of hydrogen-based steel varies across geographies and value chain configurations to a
  **similar degree as conventional steel**, and other geographically varying factors such as
  **labour costs can be as important** for relocation."
- A **modest "strategic push"** — subsidies lowering the cost of hydrogen or of capital — is
  enough to **cancel the pull effect**.
- Warning: overemphasising renewables-pull-driven relocation risks *delaying* the innovation
  needed for 2030 targets, because it defers investment while waiting for cheap-RE regions.

**This is the most important complication to the picture in the context.** If the pull is weak
and cancellable by modest EU subsidy, a model that shows Australia/Brazil HBI dominating on
energy cost alone is probably over-reading the effect; capital cost (WACC) and labour deserve
sensitivity treatment on par with electricity price.

### 2.4 Devlin, Kossen, Goldie-Jones & Yang (2023), Nature Communications

**Citation.** A. Devlin, J. Kossen, H. Goldie-Jones, A. Yang, "Global green hydrogen-based steel
opportunities surrounding high quality renewable energy and iron ore deposits",
*Nature Communications* **14**, 2578, **4 May 2023**. DOI 10.1038/s41467-023-38123-2.

- https://www.nature.com/articles/s41467-023-38123-2
- **Open PDF: https://ora.ox.ac.uk/objects/uuid:70994333-c901-4d07-893c-30c2052d4184/files/r79407z05g**

Method: optimisation of **islanded** solar+wind H2-DRI-EAF systems co-located at iron ore mine
sites in 44 regions, then a machine-learning extrapolation to **>300 iron ore deposits**
(R² = 0.96 for LCOE, 0.85 for LCOS excluding ore and labour, ±$26/t).

Numbers (USD 2020):

- **LCOS $535–972 /t** without scrap charging; by 2050 **$535–831 /t**.
- LCOH2 **$1.63–2.80 /kg** by 2050; LCOE **$16–50 /MWh**.
- Average LCOE/LCOH2 fall from **$43/MWh and $3.2/kgH2 (2030)** to **$30/MWh and $2.1/kgH2 (2050)**.
- 2021 BF-BOF opex **$621–782 /t** in the same locations (vs $428–547/t in 2020) — the comparison
  is very sensitive to met-coal prices.
- Best locations are **near the tropics** (Iran, Peru, South Africa, Chile), solar-dominant
  (80–100% solar share of RE capacity in 2050).
- Islanded systems: electrolyser oversizing **1.3–3.7×**, EAF oversizing only **1.1–1.4×**;
  ~50% of H2 stored as compressed gas; 91% of storage cost is in H2, not batteries.
  **Flexible hydrogen production matters far more than flexible steelmaking.**
- Ore sensitivity: 10% change in ore cost → **3% change in LCOS**; ore is on average **27%** of
  steel production cost, but 45% in a low-Fe case (Kazakhstan, 20% Fe, DR-pellet premium tripling
  from $40/t to $122/t ore).
- Scrap: 25% scrap → 5% cost reduction, 50% scrap → 9% in 2030; only 2% by 2050.
- **Transport is small:** global average **3% added cost for inland transport and 7% for marine**
  (FOB and CFR to Qingdao). "It is unlikely that transport of product will dictate production
  facility locations."
- Grid vs islanded: to equalise LCOS, grid electricity would have to be **$80 / $70 / $60 per MWh**
  in 2030 / 2040 / 2050 on average (lowest $62 / $54 / $46). The **largest islanded-vs-grid gaps
  are in Brazil, Chile and Australia** — i.e. in exactly the two modelled countries, islanded RE
  beats grid supply, so the model should not assume grid-connected supply there.
- Emissions: islanded systems **0.1–0.3 tCO2/t steel** (including embodied PV/wind). Grid-based
  Canada and Sweden, closely followed by **Brazil**, reach **0.2 tCO2/t steel by 2030** thanks to
  hydro/nuclear. Iran's gas grid gives **2.7 tCO2/t steel even in 2050** — worse than BF-BOF.
- Australia is disadvantaged by **high wages**; without labour costs Australia's ranking improves
  materially. Australia would need **>4× its current national grid capacity** by 2050 (60%
  technology diffusion) and a **70-fold increase in current steel production**.

Explicit recommendation: "the trade and transport of primary or intermediary products
(iron ore, green H2 or **HBI**) should be explored."

### 2.5 Devlin & Yang (2022) — the earlier HBI-trade paper

"Regional supply chains for decarbonising steel: Energy efficiency and green premium mitigation",
*Energy Conversion and Management* 254 (2022) 115268.
https://www.sciencedirect.com/science/article/pii/S0196890422000644
Open copy: https://ora.ox.ac.uk/objects/uuid:93d9b6ec-1eb1-4fc0-ae8f-27cffdaf564b/files/s6q182k972

Case studies Australia→Japan and South Africa→Europe. Finding: **exporting hydrogen-reduced iron
as HBI is more energy-efficient and cost-competitive than exporting hydrogen plus iron ore**.
For export-dominant chains, **liquefied hydrogen beats ammonia** (process flexibility matched to a
solar profile; avoids ammonia cracking before direct reduction) — but HBI beats both.

### 2.6 Bilici et al. (2024) — the scenario study of a global green-iron market

**Citation.** S. Bilici, G. Holtz, A. Jülich, R. König, Z. Li, **H. Trollip**, B. McCall, A. Tönjes,
S. S. Vishwanathan, O. Zelt, S. Lechtenböhmer, S. Kronshage, A. Meurer, "Global trade of green iron
as a game changer for a near-zero global steel industry? — A scenario-based assessment of
regionalized impacts", *Energy and Climate Change* **5** (December 2024) 100161.
DOI 10.1016/j.egycc.2024.100161.

- https://www.sciencedirect.com/science/article/pii/S2666278724000370
- SSRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4841564
- https://portal.research.lu.se/en/publications/global-trade-of-green-iron-as-a-game-changer-for-a-near-zero-glob
- Wuppertal repository: https://epub.wupperinst.org/frontdoor/index/index/year/2024/docId/8721

(Trollip is a co-author here — this is the main Trollip-linked green-iron trade paper.)

Three scenarios to 2050: **Domestic** (strict regional co-location of iron and steel),
**Max Trade** (early emergence of a global green-iron market), **Intermediate Trade** (late).

Results:
- **12–21% of global crude steel** produced from **traded** green iron in 2050.
- **15–26 Mt/a of hydrogen consumption** relocated to global "sweet spots".
- Cost savings of only **2.2–3.9% of global annual steel production costs**.

The 2–4% global saving is a *much* more modest figure than the 13–18% regional relocation savings
in Verpoort et al., because it is a global average including regions that gain nothing. For a
Germany-vs-Australia bilateral comparison the Verpoort/Seibold numbers are the right anchor;
for a claim about the world steel system, this is.

### 2.7 Seibold, Neumann, Ueckerdt & Brown (2025) — PyPSA-based, German-specific, most directly comparable

**Citation.** Toni Seibold (TU Berlin), Fabian Neumann (TU Berlin), Falko Ueckerdt (PIK),
Tom Brown (TU Berlin), "Balancing Cost Savings and Import Dependence in Germany's Industry
Transformation", **arXiv:2510.00918v1, 1 October 2025**.
https://arxiv.org/pdf/2510.00918 · code: https://github.com/toniseibold/de-import

Method: **PyPSA-DE** (European energy system, sector-coupled) **coupled to TRACE** global supply
curves for 51 countries. Precursors modelled: **HBI**, ammonia, methanol. Five scenarios:
Base (domestic, H2 self-sufficient), EH (European H2 imports), WH (world H2 imports),
EHP (European precursor imports), WHP (world H2 + precursor imports).

Results:

- German **industry consumer costs €36.8 bn/a** in Base.
  - EH (European H2): −€2.0 bn/a (**−5.4%**)
  - WH (world H2): −€3.9 bn/a (**−10.7%**)
  - **EHP (European precursors): −€4.1 bn/a (−11.2%)**
  - **WHP (world H2 + precursors): −€8.6 bn/a (−23.3%)**, down to €28.2 bn/a
- **Relocation to European partners only captures 47.7% of the full non-European benefit.**
  This is the single cleanest "Spain vs Australia" number in the literature: sourcing HBI from
  Spain/Denmark/UK rather than Australia/Brazil gives you about **half** the saving.
- **German HBI wholesale price falls by up to 20.2%** — roughly **€843/t** in the domestic Base
  scenario down to the low **€600s/t** with precursor imports. (Ammonia −28.2%, methanol −21.0%.
  HBI falls least because iron ore cost dilutes the energy share.)
- **The model prefers precursor imports over hydrogen imports** whenever both are allowed, because
  it avoids hydrogen infrastructure: European H2 pipeline build falls from **71.1 TWkm (Base) to
  35.0 TWkm (EHP)**. German H2 imports collapse from 188.1 TWh (EH) to **4.7 TWh (EHP)**.
- German H2 marginal price: ~€98.8/MWh Base → **€89.9/MWh (EH)** → **€81.0/MWh (WH)**.
- **In EHP, precursor production shifts to Spain** (solar) **and Denmark/UK** (wind); Spain hits
  its maximum technical renewable potential in the model.
- Transport assumption to watch: intra-European HBI moves by **electric truck over 2,000 km at
  ~€5.5/t** (25 t cargo, 110 kWh/100 km, €62/MWh). Non-European HBI transport is handled inside
  TRACE. HBI export volume from any one country is capped at 238 Mt.
- WACC sensitivity: raising non-European WACC from 7% to **10%** makes the **European-only (EHP)
  scenario reach 78.3%** of the global (WHP) benefit — i.e. **country risk premium is the single
  lever that most changes the Australia/Brazil vs Spain answer**.
- Precursor plant flexibility assumed conservative: min part load **80% EAF, 90% DRI**,
  no shutdowns.

### 2.8 Caiafa, de Kleijne & de Coninck (2026), One Earth — Brazil-specific

**Citation.** Clara Caiafa, Kiane de Kleijne, Heleen de Coninck (TU Eindhoven / Radboud),
"Co-locating green steel and hydrogen production in renewables-rich developing regions delivers
greater climate and economic benefits than hydrogen trade", *One Earth* **9**(2), 101614,
**1 February 2026**. DOI 10.1016/j.oneear.2026.101614.

- https://www.cell.com/one-earth/abstract/S2590-3322(26)00015-1
- Open access via https://research.tue.nl/en/publications/3907ecc5-9fd4-45ce-ab4a-0ab5fa066c33

Case: the **Ceará (Brazil) – Netherlands green hydrogen corridor**. Finding: producing green
steel locally in Ceará and shipping the steel beats shipping hydrogen on **cost, emissions and
local socioeconomic outcomes**. Note this study pushes relocation **one step deeper than HBI** —
it argues for exporting *steel*, not iron, from Brazil. That is a direct challenge to the
"HBI into an EU EAF" configuration.

### 2.9 IRENA (2025) — global trade optimisation including DRI

**Citation.** IRENA, "Analysis of the potential for green hydrogen and related commodities trade",
**June 2025**. Contributors include Arno van den Bos, Karan Kochhar, Deepti Siddhanti,
Adrian Gonzalez, Yong Chen, Patricia Wild.

- https://www.irena.org/Publications/2025/Jun/Analysis-of-the-potential-for-green-hydrogen-and-related-commodities-trade
- PDF: https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2025/Jun/IRENA_TEC_GH2_and_commodities_trade_2025.pdf

Cost-optimisation across **35 global regions** for 2050. Two scenarios: **Same WACC** and
**Differentiated WACC**.

- Share of 2050 demand met by **trade**: ammonia 30%, e-methanol 18%, pure hydrogen 14.4%,
  **DRI 14%**.
- **73–80% of total green-hydrogen-equivalent trade travels as commodities, not as hydrogen** —
  "owing to lower transport costs and higher efficiency". Same conclusion as Verpoort et al.,
  reached by a different model.
- Conversion used: **1 Mt H2-eq = 16 Mt DRI** (also 5.67 Mt NH3, 8 Mt MeOH).
- Total infrastructure investment **USD 2.49 trillion** (4.7 TW renewables, 2.1 TW electrolysers,
  0.9 TWh batteries).
- **Same WACC**: exporters are Latin America, Middle East, North Africa, sub-Saharan Africa.
  The **Middle East imports iron ore (from Australia and sub-Saharan Africa) and exports DRI**.
- **Differentiated WACC**: **Australia becomes the top DRI exporter**, followed by the USA;
  together **~80% of all green DRI exports**. The Middle East flips to self-sufficiency.
  DRI traded share rises above 14%, and the **DRI-to-iron-ore trade ratio moves from 1:4 to
  1:2.2** — more iron shipped as DRI, less as ore.
  Notably, "the portion [routed via the Middle East] to produce DRI for export to Europe is
  reduced. This portion is relocated to **Australia, which directly produces DRI and exports it
  to Europe**."
- Europe, Japan, Southeast Asia and Korea are importers under **both** WACC scenarios.
- IRENA's framing of CBAM is that it "actively incentivise[s] the production and export of clean
  hydrogen and green commodities by putting a price on embedded emissions" — no discussion of the
  scope-2 exclusion.

**The WACC point recurs.** IRENA, Seibold et al. and Verpoort et al. all find that the *cost of
capital*, not the electricity price, decides whether the exporter is Australia, the Middle East,
or Africa. For a model comparing Australia and Brazil, differential WACC is likely a
first-order lever, not a refinement.

### 2.10 SteelWatch explainer (23 April 2025) — useful compiled numbers, no CBAM analysis

https://steelwatch.org/steelwatch-explainers/steelwatch-explainer-why-green-iron-trade-will-catalyse-steel-industry-decarbonisation/

- Global steel production cost falls **2–4%** with green iron trade (citing Bilici et al.).
- **Japan ~30% cheaper** with imported H2-DR iron than domestic production.
- **Western/Central Europe and South Korea ~20% cheaper** with Canadian imports;
  Germany-vs-Australia gives similar savings.
- **Shipping green iron needs less than a third of the volume** of shipping the ore and hydrogen
  separately: a 2.5 Mtpa DRI plant needs **0.75–1 million m³** for green iron vs **3.5 million m³**
  for the feedstocks.
- Australia's green iron manufacturing potential valued at **$60–185 bn/yr**.
- **Does not discuss CBAM, scope 2, or the fossil-electricity risk at all.** Notable gap in the
  advocacy literature.

### 2.11 Related work to be aware of

- Lopez, Farfan, Breyer — "Trends in the global steel industry: Evolutionary projections and
  defossilisation pathways through power-to-steel", *J. Cleaner Production* (2022); and
  "Global demand for green hydrogen-based steel: Insights from 28 scenarios",
  *Int. J. Hydrogen Energy* (2024) — https://www.sciencedirect.com/science/article/pii/S0360319924026624
  Projects **2,809–4,371 TWh_H2** of global steel-sector hydrogen demand by 2050. Lopez et al.
  argue that if the EU imports **HBI rather than raw materials**, it achieves both lower cost and
  Global South development.
- Colen et al. (2025), "Unpacking the Renewable Pull Effect: Conditions for Green Industrial
  Relocation", *Business Strategy and the Environment* —
  https://onlinelibrary.wiley.com/doi/10.1002/bse.4301
- "Toward a Renewables-Driven Industrial Landscape: Evidence on investment decisions in the
  Chemical and Steel Sectors" (Research Square preprint) —
  https://www.researchsquare.com/article/rs-5519615/v1 — survey of 300 chemical/steel
  decision-makers; **92% expect their company to relocate facilities as it decarbonises**.
- RIFS Potsdam blog on renewables pull evidence (Dec 2023):
  https://www.rifs-potsdam.de/en/blog/2023/12/how-decarbonization-will-transform-geography-industrial-production-new-evidence
- Trollip et al., "Strategising steel sector capacities and employment in the Global South: the
  case of South Africa", *npj Clean Energy* (2026) —
  https://www.nature.com/articles/s44406-026-00020-0
- "Decarbonizing global steel production", *Nature Reviews Earth & Environment* (2026) —
  https://www.nature.com/articles/s43017-026-00786-y
- Agora Industry, "The role of green iron trade in accelerating competitive steel transformation" —
  https://www.agora-industry.org/publications/the-role-of-green-iron-trade-in-accelerating-steel-transformation
- Agora Industry, "Green iron trade: Unlocking opportunities for Brazil", **November 2025** —
  https://www.agora-industry.org/fileadmin/Projekte/2024/2024-28_IND_Green_Iron/A-IND_377_Green_iron_trade_Brazil_WEB.pdf
- Agora Industry, "The global steel industry can achieve net-zero emissions by the early 2040s" —
  https://www.agora-industry.org/news-events/the-global-steel-industry-can-achieve-net-zero-emissions-by-the-early-2040s-1-1
- WEF, "How Australia can lead the global transition to green iron", December 2025 —
  https://www.weforum.org/stories/2025/12/australia-global-transition-green-iron/

---

## 3. Synthesis of the numbers that would change a model

| Question | Best available answer | Source |
|---|---|---|
| Is HBI the right thing to trade, vs H2 or steel? | Yes for H2 (13% vs 1% relocation saving, shipped H2); ~72% of the full-relocation saving is captured at the iron stage | Verpoort 2024 |
| Corroboration | 73–80% of green-H2-equivalent trade flows as commodities not H2 | IRENA 2025 |
| Contradiction | Ceará case argues export *steel*, not iron | Caiafa 2026 |
| Full relocation saving, steel, Δ€40/MWh | 18.3% | Verpoort 2024 |
| Same, Germany-specific, system model | 23.3% of industry consumer cost; HBI price −20.2% | Seibold 2025 |
| EU-internal (Spain) vs global (Australia/Brazil) | Europe-only captures **47.7%** of the global benefit; rises to **78.3%** if non-EU WACC is 10% rather than 7% | Seibold 2025 |
| Global average saving | only 2.2–3.9% | Bilici 2024 |
| Is renewables pull real? | Weaker than previously found; cancellable by modest subsidy; labour and capital cost matter as much | **Nykvist 2025** |
| Green iron freight | €4.44–39.27/t (Lund model); ~7% of cost marine (Devlin); €5.5/t for 2,000 km e-truck intra-EU (Seibold); <⅓ the volume of shipping ore+H2 (SteelWatch) | multiple |
| CBAM scope-2 gap, quantified | **12.9 Mt CO2** of imported embedded emissions avoided if CBAM covered scope 2 | Johnson 2025 |
| Free allocation offset for EU producers | DRI, pellets and electrolytic H2 all get free allocation 2026–2033, full phase-out 2034 | Johnson 2025 |
| Resource shuffling, named for iron | Yes — Samadi footnote 9 | Samadi 2023 |

---

## 4. The CBAM indirect-emissions exclusion — legal basis, current status, quantification

### 4.1 The operative provision is Article 7(1) + Annex II, not Annex IV

Regulation (EU) 2023/956: https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R0956

- **Art. 7(1)** determines embedded emissions; **Annex II** is headed exactly *"List of goods for
  which only direct emissions are to be taken into account, pursuant to Article 7(1)"*.
- **Annex II covers all of CN Chapter 72 ("Iron and steel") except 7202 (ferro-alloys) and 7204
  (ferrous waste and scrap)**, plus listed Chapter 73 goods, all listed aluminium goods, and
  **2804 10 00 – Hydrogen**. Regulation (EU) 2025/2083 later added electricity to Annex II.
- **CN 7203 (DRI / spongy ferrous products) is inside Annex II** → DRI and HBI are charged on
  direct emissions only.
- **Agglomerated iron ore (sinter/pellets, CN 2601 12 00) is NOT in Annex II** → it *does* carry
  indirect emissions. It is the only part of the iron and steel chain that does. If a model
  routes pelletising separately from reduction, the two steps are treated differently at the
  border.
- **Art. 3(22)**: *"'indirect emissions' means emissions from the production of electricity which
  is consumed during the production processes of goods, **irrespective of the location of the
  production of the consumed electricity**."*
- **Annex IV point 4.3** sets the default-value method for indirect emissions (EU grid EF,
  country-of-origin grid EF, or CO2 EF of price-setting sources) — but it only bites on the
  non-Annex-II goods.

**Recital 19, verbatim:**

> "The CBAM should also apply to indirect emissions. … The inclusion of indirect emissions would
> further enhance the environmental effectiveness of the CBAM … **Indirect emissions should,
> however, not be taken into account initially for the goods in respect of which financial
> measures apply in the Union that compensate for indirect emissions costs incurred from
> greenhouse gas emission costs passed on in electricity prices. Those goods are identified in
> Annex II to this Regulation.** Future revisions of the EU ETS in Directive 2003/87/EC and, in
> particular, revisions of the compensation measures of the indirect costs should be appropriately
> reflected as regards the scope of application of the CBAM."

The cross-reference is to **Article 10a(6) of Directive 2003/87/EC** — Indirect Cost Compensation
under the ETS State aid guidelines (C(2020) 6400, OJ C 317/5, 25.9.2020).

**The WTO rationale is not in the regulation.** The regulation gives only the ICC-overlap reason.
The WTO argument (charging imports for something domestic producers are *compensated* for would
breach GATT Art. III:2 national treatment) is made in the secondary literature — Johnson et al.
2025 being the clearest statement.

**Review clauses that matter:**
- **Art. 30(2)(a)(i)** requires the end-of-transition report to assess *"the possibility to extend
  the scope to: (i) embedded indirect emissions in the goods listed in Annex II"*. Recital 74 adds
  "as soon as possible".
- **Art. 30(6)(a)(vi)** requires the biennial report from 1 January 2028 to assess CBAM's impact on
  *"international trade, **including resource shuffling**"*. **Resource shuffling is a named
  statutory review item.**

**Official confirmation for hydrogen** — DG TAXUD *Guidance No. 5b, Sector-specific guidance
document on hydrogen*:

> "Indirect emissions were monitored and reported during the transitional period but, for the
> definitive period, hydrogen is included in Annex II to the CBAM Regulation … so indirect
> emissions are not taken into account in the embedded emissions of hydrogen for CBAM purposes."

https://taxation-customs.ec.europa.eu/document/download/04ceca3b-466e-4ed9-a1a5-a47c2e5e5be0_en?filename=Guidance+No.+5b+-+Sector-specific+guidance+document+on+hydrogen.pdf

### 4.2 European Commission, *Technical Study on Indirect Emissions in the CBAM* (June 2026)

Commissioned by DG TAXUD Unit C.2, Specific Contract TAXUD/2024/AO-03. Contractors:
**Ricardo Energy & Environment, Öko-Institut, LBST, Trinomics**. Announced **8 June 2026**.

- Announcement: https://taxation-customs.ec.europa.eu/news/technical-study-indirect-emissions-cbam-2026-06-08_en
- Final Report (Reyna-Bensusan, Healy, Robelin, Serpi, Skribbe, Seebach, Harthan, Gailhofer, Altmann):
  https://op.europa.eu/en/publication-detail/-/publication/cbff6069-5a3a-11f1-aa6d-01aa75ed71a1/language-en
  ISBN 978-92-68-40108-8, DOI 10.2778/5807122
- Task 1 (default emission factors): https://op.europa.eu/en/publication-detail/-/publication/348b382c-5a3d-11f1-aa6d-01aa75ed71a1/language-en
- Task 2 (actual emissions claims, PPAs): https://op.europa.eu/en/publication-detail/-/publication/912753c2-5a3b-11f1-aa6d-01aa75ed71a1/language-en
- **Task 3 (scope extension)** (Reyna-Bensusan, Robelin, Serpi, Hermann, Graichen):
  https://op.europa.eu/en/publication-detail/-/publication/b3445a37-5a3a-11f1-aa6d-01aa75ed71a1/language-en
  ISBN 978-92-68-40111-8, DOI 10.2778/1051014
- Annexes 4–6: https://op.europa.eu/en/publication-detail/-/publication/cb6ffb4e-5a3b-11f1-aa6d-01aa75ed71a1/language-en

**The definitive official statement:**

> "In particular, iron and steel (**other than agglomerated iron ore CN Code 2601 12 00**),
> aluminium, and hydrogen are exempted from paying for indirect emissions under the CBAM (**while
> indirect emissions reporting remains**), because these sectors are eligible for ICC under the
> State Aid Guidelines."
>
> "By contrast, cement and fertilisers require importers to declare and pay for both direct and
> indirect emissions, and agglomerated iron ore is treated similarly because it does not qualify
> for ICC under the state aid guidelines."

Note: indirect emissions are still **reported**, just not **paid for**. The data exists.

**Task 3 §3.1, the Commission's own contractors on the gap:**

> "the CBAM's current indirect emissions only covers a portion of the sectors under the
> Regulation. This limits the CBAM's effectiveness in comprehensively addressing carbon leakage
> risks and incentivising efficiency measures in electricity consumption."
>
> "**Non-EU producers in sectors excluded from the CBAM's financial obligations for indirect
> emissions face weaker incentives to adopt cleaner electricity sources** or invest in energy
> efficiency measures. At the same time, EU producers in Member States with little or no ICC
> remain exposed to higher carbon costs from electricity consumption."

**Task 3 conclusion:**

> "options delivering immediate and full coverage of indirect emissions score highest on
> environmental integrity, but also impose the greatest adaptation challenges… **Transitional
> approaches that link CBAM expansion to a gradual ICC phase-out** emerge as a potential
> compromise"

**Task 1 conclusion — how an exporter's grid would be scored if the exclusion were lifted:**

> "grid-average electricity emission factors emerged as the most operationally defensible basis
> for third-country default values… More complex options (e.g. price-setting or marginal
> approaches) introduce methodological and implementation risks that outweighed their theoretical
> advantages."

Stakeholders "widely preferred" **IEA-based national grid emission factors**; PPAs "considered
reliable only if strict criteria ensured geographic and temporal correlation and additionality."
Five Technical Solutions assessed, from immediate 100% coverage to coverage only after full ICC
phase-out.

**The official counter-argument to answer.** Task 3 argues that extending CBAM to indirect
emissions could **discourage electrification** in coal-grid countries. Worked example: coal grid
at ~0.9 kgCO2/kWh; heat pump at COP 3.0 gives 0.30 kgCO2/kWh of useful heat vs 0.22 kgCO2/kWh from
a gas boiler (EF ~0.202). Pricing indirect emissions makes the *less efficient* option cheaper.
Also flagged: "limited producer control over grid carbon intensity", "risks of technology
lock-in", and "**resource shuffling**".

### 4.3 Nothing changed in Omnibus I, and nothing changed in December 2025

**Regulation (EU) 2025/2083** (8 October 2025, "Omnibus I" CBAM simplification), OJ L, 2025/2083:
https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=OJ:L_202502083
- 50 t/yr de minimis mass threshold (excludes ~90% of importers, retains ~99% of emissions);
  hydrogen and electricity excluded from de minimis; payment accrues 2026 but due from 2027.
- **Its only Annex II change was to ADD electricity to the direct-emissions-only list.** Steel,
  aluminium and hydrogen untouched. It also *excluded* emissions from finishing/downstream
  manufacturing from scope — a narrowing.
- ICAP summary: https://icapcarbonaction.com/en/news/eu-adopts-simplifications-cbam-rules-ahead-compliance-phase-starting-2026

**The 17 December 2025 package** (the Art. 30 deliverables):
- **COM(2025) 989 final** — proposal amending Reg. 2023/956 on extension to downstream goods and
  anti-circumvention:
  https://taxation-customs.ec.europa.eu/document/download/f270fb87-7fbe-4149-ba48-e8a568179db3_en?filename=COM_2025_989_1_EN_ACT_part1_v8.pdf
- **SWD(2025) 988 final** Impact Assessment, Parts 1–2:
  https://taxation-customs.ec.europa.eu/document/download/71d4f753-4a2d-4367-bb10-c1967cb7f28d_en?filename=SWD_2025_988_1_EN_impact_assessment_part1_v4.pdf
  https://taxation-customs.ec.europa.eu/document/download/61b0d225-8b3f-42e0-996a-58c228a4bae7_en?filename=SWD_2025_988_1_EN_impact_assessment_part2_v4.pdf
- Review report on the transitional period published 16 December 2025.

**Content on indirect emissions: none.** The IA's three problems are downstream carbon leakage,
CBAM avoidance/circumvention, and ineffective treatment of electricity imports. Indirect emissions
appear only incidentally (transitional reporting, indirect PPAs, indirect price effects).

**Current status as of September 2026: the Commission has NOT proposed including indirect
emissions for steel.** The forward commitment is a **2027 assessment** of extension to further ETS
sectors, further downstream goods, or "indirect emissions from existing CBAM sectors". EU ETS
revision expected summer 2026; CBAM review 2027.
- https://icapcarbonaction.com/en/news/eu-cbam-enters-compliance-phase-and-outlines-path-ahead
- Mayer Brown (Mizulin, Vander Schueren, Geraets, Nosowicz de Chillaz), 18 Dec 2025 — confirms no
  scope-2 guidance for iron and steel in the package:
  https://www.mayerbrown.com/en/insights/publications/2025/12/european-commission-issues-cbam-operational-rules-and-proposes-downstream-extension-of-the-cbam-scope
- https://www.fieldfisher.com/en-be/locations/belgium/insights/upcoming-cbam-amendments-under-the-second-implemen

**Commission on the record that it will not fix this** — Laura Roberts, "EU set to keep indirect
emissions out of CBAM for metals; adds pre-consumer scrap as separate product", **Fastmarkets,
7 November 2025**, reporting an **ERCST webinar of 30 October 2025**.
https://www.fastmarkets.com/insights/eu-to-keep-indirect-emissions-out-of-cbam-for-metals/
Mirror: https://eurometal.net/eu-set-to-keep-indirect-emissions-out-of-cbam-for-metals-adds-pre-consumer-scrap-as-separate-product/

**Martin Becker, deputy head of unit, European Commission:**

> "Our main concern at the moment is metals, both steel and aluminium, that are both CBAM products
> where indirect emissions are not in scope."
>
> "**There is no intention in December to make a proposal to extend the scope to indirect emissions
> for these two sectors.**"

Aluminium producers (Norsk Hydro among them) lobbied *for* keeping the exclusion.

**Political blockage.** Sandbag records a **joint declaration of Belgium, Italy, France,
Luxembourg, Romania, Slovakia and Spain** (February 2025, French Ministry of Economy) opposing
early inclusion of indirect emissions pending "further in-depth sectoral analysis" and consistency
with ICC. Sandbag's read: *"The inclusion of indirect emissions … [is] likely to be proposed at a
later stage, i.e. **after 2030**."*

**Modelling implication:** scope-2 exclusion is the base case through 2030 with high confidence.
Inclusion belongs in a scenario tied to **ICC phase-out**, not to the 2026 or 2027 review.

### 4.4 Who has quantified the gap

**(1) Johnson et al. (2025) — the only steel-specific published number: 12.9 Mt CO2/yr.**
See §1.4 above. Nature Communications 16, 9087, 13 October 2025.
Companion paper by the same group, directly on the modelled configuration:
**Li, Åhman, Algers & Nilsson, "Decarbonizing the Asian steel industries through green Hot
Briquetted Iron trade", *Resources, Conservation and Recycling* 219, 108275, June 2025.**
https://doi.org/10.1016/j.resconrec.2025.108275

**(2) Sandbag / Konrad-Adenauer-Stiftung — country-level € figures for the indirect extension.**
Adrien Assous, Chloé Barré, Duncan Woods, Meili Vanegas-Hernandez, *"The EU CBAM: A Two-Way Street
to Climate Integrity?"*, **25 August 2025**.
https://sandbag.be/2025/08/25/the-eu-cbam-a-two-way-street-to-climate-integrity/
PDF: https://sandbag.be/wp-content/uploads/Sandbag-KAS-EU-CBAM-Report-.pdf

Additional annual CBAM fees **if indirect emissions were included for all CBAM goods**
(post-free-allocation phase-out, current product scope):

| Country | Extra annual fee |
|---|---|
| China | **€718 m** |
| India | **€615 m** |
| Russia | **€535 m** |
| Türkiye | **€225 m** |

Their electricity-intensity inputs (JRC):

| Route | GJ/t steel | ≈ MWh/t |
|---|---|---|
| BF-BOF | 0.39 | 0.11 |
| **DRI-EAF** | **2.42** | **0.67** |
| Scrap-EAF | 2.07 | 0.58 |

**Important caveat: those are natural-gas DRI-EAF figures covering the EAF and DRI auxiliaries,
NOT electrolysis.** Sandbag's model does not capture the H2-DRI case at 3.5–4.5 MWh/t. Their
country numbers are therefore a **large underestimate** for an electrolytic route.

**(3) Sandbag — "Extending the CBAM to indirect emissions", Adrien Assous, July 2025.**
https://sandbag.be/2025/08/01/why-the-cbam-should-cover-indirect-emissions/
PDF: https://sandbag.be/wp-content/uploads/2025.07.30-ICC-and-CBAM-indirect-emissions.pdf

> "As electrification accelerates and carbon prices continue to rise, indirect emissions will
> account for a growing share of industrial emissions. EU producers will face more indirect costs
> due to electrification and higher carbon prices, **while importers will use more electricity to
> exploit the lack of coverage of indirect emissions in the CBAM.**"

> "hydrogen can help reduce direct emissions from industries like steelmaking, but **if it is made
> from high-emission electricity, it could emit more CO2 than steel from blast furnaces**. In
> addition, it fails to incentivise third countries to decarbonise their electricity systems."

Fiscal: ICC cost **France €600 m** and **Germany €1.6 bn** in 2023 (Commission 2024 Carbon Market
Report). Says **only 14 Member States** operate ICC schemes — note this conflicts with Johnson et
al.'s "15 out of 27" and with Sandbag's own "remaining 16 countries" elsewhere in the brief; treat
the count as **14–15 of 27**.

Recommendations: extend CBAM to indirect emissions; **reform ICC to compensate only the non-fossil
share** of electricity and only in hours when fossil generation sets the marginal price in the
installation's bidding zone (which makes ICC and CBAM non-overlapping and dissolves the legal
objection); and mandate default values, or move to *induced marginal* emissions, to block PPA
shuffling.

**(4) Bellona Europa** — advocacy, no quantification.
- Anna Pijaca / Bellona EU, "A CBAM without indirect emissions? A half-built bridge to
  decarbonisation", **23 June 2025**:
  https://eu.bellona.org/2025/06/23/a-cbam-without-indirect-emissions-a-half-built-bridge-to-decarbonisation/
  Excluding electricity emissions for aluminium, iron & steel and hydrogen means "ignoring the
  lion's share of the climate footprint"; indirect emissions are "a significant share, if not the
  majority, of the carbon footprint" of these sectors. Distinctive argument: **EU marginal
  pricing** means even renewable-supplied EU producers carry an embedded carbon cost in their
  power price, which imports do not.
- Bellona EU, "Rethinking ICC and CBAM: can indirect emissions unlock a better carbon pricing
  system?", **13 May 2026** (reporting a 28 April 2026 webinar):
  https://eu.bellona.org/2026/05/13/rethinking-icc-and-cbam-can-indirect-emissions-unlock-a-better-carbon-pricing-system/
  **Outokumpu on record:** *"excluding Scope 2 emissions remains a structural weakness in CBAM,
  particularly for electricity-intensive sectors such as steel."* **WWF:** ~**€100 bn** channelled
  through free allocation including ICC over 2013–2021.

**(5) OECD** — notes CBAM covers scope 2 only for cement and fertilisers; does not quantify.
Dechezleprêtre, Haramboure, Kögel, Lalanne & Yamano, "Carbon Border Adjustments: The potential
effects of the EU CBAM along the supply chain", **OECD STI Working Paper 2025/02**, 27 Jan 2025
(rev. 4 Jun 2025). https://doi.org/10.1787/e8c3d060-en ·
https://www.oecd.org/content/dam/oecd/en/publications/reports/2025/01/carbon-border-adjustments_b9049067/e8c3d060-en.pdf

**(6) Nobody has published the H2-DRI/HBI per-tonne arithmetic.** The following is a **computed
anchor, not a citation**:

| Electricity | Grid EF | Uncounted CO2 | @ €75/t | @ €100/t |
|---|---|---|---|---|
| 3.5 MWh/t | 0.8 tCO2/MWh | 2.80 t/t | €210/t | €280/t |
| 4.0 MWh/t | 0.9 | 3.60 t/t | **€270/t** | €360/t |
| 4.5 MWh/t | 1.0 | 4.50 t/t | €338/t | €450/t |
| 4.0 MWh/t | 0.71 (India) | 2.84 t/t | €213/t | €284/t |
| 4.0 MWh/t | 0.58 (China) | 2.32 t/t | €174/t | €232/t |

For scale: the CBAM **DRI-EAF crude-steel benchmark is 0.481 tCO2e/t** and **BF-BOF is 1.370
tCO2e/t**. A coal-grid H2-DRI plant's *uncounted* indirect emissions (~3–4 tCO2/t) exceed the
*entire counted BF-BOF benchmark* by a factor of 2–3, while its CBAM liability is near zero.
Set against the Superpower Institute's estimated green-iron cost gap of AUD 110–850/t vs
international fossil HBI, an uncounted €210–340/t is the same order of magnitude as the gap the
subsidy debate is about.

**(7) Ancillary per-tonne figures in circulation — documentation-driven, not scope-2-driven; do
not conflate.** €87.94/t difference between actual (2.0 tCO2e/t) and default (3.167 tCO2e/t China
slab) at the Q1 2026 certificate price of €75.36/tCO2e; €649/t for Indonesian stainless coil at
default values; hydrogen grey-vs-green certificate gap of **€675–900/t H2** (10.4 tCO2/t default
grey vs ~0 for RFNBO green).
https://cbamguide.com/sectors/steel/ · https://cbamguide.com/sectors/hydrogen/ ·
https://www.miningsee.eu/steel-imports-to-the-eu-under-cbam-when-the-price-of-carbon-becomes-the-price-of-evidence/

Hasanbeigi, Springer & Chobthiangtham, "The Impact of the EU CBAM on Global Steel Trade",
**Global Efficiency Intelligence, March 2025** — import charges of **$72/t (Korea) and $83/t
(India) in 2030**, rising to **$210/t and $243/t by 2034**; 20–30% of steel price by 2034.
https://www.globalefficiencyintel.com/the-impact-of-the-eu-cbam-on-global-steel-trade

### 4.5 Default values and benchmarks — can a DRI/HBI importer declare near-zero? Yes.

**Transitional period** (Commission default values, Dec 2023, per IR (EU) 2023/1773):
**CN 7203 (DRI / spongy ferrous products): 4.81 tCO2e/t direct, 0.00 tCO2e/t indirect.**
The indirect column is literally zero for the whole of Annex II.
*(Verify against DG TAXUD's own file before citing — the mirror at
https://www.scribd.com/document/728297697/CBAM-Default-Values did not render the table.)*

**Definitive period** — **Commission Implementing Regulation (EU) 2025/2621** (adopted Dec 2025,
OJ 31 Dec 2025; corrected by IR (EU) 2026/1740 effective 31 Jul 2026). Crude steel route
benchmarks:

| Route | tCO2e/t |
|---|---|
| BF-BOF | **1.370** |
| DRI-EAF | **0.481** |
| Scrap-EAF | **0.072** |

Country defaults are country-average **plus a mark-up: +10% (2026), +20% (2027), +30% (2028 on)**.
Example: China slab **3.167 tCO2e/t**. Defaults to be reassessed by December 2027.
The per-CN-code table including 7203 is in DG TAXUD's "Default values definitive period" Excel
(reposted 10 Aug 2026) — **not retrieved; get it from the Commission CBAM legislation page.**
https://cbamguide.com/sectors/steel/ · https://cbamguide.com/sectors/steel/benchmarks/ ·
https://eurometal.net/eu-commission-finalizes-cbam-benchmarks-default-values-ahead-of-january-2026-launch/

**Why near-zero is available:**
- Actual verified emissions are the default pillar of CBAM (Art. 7(2)); defaults are the fallback.
  An H2-DRI installation's *direct* emissions are genuinely near zero.
- **Hydrogen is itself an Annex II good**, so its embedded emissions are direct-only. Electrolytic
  H2 has ~zero direct emissions **regardless of the electricity used**. The exclusion therefore
  compounds down the chain: coal electricity → electrolysis → hydrogen → DRI has zero counted
  emissions at every link.
- **An HBI/DRI cargo (CN 7203) produced by electrolysis on a coal grid can legitimately be
  declared at close to 0 tCO2e/t under CBAM.** No default, no benchmark, no MRV path captures the
  grid emissions.
- Counterweights: sintered/pelletised ore as a precursor *does* carry indirect emissions; and
  RFNBO certification (additionality, temporal and geographic correlation) separates genuinely
  green H2 — but **RFNBO status is a Renewable Energy Directive concept, not a CBAM requirement**.

---

## 5. "Resource shuffling" — the critique, and it is named in the regulation

### 5.1 Origin and academic transfer

- **Fowlie, Meredith; Petersen, Claire; Reguant, Mar (2021)**, "Border Carbon Adjustments When
  Carbon Intensity Varies across Producers: Evidence from California", *AEA Papers and Proceedings*
  **111**, 401–405, May 2021. https://doi.org/10.1257/pandp.20211073
  > "**Simulations suggest significant potential for leakage via resource shuffling. Realized
  > emissions outcomes indicate that this potential has not been fully realized.**"
- **Mehling, Michael A. & Ritz, Robert A. (2023)**, "From theory to practice: determining emissions
  in traded goods under a border carbon adjustment", *Oxford Review of Economic Policy* **39**(1),
  123–133. https://doi.org/10.1093/oxrep/grac043 — **the** paper on the default-vs-actual
  trade-off:
  > "**Requiring producers to demonstrate their actual carbon intensity captures additional
  > economic benefits of carbon pricing and improves the overall legal prospects of a BCA, but adds
  > to its administrative complexity and creates risk of avoidance practices such as 'resource
  > shuffling'.**"
- Fowlie & Reguant (2018), "Challenges in the Measurement of Leakage Risk", AEA P&P 108, 124–129.
  https://doi.org/10.1257/pandp.20181087

### 5.2 The Commission's own 2021 impact assessment defines it — and lists clean electricity first

**SWD(2021) 643 final, 14.7.2021, §5.2.1.10 "Resource shuffling":**
https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:52021SC0643

> "**Resource shuffling refers to the allocation or attribution of less emissions-intensive
> materials production … towards markets with higher carbon costs, while the overall carbon
> intensity of production in the home market remains constant.** There exist three main mechanisms
> through which resource shuffling can take place:
> – **Attribution of low-carbon input factors (low-carbon electricity, low-carbon heat, biomass) to
>   imported materials.**
> – Attribution of GHG emissions of a production process to co-products …
> – **Attribution of shares of recycled material to imported or exported goods.**
>
> **Incentives for resource shuffling exist for any emissions-related policy that includes traded
> goods … where the carbon intensity of imported or exported products does not rely on default
> values only, but on actual emissions.** … exporting low-carbon products to the EU would imply
> lowering the carbon costs these importers face and therefore **undermine the carbon leakage
> protection which the CBAM provides, without leading to a decrease of global emissions.**"

**OECD 2025 taxonomy** (Dechezleprêtre et al., pp. 32–33):

> "**Within-country resource shuffling – also known as 'backfilling' (Fowlie, Petersen and Reguant,
> 2021) – captures the idea that nations may domestically consume carbon-intensive products (or
> export them to other nations) while exporting cleaner alternatives to the CBAM market (Mehling
> and Ritz, 2023).** This aspect is not covered in the model … **The results of this paper can be
> interpreted as corresponding to a hypothetical version of the CBAM in which embedded emissions
> would be based on country of origin default values and not on actual plant level carbon
> intensity.**
> **Between-country resource shuffling** implies that importers choose the country of origin of a
> specific input based on its emission intensity. The model incorporates this …
> **Downstream reshuffling** is prompted by the CBAM's limited coverage of products …"

**The December 2025 impact assessment (SWD(2025) 988 Part 1) restates the installation-level
mechanism without naming it:**

> "This risk of circumvention **cannot be addressed by the accredited verifiers since a given
> installation in a third country can produce different grades of the same product** … **it is
> possible that the installation exports to the EU highly emission-intensive products, while the
> CBAM adjustment for the importers would be based on an actual value for emissions, which is too
> low compared to the carbon footprint of the product exported.** NCAs and Customs Authorities
> conveyed their concern over this point **as well as over 70% of 142 circumvention
> stakeholders**."

Remedies proposed (Options 1 and 2): empower the Commission to **further detail CN codes** to
capture within-code composition; attach **additional conditions to the use of actual emissions**
for high-risk good/origin combinations; require **proof of place of production** (Steel Mill
Certificate) before actual emissions may be used.

### 5.3 Quantified erosion — the two hard numbers

**(1) Commission JRC-GEM-E3 modelling, SWD(2021) 643, Annex 10, Table 10-5.**
Assumption: exporters claim emission intensities **50% lower for cement and iron & steel, 80%
lower for aluminium**, per **Stede, Pauliuk, Hardadi & Neuhoff (2021), DIW Discussion Paper 1935**.

2030, Option 4 (CBAM + free-allocation phase-out), with vs without shuffling:

| Metric | Option 4 | Option 4 **with resource shuffling** |
|---|---|---|
| **Carbon leakage, iron & steel (%)** | **−24** | **0** |
| Carbon leakage, cement & lime (%) | 7 | 13 |
| **Carbon leakage, aluminium (%)** | **−89** | **+8** |
| Imports, iron & steel (% vs baseline) | −11.98 | −2.38 |
| Imports, cement & lime | −15.12 | +6.97 |
| Imports, aluminium | −4.41 | +1.75 |
| **Revenue collected at the border (€bn)** | **2.1** | **1.3** |
| Total revenue (€bn) | 9.1 | 8.2 |

**Resource shuffling takes iron and steel from CBAM *reversing* leakage (−24%) to CBAM doing
nothing at all (0%)**, and takes aluminium to actively negative protection. Border revenue falls
**38%**. The Commission conceded: "in the event that the risk of resource shuffling materialises,
the reduction in imports induced by the CBAM could be substantially limited."

**(2) Sandbag/KAS (August 2025) — the modern estimate.**

| Scenario (post-free-allocation phase-out, current scope) | Gross CBAM fees | Net cost after price pass-through |
|---|---|---|
| Business as usual | **€11.3 bn/yr** | €5.0 bn |
| **Resource shuffling** | **€7.3 bn/yr** | **€995 m** |
| Third countries adopt €50 carbon price | €7.0 bn | €715 m (~0.07% of import value) |

**Resource shuffling wipes out 35% of gross CBAM fees and 80% of net cost, with zero emissions
benefit.**

Their shuffling assumptions (the replication recipe):
- Flat steel: EAF at **64% scrap / 36% DRI** where EAF capacity allows; otherwise max 55% scrap in
  EAF with the remainder BF-BOF. (Anchored on US flat-steel practice.)
- Aluminium imports: 80% remelted scrap. Cement: 20% clinker.
- **Indirect emissions: "resource shuffling enables trade partners to cover 50% of the electricity
  use through PPAs (based on load factors achievable through renewables), meaning that only half of
  the electricity is declared with country emission intensity."**

Country-level fee reductions from shuffling: **UAE −56%**, **Taiwan −40%**; China, Türkiye,
Montenegro, Indonesia, Vietnam, Japan, South Korea also large beneficiaries; **US and Canada only
−1% to −3%**. India drops from #1 to #3 *"because India produces enough steel scrap and DRI to
export goods with lower reported emissions"*. Ukraine is nearly immune (EAF ~5% of total, already
high scrap).

### 5.4 The DRI-specific reshuffling and free-allocation interaction

Sandbag flags a **structural DRI discount** distinct from shuffling but compounding it:

> "the EU ETS grants free emission permits to factories (which may not produce finished goods)
> whereas the CBAM will charge finished goods … **DRI (direct reduced iron) manufacturing earns a
> lot of free allowances in Europe, but European long steel mills typically do not use any DRI so
> they receive few free allowances. In contrast, CBAM fees on imported long steel made from DRI
> could be massively discounted, possibly down to zero**, if the discount applied strictly mirrored
> the EU ETS allocation regime."

**Sandbag, "CBAM DRI loophole requires new free allocation reform", 13 December 2024.**
https://sandbag.be/2024/12/13/cbam-dri-loophole/
The 2024 Free Allocation Regulation reform extended the **hot metal benchmark (1.248 tCO2/t)** to
DRI, whose actual EU emissions are **0.39 tCO2/t** (JRC) — a windfall of **0.858 EUA per tonne of
DRI**. Because CBAM liability is net of the benchmark, imported long steel made from natural-gas
DRI (~0.5 tCO2/t DRI) yields a **negative CBAM value**. Proposed fix: allocate on steel output
regardless of route (flat 1.123 EUA/t; long 0.145–0.176 EUA/t).

**Sandbag, "The EU CBAM gives a boost to Algeria's iron exports", 20 February 2026.**
https://sandbag.be/2026/02/20/eu-cbam-algeria-iron-exports/
**The best quantified demonstration that CBAM actively rewards importing DRI/HBI:**

| Quantity | Value |
|---|---|
| Algerian long steel exported to EU (2024) | 350,000 t/yr |
| Algerian DRI exported to EU (2024) | 230,000 t/yr |
| Total annual CBAM net cost, all Algerian goods | €136m |
| Long steel net **cost** | €29m |
| **DRI net PROFIT** | **€54 per tonne** |
| DRI annual gain on 230 kt | €12m |
| EU BF carbon cost raises the DRI market price by | **€94/tonne** |
| Assumptions | EUA €80, pass-through 80% |

Switching 350 kt of long-steel exports to DRI turns "a €29m net loss into a €19m net gain."
**This €54/t arises from benchmark mechanics and carbon-cost pass-through, and stacks on top of
the indirect-emissions exclusion.** Algerian DRI is gas-based; the analysis never touches scope 2.
For a coal-electricity H2-DRI exporter the two effects **add**.

### 5.5 Related framings and the anti-shuffling instrument

- **PPA / book-and-claim reshuffling of electricity** — where scope 2 and resource shuffling meet.
  Sandbag: *"It is therefore **possible for third country exporters to reduce their CBAM fees by
  signing PPAs or connecting their factories to existing RES, without benefitting the climate**.
  This could be avoided if default values was the only possible option."* Their alternative is
  **"induced emissions"** — pricing on the *marginal* grid emission intensity at the hours the
  plant runs: *"if the plant functions at times of high electricity demand, it will create extra
  demand for fossil electricity generation even if it uses a dedicated RES, because that RES ends
  up not meeting the broader demand."*
  The Dec 2025 IA closed the crudest version for *imported electricity*: **virtual PPAs are
  systematically excluded**; only physical (direct and, newly, sleeved) PPAs count. But that
  discipline does not touch the electricity a third-country steel plant consumes domestically —
  and for Annex II goods that electricity is not priced at all.
- **The scrap loophole** — the dominant framing for steel. Sandbag: *"Currently considered as
  emissions-free, **metal scrap is currently the main enabler of 'resource shuffling' for steel and
  aluminium imports**."* The Dec 2025 IA devotes a sub-section to it, citing **European Aluminium's
  Ramboll study (March 2025)**, **Hydro** ("CBAM: Europe's low-carbon aluminium is threatened by a
  big loophole", ~April 2025,
  https://www.hydro.com/en/global/about-hydro/stories-by-hydro/greenwashing-via-cbam-loophole-s-threaten-european-green-products-market/)
  and **Alcoa**. COM(2025) 989 adopts **pre-consumer scrap only** as a precursor.
- **"Greenwashing of exports"** — Hydro's term for the aluminium version.
- **Standardised/default values as the anti-shuffling instrument.** Mehling & Ritz frame the
  trade-off; Sandbag operationalises it. See also "Industrial decarbonization in a fragmented
  world: Carbon pricing with border adjustments using standardized values", *Energy Policy*,
  https://www.sciencedirect.com/science/article/pii/S0301421526003393 — argues production-specific
  values create higher circumvention risk than standardised values, and that the December 2025
  proposal's *voluntary* defaults leave the shuffling risk in place because importers keep the
  right to use production-specific values.
- **Sandbag's structural fix** — "An opt-in solution to CBAM circumvention and complexity",
  **7 April 2026**: https://sandbag.be/2026/04/07/position-paper-cbam-amendment/ — let third
  countries **opt in to systematic country-average default values**. Rationale: verified emissions
  "sometimes fails to capture the wider impact of imports on the overall emissions of the country
  of origin"; systematic defaults make countries "reduce their overall emissions rather than just
  reallocating them from EU imports to other markets."
  **Revealed-preference evidence:** in the transitional data, **only ~5% of declarants chose actual
  data through Q2 2024**, and even when defaults were capped at a 20% mark-up in Q3 2024, **~50%
  still claimed defaults**. The actual-emissions route is used selectively — i.e. when it pays.
- Other Sandbag: "Strengthening the CBAM by Default", 6 Aug 2025,
  https://sandbag.be/2025/08/06/strengthening-the-cbam-by-default/ (cited in the Commission's
  Dec 2025 IA, footnote 23) · "Mind the scrap: CBAM at risk", 21 Feb 2023,
  https://sandbag.be/2023/02/21/ignoring-embedded-emissions-puts-the-cbam-at-risk/ ·
  "CBAM Extension: Closing the Emissions Gap", 22 May 2025,
  https://sandbag.be/2025/05/22/cbam-extension-closing-the-emissions-gap/
- **European Parliament written question E-001119/2025**, "Addressing resource shuffling and
  ensuring the CBAM's effectiveness in the European steel industry":
  https://www.europarl.europa.eu/doceo/document/E-10-2025-001119_EN.html — **text not retrieved**
  (EP site blocks automated fetch). The only EP-level document naming resource shuffling and steel
  together. Worth retrieving manually.
- **Energy Innovation / iGDP / IFS**, "China and the European Union's Carbon Border Adjustment
  Mechanism": https://energyinnovation.org/wp-content/uploads/China-and-the-EUs-Carbon-Border-Adjustment-Mechanism.pdf
  > "resource shuffling … results in **lower carbon intensity for export production but a largely
  > unchanged energy system** … **policy design cannot eliminate it completely**."
  Cost estimates for China: €146 m net / €174 m direct in 2026; €208 m / €484 m in 2035;
  Kuusi et al. €905 m at €25/t and €2,174 m at €60/t.
- **Carbon Market Watch**, "Proposed CBAM reforms serving industrial lobbies need climate refocus",
  **18 December 2025**: https://carbonmarketwatch.org/2025/12/18/proposed-cbam-reforms-serve-industrial-lobbies/
  · CBAM FAQ: https://carbonmarketwatch.org/2024/07/30/faq-the-eu-carbon-border-adjustment-mechanism/

### 5.6 A gap in the literature

No dedicated 2024–26 analysis ties CBAM reshuffling to **EU scrap export restrictions**. The
symmetric argument — EU scrap exported out, low-emission scrap-based product re-imported in — is
implied by the mechanism but nobody has modelled it.

## 6. Think-tank work on green iron / HBI trade — and its silence on scope 2

### 6.1 Agora Industry (Sept–Nov 2025) — the reference work, on the same modelling stack

**"The role of green iron trade in accelerating competitive steel transformation"**
Camilla Oliveira (project lead), Leandro Janke, Darlene D'Mello, Niklas Wagner, Ysanne Choksley,
Karina Marzano, Julian Somers, Zaffar Hussain. **17 September 2025.**
https://www.agora-industry.org/publications/the-role-of-green-iron-trade-in-accelerating-steel-transformation
Press: https://www.agora-industry.org/news-events/europes-steel-sector-can-cut-costs-preserve-jobs-and-accelerate-decarbonisation-through-green-iron-trade

Six country/region deep dives:
- **EU + German case study**, v1.0 Sept 2025:
  https://www.agora-industry.org/fileadmin/Projekte/2024/2024-28_IND_Green_Iron/A-IND_377_Green_iron_trade_EU_WEB.pdf
- **Brazil, Nov 2025**:
  https://www.agora-industry.org/fileadmin/Projekte/2024/2024-28_IND_Green_Iron/A-IND_377_Green_iron_trade_Brazil_WEB.pdf
- South Africa, Nov 2025: `…/A-IND_377_Green_iron_trade_South_Africa_WEB.pdf`
- Japan: `…/A-IND_377_Green_iron_trade_Japan_WEB.pdf`
- China: `…/A-IND_377_Green-iron-trade_China_WEB.pdf`
- South Korea: `…/A-IND_377_Green-iron-trade_South_Korea_WEB.pdf`

**Methodological note of direct relevance:** Agora's green iron modelling runs on the
**PyPSA-BOA model** plus the **PtX Business Opportunity Analyser v2.0** (Oeko-Institut, Agora
Energiewende and Agora Industry, 2024). Same tooling family as a PyPSA-based green-iron model —
their Brazil numbers are a like-for-like benchmark rather than a loose comparison.

**Headline numbers:**

| Finding | Value |
|---|---|
| EU steelmaking cost reduction from green HBI imports by 2040 | **12–15%** |
| — from MENA sourcing | 12% |
| — from Australia, Brazil or South Africa | 15% |
| Imported green HBI production-cost advantage | **up to 24% lower** than domestic |
| Announced EU H2-DRI capacity by 2030 | **34 Mt** |
| — at FID or under construction | **12 Mt** (≈ half of EU low-carbon iron need) |
| EU virgin iron requirement | **50–60 Mt** |
| EU crude steel production 2024 | 130 Mt (from ~152 Mt in 2021) |
| EU steel share of EU emissions | ~5–6%; **>25% of industrial emissions** |
| EU production split | 55% BF-BOF / 45% EAF; scrap could exceed 50% |
| Average EU blast furnace age | 50 years; 32 BFs have announced retirement dates |
| Current global DRI trade (2023) | **11.1 Mt**, mostly HBI — 8% of DRI output, <1% of iron production |
| Projected EU carbon price | **€132 by 2030, €194 by 2045** |
| De-risked cost of capital assumed for exporters | **4.3%** |
| Share of steel-sector employment in steelmaking + finishing | **~90%** |
| BF-BOF emissions | ~2.2 tCO2/t steel |
| H2 + high-grade ore share of H2-DRI cost | **>half** |

**Brazil deep dive (Nov 2025, co-produced with Instituto E+ Transição Energética):**

| Item | Value |
|---|---|
| **HBI production cost, Brazil 2040** | **509 USD/t HBI** |
| **HBI selling price** | **654 USD/t HBI** |
| Payback / IRR | **8 years / 15%** |
| Annualised investment | 1.3 USD bn/yr (20 yr life, 4.3% discount) |
| Iron ore export revenue today → green iron 2040 | **USD 2.9 bn/yr → USD 6.6 bn/yr** |
| At 10 Mt/yr green iron by 2040 (5% of projected 200 Mt global iron demand 2050) | **12.8 Mt CO2 avoided, ~35,500 jobs** |
| Scenario set | 5 Mt (2.5% of global iron demand) / 10 Mt (5%) / 15 Mt (7.5%); CO2 avoided 6.4 / 9.4 / 12.8 / 18.8 / 19.1 / 28.2 Mt across variants |
| Brazil iron ore production 2023 | 582 Mt — **Minas Gerais 392 Mt at ~62% Fe; Pará 175 Mt at ~65% Fe**; exports 378 Mt = USD 30.6 bn; 85–90% exported, 65–70% of that to China |
| Brazil crude steel 2023 | 32 Mt |
| DR-grade pellet cost | **207 USD2024/t** |
| DRI capex | **633 USD2024/t DRI-yr** |
| Crude steel capex, DRI / BF route | 468 / 326 USD2024/t CS-yr |
| Coking coal | 257 USD2024/t |
| Alkaline electrolyser capex / fixed opex | **657 USD2024/kWel** / 13 USD2024/kWel-yr |

Note the Pará ore at **~65% Fe** vs Minas Gerais at ~62% — relevant to a DR-grade beneficiation
assumption for a Brazilian plant, and a real advantage over Pilbara ore at 56–62% Fe.

Cost driver, verbatim: *"Green HBI production costs are mainly driven by **cost of capital in
potential exporting countries** and by hydrogen costs from high energy costs in potential importing
countries."*

2040 producer ranking: cheapest **Brazil and South Africa** (solar+wind mix, domestic DR-grade ore,
lower labour costs), then **Australia and Saudi Arabia** (lower capital costs offsetting ore
quality and logistics). Structurally high-cost: **Japan, South Korea, Germany**. China has ~2.5 Mt
H2-DRI capacity now, could be an exporter by 2030 and a net importer by 2040.

**Why iron rather than steel or hydrogen — Agora's three reasons:**
1. **Transportability**: HBI is "a solid material that can be cooled and shipped globally, much
   like iron ore," avoiding the efficiency losses of shipping liquid hydrogen or ammonia.
2. **Value-chain disaggregation**: "the DRI-EAF steelmaking routes allow the iron and steelmaking
   steps to be disaggregated" — the energy-intensive step relocates, the labour-intensive step
   (~90% of jobs, higher gross value added) stays.
3. **Infrastructure**: no new transport infrastructure beyond existing ports and rail.
Plus: *"This reduces the demand for domestic or imported H2 and associated renewable energy and
infrastructure."* Framework adapted from **Verpoort et al. (2023/2024)**.

**CBAM treatment — the notable silence.** In 24 pages of the EU deep dive, **CBAM appears exactly
once**, as a bullet: *"Integrate green iron into trade agreements, strengthen CBAM and align it
with global climate goals."* In the Brazil deep dive, once: *"Use harmonised green product
criteria, trade agreements and CBAM to align energy, industrial and climate policies."* Plus once
in each abbreviation list. **No discussion of indirect or scope-2 emissions anywhere.** The press
release does not mention CBAM at all.

**Technology scope:** the Brazil deep dive's abbreviation list includes **MOE (molten oxide
electrolysis)** and **AEL (alkaline iron electrolysis)**, so Agora's technology scope overlaps a
MOE/electrowinning model — but the cost modelling is H2-DRI-centric. Agora/Wuppertal put
**MOE market readiness at 2035 and AEL at 2040**.

**Three-phase EU strategy** (the policy spine): Phase 1 — targeted support for a first wave of
domestic H2-DRI replacing blast furnaces; Phase 2 — EU single-market value chains and intra-EU
green iron trade; Phase 3 — global green iron imports from renewable-rich regions.
Instruments: **CTIPs (Clean Trade and Investment Partnerships)**; broadening **H2Global's
double-auction to green iron**; IPCEI electrolyser funding; European Hydrogen Bank; de-risking via
ECAs/MDBs; harmonised standards; and "**Consider defining green HBI as a strategic raw material
under the Critical Raw Materials Act (CRMA)**."

Other Agora assets:
https://www.agora-industry.org/data-tools/global-steel-transformation-tracker ·
https://www.agora-industry.org/publications/low-carbon-technologies-for-the-global-steel-transformation ·
https://www.agora-industry.org/news-events/brazils-green-iron-opportunity-could-drive-jobs-trade-and-global-steel-transformation ·
https://www.agora-industry.org/international/countries-and-regions/brazil ·
https://www.agora-industry.org/publications/12-insights-on-hydrogen-brazil-edition

### 6.2 Agora + Wuppertal, "15 insights on the global steel transformation" (15 June 2023)

https://www.agora-industry.org/publications/15-insights-on-the-global-steel-transformation
PDF: https://www.agora-industry.org/fileadmin/Projekte/2021/2021-06_IND_INT_GlobalSteel/A-EW_298_GlobalSteel_Insights_WEB.pdf
Witecka, von Eitzen Toni, Somers, Reimann (Agora); Zelt, Jülich, Schneider, Lechtenböhmer
(Wuppertal).

The origin of the "trade iron not hydrogen" thesis: *"Transporting embodied hydrogen as green iron
proves significantly cheaper than shipping hydrogen directly."* Also: net-zero steel technically
feasible by the early 2040s; **DR-grade ore demand will exceed supply by 2030** absent urgent
action; 2030 low-carbon hydrogen pipelines fall short of steel demand in all scenarios; **CCS on
blast furnaces "will not play an important role"** — caps reductions at ~73% and ignores upstream
coal methane. No CBAM/scope-2 content.

### 6.3 IEEFA (Simon Nicholas) — the best iron-vs-hydrogen shipping numbers anywhere

"Australia faces growing green iron competition from overseas", **IEEFA, September 2023**, 41 pp.
PDF (works — the HTML page is Cloudflare-blocked):
https://ieefa.org/sites/default/files/2023-08/Australia%20faces%20growing%20green%20iron%20competition%20from%20overseas_Sep23.pdf

| Finding | Value |
|---|---|
| **H2 transport cost via carriers by 2030** | **US$2.5–4.5 per kg H2 delivered** |
| At **51 kg H2 per tonne of steel** | **+US$128–230 per tonne of crude steel** |
| **HBI weight advantage** | DRI/HBI is **~30% lighter than iron ore pellets** (the oxygen is gone) |
| HBI metallisation | **>90%**; "the best ore-based metallic that can be used in EAFs" |
| DRI ore grade requirement | typically **≥67% Fe**; **Pilbara commercial deposits are 56–62% Fe** |

Structure: Figure 2 (hydrogen-carrier conversion/reconversion energy losses); **Figure 3 "HBI
Export vs Iron Ore and Hydrogen Export"**; Figure 4 (transport makes imported H2 expensive even at
US$1/kg production cost); Figure 5 (H2 DRI-EAF production cost 2050, US$/t). Hubs assessed:
Australia, Brazil, Middle East, Africa.

Verbatim, from a section headed **"Export Green Iron Not Green Hydrogen"**:

> "Given the inefficiency and cost of shipping green hydrogen or green ammonia, it makes sense to
> use more of this planned production domestically to produce value-added products like green iron
> from iron ore."
>
> "Green iron – made via DRI and exported as hot briquetted iron (HBI) – is now being considered
> for import by major global steelmakers **instead of importing both iron ore and green hydrogen**."
>
> "many nations are likely to be reticent to fully offshore their steelmaking capacity and may
> strategically prefer to import green iron that could be processed into steel domestically via
> EAFs that could be powered by renewable energy."

Quotes **Kajsa Ryttberg-Wallgren (EVP, H2 Green Steel)**: *"You will be out of cost. So they want
to buy green iron, or HBI … from places in the world where it makes sense to produce it."*
**Kobe Steel** plans to import green iron from H2 Green Steel's Swedish plant rather than import
hydrogen. Cites **MRIWA**: "pathways which involve the development of intermediate iron products,
such as HBI, are the most prospective for Western Australia."

Investment signals: **POSCO** considering US$40 bn in Australia (US$28 bn green H2 + **US$12 bn
green HBI production and export**); **Nippon Steel** ~US$700 m, considering Australia *or Brazil*;
**China Baowu** looking at WA but also South America, Africa, Middle East, with an Aramco/PIF Saudi
DRI agreement. **Vulcan Green Steel**: 5 Mt H2-DRI targeting Middle East, Europe and Japan, plus a
second 5 Mt H2-ready DRI/HBI project targeting Asia and Europe. **Russia (Metalloinvest) was the
largest HBI supplier to Europe in 2021**; **LKAB** is Europe's key magnetite pellet supplier.

**No CBAM or indirect-emissions discussion.** IEEFA's steelmaking carbon intensities are the source
the Superpower Institute used for its BF-BOF / DRI-EAF / scrap-EAF figures.

Further IEEFA items identified but not retrievable (Cloudflare on HTML; PDFs at
`ieefa.org/sites/default/files/` do fetch if you have the exact path):
- https://ieefa.org/resources/which-dri-solutions-can-deliver-green-iron-australia
- https://ieefa.org/august-2025-australian-green-iron-tracker
- https://ieefa.org/articles/australia-must-act-quickly-overseas-competition-green-iron-grows
- https://ieefa.org/articles/south-australia-leads-nation-towards-green-iron-amid-growing-global-competition
- https://ieefa.org/resources/australias-iron-ore-sector-crossroads-business-usual-or-time-embrace-green-iron
- https://ieefa.org/resources/can-australia-keep-pace-evolving-green-iron-market
- https://ieefa.org/resources/2025-australia-needs-accept-its-iron-ore-and-coal-markets-are-changing-permanently-and

Secondary reporting: Australian green iron export market **~AU$300 bn**; **only two pilot-scale
DRI/HBI projects under construction, no large-scale project at FID**.

### 6.4 The Superpower Institute — the most detailed cost modelling, and the most CBAM-engaged

Ross Garnaut / Rod Sims institute. **"A Green Iron Plan for Australia: Securing prosperity in a
decarbonising world", May 2025.** Contributors include **Ingrid Burfurd** (Carbon Pricing and
Policy Lead); modelling in partnership with **Bivios**. Underlying analysis: Finighan,
*The New Energy Trade*.

- Landing: https://www.superpowerinstitute.com.au/work/green-iron-plan
- **Full PDF (13 MB):** https://www.superpowerinstitute.com.au/resource/file-f8963eebb1fe8a9d1e0ae4a15f6c9b994ab7eccc-pdf/TSI_A-Green-Iron-Plan-for-Australia_May-2025-06.pdf
- Exec summary: https://www.superpowerinstitute.com.au/resource/file-0c29f1c4d244e8986808ec35cb651d38a2c59fb0-pdf/TSI_Executive-Summary_A-Green-Iron-Plan-for-Australia_May-2025.pdf
- APO record: https://apo.org.au/node/330769
- Launch speech, Asst Minister Andrew Leigh: https://ministers.treasury.gov.au/ministers/andrew-leigh-2025/speeches/address-launch-superpower-institutes-green-iron-plan-australia

> "the future energy trade will not be dominated by fossil fuels, but by trade in goods that embody
> clean energy. Energy-intensive industries will migrate to regions where cheap renewable energy
> exceeds domestic needs. Australia is one of those rare regions."
>
> "**Australia's comparative advantage is stronger in iron-making than steel-making, because it is
> the more energy-intensive process.**"

All figures AUD 2024 unless noted:

| Item | Value |
|---|---|
| Green iron export potential by 2060 | **$386 bn/yr** (vs ~$120 bn/yr iron ore today) |
| Emissions abatement potential | **~4% of global emissions** — >3× Australia's domestic emissions |
| Coal + gas exports at risk | ~$120 bn/yr (coal ~$70 bn, LNG ~$50 bn) |
| **Carbon price used (EU ETS 2030 forecast)** | **$155/tCO2** (range $110–225) |
| **HBI traded price** | **$345–712/t**, 5-yr weighted average **$554/t** |
| **Pig iron traded price** | **$690–924/t**, 5-yr weighted average **$779/t** |
| Illustrative BF pig iron price used | $400/t; BF-BOF hot-rolled band ~$640/t (USD440 @ 1.45) |
| Australian gas-based DRI modelled cost | ~**$570/t** |
| **Cost gap, green vs Australian gas DRI** | **~$100 to >$800/t** |
| **Cost gap vs international fossil HBI ($554)** | **$110–850/t** |
| Cost gap vs pig iron at $400 | **$270 to >$1,000/t** |
| Emissions intensities used | **2.0 tCO2e/t pig iron; 1.1 tCO2e/t fossil DRI; 0.518 tCO2e/t Australian gas HBI** |
| Carbon price effect at $155/t | pig iron **+>$300/t**; international fossil DRI **+~$170/t**; Australian gas HBI **+~$80/t** |
| **Recommended green iron production tax credit** | **≥$170/t of green iron in 2030** (inclusive of the Hydrogen Production Tax Incentive) — "would have a very similar effect to a carbon price" |
| Recommended capital support | up to **30% of investment cost** (15% tax benefit + 15% grants), drawing on the announced **$1 bn green iron investment fund** |
| Renewables capex | Pilbara inflexible ~**$27 bn**; Eyre Peninsula flexible ~**$10 bn** |
| Electricity for one plant | **>10 TWh** ≈ two-thirds of South Australia's 15.7 TWh (2024) generation |
| Grid-connection constraint effect | raises cost from ~**$1,000/t to >$1,200/t** |
| **Capex share of green iron cost** | **>60%** |
| FOAK cost penalty (Eyre Peninsula) | **+$62/t** (~10%) |
| Water | 0.1–0.14% of total cost — negligible |
| Eyre Peninsula advantage | 75% of South Australian power already renewable |

**Sites modelled:** Pilbara (NW WA), Geraldton (mid-west WA), Kwinana (SW WA), Eyre Peninsula (SA),
Gladstone (QLD).

**Counter-intuitive finding worth carrying into a model:** *"the Pilbara is unlikely to be one of
Australia's lower-cost locations… It may make economic sense to **ship ore from the Pilbara to
other locations in Australia** where green iron can be produced more cheaply."*

**Flexible (ramping) vs inflexible (continuous) technology is the single biggest cost lever**, and
a grid connection helps by letting the plant arbitrage — the same dispatch question a PyPSA model
resolves endogenously.

**Third-party cost-gap comparators — a natural sensitivity range.** European industry estimates:
~$760/t conventional BF-BOF steel, ~$950/t grey HBI-BOF steel (+$190), **~$1,350/t green-DRI-BOF
steel (+$590)**. TSI: *"Our production cost gap is larger than industry consensus of $0 to $150 per
tonne of green HBI."* Required carbon price for competitiveness: **TSI AUD 155/t (2030 EU forecast)
· SteelConsult AUD ~240/t (by 2035) · CRU AUD ~450/t (2030)**.

**Chapter 6.1 is entirely about the EU CBAM** ("The EU CBAM demonstrates how a carbon price creates
demand for green iron"):

> "with the CBAM helping to create a new market for green iron, **EU companies will need to import
> up to 13 million tonnes of green iron by 2030, and up to 18 million tonnes by 2045.**"

Also: the Commission expects ~30% of EU primary steel decarbonised with renewable hydrogen by 2030;
the CBAM price ramps to the full ETS price from 2034; Australia exported **<AUD 30 m** of iron ore
to the EU in 2023, so an EU market for Australian green iron would be new rather than substitution.

**Certification (§6.1.1):** Australia's **Guarantee of Origin** scheme (legislated late 2024, in
force late 2025) — REGO + PGO certificates — is being extended to green iron, steel, aluminium and
liquid fuels, with **explicit alignment to EU CBAM** as a design goal; accreditors able to report
embedded carbon from 2026. Recommendation: align GO with CBAM "at the earliest" opportunity.

**The critical gap.** The report **never mentions the indirect-emissions exclusion** and treats
CBAM as an unqualified demand signal. It applies **$155/tCO2 to the full carbon intensity** of
competing iron — a price CBAM does not levy on iron's electricity emissions. And because Pilbara
sits on a gas/coal-heavy grid while the Eyre Peninsula is 75% renewable, **the scope-2 exclusion
means CBAM cannot distinguish those two Australian sitings at the EU border** — which undercuts
TSI's own siting analysis. This is the sharpest available wedge.

Three market failures and the policy asks: (1) unpriced emissions → the **≥$170/t production tax
credit**; (2) under-provision of common-user infrastructure (transmission, pipelines, storage,
ports) → public investment; (3) innovation spillovers / early-mover risk → **capital support up to
30%**. TSI explicitly ranks an international carbon price with border adjustments as first-best and
subsidies as "a second-best option" that "simulate the effect of a carbon price."

Also cited: **ArcelorMittal and H2 Green Steel** — the two largest investors in European green iron
plants — "have both remarked that **Europe will not be able to produce most of its own green iron**."
Related: https://www.superpowerinstitute.com.au/the-opportunity ·
https://www.superpowerinstitute.com.au/news/building-blocks-for-a-green-iron-industry

### 6.5 Grattan Institute — a correction, and a genuine disagreement

**There is no Grattan report called "Go for green iron."** What exists:

- **"Start with steel: A practical plan to support carbon workers and cut emissions",
  Grattan Institute, May/June 2020.** https://grattan.edu.au/report/start-with-steel/ ·
  PDF: https://grattan.edu.au/wp-content/uploads/2020/05/2020-06-Start-with-steel.pdf ·
  https://grattan.edu.au/news/australians-want-industry-and-theyd-like-it-green-steel-is-the-place-to-start/
  Assesses aviation fuel, ammonia and steel; picks **green steel** (not iron). Capturing
  **~6.5–7% of the global steel market** → ~**$65 bn/yr** export revenue and **25,000
  manufacturing jobs** in Queensland and NSW. Australia is 38% of iron ore production but 0.3% of
  traded steel.
- "Green metals: Delivering Australia's opportunity", consultation paper, **July 2024**:
  https://grattan.edu.au/wp-content/uploads/2024/07/Green-metals-consultation-paper-2024.pdf
- Tony Wood & Alison Reeve, "How to forge a Future Made in Australia":
  https://www.aph.gov.au/DocumentStore.ashx?id=4cfea496-bc62-49fc-861a-67ca4e00aff1&subId=760974

**Grattan argues for green *steel* export; the Superpower Institute argues for green *iron*
export.** That is a genuine, citable disagreement between two Australian think tanks on exactly the
question being modelled. Related figure from the Grattan orbit: converting **40% of Australian iron
ore to DRI** with renewable power and green hydrogen could roughly double export revenues from
AUD 138 bn to **AUD 250 bn**.

### 6.6 E3G, RMI, LeadIT, T&E, Climate Strategies

- **E3G Steel Policy Scorecard 2025**, Aleksandra Waliszewska, Laith Whitwham, Johanna Lehne &
  Mark Hagen, Briefing Paper, **July 2025** (data cut-off 15 May 2025; G7 plus key non-G7
  producers; methodology unchanged since 2023 for comparability).
  https://www.e3g.org/publications/e3g-steel-policy-scorecard-2025-maintaining-the-momentum/ ·
  PDF: https://www.e3g.org/wp-content/uploads/E3G-Briefing-Steel-Policy-Scorecard-2025.pdf
  Recommendation #5, "Optimise: Build green iron corridors":
  > "The G7 should use trade and industrial policy to support efficient, low-emissions global supply
  > chains. That includes **enabling green iron production in renewable energy-rich countries and
  > facilitating offtake agreements to supply downstream steelmaking hubs**. Supporting this
  > international division of labour – while ensuring benefits for both producer and buyer nations
  > – can reduce costs, cut emissions and enhance global competitiveness."
  On CBAM, only MRV capacity-building for developing countries. **No indirect/scope-2 discussion,
  no cost numbers.** Headline: ambition up, delivery insufficient; near-zero steel projects
  stalling or cancelled in most G7 countries. Companion: **"The State of the European Steel
  Transition (2025)"**, E3G / CAN Europe / Bellona Foundation — the source of Agora's "32 blast
  furnaces with announced retirement dates".
- **RMI, "Brazil's Green Iron Opportunity"** — Rachel Wilmoth, Ariane DesRosiers, Thanh Ha,
  Chathurika Gamage, 2025. https://rmi.org/resources/brazils-green-iron-opportunity/
  Renewable hydrogen is **up to 50% of the final cost** of H2-DRI steel; Brazilian production
  potentially **"65% cheaper than other steelmaking countries."** Full report form-gated.
  **No CBAM or scope-2 discussion.**
- **LeadIT / SEI, "Demands for renewable hydrogen and electricity to drive the EU's green iron and
  steel transition"**, Leadership Group for Industry Transition Secretariat, Stockholm, 2025.
  https://www.sei.org/publications/demands-renewable-hydrogen-electricity-iron-steel-transition/
  (Cloudflare-blocked; not retrieved) ·
  press: https://www.sei.org/about-sei/press-room/hydrogen-and-renewable-electricity-for-the-eu-steel-and-iron-transition/ ·
  https://www.sei.org/projects/leadership-group-for-industry-transition-leadit/
  **Highest-value unretrieved item** — it sizes the "how much renewable power would the EU need if
  it *didn't* import iron" counterfactual.
- **Transport & Environment**: nothing on CBAM/steel/green iron. Current output is critical raw
  materials and circularity —
  https://www.transportenvironment.org/articles/europes-waste-creates-value-elsewhere (31 Jul 2026) ·
  https://www.transportenvironment.org/articles/eu-crm-centre-how-to-secure-critical-minerals (4 Sep 2026).
  **Not a player on this topic.**
- **Climate Strategies**: https://climatestrategies.org/publications/ — nothing found on
  CBAM/steel/green iron; current focus is just transitions and Asia energy.

### 6.7 Brazilian organisations

- **Instituto E+ Transição Energética** is a **co-author** of the Agora Brazil deep dive — that is
  the substantive Brazilian institutional contribution found.
- **Instituto Talanoa** (https://institutotalanoa.org/en/), **CEBRI**, and **Instituto Clima e
  Sociedade / ICS** (https://climaesociedade.org/): **no green iron or CBAM publication found**.
  Talanoa is ICS-funded and climate-policy-focused but not evidently active on this file.
  A Portuguese-language search ("ferro verde", "MACF" — *Mecanismo de Ajuste de Carbono na
  Fronteira*) is the obvious next step.
- Industry: https://www.braziliron.com.br/cbpm-and-brazil-iron-strengthen-green-iron-innovation-with-rwe-in-germany/

### 6.8 Others

- **Energy Transitions Commission**, "Unlocking the First Wave of Breakthrough Steel Investments",
  April 2023 — the source of the HBI mass-loss figure IEEFA cites.
- **Climate Energy Finance**, "Green Metal Statecraft: Policy, Investment and Technology Trends in
  the Green Iron Evolution", **30 April 2026** — Australian, very recent, not opened:
  https://climateenergyfinance.org/wp-content/uploads/2026/04/CEF_Green-Metal-Statecraft_-Policy-Investment-and-Technology-Trends-in-the-Green-Iron-Evolution.pdf
- **SteelWatch**, "Steel decarbonisation in 2025: stagnant but far from static":
  https://steelwatch.org/commentary/steel-decarbonisation-in-2025-stagnant-but-far-from-static/ —
  2025 dominated by tariffs, cost pressure and securitisation; **high energy prices delayed H-DRI
  commercialisation in Western Europe.**
- **WEF**, "Unlocking Asia-Pacific as a First Mover: Australia's Green Iron Opportunity":
  https://www.weforum.org/publications/unlocking-asia-pacific-as-a-first-mover-australia-s-green-iron-opportunity/ ·
  https://www.weforum.org/stories/2025/12/australia-global-transition-green-iron/ (Dec 2025)
- **CSIRO**, "Australia's green metals gambit: the technologies to decarbonise steel using Pilbara
  iron ores", **March 2026**: https://www.csiro.au/en/news/All/Articles/2026/March/Green-metals —
  relevant to the low-grade-ore constraint on Australian cases.
- **ODI**, "EU climate, trade and industrial policy: navigating CBAM implementation and steel
  decarbonisation": https://odi.org/en/publications/eu-climate-trade-and-industrial-policy-navigating-cbam-implementation-and-steel-decarbonisation/
  (403; not retrieved) — likely relevant on developing-country exporter exposure.
- **Wood Mackenzie**, "How will the EU's CBAM impact global iron and steel?":
  https://www.woodmac.com/news/opinion/how-will-the-eus-cbam-impact-global-iron--steel/
- **EUROFER** (incumbent-industry counterpart): "The CBAM must be fixed and launched urgently":
  https://www.eurofer.eu/publications/position-papers/the-cbam-must-be-fixed-and-launched-urgently
- **Norsk Hydro** via Fastmarkets: "System flaws and loopholes threaten CBAM's decarbonization
  objectives": https://www.fastmarkets.com/insights/system-flaws-and-loopholes-threaten-cbams-decarbonization-objectives-hydro-says/
  — but Hydro **welcomed** the indirect-emissions exclusion. The loophole has organised defenders.
- **Aluminium comparator with hard numbers** (trade-compliance source, indicative): a coal-grid
  smelter can carry **12–16 tCO2e/t** total, while CBAM prices only the **1.5–2.1 tCO2e/t** from
  anode consumption and PFCs. **CBAM prices roughly one eighth of the real footprint in the
  coal-powered case** — a direct analogue for the H2-DRI framing.

---

## 7. Additional CBAM detail worth recording

### 7.1 Article 7(2) — the default-value fallback, and its telling phrasing

> "Embedded emissions in goods other than electricity shall be determined based on the **actual
> emissions** in accordance with the methods set out in points 2 and 3 of Annex IV. Where the actual
> emissions cannot be adequately determined, **as well as in the case of indirect emissions**, the
> embedded emissions shall be determined by reference to **default values** in accordance with the
> methods set out in point 4.1 of Annex IV."

Actual emissions are the default pillar; defaults are the fallback. For Annex II goods the
electricity term is removed from the charge base entirely, so an actual-emissions declaration for
CN 7203 produced by electrolysis on a coal grid legitimately lands near zero.

### 7.2 Annex II precise scope

Annex II covers, for iron and steel: **"72 – Iron and steel" except 7202 (ferro-alloys) and 7204
(ferrous waste and scrap)**, plus headings **7301–7311, 7318, 7326**; aluminium **7601–7616**; and
under Chemicals, **2804 10 00 – Hydrogen**. Scrap (7204) is excluded because it is outside CBAM
scope altogether — the separate scrap loophole.

### 7.3 The December 2025 legislative package in full

- **Review report on the transitional period**, 16 Dec 2025 — COM(2025) 783 final:
  https://eur-lex.europa.eu/resource.html?uri=cellar%3A05f0b7f5-da86-11f0-8da2-01aa75ed71a1.0001.02%2FDOC_1&format=PDF
- **COM(2025) 989 final**, 17 Dec 2025 — extension to downstream goods and anti-circumvention.
- **SWD(2025) 988 final** Parts 1–2 — Impact Assessment.

What the proposal **does**: ~180 steel- and aluminium-intensive downstream CN codes from
**1 January 2028** (chargeable only on precursor emissions the ETS would cover); **pre-consumer
steel and aluminium scrap becomes a CBAM precursor**; an empowerment to sub-divide CN codes to
capture composition; an empowerment to demand proof of place of production (**Steel Mill
Certificate**) before actual emissions may be used; clarified electricity rules (**virtual PPAs
systematically excluded**, indirect/sleeved physical PPAs allowed, congestion/nomination conditions
relaxed).

What it **does not** do: indirect emissions are outside the scope of the impact assessment
entirely. "Indirect" appears in Part 1 only in relation to *indirect PPAs*, *indirect price effects*
and *indirect exposure*. No scope-2 policy option, no modelling, no quantification.

Commentary: https://www.akingump.com/en/insights/alerts/eu-carbon-border-adjustment-mechanism-financial-obligations-commence-amid-proposed-scope-expansion-to-include-new-downstream-products ·
https://www.insideenergyandenvironment.com/2026/02/eu-cbam-proposes-expansion-to-complex-metal-products/ ·
https://sustainablefutures.linklaters.com/post/102me4c/eu-commission-proposes-to-extend-cbam-scope-and-adopts-implementing-legislation ·
https://www.pwc.com/mt/en/publications/tax-legal/regulation-eu-2025-2083-simplifying-the-carbon-border-adjustment-mechanism-cbam.html

### 7.4 Source of the Commission's 50%/80% reshuffling assumptions

**Stede, J., Pauliuk, S., Hardadi, G., Neuhoff, K. (2021), "Carbon pricing of basic materials:
Incentives and risks for the value chain and consumers", DIW Discussion Paper No. 1935:**

> "recent academic literature focusing on the EU approximates the scale of potential risk from
> resource shuffling from a CBAM at **around 50% for steel and 80% for aluminium** — the latter
> driven by the higher opportunities to source or attribute the production of aluminium to clean
> electricity."

Note what that says: **aluminium's higher reshuffling risk is precisely because its emissions are
electricity emissions.** Electrolytic ironmaking (MOE, aqueous electrowinning) has aluminium's
emissions profile, not steel's — so the 80% figure is arguably the right analogue for those routes,
not the 50%.

### 7.5 OECD context numbers for scale

Dechezleprêtre et al., OECD STI Working Paper 2025/02: CBAM covers **0.37% of global goods-and-
services trade value**, **3% of EU non-EU imports** (2022); emissions embedded in CBAM imports are
**0.31% of global GHG**; **USD 132 bn** of covered imports of which **USD 78 bn basic metals**.
CBAM is only **8%** of the revenue raised by the modelled ETS package (free-allocation removal 49%,
ETS price rise 43%, at €80/tCO2).

Condensed policy version, OECD March 2025: *"**Within-country resource shuffling**—where
emissions-intensive CBAM goods are consumed domestically while cleaner goods are exported to the
EU—could diminish the policy's effectiveness in reducing carbon leakage and global emissions."*

### 7.6 Sandbag scrap-substitution ceilings (useful if modelling the shuffling counterfactual)

Aluminium ingot is already made at primary quality from **75% post-consumer scrap**; steel long
products **100% scrap** via EAF; steel flat products **up to 80%**. Ferrous scrap traded around
**€300–400/t** (LME, end-July 2025), stainless up to **€1,500/t**.

Commission IA 2025 fn.36: *"Since pre-consumer aluminium and pre-consumer steel scrap are assigned
zero-emissions, imported goods using pre-consumer aluminium and pre-consumer steel scrap as input
material are subject to a lower carbon price compared to goods produced in the EU, thus weakening
the effectiveness of the CBAM."*

### 7.7 Further academic

- **Berahab**, "The European Union's CBAM: averting emissions leakage or promoting the diffusion of
  carbon pricing?", *Journal of Environmental Policy & Planning*, 2025,
  DOI 10.1080/1523908X.2025.2591794. Notes the December 2025 proposal *"facilitates the use of
  voluntary default values … but importers **retain the right to use production-specific values and
  resource shuffling risks therefore remain**."* (Paywalled.)
- **"Early signs that the EU carbon border adjustment mechanism is reshaping EU–India steel trade"**,
  *Nature Climate Change*, 2026: https://www.nature.com/articles/s41558-026-02607-y — firm-level
  empirical evidence of CBAM trade effects. Not retrieved.
- **"Decarbonizing global steel production"**, *Nature Reviews Earth & Environment*, 2026:
  https://www.nature.com/articles/s43017-026-00786-y
- **Gielen et al. (2020)**, "Renewables-based decarbonization and relocation of iron and steel
  making: A case study", *Journal of Industrial Ecology*:
  https://onlinelibrary.wiley.com/doi/abs/10.1111/jiec.12997 — the original relocation case study.

### 7.8 Retrieval gaps to close manually

- **EP written question E-001119/2025**, "Addressing resource shuffling and ensuring the CBAM's
  effectiveness in the European steel industry" — the EP site blocks automated fetch. The only
  EP-level document naming resource shuffling and steel together.
- **DG TAXUD "Default values definitive period" Excel** — for the per-CN-code 7203 row under
  IR (EU) 2025/2621. TAXUD rate-limits automated access.
- **LeadIT/SEI** EU hydrogen-and-electricity demand report (Cloudflare).
- **ODI** CBAM/steel paper (403).
- Agora deep dives for Japan, China, South Korea, South Africa.
- **Nykvist et al. (2025)** *Applied Energy* full text — CC-BY but ScienceDirect blocks fetch; the
  numbers behind "weaker than previously found" are worth having.

## 8. Australian and Brazilian policy context — does either price an export-oriented green iron plant?

**Short answer: no, in both countries — and in Australia the plant is actually *paid*.**

Access note: `dcceew.gov.au` and `industry.gov.au` were largely unreachable to automated fetch;
some points below rest on secondary sources and are flagged.

### 8.1 Australia — the Safeguard Mechanism covers scope 1 only, and credits a new green iron plant

**Confirmed.** The Clean Energy Regulator states covered emissions are "a facility's **scope 1
emissions**", with carve-outs for legacy landfill waste, Greater Sunrise, and grid-connected
generators in sectoral-baseline years (https://cer.gov.au/schemes/safeguard-mechanism, updated
4 Sep 2026). The statutory hook is **s 7 of the Safeguard Rule**
(https://www.legislation.gov.au/F2015L01637), which carves covered emissions out of the universe of
*scope 1* emissions. **Scope 2 and scope 3 are entirely outside the mechanism.**
(Corroborated by ICAP: https://icapcarbonaction.com/en/ets/australia-safeguard-mechanism)

**Electricity is handled by a dead letter.** Grid-connected generators sit under a **collective
sectoral baseline of 198 MtCO2e/yr** (NEM, SWIS, NWIS, Darwin–Katherine, Mt Isa–Cloncurry), set on
2009-10 to 2013-14 emissions. Individual generators are covered only if the sector breaches it.
**It has never been breached.** So Australian electricity carries no Safeguard carbon cost, and an
electrolyser's or EAF's power consumption carries none either.

**Coverage threshold:** >100,000 tCO2e/yr covered (scope 1) emissions — Safeguard Rule s 8.

**Decline rate — 4.9%/yr compounding to 30 June 2030** (s 31), expressed as an emissions reduction
contribution multiplier:

| FY | ERC |
|---|---|
| 2023–24 | 0.951 |
| 2025–26 | 0.853 |
| 2026–27 | 0.804 |
| 2029–30 | 0.657 |
| 2030–31 on | 0.624 (**3.285%/yr — indicative only**) |

The post-2030 trajectory is **not legislated**; DCCEEW sets it in five-year blocks. Genuine open
risk in either direction beyond FY2030.

**New facilities — the international best practice benchmark.** Safeguard Rule **s 29(2)**: a
facility is *new* if it has no historical or transitional production variables — operationally,
all commercial production began **on or after 1 July 2023**. Its baseline is

> Baseline = ERC × Σ (best-practice emissions intensity × production quantity)

with fallback to the industry-average default only where Schedule 1 specifies no best-practice
number. **Mandatory, not elective**, and it applies equally to a new production variable added at
an existing site.

**Schedule 1 s 39 defines "primary iron"** as "physical and chemical processing of iron feed
material into a **crude iron product suitable for export from the facility**", with the statutory
example: "**Pig iron, hot briquetted iron, direct reduced iron and cast iron** are each a crude iron
product…" and "The production of crude iron products from iron ore pellets using **direct
reduction**." So HBI is named in the rule.

| Production variable | Default (tCO2e/t) | **Best practice (tCO2e/t)** |
|---|---|---|
| **Primary iron (s 39)** | 2.08 | **1.77** |
| Iron ore pellets (s 40) | 0.0526 | 0.0501 |
| Primary steel (s 41) | 2.07 | none specified |

*Sourced from Safeguard Rule compilation No. 13 (in force 31 Aug 2024). A later Production
Variables Update Rule could have revised 1.77 — worth one check at legislation.gov.au/F2015L01637.*

**What this means for a 2 Mt/yr HBI plant:**
- FY2026–27 baseline = 0.804 × 1.77 = **1.423 tCO2e/t HBI** → ~2.85 MtCO2e
- FY2029–30 baseline = 0.657 × 1.77 = **1.163 tCO2e/t HBI**

Against actual scope 1:

| Route | Scope 1 (tCO2/t HBI) | Outcome |
|---|---|---|
| **H2-DRI** | ~0.03–0.10 (flux calcination, carburisation carbon, residual firing) | 60–200 ktCO2e/yr. If >100 kt → covered and **generates ~2.6–2.8 m SMCs/yr**. If <100 kt → outside the scheme, no cost and **no credit** |
| **MOE** | ~0 with inert anodes | Almost certainly **below 100 kt → outside the scheme entirely** |
| **Aqueous electrowinning** | ~0 | Same |
| **NG-DRI** | ~0.50–0.70 | Covered, still ~0.7–0.9 t/t **below** baseline → ~1.4–1.8 m SMCs/yr |

**The benchmark was calibrated on blast-furnace-era performance, so every DRI-route new entrant —
hydrogen or gas — clears it by a wide margin and is credited rather than charged.** A gas or
blended-gas transitional reductant is not penalised in Australia this decade; on a continued
3.285%/yr decline an NG-DRI plant's baseline only crosses ~0.60 t/t around FY2041.

**Modelling cautions.** SMC value is **not** the cost-containment price — SMCs trade against the
ACCU market, so treat SMC revenue as an upside option bounded above by the ceiling and plausibly
well below it. Note also **s 10(1)**: baselines below 100,000 t are rounded up to 100,000, but that
floor **does not apply** when calculating SMC entitlement. And note the perverse cliff: a plant
*below* 100 kt gets nothing at all, so MOE and aqueous electrowinning likely fall out of the scheme
entirely and forgo the credit.

**Cost containment measure (the price ceiling)**
(https://cer.gov.au/schemes/safeguard-mechanism/managing-excess-emissions/cost-containment-measure,
updated 3 Aug 2026): 2023–24 **A$75.00/ACCU**; **2026–27 A$87.72/ACCU**; indexed **CPI + 2%**.
Available only to facilities in excess who confirm they cannot buy ACCUs cheaper; must be
surrendered immediately, not tradeable or bankable. FY2026–27 window 16 Feb – 4 Mar 2027.

**TEBA (trade-exposed baseline-adjusted)** exists and green iron passes the product test —
**Schedule 2 item 15, "Tonnes of primary iron"**, a *manufacturing* production variable. Cost test:
compliance cost > **3% of EBIT** with full relief at **10%**; decline-rate floor **1%** for
manufacturing (2% otherwise); three-year relief periods. **Irrelevant to a green iron plant** — TEBA
slows baseline decline, and a plant 90%+ below its baseline has no cost impact to measure. It
matters to incumbent Australian steelmakers as comparators.

**2026–27 review:** the CER states (4 Sep 2026) that "the department will review the Safeguard
Mechanism's policy settings in 2026-27." It is the designated vehicle for post-2030 decline rates
*and* the CBAM decision. Commencement status, ToR and consultation paper **not verified** — the
biggest open item on the Australian side.

**Bottom line: an export-oriented Australian green iron plant faces no effective carbon price.
Its carbon "price" is a negative number.**

### 8.2 Australia — Green Iron Investment Fund: A$1 bn, capex grants, no awards confirmed

**Announced 20 February 2025** by the PM and Minister for Industry
(https://www.pm.gov.au/media/albanese-government-building-australias-green-iron-future). Split:
**A$500 m National Development Stream** (open nationally) and **up to A$500 m earmarked for
Whyalla Steelworks** under new ownership.

**Form of support: capital grants, not a production credit.** From business.gov.au
(https://business.gov.au/grants-and-programs/green-iron-investment-fund-national-development-stream):

| Parameter | Value |
|---|---|
| Grant share | **up to 25% of eligible project expenditure** |
| Applicant co-contribution | ≥75% of total eligible expenditure |
| Cap across *all* government sources | 65% |
| **Minimum production capacity** | **1,000,000 t/yr** |
| **Minimum TRL at application** | **TRL 7** |
| Must be producing and selling green iron by | **31 March 2031** |
| FID required within | 18 months of grant agreement execution |

Eligible = "reduction of iron ore into a concentrated iron metal" using renewable hydrogen,
renewable energy, **or natural gas with a documented pathway to renewables** (assessed by a
Technical Assessment Panel). Ineligible: input costs (ore, pellets, reducing agents), mining,
fossil extraction infrastructure. Stackable with NRF debt/equity, state grants and Export Finance
Australia, subject to the 65% cap.

*Date discrepancy: business.gov.au says applications closed **17 February 2026**; DISR's news page
and secondary guides say **28 October 2025 → 16 January 2026**. Probably a one-month extension.*

**Awards: none confirmed.** Assessment scoped at ~9 weeks with announcement ~12 weeks (i.e. ~May
2026). As at **12 May 2026** the Superpower Institute stated that Future Made in Australia funding
"**has not flowed to projects on the ground**"
(https://www.superpowerinstitute.com.au/news/federal-budget-response-budget-shows-progress-on-fairness-but-gaps-remain).

**Whyalla** entered administration February 2025 after the SA government forced it there over
GFG/Sanjeev Gupta's unpaid royalties (A$18.5 m). A **A$2.4 bn rescue package** followed (A$1.9 bn
for infrastructure under a new owner). Sale run by KordaMentha/333 Capital; 70+ parties, five
shortlisted, binding bids pending in early 2026. The A$500 m is contingent on a buyer.

**NeoSmelt** (BlueScope + BHP + Rio Tinto + Woodside + Mitsui) — WA pilot, **30,000–40,000 t/yr
molten iron, operations from 2029**. A *pilot* an order of magnitude below the GIIF 1 Mt/yr floor,
so not a National Development Stream candidate.

### 8.3 Australia — Hydrogen Headstart, and the production tax credits that exclude green iron

**Hydrogen Headstart** pays a **Headstart Production Credit** **per unit of production over a
10-year operating period**, once the plant is running. **No published A$/kg strike price** — the
rate is contracted per project. Do not model it as a headline number.

*Round 1* (A$2 bn envelope):
- **Murchison Green Hydrogen** (Copenhagen Infrastructure Partners), WA — **A$814 m, 20 March
  2025**. 1,500 MW electrolysis, 3,600 t/day ammonia, ~1.2 GW solar + 1.7 GW wind.
  https://arena.gov.au/news/murchison-green-hydrogen-project-given-a-headstart/
- **Hunter Valley Hydrogen Hub** (Orica), NSW — **A$432 m, July 2025**.
  https://arena.gov.au/funding/hydrogen-headstart/
- Total committed ≈ A$1.25 bn of A$2 bn.

*Round 2* — announced at **A$2 bn**, **cut to ~A$1 bn in the 12 May 2026 Budget**. Garnaut: "The
halving of Hydrogen Headstart appropriations is disappointing." **Shortlist announced 13 May 2026 —
seven projects, ~2.18 GW, full applications due early September 2026**
(https://arena.gov.au/news/arena-shortlists-major-projects-to-scale-australias-renewable-hydrogen-industry/):

| Project | State | MW | Offtake |
|---|---|---|---|
| Bell Bay Powerfuels | TAS | 300 | Methanol |
| European Energy — SE Qld Power-to-X | QLD | 150 | Methanol |
| HAMR Energy — Portland Renewable Fuels | VIC | 220 | Methanol + SAF |
| HIF Tasmania e-Fuel Facility | TAS | 140 | Methanol |
| Murchison Stage 1B | WA | 500 | Ammonia |
| Perdaman Helios (Karratha) | WA | 750 | Urea |
| Summit Hydro — Gladstone | QLD | 120 | Alumina |

**Not one green iron project on the shortlist.** Every awardee and shortlistee to date is ammonia,
methanol, SAF, urea or alumina.

**Production tax credits — hydrogen yes, green iron never legislated.** The *Future Made in
Australia (Production Tax Credit and Other Measures) Bill 2024* passed the Senate **11 February
2025**; combined value **A$13.7 bn over 10 years**.

**Hydrogen Production Tax Incentive**
(https://cer.gov.au/schemes/guarantee-origin-scheme/product-guarantee-origin/hydrogen-production-tax-incentive,
updated 25 Aug 2026):
- **A$2.00 per kg** of eligible low-emissions hydrogen
- Emissions threshold **≤ 0.6 kgCO2e per kg H2**
- Income years **1 July 2027 – 30 June 2040**, max **10 years per project**
- **Refundable** tax offset
- Requires registration under the **Product Guarantee of Origin** scheme
- 2026 additions: a **Grid Matching Requirements Instrument** (how grid-connected facilities must
  source renewables); Community Benefit Principles still under development

**Critical Minerals Production Tax Incentive:** 10% of eligible processing/refining costs, up to
10 years, same window. **Iron is not on Australia's Critical Minerals List or Strategic Materials
List** — it is a bulk commodity, so green iron does not qualify. *(High confidence but not verified
against the DISR list page.)*

**A green iron production credit has never been announced, introduced or legislated.** It exists
only as a think-tank recommendation. The 2026–27 Budget contained none. Superpower Institute: *"For
green iron specifically, the missing piece is a production credit… Without it, Australian green
iron projects will continue to struggle to reach final investment decision."*

**Guarantee of Origin — green metals announced, not yet live.** The Product GO currently covers
**only electrolytic hydrogen**. CER (24 Aug 2026): *"In 2026, it will expand to include… **green
metals**"* — listed under "Upcoming products", not in force
(https://cer.gov.au/schemes/guarantee-origin-scheme/product-guarantee-origin). Framework: FMIA
(Guarantee of Origin) Act 2024 (https://www.legislation.gov.au/C2024A00112) + Rules 2025 +
Methodology Determination 2025. **The GO methodology covers scope 1 + scope 2 + scope 3** — a much
wider boundary than the Safeguard's scope-1-only liability, and the instrument that would
substantiate an emissions claim to an EU CBAM declarant.

### 8.4 Australia — the Carbon Leakage Review: a CBAM for cement first, and no export rebates

Led by **Prof. Frank Jotzo** (ANU). Commenced 1 July 2023; two consultations (Nov–Dec 2023,
Nov–Dec 2024). **Final report dated February 2025 but not publicly released until February 2026** —
a ~12-month delay.
- https://www.dcceew.gov.au/sites/default/files/documents/carbon-leakage-review-final-report.pdf
- https://www.dcceew.gov.au/climate-change/emissions-reduction/review-carbon-leakage
- https://www.pwc.com.au/tax/tax-alerts/final-report-from-the-australian-carbon-leakage-review.html

Recommendations (p.19):
1. (a) BCA for **cement and clinker**. (b) BCA **considered for lime; hydrogen, ammonia and
   derivatives; steel and iron; glass — subject to further assessment**. **(c) A BCA providing
   rebates for exports should not be pursued** (trade-law and target-integrity grounds).
2. (b) Liability assessed on **scope 1 emissions above benchmarks set in line with Safeguard
   baselines**, netting the **explicit** carbon price paid in the origin country. (c) TEBA removed
   for a commodity once a BCA is fully implemented for it.

Overall: *"existing policy measures mitigate leakage risks in the short to medium term."*

On HBI specifically: *"emerging low emissions iron production could be at risk of leakage from
imports of conventional iron, such as **hot briquetted iron which is globally traded**… Industry
feedback highlighted the importance of any border carbon adjustment for steel including
**precursors such as iron**."* Steel and iron were deferred because *"complex and varied production
structures… could not proceed as quickly as for other commodities such as clinker and cement."*

**Status: no Australian CBAM. No bill. No government decision.** DCCEEW: "The Review's
recommendations will be considered in the 2026-27 review of the Safeguard Mechanism."

**Critical: if Australia ever builds a CBAM, an Australian green iron exporter gets no border
rebate** — recommendation 1(c) explicitly rules it out. An Australian BCA would be a cost on
imports *into* Australia and would do nothing for HBI leaving for the EU.

### 8.5 Australia — Garnaut / Superpower Institute proposals and export targets

*A Green Iron Plan for Australia* (Superpower Institute + Bivios), **May 2025** — see §6.4 for the
full cost modelling. The two recommendations most relevant here:

**Recommendation 1:** a **stackable green iron production tax credit worth at least $170 per tonne
of green iron in 2030, inclusive of the HPTI**, rising to maintain equivalence with the EU carbon
price.

**Recommendation 2 — the one that bears on MOE and electrowinning:** *"Some nascent green iron
production technologies do not use hydrogen, but may use significant amounts of renewable energy
dedicated to iron-making. Here, the HPTI does not help close the cost gap… The government should
provide support that simulates the effect of a carbon price for **non-hydrogen-based green iron
technologies**… worth at least $170 per tonne."* **MOE and aqueous electrowinning get nothing from
the A$2/kg HPTI, and TSI names this as an explicit policy gap.**

Others: capital support at **15% of capex for up to three projects, +15% for first-of-a-kind,
capped at A$500 m/project**, for plants ≥0.5 Mt/yr (Rec 5 — note this *exceeds* the GIIF's actual
25% / A$500 m-total design); a **green hydrogen certificate scheme** letting NG-based iron count as
green by surrendering certificates (Rec 4); and **shaping the GO scheme to be CBAM-compatible**
(Rec 6).

**Export figures in circulation:**

| Figure | Source |
|---|---|
| Green iron could generate **up to A$386 bn annually by 2060** (vs ~A$120 bn/yr iron ore today) | TSI, May 2025 |
| Green iron worth **up to A$400 bn** in export revenue | TSI, 12 May 2026 |
| Australian green iron could abate **~4% of global emissions** — >3× Australia's domestic total | TSI, May 2025 |
| Global green iron demand **~852 Mt by 2050** (net zero scenario) | DISR, 7 Mar 2025 |
| Green iron + steel could add **up to A$96 bn/yr to exports by 2040** | Accenture (2023), cited by DISR |
| One **4.8 Mt/yr** plant → A$85 bn GDP, A$2.4 bn real income/yr, 1,540 FTE | MRIWA, cited by DISR |
| Australia produced **5.4 Mt of steel in 2023** | Carbon Leakage Review |

**Carbon pricing proposal — renamed and repriced.** The Carbon Solutions Levy framing is superseded
by **"The Case for Pricing Pollution"** (Finighan & Burfurd, **29 January 2026**),
https://www.superpowerinstitute.com.au/work/the-case-for-pricing-pollution:
- **Polluter Pays Levy** — upstream on carbon embedded in fossil fuels extracted or imported;
  ~140 sites / <60 companies; covering **>80% of Australian emissions** (vs 30% under the
  Safeguard). **Starts 2026 at A$17/t, rising to meet the EU carbon price in 2034** and tracking it
  thereafter. Accompanied by an **EU-CBAM-style border levy**. Revenue **A$22.6 bn/yr average
  2026–2050**.
- **Fair Share Levy** — resource-rent reform on gas. Together **A$35.6 bn/yr average to 2050**.
- **Status: proposal. Never legislated, never adopted.** TSI's 12 May 2026 budget response notes
  "the absence of progress." *(The 2024 Carbon Solutions Levy figures, ~A$90/t and ~A$100 bn/yr,
  could not be re-verified — treat as superseded.)*

TSI's critique of the Safeguard mirrors the §8.1 arithmetic: *"Emissions intensity baselines are
centrally determined… It is not possible for agencies to specify these accurately and fairly for
all industries, resulting in **wealth transfers between sectors and distorting investment
decisions**."*

**Targets and what was cut.** Legislated: 43% below 2005 by 2030, net zero 2050 (Climate Change Act
2022, https://www.legislation.gov.au/C2022A00037). **2035 NDC: 62–70% below 2005**, submitted
September 2025 — a range, and not in the Act. **Net Zero Plan** released September 2025 with six
sector plans including Industry and Resources *(contents unverified)*. **Cut:** Hydrogen Headstart
Round 2 halved (12 May 2026 Budget); Energy Bill Relief Fund not extended beyond 2025.
**New in the 2026–27 Budget:** A$1.1 bn Cleaner Fuels Program; two-year loss carry-back; start-up
loss cash-out. **No carbon pricing measure, no green iron production credit.** Other funding:
Future Made in Australia **A$22.7 bn** (green metals in the Net Zero Transformation Stream);
**FMIA Innovation Fund A$750 m for green metals**; National Reconstruction Fund **up to A$3 bn of
A$15 bn** for renewables/low-emissions.

### 8.6 Brazil — SBCE is legislated, but nothing prices before ~2030–2031

**Lei nº 15.042, signed 11 December 2024**, published 12 December, in force on publication
(https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2024/lei/l15042.htm). Creates **CBE**
(allowance, 1 tCO2e) and **CRVE** (offset credit).

**The two thresholds — Art. 30, and the reading is the reverse of the intuitive one:**
- **>10,000 tCO2e/yr → monitoring plan + emissions reporting only (MRV). No surrender.**
- **>25,000 tCO2e/yr → MRV *plus* obligation III, the periodic reconciliation, i.e. actual
  allowance surrender.**

Art. 30 §1 lets the *órgão gestor* **raise** either threshold. Art. 30 §2 limits obligations to
activities with consolidated MRV methodologies.

**Scope 2 is NOT covered.** The law regulates installations and sources that emit (Art. 1 §1) — a
point-source, direct-emissions design. "Emissões indiretas" appears once in the whole law
(Art. 1 §3), and only to *exclude* agricultural input emissions. IETA's August 2026 business brief
is explicit: *"Scope 2 emissions may be **registered but not capped** to support competitive
benchmarking."*
(https://www.ieta.org/uploads/wp-content/2026/07/18August-IETABR_BusinessBrief_Brazil-SBCE_2026_v3.pdf)
**Agriculture explicitly excluded** (Art. 1 §2).

**Phase timeline (Art. 50 gives durations, not dates):**

| Phase | Content | Duration | IETA calendarisation |
|---|---|---|---|
| I | Regulation issued | 12 mo, extendable 12 (→ 11 Dec 2026) | 2025–2026 |
| II | Operators stand up MRV | 1 yr | 2027 |
| III | Monitoring plan + reporting only | 2 yrs | 2028–2029 |
| IV | First National Allocation Plan, **free** CBE distribution, market opens | PNA ≥12 mo | **first NAP ~2030** |
| V | Full implementation, **first CBE auction** | — | **~2031** |

IETA: *"The first binding compliance obligations are expected from **2030 onwards**, following
publication of the first NAP."* Finance Ministry / World Bank roadmap:
https://www.gov.br/fazenda/pt-br/central-de-conteudo/publicacoes/guias-e-manuais/2024/241209-crtlh-implementacao-sbce-v4.pdf
ICAP factsheet: https://icapcarbonaction.com/en/ets/brazilian-greenhouse-gas-emissions-trading-system

**Regulator: no permanent one exists yet.** Decreto 12.677 of 15 October 2025
(https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2025/decreto/d12677.htm) created the
**Secretaria Extraordinária do Mercado de Carbono (SEMC)** in the Ministry of Finance, exercising
the managing-body functions **on an interim basis** with a mandate to **end-2026**; transition to a
permanent body expected 2027.

**Iron and steel ARE covered — and integrated steel is in the first MRV cohort.** The July 2026
draft *Portaria* phases 17 sectors into the PB-MRV programme:
- **Cohort 1, from 2027:** **integrated steel (siderurgia integrada)**, cement, primary
  aluminium/alumina, oil & gas E&P, refining, pulp & paper, air transport
- **Cohort 2, from 2029:** **semi-integrated steel**, chemicals, power, mining, glass, ceramics,
  waste, recycled aluminium, food & beverage
- **Cohort 3, from 2031:** road, rail, waterway transport

Each cohort runs 4 years: Y1 monitoring plans, Y2–3 data collection, Y4 assessment against
remaining criteria (including **international trade exposure**) and NAP development. **Cohort 1
data therefore feeds a NAP no earlier than ~2030.** The consultation text states that inclusion in
MRV *"does not by itself imply immediate submission to emission limits or reconciliation
obligations."*

**Definitional gap to resolve:** the cohorts split "integrated" vs "semi-integrated" steel.
**A merchant HBI/DRI plant that does not melt steel fits neither label cleanly.** Nothing published
resolves whether it lands in cohort 1, cohort 2, or under "mining".

**2025–26 implementation trail:**
- **24 March 2026** — CTCP (permanent technical advisory committee) installed:
  https://www.gov.br/fazenda/pt-br/assuntos/noticias/2026/marco/ministerio-da-fazenda-instala-comite-da-sociedade-civil-e-governos-para-implementacao-do-mercado-regulado-de-carbono-no-brasil
- **May 2026** — CTCP Resolutions 1–4/2026 create three working groups (Financeiro, MRV,
  Metodologias):
  https://www.tauilchequer.com.br/pt/insights/publications/2026/05/progress-in-the-implementation-of-the-brazilian-emissions-trading-system-ctcp-resolutions-and-sectoral-coverage-proposal
- **19 May 2026** — SEMC presents the sectoral-coverage methodology (MRV feasibility 50% / market
  structure 30% / emissions intensity 20%):
  https://www.gov.br/fazenda/pt-br/assuntos/noticias/2026/maio/ministerio-da-fazenda-apresenta-proposta-preliminar-de-cobertura-setorial-do-mercado-regulado-de-carbono
- **28 July – 28 August 2026** — two public consultations (sectoral coverage; MRV timetable draft
  Portaria). Final Portaria expected "still in 2026":
  https://www.tauilchequer.com.br/pt/insights/publications/2026/08/extraordinary-secretariat-for-the-carbon-market-submits-to-public-consultation-a-draft-ordinance-on-the-timeline-for-the-implementation-of-the-mrv-obligations-under-the-sbce
- Government target: full infralegal framework **by December 2026**:
  https://agenciabrasil.ebc.com.br/economia/noticia/2025-11/governo-quer-regulamentar-mercado-de-carbono-ate-fim-de-2026
- **February 2026** — Brazil joined ICAP as 36th member
- **29 May 2026** — **the STF struck down Art. 56 (ADI 7795)**, which had forced insurers and
  pension funds to put ≥0.5% of technical reserves into carbon assets. An announced demand channel,
  killed:
  https://noticias.stf.jus.br/postsnoticias/stf-invalida-regra-que-obrigava-seguradoras-a-aplicar-recursos-em-creditos-de-carbono/
- **No CRVE offset methodologies have been accredited.** The methodologies working group only
  started in May 2026.

### 8.7 Brazil — would an export-oriented green iron plant face a carbon price? No, on three grounds

1. **Timing.** No surrender obligation before the first NAP (~2030); first auction ~2031. Phases 1–3
   (2025–2029) carry no compliance obligations.
2. **Threshold.** Surrender attaches only above **25,000 tCO2e/yr scope 1**:
   - **H2-DRI:** residual scope 1 is carburisation carbon (~15–25 kg C/t to reach 1.5–2.5% C in HBI)
     plus back-up firing. At 2 Mt/yr that is ~100–200 ktCO2/yr — **over the threshold**.
     Carburising with bio-methane or charcoal makes it biogenic and could plausibly drop the plant
     **under 25,000 t**.
   - **MOE:** carbon-free anodes → essentially zero process CO2. **Likely below both thresholds.**
   - **Aqueous electrowinning:** same; makes iron plate/powder, so briquetting adds only handling
     emissions.
   Worth modelling explicitly — under 25,000 t you never surrender; under 10,000 t you never even
   report.
3. **Scope 2 is out.** Grid electricity — the dominant emission source for all three routes — is not
   capped.

**No operating sub-national ETS or state carbon tax** found *(not an exhaustive 26-state survey —
read as "no evidence of one")*. The voluntary market is a revenue opportunity, not a cost, but no
CRVE methodologies are accredited.

**Article 6:** an MMA draft CIM Resolution (consultation to 6 August 2026) sets **100 MtCO2e for
2031–2035, of which up to 50 MtCO2e authorisable as ITMOs**, with a 50% per-project ceiling. MoUs
with Switzerland and Singapore. **Selling ITMOs and claiming low-carbon HBI for CBAM would be
double counting.**

**Net: zero domestic carbon cost through at least 2030–2031 — but also zero domestic carbon
revenue.** Unlike Australia, where the Safeguard *pays* in SMCs, Brazil gives a green iron plant no
way to monetise its low emissions domestically. The value has to come from the EU side.

### 8.8 Brazil — green hydrogen policy: steel is a named priority, but the money moved two years right

**Lei nº 14.948, 2 August 2024** — the low-carbon hydrogen legal framework
(https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2024/lei/l14948.htm).
- **"Low-carbon hydrogen" = lifecycle GHG ≤ 7 kgCO2e/kgH2**, frozen until 31 Dec 2030
  (Art. 3 XII, §1). **Far looser than the EU RFNBO threshold of 3.38 kgCO2e/kgH2 — a Brazilian
  certificate does not automatically satisfy EU rules.**
- Creates the **SBCH2** certification system.
- **Decreto nº 13.096, 12 August 2026**
  (https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2026/decreto/d13096.htm) regulates the whole
  framework — PNH2 governance, SBCH2, REHIDRO and the PHBC auction procedure. Promised since early
  2026 and repeatedly slipped; this one is real.

**REHIDRO** (Arts. 26–29): **suspension of PIS/PASEP and COFINS** on capital goods, construction
materials and services incorporated into the project (REIDI-style, Law 11.488/2007). **A capex tax
relief, not an operating subsidy.** Runs **5 years from 1 January 2025** (→ 31 Dec 2029). Carries
**minimum local content** and **minimum R&D and domestic energy-transition spend** requirements
(Art. 26 §2; Art. 27 §6). **Art. 27 §5: firms in Export Processing Zones (ZPEs) may use REHIDRO
without losing ZPE benefits** — directly relevant, since Pecém is a ZPE.

**PHBC — R$18.3 billion, and the dates moved.** The credits were vetoed out of Law 14.948 and
re-enacted in **Lei nº 14.990, 27 September 2024**
(https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2024/lei/l14990.htm), then **shifted two years
later by Lei nº 15.269, 24 November 2025**
(https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2025/lei/l15269.htm):

| Original (2024 law) | Amended (Lei 15.269/2025) | Amount |
|---|---|---|
| 2028 | **2030** | R$ 1.7 bn |
| 2029 | **2031** | R$ 2.9 bn |
| 2030 | **2032** | R$ 4.2 bn |
| 2031 | **2033** | R$ 4.5 bn |
| 2032 | **2034** | R$ 5.0 bn |
| | **Total** | **R$ 18.3 bn** (~USD 3.3–3.4 bn) |

**If a model assumed PHBC credits from 2028, that is now wrong — nothing flows before 2030.**
Decreto 13.096 Art. 52 sole paragraph confirms credits may only be granted **2030–2034**.

Mechanism: a **CSLL credit**, compensable or **refundable in cash** within 12 months. Value = up to
**100% of the gap between the estimated price of low-carbon hydrogen and its substitute**, with the
percentage possibly **inversely proportional to GHG intensity**. Awarded by **competitive auction**
run by the Finance Ministry, minimum criterion **lowest credit per unit of product**. Winners post
guarantees up to 10% of total credit; 20% penalty for non-implementation. Unused annual amounts
**roll forward**.

**Green steel explicitly qualifies — Decreto 13.096 Art. 53:**
- §1(II): eligible bidders include *"companies or consortia that **use hydrogen for their own
  consumption** in an energy or non-energy industrial process."*
- §2: own-consumption industrial projects are ***"considerados prioritários"*** in fertilisers;
  **siderúrgico (steel)**; cement; chemical; petrochemical; heavy transport.
- §3: credits awarded **preferentially** to own-consumption industrial projects.
- Art. 4 §11(II) also prioritises lowest GHG intensity and **greater deepening of the national value
  chain**.

Link to REHIDRO (Law 14.990 Art. 4 §9): producers must be REHIDRO beneficiaries; **buyers must
purchase H2 produced by a REHIDRO beneficiary.**

**Domestic-use requirement:** the hydrogen must be **produced in Brazil**, but there is **no
explicit export prohibition**. The incentive structure strongly favours domestic own-consumption.
**Practically: an H2-DRI plant consuming its own hydrogen in Brazil and exporting HBI is the
best-favoured configuration**; exporting the hydrogen as ammonia is the worst.

**MOE and aqueous electrowinning use no hydrogen and are therefore ineligible for both REHIDRO and
PHBC.** A material asymmetry between the three routes — and it mirrors the Australian one (no HPTI
either). **In both countries the non-hydrogen routes are policy orphans.**

**Auction status: announced, regulated, none held.** **Zero PHBC auctions** to date. On **17 June
2026** the Finance Ministry indicated the first auction **will happen in 2027**, with consultation
on the *edital* by **30 November 2026** and edital publication **January 2027**
(https://eixos.com.br/newsletters/dialogos-da-transicao/primeiro-leilao-de-hidrogenio-fica-para-2027-indica-fazenda/ ·
https://portalradarenergia.com.br/noticias/leilao-de-hidrogenio-de-baixo-carbono-fica-para-2027-indica-fazenda-d57094e1).
It had been targeted for 2026. **Sequencing risk: edital early 2027, first credits payable 2030,
and the auction falls in the next presidential term.**

**Other Brazilian instruments:**
- **Fundo Clima** via BNDES — **2026 reimbursable budget R$ 27.5 bn, a record**, approved ~March
  2026. Lines include *Indústria Verde* and *Transição Energética*. *(Current interest rate not
  confirmed.)* https://www.bndes.gov.br/wps/portal/site/home/financiamento/produto/fundo-clima
- **BNDES Climate Call** — R$5 bn call; **up to R$4.3 bn selected January 2026**.
- **PATEN** (Lei nº 15.103, 22 January 2025,
  https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2025/lei/l15103.htm) — lets companies use
  federal tax-debt credits as collateral for energy-transition financing. *Unproven throughput.*
- **Nova Indústria Brasil Mission 5** — R$468 bn headline, plus R$140 bn via BNDES/Finep to Dec
  2026. **These are credit envelopes, not grants, and overlap heavily with Fundo Clima. Do not add
  them to PHBC.**
- *Mover* is automotive R&D — not relevant to iron and steel.

### 8.9 Brazil — Vale's briquette and the project pipeline: nothing at FID

**Critical distinction:** Vale's ***briquete verde*** is an **agglomerated iron ore feedstock** — a
low-temperature binder-bound pellet/sinter substitute claiming up to 10% CO2 reduction for
customers. **It is NOT metallised, and it is not HBI.** Do not conflate it. World's first briquette
plant inaugurated at **Tubarão, Vitória (ES), December 2023**
(https://saladeimprensa.vale.com/w/revolution-in-the-global-steel-industry-vale-inaugurates-the-worlds-first-briquette-plant-in-vitoria-brazil).

- **30 April 2026:** Vale announced **R$12 bn in Espírito Santo through 2030**, including briquette
  expansion and green hydrogen integration. No capacity target or dedicated briquette budget
  disclosed.
- **16 June 2026:** Vale signalled up to **R$13 bn for decarbonisation**.

**Mega Hubs (Middle East) — announced, then slowed:**
- **Saudi Arabia:** preliminary agreement January 2025, potential up to **12 Mt/yr DRI**
- **Oman (Duqm SEZ):** **USD 5 bn**, pellets + HBI, **FID expected 2026**, completion 2029
  (reported 14 Jan 2026)
- **4 August 2026 — CEO Gustavo Pimenta said the Middle East mega hubs have been decelerated
  because of the Iran war.** The first industrial centre had been projected for 2027; now "still
  being evaluated and discussed, but they have slowed down."
  https://www.infomoney.com.br/mercados/vale-desacelera-projetos-de-mega-hubs-no-oriente-medio-diante-da-guerra-no-ira/
- **Status: MOUs and preliminary agreements. No FID on any Mega Hub.**

**Brazil Mega Hub:** Vale + **Porto do Açu** (São João da Barra, RJ) — **MOU signed 26 September
2023** to *study* an HBI hub, starting on natural gas with future conversion to green hydrogen. One
report cited USD 1 bn. Explicitly an MOU.

**Other producers:**
- **ArcelorMittal Tubarão:** 28 August 2026 — up to **R$5 bn for a new rolling mill**. **Downstream
  rolling capacity, not ironmaking decarbonisation.** Green-H2 activity remains a 2023 EDP
  study/pilot MOU.
- **Aço Verde do Brasil (Marabá, PA):** operating charcoal-based integrated steel; R$100 m
  industrial-gases plant Jan 2026. **AVB has publicly rejected the hydrogen route**, favouring
  biomass/charcoal.
- **Ternium Brasil:** no green DRI/HBI announcement found. Do not assume one.

**Port hubs:**
- **Pecém (Ceará), a ZPE** — the flagship. Casa dos Ventos/TotalEnergies ~USD 8.4 bn, FID repeatedly
  slipped; as of 30 April 2026 the USD 5 bn commitment is "to be formalised by 2027". Fortescue
  ~R$20 bn "advanced". **Red flag: Casa dos Ventos' R$100 bn Ceará pipeline reportedly lacks grid
  connection capacity**
  (https://neofeed.com.br/negocios/projetos-de-r-100-bilhoes-da-casa-dos-ventos-em-hub-de-hidrogenio-verde-no-ceara-ficam-sem-energia/).
  Note this is the same corridor as the Caiafa et al. (2026) One Earth study in §2.8.
- **Suape (Pernambuco):** materially less advanced — 2025–26 news is an R&D hub, not industrial
  scale.

**No FID has been taken on any green iron / H2-DRI / HBI-for-export project in Brazil.** Everything
is MOU, feasibility study, or a capital-allocation headline with no disclosed project boundary.
**A model that uses the announced pipeline as a supply curve will overstate near-term deployment
substantially.**

### 8.10 Brazil's grid emission factor — two official numbers that differ by 4–5×

**Do not mix them.**

**Ember / Our World in Data** (production-based, all generation, attributes CO2 to bioenergy):

| Year | gCO2/kWh |
|---|---|
| 2019 | 141.9 |
| 2020 | 133.4 |
| **2021 (drought)** | **168.1** |
| 2022 | 103.4 |
| 2023 | 96.3 |
| 2024 | 106.1 |
| 2025 | 110.0 |

https://ourworldindata.org/grapher/carbon-intensity-electricity ·
https://ember-energy.org/countries-and-regions/brazil/

**MCTI "fator médio de emissão" for the SIN** — the official factor for Brazilian corporate
inventories, published monthly: 2025 Jan–Apr = **0.0237, 0.0248, 0.0215, 0.0289 tCO2/MWh**, i.e.
**~21–29 gCO2/kWh**.

The gap arises because (a) MCTI covers only the SIN, excluding isolated systems; (b) MCTI counts
fossil CO2 only, while Ember attributes CO2 to bioenergy (Brazil burns a lot of bagasse); (c) from
January 2025 ONS **widened the calculation base** to include biomass thermal and shared-connection
solar/wind *conjuntos* (https://www.gov.br/mcti/pt-br/acompanhe-o-mcti/cgcl/paginas/NT_FE_jun25.pdf).

**Recommendation:** use MCTI for anything aligning with a Brazilian inventory or certification;
Ember/OWID for cross-country comparison with EU grids; state which. For a specific plant the better
option is a **contracted/PPA-specific factor** — Brazil's free market (ACL) lets a plant contract
dedicated wind and solar, and any credible EU-facing claim needs that, not a residual mix.

**Hydrology is the tail risk.** 2021 was 168 gCO2/kWh — ~74% above 2023 — on thermal dispatch
covering depleted reservoirs. Monthly MCTI factors already range ±15% within a *normal* year.
**A dry year hits scope 2 and the power price simultaneously** — the correlation matters if both
are modelled.

### 8.11 Brazil and CBAM — the two gaps stack

- **Definitive regime started 1 January 2026.**
- **HBI/DRI is in scope: CN 7203.**
- **Indirect emissions are NOT priced for iron and steel** — Art. 7(1) + Annex II (see §4.1).

**This stacks with the SBCE asymmetry.** Brazilian green iron's dominant emission source —
electricity — is priced by **neither** SBCE (scope 1 only) **nor** CBAM (Annex II, direct only). An
H2-DRI plant on the Brazilian grid at ~100 gCO2/kWh consuming ~3.5–4 MWh/t HBI carries **~0.35–0.4
tCO2/t of scope 2 that is invisible to both systems**, while its EU competitor's grid power *is*
priced through the EU ETS on the generator and passed into the power price. A naive scope-1-only
comparison flatters the Brazilian route.

**Other CBAM mechanics worth carrying into a cash-flow model:**
- **Cash-flow timing:** **no CBAM certificates are sold during 2026 — sales open 1 February 2027.**
  The first declaration, covering calendar 2026, is due **30 September 2027**. So 2026 imports
  create a liability settled in late 2027, not a 2026 cash cost.
  https://www.dehst.de/EN/Topics/CBAM/CBAM-definitive-regime-2026/cbam-definitive-regime-2026_node.html
- **Default values carry punitive mark-ups: 10% in 2026, 20% in 2027, 30% from 2028.** **Verified
  actual emissions are strongly worth obtaining** — a genuinely low-carbon plant loses its entire
  advantage if it falls back to a default. Requires an accredited third-party verifier.
- **Omnibus, Regulation (EU) 2025/2083:** replaced the €150/consignment de minimis with a
  **50 tonne net-mass annual threshold** per importer, and cut the certificate-holding requirement
  from 80% to 50% of cumulative embedded emissions. Irrelevant to a bulk HBI cargo.
- **No Brazilian carbon price to deduct.** CBAM allows deducting a carbon price effectively paid at
  origin; Brazil prices nothing before ~2030–2031, so **the deduction is zero** — but a green iron
  plant's *direct* emissions are near zero anyway. The far bigger consequence is for Brazil's
  **conventional** exporters (Ternium, ArcelorMittal Tubarão, Gerdau semis), who pay full CBAM with
  no origin offset. That asymmetry drives the whole Brazilian green-iron narrative, and it is what
  the SBCE's own "international trade exposure" criterion is eventually meant to answer.

### 8.12 The asymmetry table

| | Scope 1 priced? | Scope 2 priced? |
|---|---|---|
| **Australia (Safeguard)** | Only above 100 kt/yr, and only above a 1.77 t/t benchmark → **credited, not charged** | **No** — sectoral baseline never binds |
| **Brazil (SBCE)** | Not before ~2030–2031, then only above 25 kt/yr | **No** — registered, not capped |
| **EU (CBAM) on imported HBI** | Yes, from 1 Jan 2026 (settled Sept 2027) | **No** — iron & steel is in Annex II |
| **EU EAF competitor** | EU ETS, free allocation phasing out 2026–2034 | **Yes**, through power prices |

Five consequences:

1. **Australia pays you; Brazil doesn't.** An Australian plant above 100 ktCO2e generates SMCs
   against a blast-furnace-era 1.77 t/t benchmark — ~2.6–2.8 m SMCs/yr for a 2 Mt/yr H2-DRI plant,
   ~1.4–1.8 m/yr even for NG-DRI. Bounded above by A$87.72/t (2026–27) and realistically well below.
   Brazil offers no equivalent.
2. **Both countries subsidise hydrogen and orphan the electrolytic routes.** Australia's A$2/kg HPTI
   and Brazil's PHBC/REHIDRO are both hydrogen-gated. **MOE and aqueous electrowinning qualify for
   neither.**
3. **Natural gas as a transitional reductant is not penalised anywhere relevant this decade.**
   Australia's GIIF explicitly permits NG-DRI with a documented pathway to renewables; the Safeguard
   credits it against 1.77 t/t; Brazil prices nothing. Only EU CBAM charges it, on direct emissions.
4. **Neither country has money on the ground.** GIIF closed early 2026 with no confirmed awards;
   Brazil's first PHBC auction is now 2027 with credits from 2030; no Brazilian green iron project
   has reached FID; Vale's Mega Hubs slowed in August 2026.
5. **Announced then cut or moved:**

| | What happened |
|---|---|
| Hydrogen Headstart Round 2 (AU) | **Halved**, A$2 bn → ~A$1 bn, 12 May 2026 Budget |
| Brazil PHBC R$18.3 bn | **Pushed two years**, 2028–2032 → **2030–2034** (Lei 15.269, 24 Nov 2025) |
| First PHBC auction (BR) | Targeted 2026 → **slipped to 2027** |
| SBCE Art. 56 insurer mandate (BR) | **Struck down by the STF**, 29 May 2026 (ADI 7795) |
| Vale Mega Hubs | **Decelerated**, 4 Aug 2026, Iran war |
| Energy Bill Relief Fund (AU) | **Not extended** beyond 2025 |
| Australian CBAM | **Cement first**; steel/iron deferred; no bill, no decision; **export rebates explicitly ruled out** |
| Green iron production credit (AU) | **Never announced or legislated** |

### 8.13 Open items worth one confirmation each

1. **GIIF closing date** — business.gov.au says 17 Feb 2026; DISR says 16 Jan 2026.
2. **Whether GIIF has made any awards** since ~May 2026. None found.
3. **Whether Safeguard Rule compilation No. 13 (31 Aug 2024) is still current** — a later Production
   Variables Update Rule could have revised the 1.77 t/t primary iron benchmark.
4. **Whether the 2026–27 Safeguard review has formally commenced** (ToR, consultation paper). The
   vehicle for post-2030 decline rates *and* the CBAM decision — the highest-value follow-up on the
   Australian side.
5. **Where a merchant HBI plant sits in Brazil's SBCE cohorts** — "integrated" vs "semi-integrated"
   steel doesn't cover a plant that never melts.
6. **Currency of the Superpower Institute's $170/t and $155/t figures** — read as AUD in context,
   not stated explicitly.
7. Australian Industry and Resources sector emissions reduction plan contents; the National Green
   Iron Statement; the 2026–27 Budget items checked against Budget Paper No. 2 rather than TSI's
   response.
8. Whether iron is genuinely absent from Australia's Critical Minerals and Strategic Materials
   Lists (high confidence, unverified against the DISR page).


## 9. What contradicts or complicates the modelling premise

Ordered by how much each would change a model's answer.

**1. The renewables pull may be much weaker than the headline literature says.**
Nykvist, Gong, Algers & Åhman (*Applied Energy* 395, 2025) harmonise assumptions across five value-
chain configurations and find the pull "sensitive to assumptions and weaker than previously found";
hydrogen-based steel cost varies across geographies "to a similar degree as conventional steel";
labour costs matter as much; and a modest subsidy on hydrogen or capital cancels the pull. Colen et
al. (*Business Strategy and the Environment*, 2025) reach a compatible conclusion from interviews:
renewables-pull factors dominate siting decisions **only when counteracting regulatory subsidies
are limited** — and EU support for domestic H2-DRI (IPCEI, European Hydrogen Bank, Agora's Phase 1)
is exactly such a subsidy. Any model showing Australia/Brazil dominating on electricity price alone
is probably over-reading the effect.

**2. Cost of capital, not electricity price, decides who exports.**
IRENA (2025) flips the top DRI exporter from the Middle East to **Australia (80% of exports with
the USA)** purely by moving from a Same-WACC to a Differentiated-WACC scenario. Seibold et al.
(2025) find that raising non-European WACC from 7% to 10% takes the European-only scenario from
47.7% to **78.3%** of the full global benefit. Agora attributes green HBI cost primarily to "cost of
capital in potential exporting countries" and assumes a **de-risked 4.3%**. The Superpower
Institute finds **capex is >60% of green iron cost**. If WACC is not a first-class sensitivity axis,
the Australia/Brazil-vs-Spain answer is not being tested.

**3. The EU producer gets free allocation the exporter does not.**
From 2026 the revised ETS benchmarks extend free allocation to **pelletising (sintered ore
benchmark), DRI reactors and sponge iron (hot metal benchmark), and electrolysis (hydrogen
benchmark)**, phasing out fully only in 2034 (Johnson et al. 2025). Sandbag quantifies the DRI
windfall: the **1.248 tCO2/t hot metal benchmark** applied to DRI whose real EU intensity is
**0.39 tCO2/t** = **0.858 EUA/t**. A model that gives the exporter a CBAM-free ride but does not
give the EU plant its free allocation has the asymmetry pointing only one way.

**4. CBAM does not merely fail to penalise imported DRI — it pays for it.**
Sandbag's Algeria case (Feb 2026): **+€54/t net profit** on DRI imported into the EU, at EUA €80
and 80% pass-through, because the benchmark net-out plus EU carbon-cost pass-through
(**+€94/t on the DRI market price**) exceeds the certificate cost. If a model applies zero CBAM cost
to exported HBI, that is *conservative*, not neutral.

**5. Both regimes miss the same emissions.** Australia's Safeguard Mechanism covers
**scope 1 only**; CBAM Annex II excludes **scope 2** for iron, steel and hydrogen. An Australian
green-iron exporter's electricity emissions are priced by **neither** jurisdiction. This is a
regulatory void, not a single-policy gap.

**6. The efficient trade depth may be steel, not iron.**
Caiafa, de Kleijne & de Coninck (*One Earth*, Feb 2026) find that producing green **steel** in
Ceará and shipping the steel beats shipping hydrogen on cost, emissions *and* local socioeconomic
outcomes. Grattan Institute argues for exporting green **steel** from Australia, against the
Superpower Institute's green **iron**. Verpoort et al. put the marginal saving from iron→steel at
only 5 percentage points (13% → 18%), but that is not zero, and the political-economy case for
stopping at iron (~90% of steel jobs are downstream) is an *importer* argument, not an exporter one.

**7. The global saving is small even where the bilateral saving is large.**
Bilici et al. (2024) find only **2.2–3.9%** of global steel production cost saved by green iron
trade, with 12–21% of global crude steel from traded iron in 2050. Verpoort's 18% and Seibold's
23.3% are Germany-specific; do not present them as a world figure.

**8. Ore grade may bind before energy cost does.**
DRI typically needs **≥67% Fe**; most Pilbara deposits are **56–62% Fe** (IEEFA). Agora puts Pará
ore at ~65% and Minas Gerais at ~62%. Devlin et al. show ore is 27% of cost on average but 45% in a
low-Fe case, with the DR-pellet premium tripling from $40/t to $122/t. Agora/Wuppertal warn
DR-grade ore demand will exceed supply by 2030. On this axis Brazil beats Australia, which is the
opposite of the electricity-cost ranking in some studies.

**9. The Pilbara is probably the wrong Australian site.**
The Superpower Institute: "the Pilbara is unlikely to be one of Australia's lower-cost locations…
It may make economic sense to ship ore from the Pilbara to other locations in Australia." Eyre
Peninsula (SA, 75% renewable) and Geraldton dominate. And a grid connection *helps* by enabling
arbitrage — a constrained connection raises cost from ~AUD 1,000/t to >AUD 1,200/t.

**10. Flexibility of the hydrogen step, not the iron step, drives cost.**
Devlin et al.: electrolyser oversizing 1.3–3.7× vs EAF oversizing only 1.1–1.4×; 91% of storage
cost is compressed hydrogen, not batteries; ~50% of H2 is stored. Seibold et al. by contrast assume
conservative **80% minimum part load for EAF and 90% for DRI, with no shutdowns**. These two
assumptions produce very different answers; which one is used should be stated.

**11. MOE and aqueous electrowinning are dated later than H2-DRI by the reference literature.**
Agora/Wuppertal: **MOE not market-ready before 2035, AEL not before 2040**. Also, an all-electric
route has aluminium's emissions profile, not steel's — and the DIW figure the Commission uses puts
reshuffling risk at **80% for aluminium versus 50% for steel**, "driven by the higher opportunities
to source or attribute the production of aluminium to clean electricity." The 80% number is
arguably the right analogue for MOE and electrowinning.

**12. Resource shuffling could erase the whole CBAM effect for steel.**
The Commission's own 2021 modelling: iron and steel leakage protection goes from **−24% to 0%** with
reshuffling; border revenue from €2.1 bn to €1.3 bn. Sandbag (2025): gross CBAM fees
**€11.3 bn → €7.3 bn**, net cost **€5.0 bn → €0.995 bn**. And Samadi et al. (2023) footnote 9 names
the exact iron case — a renewables-rich country routing existing hydro to EU-bound exports while
using fossil power domestically. Brazil, with an ~85–90% hydro/renewable grid, is the textbook
instance.

**13. Freight assumptions vary by an order of magnitude across the literature.**
Johnson et al.: **€4.44–39.27/t** by distance. Devlin et al.: ~**7% of cost** marine, 3% inland.
Seibold et al.: **€5.5/t** for 2,000 km by electric truck intra-EU. IEEFA: HBI is **~30% lighter**
than the equivalent pellets, and hydrogen shipping would add **US$128–230/t crude steel**.
SteelWatch: green iron needs **less than a third of the volume** of shipping ore and hydrogen
separately. Pick one and say which.

**14. Australia's Safeguard Mechanism pays a new green iron plant — a revenue line most models omit.**
A 2 Mt/yr H2-DRI plant sits ~1.3–1.4 tCO2e/t *below* the best-practice "primary iron" baseline of
1.77 tCO2e/t and generates roughly **2.6–2.8 million SMCs/yr**, bounded above by the cost
containment price of **A$87.72/t (FY2026-27, CPI+2%)** and realistically well below it. Even NG-DRI
earns ~1.4–1.8 m SMCs/yr. But note the cliff: MOE and aqueous electrowinning have near-zero scope 1,
so they likely fall **below the 100,000 tCO2e coverage threshold** and get nothing at all. The
Safeguard rewards a *slightly dirty* green iron plant and ignores a perfectly clean one.

**15. MOE and aqueous electrowinning are policy orphans in both exporting countries.**
Australia's A$2/kg Hydrogen Production Tax Incentive and Brazil's PHBC and REHIDRO are all
hydrogen-gated. Neither non-hydrogen route qualifies for any of them. The Superpower Institute
names this explicitly (Recommendation 2) and asks for a A$170/t credit extended to non-hydrogen
routes — a proposal, not policy. If a model gives MOE the same subsidy treatment as H2-DRI, it is
wrong in both jurisdictions.

**16. The Brazilian subsidy schedule moved two years right, and no auction has been held.**
Lei 15.269 (24 Nov 2025) pushed the R$18.3 bn PHBC credits from 2028–2032 to **2030–2034**, and the
first auction slipped from 2026 to **2027** (edital January 2027). Any model with PHBC revenue
before 2030 is out of date.

**17. If Australia builds its own CBAM, the exporter gets no rebate.**
The Carbon Leakage Review's Recommendation 1(c) explicitly rules out export rebates on trade-law and
target-integrity grounds. An Australian BCA would tax imports *into* Australia and do nothing for
HBI leaving for the EU.

**18. Brazil's grid emission factor differs by 4–5× between the two official sources.**
Ember/OWID gives ~96–110 gCO2/kWh recent years; MCTI's official SIN factor gives **~21–29
gCO2/kWh**. And 2021 hit **168 gCO2/kWh** on drought-driven thermal dispatch — a ~74% swing that
moves scope 2 and the power price together. Which factor is used, and whether hydrology risk is
correlated across the two channels, materially changes any Brazilian scope-2 scenario.

**19. Nothing is at FID.** No large-scale Australian green iron project has taken FID; the GIIF
closed in early 2026 with no confirmed awards; no Brazilian green iron project has reached FID; and
Vale decelerated its Middle East Mega Hubs in August 2026. Using the announced pipeline as a supply
curve will substantially overstate near-term deployment.
