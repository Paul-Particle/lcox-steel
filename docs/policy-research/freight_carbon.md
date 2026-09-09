# Carbon cost and emissions accounting of freight for imported green iron (HBI)

Research date: 2026-09-07. All prices in EUR unless stated. EUA spot on the day of writing: **€84.98/tCO2**.

Scope of the question: bulk-carrier sea freight Australia→EU (~19,500 km) and Brazil→EU (~9,500 km),
plus rail 300–1,500 km inside the EU. The model currently has freight tariffs but no carbon cost on freight.

---

## 0. Headline answer — numbers to put in the model

Per tonne of HBI delivered, EU ETS maritime cost at **€85/tEUA**, 2026 rules (100% phase-in, 50% of
extra-EEA voyage emissions in scope), **round-trip basis** (laden inbound leg + ballast outbound leg +
at-berth, which is what a carrier actually has to recover from the one cargo it carries):

| Route | In-scope CO2e per t cargo | ETS cost @ €85 | ETS cost @ €100 | ETS cost @ €145 (2030 base case) |
|---|---|---|---|---|
| Australia → Rotterdam, 19,500 km, capesize | **26.2 kg/t** | **€2.23/t** | €2.62/t | €3.80/t |
| Brazil → Rotterdam, 9,500 km, capesize | **12.9 kg/t** | **€1.10/t** | €1.29/t | €1.87/t |

Laden-leg-only variant (if you prefer not to charge the ballast return to the cargo):
Australia **13.9 kg/t → €1.18/t**; Brazil **6.8 kg/t → €0.57/t**.

CH4+N2O added to the ETS from 2026 contribute ~**1.6%** on liquid fossil fuel (already included above;
material only for LNG-fuelled ships).

Add **FuelEU Maritime** at roughly **+22%** of the ETS number in 2025–2029, rising to ~**+58%** from 2030
and ~**+134%** from 2035 (derivation in §2). So an all-in 2026 EU regulatory freight carbon cost of about
**€2.7/t Australia** and **€1.35/t Brazil**.

**Rail inside the EU is not separately chargeable under ETS1** if the traction is electric — the carbon
cost is already inside the electricity price the railway pays. If you want an explicit emissions number,
use GLEC v3.2 EU electric rail **10.8 gCO2e/tkm WTW**: 1,500 km = 16.2 kg/t → €1.38/t at €85; 300 km =
3.2 kg/t → €0.28/t. Diesel traction: **30.7 gCO2e/tkm WTW** → 46 kg/t over 1,500 km → €3.9/t.
**Do not add the electric-rail figure on top of an electricity price that already carries an ETS-inclusive
tariff** — that is double counting.

Context for materiality: capesize Brazil→Rotterdam freight is order **$15–25/t** and Australia→Rotterdam
order **$25–35/t**; green HBI production cost is order **€400–600/t**. So maritime ETS is **~0.4–0.6% of
delivered HBI cost** and **~7–9% of the ocean freight bill**. It is a real line item but it does not
change the export case (§7).

---

## 1. EU ETS maritime — precise scope from 2026

Legal basis: Directive (EU) 2023/959 amending Directive 2003/87/EC; MRV Regulation (EU) 2015/757 as amended.

### Ships covered
- **Cargo and passenger ships ≥5,000 GT** — in ETS since 1 January 2024.
- **Offshore ships ≥5,000 GT** — from 2027.
- **Offshore ships and general cargo ships 400–5,000 GT** — MRV reporting from 2025, ETS scope decision from 2028.
- Flag-neutral: applies regardless of flag or country of incorporation.
- A capesize bulk carrier (~180,000 dwt, ~93,000 GT) is comfortably in scope; **bulk carriers are treated
  identically to container ships** — there is no segment-specific rule. The only segment carve-outs are for
  passenger/ferry services (island derogation) and the transhipment-port anti-evasion rule, which is
  container-specific in effect.

### Voyages covered
- **100%** of emissions on voyages between two EEA ports.
- **100%** of emissions **at berth** in an EEA port (and manoeuvring within the port).
- **50%** of emissions on a voyage **departing from** an EEA port to a non-EEA port, **or arriving at** an
  EEA port from a non-EEA port.
- Important for a round-trip bulk trade: **the ballast leg out of the EU is a separate voyage and is also
  50% in scope.** So a carrier serving Brazil/Australia→Rotterdam pays on 50% of the laden inbound leg,
  100% of the berth stay, and 50% of the ballast outbound leg.

### Phase-in
| Emissions year | Share of verified in-scope emissions to surrender | Surrender deadline |
|---|---|---|
| 2024 | 40% | 30 Sep 2025 |
| 2025 | 70% | 30 Sep 2026 |
| **2026 onward** | **100%** | 30 Sep 2027 |

### Gases
- CO2 from 2024.
- **CH4 and N2O added from 1 January 2026** (already in MRV reporting from 2024). The EU cap was increased
  to reflect the added coverage.
- Practical uplift on a VLSFO/HFO-burning ship: N2O default ~0.00018 g/g fuel (GWP100 273) and CH4 slip
  ~0.00005 g/g fuel (GWP100 28) against a CO2 factor of 3.114–3.151 g/g fuel →
  (0.00018×273 + 0.00005×28)/3.151 ≈ **+1.6%**. Negligible for the model. For LNG dual-fuel ships methane
  slip can add 10–25% and is the reason the provision exists.

### Liable party and pass-through
- The liable entity is the **"shipping company"** = the shipowner, or whoever has assumed ISM Code
  responsibility for the ship's operation — in practice the **Document of Compliance (DoC) holder**
  (owner or technical manager), or a bareboat charterer where that applies. Liability cannot be avoided
  by flag, country of incorporation, manager or charterer.
- The Directive contains an **explicit right of recourse**: where responsibility for fuel purchase or for
  the operation of the ship lies with another entity by contract, the shipping company is entitled to
  **reimbursement from that entity** for the cost of surrendering allowances. This is what makes
  commercial pass-through legally underpinned rather than merely a market practice.
- In practice:
  - **Time charter**: BIMCO ETS Emission Scheme Surcharge Clause — the charterer supplies or pays for the
    EUAs, because the charterer controls speed and routing.
  - **Voyage charter / COA**: the owner prices ETS into the freight rate, usually as a separately quoted
    **ETS surcharge per tonne of cargo**, so it is auditable.
  - Container lines publish a per-container surcharge (Maersk EMS/ESS; Asia–N. Europe dry 40ft averaging
    ~$168 in 2026, 6–7% of base freight). **Dry bulk quotes it per tonne of cargo** — exactly the form
    your model wants.
  - Pass-through completeness is market-dependent: reported ~80–100% in strong markets, 40–60% absorbed by
    owners in weak markets. For a long-run cost model, assume **100% pass-through to the cargo**.
- T&E has documented container lines **over-recovering** ETS surcharges versus their actual liability, so
  a published surcharge is not a clean proxy for the underlying cost — compute it bottom-up.

### Exemptions
- **Outermost regions**: no surrender for voyages between an outermost-region port and the mainland of the
  same Member State, until 31 Dec 2030.
- **Small-island derogation**: passenger/ferry ships on routes to islands under 200,000 inhabitants with no
  road/rail link, until 2030.
- **Ice-class ships**: may surrender **5% fewer** allowances than verified emissions, until 2030.
- **Transnational public service obligation routes**: exempt to 2030.
- **Transhipment anti-evasion**: calls at designated neighbouring container transhipment ports
  (Tanger Med, East Port Said) do not break a voyage. Irrelevant to a direct HBI bulk trade.
- **No exemption applies to a direct extra-EEA dry bulk voyage.**

Sources:
- https://climate.ec.europa.eu/areas-action/transport-decarbonisation/reducing-emissions-shipping-sector/faq-maritime-transport-eu-emissions-trading-system-ets_en
- https://climate.ec.europa.eu/areas-action/transport-decarbonisation/reducing-emissions-shipping-sector_en
- https://www.emsa.europa.eu/reducing-emissions/extension-ets.html
- https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32023L0959
- https://www.bimco.org/news-insights/trending-topics/eu-ets/
- https://britanniapandi.com/2025/06/regulatory-overview-of-european-union-emissions-trading-system-eu-ets/
- https://www.qaship.net/eu-ets-for-shipowners/
- https://www.sustainable-ships.org/rules-regulations/eu-ets
- https://www.lloydslist.com/LL1148655/Maersk-and-MSC-overcharging-cargo-owners-for-EU-ETS-says-TE
- https://www.maersk.com/news/articles/2025/12/01/emissions-surcharge-ems-ess
- https://safety4sea.com/cm-new-eu-ets-rules-take-effect-from-january-2026-what-you-should-know/

---

## 2. FuelEU Maritime — Regulation (EU) 2023/1805

### What it requires
A **well-to-wake GHG intensity limit** on the energy used on board, in gCO2e/MJ, covering CO2, CH4 and N2O
across the whole fuel chain (WtT + TtW). Applies to the same ships (**≥5,000 GT** carrying cargo or
passengers commercially) with the **same geographic scoping as ETS** (100% intra-EEA and at berth, 50%
extra-EEA). It is a **fuel standard, not a carbon price** — it bites only on the gap to the target.

Baseline: **91.16 gCO2e/MJ** (2020 fleet average).

### Trajectory
| Period | Reduction vs baseline | Target GHG intensity (gCO2e/MJ) |
|---|---|---|
| 2025–2029 | −2% | 89.34 |
| 2030–2034 | −6% | 85.69 |
| 2035–2039 | −14.5% | 77.94 |
| 2040–2044 | −31% | 62.90 |
| 2045–2049 | −62% | 34.64 |
| 2050+ | −80% | 18.23 |

### Penalty mechanism
Annex IV: `Penalty (EUR) = (compliance deficit in gCO2e) / (actual GHG intensity × 41,000) × 2,400`.
The 2,400 is **EUR per tonne of VLSFO-equivalent energy** (41 GJ/t), ≈ **€0.058/MJ** of non-compliant energy.
Multiplier of +10% per consecutive year of deficit. Flexibility: **banking, borrowing (up to 2%), and
pooling** across ships, plus a wind-propulsion reward factor. In practice pooling and biofuel blending
clear well below the €2,400 headline, so the penalty is an economic ceiling, not the expected cost.

### Material per-tonne cost on a bulk trade
For a ship on 100% conventional fuel the penalty simplifies to a clean per-tonne-of-fuel figure:

`€/t fuel = (GHG_actual − GHG_target)/GHG_actual × 2,400`

Using VLSFO WtW ≈ 91.6 gCO2e/MJ (HFO WtT default 13.5 gCO2e/MJ plus TtW):
- **2025–2029**: (91.6−89.34)/91.6 × 2,400 = **€59/t fuel**
- **2030–2034**: (91.6−85.69)/91.6 × 2,400 = **€155/t fuel**
- **2035–2039**: (91.6−77.94)/91.6 × 2,400 = **€358/t fuel**

Compare EU ETS at €85/tEUA on the same fuel: 3.151 tCO2/t fuel × €85 = **€268/t fuel** (before the 50%
extra-EEA halving, which applies equally to both regimes).

So **FuelEU ≈ 22% of the ETS cost in 2025–2029, ~58% in 2030–2034, ~134% from 2035** — a real but
secondary add-on now, which **overtakes ETS in the mid-2030s**.

Applied to the routes: Australia round trip in-scope fuel per tonne of cargo ≈ 8.1 kg fuel/t cargo
→ 2026 FuelEU = 8.1 kg × €59/t = **€0.48/t of HBI**; Brazil ≈ 4.0 kg → **€0.24/t**.

Reported real-world magnitude: the first compliance cycle showed penalty exposure ≈ **10% of the annual
EU-scope fuel bill** for a typical conventionally fuelled EU-trading vessel; six figures in EUR for a
mid-size bulker.

Sources:
- https://ww2.eagle.org/en/rules-and-resources/regulatory-updates/fueleu-maritime.html
- https://ww2.eagle.org/en/rules-and-resources/regulatory-updates/fueleu-maritime/fuel-eu-faqs.html
- https://www.bettersea.tech/post/guide-how-to-calculate-fueleu-maritime-penalties
- https://www.bettersea.tech/post/guide-how-to-calculate-ghg-intensity-under-fueleu-maritime
- https://www.qaship.net/fueleu-maritime-compliance-pooling-penalties/
- https://www.sustainable-ships.org/rules-regulations/fueleu
- https://www.intercargo.org/wp-content/uploads/2025/05/2025-May-ESSF-SAPS-WS1-FuelEU-calculation-methodologies.pdf
- https://carboneer.earth/en/2025/08/fueleu-maritime-compliance-2025/

---

## 3. IMO Net-Zero Framework — status, timeline, prices, ETS interaction

### Status (as of September 2026) — NOT YET ADOPTED
- Approved at **MEPC 83, April 2025**.
- The **extraordinary MEPC session of October 2025 adjourned without adopting it**. The motion to delay,
  introduced by Singapore and formally submitted by Saudi Arabia, passed **57 in favour / 49 against /
  21 abstentions**, under heavy US pressure including threatened sanctions and tariffs.
- Talks **resume October 2026**. Nothing is in force. Commentators describe the plausible 2026 outcomes as
  "Survival, Coma, Death, Zombie, Mutant".
- If adopted in Oct 2026 on the original design, entry into force follows a 16-month tacit-acceptance
  period → **in force ~2028, first compliance year ~2028–2029**. Treat as a scenario, not a base case.

### Mechanism
A **global fuel standard** (GHG Fuel Intensity, gCO2e/MJ, well-to-wake) with **two targets** and a
**credit / remedial-unit trading scheme**. Reference value **93.3 gCO2e/MJ** (2008). Applies to
**ships ≥5,000 GT** on international voyages — ~85% of international shipping CO2. Purely domestic voyages
within a flag state's waters are excluded.

| Year | Base Target (Tier 2) | Direct Compliance Target (Tier 1) |
|---|---|---|
| 2028 | −4% | −17% |
| 2029 | −6% | −19% |
| 2030 | −8% | −21% |
| 2031 | −12.4% | −25.4% |
| 2032 | −16.8% | −29.8% |
| 2033 | −21.2% | −34.2% |
| 2034 | −25.6% | −38.6% |
| 2035 | −30% | −43% |

### Prices
- **Tier 1** (emissions between the Direct Compliance Target and the Base Target): **US$100/tCO2e**.
- **Tier 2** (emissions above the Base Target): **US$380/tCO2e**.
- Both fixed for **2028–2030**, reviewable thereafter. Both assessed on a **well-to-wake** basis.
- Ships that meet both targets pay **nothing** — unlike the EU ETS, which prices every tonne.
- Order of magnitude for a ship on plain VLSFO in 2028: Tier 2 on the 4% above base (~3.7 gCO2e/MJ) at $380
  plus Tier 1 on the 13% band (~12.1 gCO2e/MJ) at $100. Per tonne VLSFO (41 GJ): Tier 2 ≈ 0.153 tCO2e ×
  $380 ≈ $58; Tier 1 ≈ 0.496 tCO2e × $100 ≈ $50 → ~**$108/t fuel** in 2028, escalating steeply. Same order
  as FuelEU, roughly a third of ETS.

### Interaction with EU ETS — double charging
- **Article 3gg of Directive 2003/87/EC** is the review clause: if the IMO adopts a global market-based
  measure, the Commission must **review the Directive within 18 months of that adoption and before the
  global measure becomes operational**, assessing ambition against the Paris Agreement, overall
  environmental integrity, and coherence — explicitly to **avoid the risk of double payments while
  preserving the environmental integrity of the EU ETS**.
- So there is a *mandate* to adjust, but **no pre-agreed formula**. Legal commentary is clear that on the
  current texts a non-compliant ship would in theory be penalised under both regimes simultaneously.
- Article 3gg also contains the opposite lever: if **no comparable global measure is adopted by 2028**, the
  Commission is to consider **extending the surrender obligation beyond 50% of extra-EU voyage emissions**.
  Given the October 2025 adjournment this is a live risk — model a scenario where the Australia and Brazil
  routes move from 50% to 100% coverage, which simply **doubles** the §5 numbers.

Sources:
- https://www.imo.org/en/mediacentre/pressbriefings/pages/imo-net-zero-shipping-talks-to-resume-in-2026.aspx
- https://www.velaw.com/insights/imo-postpones-adoption-of-net-zero-framework/
- https://en.wikipedia.org/wiki/IMO_Net-Zero_Framework
- https://www.kslaw.com/insights/articles/imos-net-zero-framework-and-the-global-carbon-price-on-shipping
- https://globalmaritimeforum.org/news/a-guide-to-the-imos-net-zero-framework/
- https://www.opportunitygreen.org/factsheet-imo-net-zero-framework
- https://www.transportenvironment.org/uploads/files/Impact-of-the-IMOs-draft-Net-Zero-Framework-April-2025.pdf
- https://opportunitygreen.org/shipping/briefings/imo-net-zero-framework-adjournment-oct-2025-implications/
- https://www.hilldickinson.com/our-view/articles/imo-s-net-zero-framework-adoption-adjourned-and-consensus-building-continues/
- https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2026.1910244/full
- https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:02003L0087-20240301
- https://www.zerocarbonshipping.com/imo-net-zero-framework
- https://www.ricardo.com/en/news-and-insights/industry-insights/charting-a-course-to-decarbonisation-implications-of-the-imos-net-zero-framework

---

## 4. Emission factors

### 4.1 Dry bulk shipping, gCO2e per tonne-km

| Source | Vessel | Factor | Basis |
|---|---|---|---|
| GLEC Framework v3.2 (ISO 14083 aligned) | Bulk carrier **>200,000 dwt** (VLOC) | **3.1 gCO2e/tkm** | WTW, incl. load factor, empty running, +15% weather DAF |
| GLEC v3.2 | Bulk carrier **60,000–99,999 dwt** (panamax/kamsarmax) | **5.2 gCO2e/tkm** | WTW, same |
| GLEC v3.2 | Bulk carrier **0–9,999 dwt** | 31.2 gCO2e/tkm | WTW |
| GLEC 2019 (via Climatiq) | Bulk carrier **>100 dwkt, average load, MGO** | 2.6 gCO2e/tkm | **TTW** |
| CN Carbon Calculator, citing 4th IMO GHG Study 2020 | Bulk carrier (generic) | 3.94 gCO2e/tkm | — |
| 4th IMO GHG Study 2020 (range within the bulk carrier class) | — | **2.2 – 37.6 gCO2e/tkm** | — |
| NZ MfE 2020 (via Climatiq) | Bulk carrier 100,000–199,999 dwt | 3.04 gCO2e/tkm | fuel combustion |

Interpolating GLEC between the 60–100k band (5.2) and the >200k band (3.1), a **capesize (~180,000 dwt)
sits at ~3.4–3.8 gCO2e/tkm WTW**, i.e. **~2.8–3.1 gCO2/tkm TTW** once the ~20% WtT share is removed.

**Recommendation for a long-haul iron ore / HBI trade:**
- For **ETS and FuelEU costing** (which price TtW combustion), use **2.9 gCO2/tkm TTW, capesize,
  round-trip-inclusive** — the value implied by the bottom-up voyage model in §5, agreeing with GLEC to
  within ~10%.
- For **scope-3 / LCA reporting** (GHG Protocol, ISO 14083), use **3.5 gCO2e/tkm WTW** for capesize; use
  **5.2 gCO2e/tkm WTW** if the trade moves on panamax tonnage (relevant if HBI parcels are smaller than a
  full capesize stem).
- Note that HBI is **weight-limited, not volume-limited** in a bulk carrier (bulk density ~2.4–2.8 t/m3
  versus a capesize cubic of ~1.1 m3/t), so full deadweight is loadable and no stowage penalty applies —
  unlike DRI fines. This is one of the reasons HBI is the exportable form.

### 4.2 Rail freight in Europe, gCO2e per tonne-km

| Source | Mode | Factor | Basis |
|---|---|---|---|
| GLEC Framework v3.2 | **EU rail freight, electric traction** | **10.8 gCO2e/tkm** | WTW |
| GLEC Framework v3.2 | **EU rail freight, diesel traction** | **30.7 gCO2e/tkm** | WTW |
| GLEC v3.1 / UIC Railway Handbook | EU mix, weighted at 62% electrified | ~18.4 gCO2e/tkm | WTW |
| EEA, EU-27, 2018 | Rail freight average (road ~137 for comparison) | ~24 gCO2e/tkm | WTW |
| UBA (Germany), 2018 | Average freight train | 18 gCO2/tkm | — |

**Recommendation:** use **10.8 gCO2e/tkm** for an electrified EU corridor (which is what a port→EAF unit
train on a Rotterdam/Antwerp/Hamburg hinterland route will be) and **30.7 gCO2e/tkm** for a diesel leg.
A dedicated heavy bulk unit train will beat the EU average, so 10.8 is if anything conservative on the
high side for electric.

**Important — do not price EU electric rail carbon separately.** The traction electricity generator is an
ETS1 installation, so the ETS cost is already inside the traction power tariff. Charging it again in the
model double-counts against the electricity price. Diesel rail traction is **not** in ETS1; whether it
lands in **ETS2** (which covers fuels for buildings, road transport and certain additional sectors from
2027, with an active political discussion about deferral to 2028) I could not confirm within the search
budget — see open items in §8.

Sources:
- https://greencalculus.com/data/business-travel-freight-logistics-emission-factors/
- https://smart-freight-centre-media.s3.amazonaws.com/documents/2019_GLEC_Framework_July_2022.pdf
- https://www.smartfreightcentre.org/en/our-programs/emissions-accounting/global-logistics-emissions-council/glec-faq/
- https://ecofreight.co/blog/glec-v3-2-explained
- https://www.climatiq.io/data/source/glec
- https://www.climatiq.io/data/emission-factor/050ca9a9-5531-4a9b-8645-d9dd8d0764e6
- https://www.climatiq.io/data/emission-factor/009054a1-6e08-4d62-b943-2fa9322b0e9d
- https://greenvoyage2050.imo.org/wp-content/uploads/2021/07/Fourth-IMO-GHG-Study-2020-Full-report-and-annexes_compressed.pdf
- https://www.imo.org/en/ourwork/environment/pages/fourth-imo-greenhouse-gas-study-2020.aspx
- https://lune.co/blog/what-is-carbon-intensity
- https://www.cn.ca/repository/popups/ghg/Carbon-Calculator-Emission-Factors
- https://www.eea.europa.eu/en/analysis/maps-and-charts/ghg-efficiency-of-different-transport-figures
- https://www.eea.europa.eu/en/analysis/publications/rail-and-waterborne-best-for-low-carbon-motorised-transport
- https://cedelft.eu/wp-content/uploads/sites/2/2021/03/CE_Delft_190325_STREAM_Freight_Transport_2020_FINAL.pdf
- https://cedelft.eu/wp-content/uploads/sites/2/2021/05/CE_Delft_200258_Methodology_GHG_Efficiency_Transport_Modes.pdf
- https://www.cleanenergywire.org/factsheets/rail-cargo-emissions-germany

---

## 5. EUA prices and the per-tonne ETS calculation

### 5.1 EUA price
- **Spot, 7 September 2026: €84.98/tCO2** (+0.95% on the day, +3.3% on the month, +10.1% year on year;
  highest since July 2026). Source: Trading Economics EU Carbon Permits.
- 2026 trading range so far: **€74–86.60**; the Dec-26 ICE benchmark hit a six-month high of €86.60 on
  22 July 2026. EUAs traded €68–76 through 2025. All-time high €105.73 (Feb 2023).
- Analyst forecasts:
  - Survey of 10 analysts: **€92.65/t average 2026**, **€107.29/t 2027**.
  - A more bearish revision: **€82/t average 2026**.
  - **Citi / Bloomberg carbon model base case: €145/t by 2030.** Another long-term baseline: **€138/t by 2030.**
- Drivers of the upward path: cap tightening, the end of front-loaded auctioning, CBAM full phase-in with
  the parallel phase-out of free allocation.

**Recommended model inputs:** €85/t for 2026, ramping linearly to **€145/t by 2030**
(≈ €100 in 2027, €115 in 2028, €130 in 2029). Sensitivity band €80–160.

### 5.2 Voyage assumptions (bottom-up)

| Parameter | Value | Basis |
|---|---|---|
| Vessel | Capesize, 180,000 dwt | standard for long-haul iron/HBI |
| Cargo intake | 170,000 t | dwt less bunkers, stores, constants |
| Laden speed | 12.5 kn | Baltic Exchange guidance: 12 kn laden, 13 kn ballast |
| Laden consumption | 42 t/day VLSFO | Baltic Exchange ~44 t/day; 2016-built 182,620 dwt case: 41.5 t/day laden, 38.1 t/day ballast |
| Ballast speed / consumption | 13 kn / 38 t/day | as above |
| At berth, Rotterdam | 3 days at 5 t/day = 15 t | grab discharge, aux + cargo gear |
| CO2 factor, VLSFO | **3.151 tCO2/t fuel** | EU MRV / IMO (HFO 3.114, MGO 3.206) |
| CH4+N2O uplift from 2026 | ×1.016 | N2O 0.00018 g/g × GWP 273; CH4 slip 0.00005 g/g × GWP 28 |
| EUA price | €85/t | spot 7 Sep 2026 |
| Phase-in | 100% | 2026 onward |

### 5.3 (a) Australia → Rotterdam, 19,500 km

```
Distance          19,500 km / 1.852         = 10,529 nm
Laden days        10,529 / 12.5 kn / 24     = 35.1 days
Laden fuel        35.1 x 42                 = 1,474 t VLSFO
Laden CO2         1,474 x 3.151             = 4,644 tCO2
Ballast days      10,529 / 13 kn / 24       = 33.8 days
Ballast fuel      33.8 x 38                 = 1,284 t VLSFO
Ballast CO2       1,284 x 3.151             = 4,046 tCO2
Berth CO2         15 t x 3.151              =    47 tCO2

In-scope (2026 rules):
  laden leg   50% x 4,644                   = 2,322 tCO2
  ballast leg 50% x 4,046                   = 2,023 tCO2
  at berth   100% x    47                   =    47 tCO2
  subtotal                                  = 4,392 tCO2
  + CH4/N2O  x 1.016                        = 4,462 tCO2e

Per tonne of cargo   4,462 / 170,000        = 0.0262 tCO2e/t = 26.2 kg/t
ETS cost @ EUR 85    26.2 kg x 0.085        = EUR 2.23 per tonne of HBI
```
Laden-leg-only variant: 2,322 × 1.016 / 170,000 = 13.9 kg/t → **€1.18/t**.

Cross-check by the factor method: 2.9 gCO2/tkm TTW (capesize, ballast-inclusive) × 19,500 km = 56.6 kg/t
round trip; 50% in scope = **28.3 kg/t** versus 26.2 bottom-up. Agreement within 8%.

### 5.4 (b) Brazil → Rotterdam, 9,500 km

```
Distance          9,500 km / 1.852          = 5,130 nm
Laden days        5,130 / 12.5 / 24         = 17.1 days
Laden fuel        17.1 x 42                 =   718 t VLSFO
Laden CO2         718 x 3.151               = 2,262 tCO2
Ballast days      5,130 / 13 / 24           = 16.4 days
Ballast fuel      16.4 x 38                 =   624 t VLSFO
Ballast CO2       624 x 3.151               = 1,966 tCO2
Berth CO2                                   =    47 tCO2

In-scope: 0.5 x 2,262 + 0.5 x 1,966 + 47    = 2,161 tCO2
  + CH4/N2O x 1.016                         = 2,196 tCO2e

Per tonne of cargo  2,196 / 170,000         = 0.0129 tCO2e/t = 12.9 kg/t
ETS cost @ EUR 85   12.9 kg x 0.085         = EUR 1.10 per tonne of HBI
```
Laden-leg-only variant: 1,131 × 1.016 / 170,000 = 6.8 kg/t → **€0.57/t**.

Factor cross-check: 2.9 × 9,500 = 27.6 kg/t round trip; 50% = **13.8 kg/t** versus 12.9. Within 7%.

### 5.5 Price sensitivity

| EUA price | Australia €/t HBI | Brazil €/t HBI | Delta (AU − BR) |
|---|---|---|---|
| €80 | 2.10 | 1.03 | 1.07 |
| **€85 (spot)** | **2.23** | **1.10** | **1.13** |
| €100 | 2.62 | 1.29 | 1.33 |
| €120 | 3.14 | 1.55 | 1.59 |
| €145 (2030 base) | 3.80 | 1.87 | 1.93 |
| €160 | 4.19 | 2.06 | 2.13 |
| €145, 100% coverage scenario (Art. 3gg fallback) | 7.60 | 3.74 | 3.86 |

**The Australia-vs-Brazil ETS differential is only ~€1.1/t at today's carbon price and ~€1.9/t at €145.**
The ocean freight differential itself is order **$10/t**. Carbon is a second-order term in the
Australia-vs-Brazil comparison.

### 5.6 Independent cross-checks on the ETS numbers
- **CRU Group**: "At current carbon prices of $60/t, ETS-related freight costs would add around **$1–3/t or
  0.5–1.5% to delivered iron ore prices**", modelling Australia→Rotterdam capesize, Brazil panamax and
  Russian handymax routes at ~0.07 tCO2e per tonne of cargo round trip. My Australia figure at ~$60 carbon
  would be ~$1.6/t — inside their range.
- **Widely reproduced worked example**: Capesize, Ponta da Madeira (Brazil) → Rotterdam, 4,100 nm, 14 kn,
  ~12 days, **62 t/day**, **2,300 tCO2** laden leg; 50% in scope = 1,150 t; at the early phase-in rate the
  cited cost was €19,550 at €85/t. Scaled to 100% phase-in that is 1,150 t × €85 = €97,750 over ~170,000 t
  = **€0.57/t** — exactly my laden-only Brazil figure.
- **shipfinex**: 82,000 dwt panamax on iron ore Australia/Brazil→Europe with 50–60% EU exposure,
  **24,255 tCO2 covered annually**, **€1,940,400/yr at €80/t**. Over ~6 laden EU voyages a year at 76,000 t
  that is ~€4.3/t — higher, as expected, because a panamax is ~40% less efficient per tonne-km than a
  capesize. **If HBI parcels ship panamax-size, scale the §5 numbers up by ~1.5×.**
- **Dry bulk market estimate**: capesize Brazil–Netherlands ETS surcharge **$0.30/t in 2024 (40% phase-in)
  rising to $0.74/t in 2026 (100% phase-in)**. My Brazil laden-only figure of €0.57 ≈ $0.62 is consistent;
  their $0.74 sits between my laden-only and round-trip figures.

All four independent checks land in the **$0.5–$4/t** band and my central numbers sit in the middle of it.

Sources:
- https://tradingeconomics.com/commodity/carbon
- https://carboncredits.com/carbon-prices-today/
- https://gmk.center/en/news/the-price-of-carbon-in-the-eu-in-2026-will-be-e92-6-t-analysts/
- https://www.homaio.com/post/2030-eua-price-predictions-expert-analysis-of-3-scenarios
- https://www.abnamro.com/research/en/our-research/carbon-market-strategist-carbon-prices-heat-up-in-2026
- https://www.materiaintel.com/insight/eu_ets_review_eua_price_forecast
- https://www.crugroup.com/en/communities/thought-leadership/sustainability/proposed-revisions-to-the-eus-ets-impact-on-maritime-transport-costs/
- https://safety4sea.com/prepare-for-higher-shipping-costs-but-the-eu-ets-should-be-a-manageable-change/
- https://www.shipfinex.com/blog/eu-ets-for-shipping
- https://shipandbunker.com/news/emea/779990-new-guidance-on-capesize-fuel-use-from-baltic-exchange
- https://www.lloydslist.com/LL1146970/EU-ETS-explainer-How-will-the-worlds-first-carbon-tax-on-shipping-work
- https://searoutes.com/2026/05/28/eu-ets-surcharges-shippers-audit-2026/

---

## 6. How freight is treated in each accounting framework

### (a) CBAM — **transport is OUTSIDE embedded emissions**
- Regulation (EU) 2023/956 Art. 3(4) and Annex IV define embedded emissions as **direct emissions from the
  production processes of goods** (including heat/cooling consumed in production) plus, for some goods,
  **indirect emissions from electricity consumed in production**. Both definitions are bounded by
  "**during the production processes**".
- There is **no transport term**. Emissions from moving the finished good from the third-country plant to
  the EU port are not embedded emissions and generate no CBAM certificate obligation. Precursor transport
  is likewise excluded.
- Consequence for your model: **CBAM does not penalise distance.** The maritime ETS does — but the maritime
  ETS charge sits on the shipowner's side of the boundary, so the two never interact and there is no double
  count between CBAM and maritime ETS on the same tonne.
- **CN code caveat worth verifying**: DRI is CN **7203** ("ferrous products obtained by direct reduction of
  iron ore and other spongy ferrous products"), which is on the CBAM Annex I list. Sources disagree on
  whether **hot-briquetted iron** classifies under 7203.10 or 7204.10 — trade data services show HBI
  imports/exports declared under 7203, while at least one guidance site asserts 7204.10. This matters
  because 7204 is ferrous waste and scrap. **Check the actual customs classification for your HBI before
  relying on the CBAM treatment**, since it determines whether the import carries a certificate obligation
  at all. There is a 50 t/year mass de minimis on all Chapter 72/73 CBAM goods.
- For calibration of the competing carbon costs: EU hot metal carries a net carbon cost of roughly
  **€30/t of metallic in 2026**, while imported DRI declared on default values rises from about
  **€50 to €80/t** — one to two orders of magnitude above the €1–2/t freight ETS number.

Sources:
- https://eur-lex.europa.eu/eli/reg/2023/956/oj/eng
- https://eur-lex.europa.eu/legal-content/EN/TXT/PDF/?uri=CELEX:32023R0956
- https://taxation-customs.ec.europa.eu/document/download/29b9eec7-1a4b-4eb6-ab85-96a0c9e35fd0_en?filename=Guidance+No.+3+-+CBAM+methods+for+the+calculation+of+emissions+embedded+in+goods.pdf
- https://www.carbonchain.com/cbam/eu-cbam-cn-codes
- https://cbamguide.com/sectors/steel/cn-codes/
- https://www.zauba.com/import-HOT+BRIQUETTED+IRON/hs-code-7203-hs-code.html
- https://www.metallics.org/about-metallics/hbi/
- https://www.fastmarkets.com/insights/the-free-allocation-gap-keeping-eu-pig-iron-carbon-cost-competitive-against-imported-low-emissions-dri/

### (b) RFNBO methodology, Delegated Regulation (EU) 2023/1185 — **transport is INSIDE, via e_td**
- The DR supplements Directive (EU) 2018/2001, setting the minimum GHG saving threshold for recycled carbon
  fuels and the methodology for assessing GHG savings from RFNBOs — **70% saving** against a fossil
  comparator of **94 gCO2e/MJ**.
- Life-cycle emissions are `e = e_i + e_p + e_td + e_u − e_ccs …`, and the **e_td term explicitly covers
  emissions from transport and distribution**, including "emissions from the storage and distribution of
  the finished fuels", plus transport of inputs and intermediates. All stages — electricity generation,
  conversion, compression, transport and storage — are in the boundary.
- Relevant to your model **only in that the hydrogen or e-fuel feeding the DRI shaft must carry its own
  e_td**. It does **not** attach to the HBI itself: HBI is not a transport fuel and is outside the RFNBO
  methodology's product scope. If your Australian or Brazilian plant imports or moves hydrogen or ammonia,
  that movement is inside e_td and must still clear the 70% threshold; the HBI ocean leg is not.
- Worth flagging the asymmetry: **RFNBO rules charge fuel transport, CBAM does not charge product
  transport.** The EU's two regimes are inconsistent on this point.

Sources:
- https://eur-lex.europa.eu/eli/reg_del/2023/1185/oj/eng
- https://www.efta.int/eea-lex/32023r1185
- https://rsb.org/wp-content/uploads/2025/11/rsb-std-11-103-v.1.0-eu-standard-for-rfnbos-and-rcfs.pdf
- https://www.ecologic.eu/sites/default/files/publication/2024/60022-FAQ-EU-requirements-green-hydrogen-and-PtX.pdf
- https://ptx-hub.org/wp-content/uploads/2023/04/International-PtX-Hub_EU-Requirements-for-green-hydrogen-and-PtX.pdf

### (c) worldsteel CO2 Data Collection / ISO 14404 — **transport is OUTSIDE**
- ISO 14404-1/-2/-3/-4 (CO2 emission intensity from iron and steel production; -4 covers the
  DRI/scrap-based EAF route) is built on worldsteel's **CO2 Emissions Data Collection User Guide** (v11),
  collecting site-level data from 200+ facilities annually since 2007. ISO 14404 is based on the worldsteel
  methodology but is not identical to it.
- The system boundary is the **steel plant site**, extended to selected upstream carriers.
  **Indirect emissions from transport of raw materials are not included** in the worldsteel CO2 data
  collection methodology.
- So an EU EAF reporting under worldsteel/ISO 14404 shows **no penalty for importing HBI from Australia
  rather than from a nearby EU plant** — the ocean and rail legs simply do not appear in its intensity
  figure. This is a real reporting blind spot for your comparison and worth stating explicitly if the model
  is used to argue on the basis of reported steel intensities.

Sources:
- https://worldsteel.org/wp-content/uploads/CO2_User_Guide_V11.pdf
- https://cdn.standards.iteh.ai/samples/77622/8771f743065640c7976465afd639c323/ISO-14404-4-2020.pdf
- https://www.iso.org/standard/57298.html
- https://iea.blob.core.windows.net/assets/8f6568aa-1dd8-4578-bc61-24ceba4a07dd/EmissionsMeasurementandDataCollectionforaNetZeroSteelIndustry.pdf
- https://worldsteel.org/climate-action/climate-change-and-the-production-of-iron-and-steel/

### (d) GHG Protocol — **Scope 3, Category 4 for the steelmaker**
- Inbound freight of purchased HBI falls in **Scope 3 Category 4, Upstream Transportation and Distribution**.
- The Corporate Value Chain (Scope 3) Standard is explicit: Category 4 covers transportation and
  distribution of products purchased in the reporting year between the company's **tier-1 suppliers and its
  own operations**, in vehicles and warehouses not owned or controlled by the reporting company — and
  **"inbound transportation of purchased goods should always be reported under Category 4, even if the
  reporting company does not explicitly pay for it"** (i.e. even on CIF/CFR terms where the seller pays freight).
- It is **not** Scope 1 (the steelmaker owns no ships or locomotives) and **not** Scope 2.
- Overlap with Category 1: cradle-to-gate emissions of the HBI itself are Category 1 (Purchased Goods and
  Services). If the supplier's cradle-to-gate factor already includes delivery, avoid counting the ocean
  leg in both Category 1 and Category 4.
- The rail leg from EU port to EAF is also Category 4 if the EAF operator does not own the locomotives.
- Methodologically, GHG Protocol Cat. 4 is calculated exactly as ISO 14083 / GLEC prescribes — so the
  **WTW** factors in §4 (3.5 gCO2e/tkm capesize, 10.8 gCO2e/tkm EU electric rail) are the right ones for
  disclosure, while the **TTW** factors are the right ones for ETS cost.

Sources:
- https://ghgprotocol.org/sites/default/files/2022-12/Chapter4.pdf
- https://ghgprotocol.org/sites/default/files/standards/Scope3_Calculation_Guidance_0.pdf
- https://plana.earth/glossary/scope-3-category-4
- https://www.climatiq.io/docs/guides/understanding/freightv3-ISO14083

### Summary table

| Framework | Ocean leg AU/BR → EU port | EU rail leg port → EAF | Who bears it |
|---|---|---|---|
| **EU ETS maritime** | **IN**, 50% of voyage emissions (each leg) + 100% at berth | out of ETS1 scope directly; electric-rail carbon is inside the power price | shipowner / DoC holder, with statutory right of recourse → passed to cargo |
| **FuelEU Maritime** | **IN**, same 50%/100% scoping, WtW intensity basis | n/a | shipowner, passed via charter/freight |
| **IMO Net-Zero Framework** | would be **IN**, 100% of the international voyage, WtW — but **not adopted**; talks resume Oct 2026 | n/a | shipowner |
| **CBAM embedded emissions** | **OUT** — production-process boundary only | **OUT** | n/a |
| **RFNBO DR (EU) 2023/1185** | **OUT** for the HBI; **IN** via `e_td` for any hydrogen/ammonia moved | same | fuel producer / certifier |
| **worldsteel / ISO 14404** | **OUT** — raw material transport not in the boundary | **OUT** | n/a |
| **GHG Protocol** | **IN — Scope 3 Category 4** | **IN — Scope 3 Category 4** | the EU steelmaker |

**The key structural point**: the frameworks that *price* carbon (ETS, FuelEU) charge the freight leg, but
only at 50% and only to the shipowner; the frameworks that *measure* the product's carbon for regulatory
purposes (CBAM, worldsteel/ISO 14404) exclude it entirely; and the framework that fully counts it
(GHG Protocol Scope 3.4) attaches no price. Distance is therefore **priced weakly and disclosed
inconsistently**. Your model is right to add a freight carbon cost, but it should be added as a **cost**,
not as an addition to the HBI's declared embedded emissions.

---

## 7. Published analysis: does maritime carbon cost erode the green iron export case?

**The consistent finding across every source I found is: no, not materially.**

- **CRU Group**: ETS on freight adds **$1–3/t, or 0.5–1.5% of delivered iron ore prices**, at $60/t carbon.
  Their framing is explicitly that this is a manageable change, not a structural shift.
- **OECD, "Green Iron opportunities in Australia"**: **shipping costs are a relatively small component of
  total production costs and would not disadvantage Australian green iron producers.** The OECD's actual
  competitiveness concerns are elsewhere: (i) green iron is **not yet cost-competitive** because of the
  energy intensity of production; (ii) **carbon pricing among trading partners is currently insufficient to
  provide the price signal** needed to justify investment; and (iii) **DRI/HBI is more expensive to ship
  per unit of contained Fe than iron ore is per unit of ore** — though because you move ~90%-Fe metallic
  instead of ~62%-Fe ore, the freight per tonne of *contained Fe delivered* improves. The OECD notes that
  for **markets far from Australia such as Europe these cost factors weigh heavier and favour MENA
  competitors** — but the driver is the base freight rate and energy cost, not the ETS.
- **SteelWatch**: shipping green iron as HBI needs **3.5–4.7× less transport capacity** than shipping the
  hydrogen and the ore separately — a **~75% saving in transport capacity**. This is the strongest freight
  argument in favour of the export model and it dwarfs the ETS effect.
- **UWA / Curtin work on an IMO-style global levy applied to Port Hedland–Shanghai iron ore**: the finding
  is about **fuel switching**, not about killing the trade — blue ammonia may become the cheaper bunker
  fuel around 2030, and green ammonia could compete with HFO in the 2030s even without subsidies. A carbon
  price on shipping changes which fuel the ship burns long before it changes whether the cargo moves.
  The companion journal paper on the WA–East Asia green shipping corridor reaches the same conclusion.
- **Fastmarkets (2026)**: the competitive battleground for imported low-emission DRI/HBI into Europe is the
  **free-allocation gap** — EU hot metal at ~**€30/t of metallic** net carbon cost in 2026 versus imported
  DRI on default values rising from **€50 to €80/t**. That carbon asymmetry is **20–50× larger than the
  €1–2/t freight ETS charge**, and it runs in the *opposite* direction from what freight distance would
  suggest.
- **Wang et al., South Australian Green Iron Supply Chain Study**, and IEEFA's "Costing Australia's green
  iron export ambitions" reach the same conclusion from the cost-stack side: the binding constraints are
  **renewable electricity cost, electrolyser capex, and the AU$170bn/yr capital requirement** (against a
  projected AU$96bn/yr export revenue by 2040), not freight.

**Bottom line for your model.** Maritime ETS costs **€2.2/t on the Australian route and €1.1/t on the
Brazilian route** at today's carbon price — about **0.4–0.5% of a €450–500/t green HBI cost** and about
**€1.1/t of Australia-vs-Brazil disadvantage**. Even at €145/t EUA with a hypothetical extension to 100%
voyage coverage, the Australia figure reaches €7.6/t, still **under 2%** of delivered cost. Include it for
completeness and because it is directional (it always penalises the longer haul), but **it should not be
expected to flip any siting decision in the model**. The variables that will flip siting are renewable
electricity LCOE, the CBAM free-allocation glidepath, and the base freight rate — the last of which is
roughly **10× larger and far more volatile** than the carbon component.

Sources:
- https://www.crugroup.com/en/communities/thought-leadership/sustainability/proposed-revisions-to-the-eus-ets-impact-on-maritime-transport-costs/
- https://www.oecd.org/en/publications/green-iron-opportunities-in-australia_bbd1e2b8-en/full-report/component-6.html
- https://steelwatch.org/press-releases/importing-green-iron-saves-75-of-transport-capacity-needs-compared-to-separate-hydrogen-and-iron-ore-shipments/
- https://www.uwa.edu.au/news/article/2025/june/a-carbon-levy-on-global-shipping-promises-to-slash-emissions-we-calculated-what-that-means-for-australias-biggest-export
- https://www.sciencedirect.com/science/article/abs/pii/S0306261925001953
- https://www.fastmarkets.com/insights/the-free-allocation-gap-keeping-eu-pig-iron-carbon-cost-competitive-against-imported-low-emissions-dri/
- https://www.fastmarkets.com/insights/global-green-steel-markets-in-2026-regulation-costs-and-regional-divergence/
- https://ieefa.org/articles/costing-australias-green-iron-export-ambitions
- https://ieefa.org/resources/australia-faces-growing-green-iron-competition-overseas
- https://www.energymining.sa.gov.au/__data/assets/pdf_file/0011/1006985/66ed5bccec444ce135c4d0697da4f4450da87002.pdf
- https://www.nature.com/articles/s41467-025-64440-9

---

## 8. Suggested model implementation

```python
# per tonne of HBI delivered
sea_co2_intensity_ttw  = 2.9e-6   # tCO2 per tonne-km, capesize, ballast-inclusive
ets_scope_extra_eea    = 0.50     # 50% of extra-EEA voyage emissions
ch4_n2o_uplift         = 1.016    # CH4 + N2O added to the ETS from 2026
eua_price              = 85.0     # EUR/tCO2 in 2026; ramp to 145 by 2030
fueleu_multiplier      = 1.22     # 2025-2029; 1.58 from 2030; 2.34 from 2035

sea_carbon_cost = (distance_km * sea_co2_intensity_ttw
                   * ets_scope_extra_eea * ch4_n2o_uplift
                   * eua_price * fueleu_multiplier)

# Australia 19,500 km -> EUR 2.72/t ; Brazil 9,500 km -> EUR 1.33/t
```

Rail: **do not add a separate carbon charge for electric traction** — it is already in the traction
electricity price. If the model prices rail energy separately and wants an explicit emissions number for
reporting, use **10.8 gCO2e/tkm** electric / **30.7 gCO2e/tkm** diesel, WTW.

### Open items to check
1. **HBI customs classification (CN 7203 vs 7204)** — determines whether the import is a CBAM good at all.
2. **ETS2 coverage of rail diesel**, and whether ETS2 starts in 2027 or is deferred to 2028. The web-search
   budget ran out before I could confirm; only relevant if a diesel rail leg is modelled explicitly.
3. **Panamax vs capesize parcel size** — if HBI ships in panamax lots, multiply the sea carbon cost by ~1.5.
4. **Article 3gg fallback**: if no comparable IMO measure is adopted by 2028, extra-EU coverage may rise
   above 50%. Worth a scenario switch (`ets_scope_extra_eea = 1.0`) given the October 2025 adjournment.
5. **Load-factor assumption**: my numbers assume 170,000 t of cargo on a 180,000 dwt capesize. If HBI
   parcels are smaller or the trade uses part cargoes, the per-tonne figures scale inversely with intake.
