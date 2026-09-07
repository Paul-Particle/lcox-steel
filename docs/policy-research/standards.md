# Voluntary and emerging-regulatory standards for low-carbon / green steel and iron

Research compiled 7 September 2026. Focus: system boundary, treatment of scope 2
(electricity), thresholds, scrap, transport, adoption — with the specific question of
whether an exported HBI made on a coal-heavy grid in Australia or Brazil, melted in an
EU EAF, is visible to any of these frameworks.

**Framing.** The EU CBAM excludes indirect (electricity) emissions for iron, steel and
hydrogen, so the *regulatory* border regime is blind to the grid the iron was made on.
Everything below is the voluntary / labelling layer that could in principle see it.

---

## 0. The landscape in one paragraph

Almost every credible framework in this space descends from a single 2022 IEA proposal:
a **sliding scale** of allowed CO2e intensity running from ~0.4 t/t crude steel at 0%
scrap down to ~0.05 t/t at 100% scrap, on a **cradle-to-gate** boundary that includes
**Scope 1 + Scope 2 + partial upstream Scope 3**. ResponsibleSteel, LESS, the First
Movers Coalition, IDDI and China's CISA C2F all adopt that shape. The consequence for
this project is structural: because the boundary is cradle-to-gate at the *steelmaking
site*, **purchased HBI enters the EU EAF's account as Scope 3.1 with the HBI producer's
own Scope 1 + 2 + 3 embedded in it**. The grid that made the iron is therefore *not*
invisible to these standards, unlike under CBAM.

---

## 1. ResponsibleSteel International Standard v2 — Decarbonisation Progress Levels

Primary sources:
- ResponsibleSteel, *Fundamentals for GHG Emissions Accounting and Classification*,
  Version 1.0 (published 18 June 2025).
  PDF: https://cdn.prod.website-files.com/653ed7060b01292cd4518d0e/68a72ce90dc0571ac8e65679_ResponsibleSteel%20Fundamentals%20for%20GHG%20Emissions%20Accounting%20%26%20Classification%20-%20Clean%20-%20V1.0.pdf
  Landing page: https://www.responsiblesteel.org/resources/fundamentals-for-ghg-emissions-accounting-and-classification
- https://www.responsiblesteel.org/emissions-accounting-classification
- https://www.responsiblesteel.org/steel-decarbonisation-scale
- https://www.responsiblesteel.org/resources/the-decarbonisation-progress-levels
- ResponsibleSteel International Production Standard v2.1 / v2.1.1
- https://www.responsiblesteel.org/become-certified

### (a) System boundary

Two *distinct and complementary* metrics:

1. **Site-level "ResponsibleSteel crude steel GHG emissions intensity"** — the metric that
   drives the Decarbonisation Progress Level. Boundary is **cradle-to-crude-steel**,
   measured at the *site*. Crude steel production is "measured at the point that
   continuous casting or ingot casting has been [completed]". Hot rolling and all
   downstream processing are **excluded**.
   - `intensity = total GHG emissions (t CO2e) for the previous year of operation
     / saleable tonnes of crude steel produced in the previous year`.
2. **Product-level Product Carbon Footprint (PCF) / EPD** — required at Progress level for
   any steel *marketed or sold* as ResponsibleSteel certified. Boundary is **cradle-to-gate
   including downstream processing**, "Inclusive, as a minimum, of the 'cradle-to-gate'
   emissions associated with extracting and processing raw materials, transportation of
   materials to site, and product manufacturing (i.e. life cycle Modules A1-A3)."
   Module D (recycling) may be added but "must be reported separately".

Gases: **CO2, N2O, CH4, HFCs, PFCs, SF6, NF3** — full CO2e, not just CO2.

**Does it cover iron/DRI separately from crude steel?** Not as a graded metric.
"The aspects of Principle 10 contained in this document apply to operational steelmaking
sites. The methodology for the calculation of crude steel GHG emissions intensity at site
level (Section 1) and classification system ... (Section 2) apply to sites producing carbon
steels, which contain less than 8% alloys by mass." A *stand-alone* DRI/HBI plant can hold
**Core Site Certification** (a site "may consist of stand-alone facilities for the production
of steelmaking raw materials"); ArcelorMittal Hamburg (Europe's only DRI-EAF site) holds Core
Site Certification. But there is **no Progress Level scale for merchant iron**. Merchant HBI
is handled instead as an *input* to whoever melts it (see (f) below).

High-alloy / stainless (>= 8% alloy) sites: "an adjusted emissions accounting methodology and
unique Decarbonisation Progress Levels will be developed in the future".

### (b) Scope 2 — counted, and MARKET-BASED IS ALLOWED

Decisive text. **Criterion 10.4.4, Energy-related indirect (Scope 2) GHG emissions:**

> a) Imported electricity:
> - GHG emissions are quantified in accordance with the requirements of **ISO 14064-1:2018
>   Annex E.2 Treatment of imported electricity, using the emission factor that best
>   characterises the pertinent grid, i.e. dedicated transmission line, local, regional or
>   national grid-average emission factor.**
> - Grid-average emission factors are from the emissions year being reported, if available,
>   or the most recent year, if not. Grid-average emission factors for imported consumed
>   electricity are based on the average consumption mix of the grid from which the
>   electricity is consumed.
> - **Determining energy-related indirect (Scope 2) GHG emissions may be based on the use of
>   renewable energy certificates, power purchase agreements, virtual power purchase
>   agreements, or green tariffs paid in relation to the site's sourcing of electricity, where
>   these meet the requirements of ISO 14064-1:2018 E.2.2 Additional information.**
> - Imported electricity that is used upstream of crude steel production at the site, and that
>   has been generated from process gases or waste energy while producing crude steel at the
>   site, is excluded ...

So: **default is location-based (grid-average), but a market-based claim via REC / PPA / vPPA /
green tariff is expressly permitted**, subject only to ISO 14064-1 E.2.2. There is **no
hourly-matching, deliverability or additionality test**. A PPA can make grid electricity count
as zero.

The counterweight is **disclosure, not prohibition**. Criterion 10.7.1.b requires public
disclosure of:

> (iii) Whether the determination includes the purchase of renewable energy certificates or
> similar mechanisms, such as power purchase agreements, virtual power purchase agreements,
> or green tariffs paid in relation to the sourcing of the site's electricity, and, if so, **a
> description of the source and quantity of certificates or agreements**

**Carbon offsets are banned** for all three scopes:
"(10.4.3.b) GHG offsets are not recognised for the purpose of determining the site's GHG
emissions intensity, in relation to its direct (Scope 1), energy indirect (Scope 2) or
upstream indirect (Scope 3) GHG emissions." RECs/PPAs are treated as *energy attribute*
instruments, categorically distinct from offsets.

Recognised measurement standards: "The GHG Protocol and EN 19694 (parts as applicable) ...
ISO 14404 (parts as applicable) for the measurement of CO2 emissions by steelmaking sites."

### (c) Thresholds — Decarbonisation Progress Levels (DPL 1-4)

Linear sliding scale `y <= b - m(x)`, `y` = t CO2e / t crude steel, `x` = scrap share of
metallic input (fraction).

| Level | b (0% scrap) | m (gradient) | value at 100% scrap |
|---|---|---|---|
| **DPL 1** (better than industry average, "Basic Threshold") | 2.80 | 2.30 | 0.50 |
| **DPL 2** | 2.00 | 1.65 | 0.35 |
| **DPL 3** | 1.20 | 1.00 | 0.20 |
| **DPL 4 ("near zero")** | 0.40 | 0.35 | 0.05 |

DPL 4 "has been aligned with the International Energy Agency (IEA)'s proposed threshold for
'near-zero emission production' of steel"; DPL 2/3 with IEA's proposed intermediate classes.
Worked example: a site at 20% scrap needs <= 2.80 - 2.30(0.20) = **2.34 t CO2e/t crude steel**
for DPL 1.

Progress Levels "will be reviewed on a five-yearly basis"; revisions apply after a 2-year
transition. Levels were specified excluding high-alloy sites, so are based on global
performance for carbon steel.

### (d) Scrap

Scrap share is measured as a share of **metallic input**, not of output:

```
Scrap share of metallic inputs (%)
  = SUM_s(f_met,s * Q_in,s) / [ SUM_p(f_met,p * Q_in,p) + SUM_s(f_met,s * Q_in,s) ]
```

- `Q_in` = tonnes of material input into steelmaking
- `f_met` = metallic fraction; absent primary data, **98% for scrap, 94% for DRI/HBI, 94% for
  pig iron** (as defined in the SBTi Steel Sector Guidance)
- `s` = secondary metallics: home scrap, manufacturing scrap, end-of-life scrap —
  **excluding internal scrap** (crude steel returned within the same BOF/EAF unit)
- `p` = primary metallics: **pig iron, DRI, HBI**, ferro-alloys, non-ferrous metals

Scrap carries **zero embodied burden** (cut-off): "The use of primary data is not applicable in
the case of scrap and post-consumer reclaimed material, for which the default embodied GHG
value of zero always applies." But scrap *collection and transport* IS counted: "In the case of
scrap and other recycled or reclaimed materials the GHG emissions associated with
transportation should be estimated from the commercial collection point to the ResponsibleSteel
certified site gate."

**No end-of-life recycling credit.** (JRC: none of the six frameworks it compared allow EoL
credits; IDDI estimates this is worth up to 14% of reported emissions.)

### (e) Transport / freight — INCLUDED, and load-bearing for imported HBI

Scope 3.4 upstream transport is inside the boundary:

> (10.4.5.c, d) It is the responsibility of the purchaser to ensure that an estimate for the GHG
> emissions associated with **transportation of the input material up to the point of delivery**
> has been provided in accordance with the point of delivery specified in the purchase contract
> (e.g. free on rail at mine gate, free on board, or including carriage, insurance and freight).
> **The purchaser is responsible for determining any additional estimated GHG emissions
> associated with further carriage of the material by the purchaser.**
> Estimates should consider the transportation distance, mass of material and the mode of
> transportation (road, rail, ship) and the related carrier type.

Outbound transport to the steel customer is **not** in the site metric (outside cradle-to-gate).

### (f) THE KEY MECHANISM FOR IMPORTED HBI — Criterion 10.4.5

An EU EAF charging imported HBI must count that HBI's full **cradle-to-gate Scope 1 + Scope 2 +
Scope 3** in its own Scope 3.1:

> a) Determining the site's upstream indirect (Scope 3) GHG emissions includes the direct
> (Scope 1), energy-related indirect (Scope 2), and upstream indirect (Scope 3) GHG emissions
> from 'cradle-to-gate' for the following input materials, if applicable:
> - **Ferrous-containing materials: cold iron, direct reduced iron (DRI), granulated pig iron
>   (GPI), hot briquetted iron (HBI), iron ore, pellets, scrap, sinter, steel slab**

Primary data is *mandatory where available*, with an anti-gaming rule:

> (10.4.5.b) ... when the steelmaker has received primary data from a supplier for the embodied
> GHG value for the supplied input material the steelmaker **must** use these data ... and **may
> not use the default embodied GHG value for the material even if the default value is lower**.
> If a steelmaker has primary data provided by some but not all suppliers, primary data must be
> used for the proportion of the material for which primary data is available ...

**ResponsibleSteel default embodied GHG values (Annex 5, Table A1)** — used when no supplier
primary data exists (all from "CRU methodology for ResponsibleSteel" unless noted):

| Input | t CO2e / t |
|---|---|
| **Hot briquetted iron (HBI)** | **1.219** |
| **DRI, gas-based** | **1.219** |
| **DRI, coal-based** | **2.623** |
| Cold iron, generic | 2.623 |
| Cold iron, charcoal based | 2.350 |
| Granulated pig iron (GPI) | 2.623 |
| Iron ore | 0.025 (2024.1 GaBi / Sphera) |
| Pellets | 0.235 |
| Sinter | 0.365 |
| **Scrap** | **0.000** |
| Steel slab, BOF | 2.460 (= RS DPL 1 at 15% scrap) |
| Steel slab, EAF | 0.620 (= RS DPL 1 at 95% scrap) |
| Non-ferrous metal & ferro-alloy replacement value | 2.186 |

Note the asymmetry that punishes an unverified importer: the **default for HBI is 1.219
t CO2e/t**, so an EU EAF importing HBI without supplier data is charged ~1.2 t CO2e per tonne
of HBI regardless of how the iron was actually made. A genuinely green Australian/Brazilian HBI
producer must supply verified primary data to beat that default — and, crucially, that primary
data *includes the HBI plant's own Scope 2*.

Also: "Steel (non-scrap) — If steel other than scrap is imported to the site as an input for
production of crude steel at the site, and if primary data for its upstream emissions is not
available, it is assigned a default upstream emission factor as for other ferrous input
materials as listed in Table A1. If steel is imported to the site for further downstream
processing, the upstream emissions associated with its production are not included when
determining the crude steel GHG emissions intensity for the site."

### (g) Status and adoption

- **International Standard v2.0** published 2022; **International Production Standard v2.1**
  (2024); **v2.1.1** current. *Fundamentals for GHG Emissions Accounting and Classification*
  v1.0 published **18 June 2025**.
- ~**90 ResponsibleSteel-certified sites** worldwide; certified sites account for the equivalent
  of **7.3% of the steel sector's CO2e emissions**.
- **Issued certificates directory** (https://www.responsiblesteel.org/issued-certificates-and-ongoing-audits):
  **47 issued certificates**, of which **45 are Core Site Certification** and only **2 are Certified
  Steel — both Big River Steel** (the EAF facility at 57% scrap that LESS cites as the first
  ResponsibleSteel steel certification). Regional split: Europe 23, South America 10, Asia 8, North
  America 3, Oceania 2, Middle East 1.
  - **Brazil: 10 sites** — ArcelorMittal Brasil (Monlevade, Pecém, Sul Fluminense, Piracicaba,
    Sabará/São Paulo, Tubarão, Vega) and Aperam Inox América do Sul.
  - **Australia: BlueScope ASP Manufacturing Cluster.**
  - **No Decarbonisation Progress Level designations appear in the public certificate listings.**

  This matters for the export case: **both Brazil and Australia already have ResponsibleSteel Core
  Site certificates in place**, so the certification infrastructure exists in both origin countries —
  but the graded, emissions-intensity part of the standard (Certified Steel + a Progress Level) has
  been awarded exactly twice worldwide, and never in either country.
- Two tiers: **Core Site Certification** (ESG, >300 requirements) and **Certified Steel** (adds
  Decarbonisation Progress Levels + responsible sourcing).
- **14 November 2025 (around COP30, Belém)**: agreements with **LESS** and with **CISA** (China
  Iron and Steel Association) to build *conversion tools* under ResponsibleSteel's **"Framework
  for Credible Interoperability"**, so a site can claim equivalency across schemes. The three
  organisations' members represent **~60% of world steel production**.
  https://www.responsiblesteel.org/news/landmark-agreements-link-majority-of-worlds-steel-production-under-global-and-regional-standards-for-low-emission-steel
- Caveat in the Fundamentals document itself: "ResponsibleSteel certification is not currently
  available in relation to conformity with the requirements of this document ... use of this
  methodology does not entitle the user to make any claims of compliance with the Production
  Standard, claims of equivalency or claims of Certified Steel."
- ResponsibleSteel Progress Report 2026:
  https://responsible-steel.shorthandstories.com/responsiblesteel-progress-report-2026/index.html

---

## 2. Low Emission Steel Standard (LESS) and the joint RS-LESS Steel Decarbonisation Scale

Primary sources:
- LESS Rulebook v1.1 (March 2025) — https://lowemissionsteelstandard.org/downloads
  (`Rulebook_LESS_v1.1.pdf`)
- LESS, *Technical Recommendations for a European Steel Classification System*, 23 July 2026
  (`20260723_LESS_Technical Recommendations for a European Steel Classification System.pdf`)
- LESS, *Public Statement Sliding Scale*, 2 December 2025
- LESS, *Public Statement ESPR Delegated Act*, 12 May 2026
- ResponsibleSteel & LESS, *The Steel Decarbonisation Scale: A briefing for policy makers in the
  EU designing a label for low-emission steel*, June 2025 —
  https://lowemissionsteelstandard.org/files/less/home/insights/RS-LESS_The%20Steel%20Decarbonisation%20Scale_final.pdf
- https://lowemissionsteelstandard.org/ ; https://lowemissionsteelstandard.org/certification-system
- https://www.thyssenkrupp-steel.com/en/company/sustainability/environmental-protection/less/less-low-emission-steel-standard.html

LESS was developed by **WV Stahl** (German Steel Federation) out of the BMWK stakeholder process
"Leitmärkte für klimafreundliche Grundstoffe" (Guidehouse / Fraunhofer ISI / Wuppertal Institute,
2023). **LESS aisbl** (Belgian non-profit) was established **October 2024** as system owner.
Launched April 2024; Rulebook v1.0 April 2024, v1.1 March 2025.

### (a) System boundary

**Cradle-to-hot-rolled-steel** — deliberately wider than IEA/ResponsibleSteel:

- "the present proposal takes up the **IEA's six-stage model** and the normative approach, but at
  the same time **extends the scope of the balance sheet to include steel refining, casting and
  hot rolling as the first further processing stage, including initial reheating (1st heat)**."
- "Only the initial pass in the hot rolling mill (first heat) is relevant in determining the
  emissions ... Further passes should not be taken into account."
- "Processing and finishing steps downstream from the hot rolling mill are not taken into
  account, including heat treatment or the production of cold drawn steel."
- Non-production facilities (admin buildings, research) excluded.
- Declared unit: **hot-rolled steel** — "The product as defined by the CS is hot-rolled steel."
- Scope 3 in: **3.1** (raw materials incl. scrap, ore, ferro-alloys, refractories, consumables),
  **3.3** (energy-related), **3.4** (transport of those materials to the production site),
  **3.8** (leased plant). Explicitly out: 3.2 capital goods, 3.5 waste disposal, 3.6 business
  travel, 3.7 commuting, and all downstream 3.9-3.15.
- Gases: **Scope 1 and Scope 2 are CO2 only**, "in line with the EU ETS". Scope 3 covers all GHGs
  (CO2, CH4, N2O, HFCs partial, SF6). A real divergence from ResponsibleSteel (full CO2e).
- **Process gases are NOT in the boundary** — the JRC flags this as the single largest source of
  non-comparability with the other frameworks, since process-gas combustion is "10-40% of total
  emissions" (Bonaplata et al., 2021).

Iron/DRI is not classified separately — the classified product is hot-rolled steel. But DRI-EAF
and **DRI-SAF** are both modelled as explicit reference routes: "The model was used to derive two
variants in relation to direct reduction plants (DR plants); DRI EAF (electric arc furnace) and
DRI SAF (submerged arc furnace), which here too ensures openness to [technology]".

### (b) Scope 2 — counted, CO2 only, and MARKET-BASED IS EXPRESSLY PERMITTED

**LESS Rulebook v1.1:**

> **If electricity is purchased by the operator, the determination of Scope 2 emissions requires
> the real emission factor from the actual electricity mix used by the plant operator to be
> applied.**

followed immediately by:

> Also utilizing **accounting-related allocations relating to the use of energy sources and the
> resulting emissions is expressly permitted**, provided that **guarantees of origin from
> officially recognised registers** are referenced. Currently and with significance for the CS,
> these are: **Guarantees of origin for green electricity purchased on an accounting-related
> basis and not physically delivered (German Environment Agency Register)**, accounting-related
> procurement of biogas (German Biogas Register, German Energy Agency). It may also be
> permissible in the course of reviews to include additional guarantees of origin as necessary
> in relation to systems currently under construction (such as hydrogen), provided that this is
> also technically expedient.

Two consequences. First, **market-based accounting is allowed and green electricity need not be
physically delivered** ("not physically delivered" is explicit). Second, the *registers named are
German*; "officially recognised registers" is a European construct, and LESS's declared
geographical area is the **EU** (JRC comparison table). A non-EU HBI producer has no named
register to point at under the current rulebook.

**LESS Rulebook v1.1, Annex VI "Standard Emission Factors" carries an explicit renewable-electricity
row** — the clearest statement anywhere in the survey that a market-based zero is intended:

| Standard factor | Unit | col. 1 | col. 2 | col. 3 |
|---|---|---|---|---|
| **Electricity** | kWh | **0.4280** | 0.0570 | 0.0285 |
| **Electricity (renewable)** | kWh | **0.0000** | 0.0262 | 0.0131 |
| Natural gas | GJ | 56.0000 | 11.3000 | 5.6500 |
| Oxygen (renewable) | m3 | — | 0.0002 | 0.0001 |
| Graphite electrodes | kg | 3.6337 | — | — (DEHSt electrode burnup) |
| COG gas | GJ | 44.4 | — | — (Reg. (EU) 601/2012) |
| BF gas | GJ | 260.0 | — | — (Reg. (EU) 601/2012) |
| OSBF offgas | GJ | 165.0 | — | — ("Estimation based on the expected gas composition / amortized by use in the DRI") |
| BOF gas | GJ | 182.0 | — | — (Reg. (EU) 601/2012) |
| Steelmaking slag | kg dry | 0.0060 | — | — |
| Granulated BF slag | kg dry | 0.0010 | 0.1000 | — |
| DRI | kg | *not publicly available* | — | — (Duarte et al., "Decarbonising the steelmaking industry", *Steel Times International*, 2021) |
| Hot metal, hydrogen, oxygen, ferro-alloys, anthracite | — | *not publicly available* | — | — (ISO 19694-2) |
| Aluminium, ferrosilicon, olivine | — | — | — | — (**Cut-off**) |

Column 1 reads as the Scope 1 / Scope 2 factor (56 kg CO2/GJ is the standard natural-gas combustion
factor; 0.428 kg CO2/kWh is a German grid factor); columns 2 and 3 are the Scope 3.3 upstream factor
and — *inferred*, since they are in an exact 2:1 ratio in every row — the same factor with the "~50%
reduction in Scope 3 upstream emissions" the rulebook applies when deriving the near-zero reference
plant.

The substantive point: **renewable electricity gets a Scope 2 factor of exactly 0.0000, but still
carries 0.0262 kg CO2e/kWh of upstream Scope 3.3 burden** — i.e. LESS zeroes the operational term but
does not pretend renewable generation is free of lifecycle emissions. Grid electricity is 0.4280
Scope 2 + 0.0570 Scope 3.3. The whole gap between them, ~0.40 kg CO2e/kWh, is what a Guarantee of
Origin buys.

Note also that **the DRI standard factor's value is "not publicly available"** — LESS publishes the
*source* (Duarte et al. 2021) but not the number, in the public version of the annex.

LESS is candid that the electricity factor drives everything:
> "The emission factor of the electricity used in the secondary route plays a decisive role in
> the classification system. This is also the largest influencing factor upon which producers can
> make a [difference]."

and that a scrap-EAF on average grid power reaches a level "but not level A. To reach this higher
level, additional efforts such as the use of hydrogen or biogenic [carbon are needed]."

Reference cases in the rulebook: integrated route ~26%/36%/38% Scope 1/2/3 split; scrap EAF "the
proportion of Scope 2 emissions at 42% and Scope 3 emissions at 38% is very high"; one reference
case notes "no electricity is purchased and the Scope 2 emissions are equal to zero."

**LESS's July 2026 recommendation to the European Commission** (Technical Recommendations, §2) —
this is what LESS wants the EU voluntary label to be:

> For producers within the European Union, Scope 1 emissions should be based on verified CO2
> emissions reported under the EU Emissions Trading System (EU ETS). **For producers outside the
> European Union, equivalent Carbon Border Adjustment Mechanism (CBAM) data should be used.**
>
> Scope 2 emissions should be derived from CO2 verified emission factors of the electricity
> supplied to the production site. **Market-based instruments, such as Power Purchase Agreements
> (PPAs) and Guarantees of Origin (GoOs), should be recognised for the determination of Scope 2
> emissions, provided that they are substantiated by officially recognised registries and robust
> verification systems. Appropriate eligibility rules should ensure that market-based instruments
> maintain the environmental integrity and credibility of the classification methodology.**

Significant for this project: LESS proposes to bolt a **Scope 2 term onto CBAM Scope 1 data for
non-EU producers** — i.e. to close exactly the gap CBAM leaves open.

Also from the same document, the full category list:

> - Scope 1: Direct process emissions at the production site.
> - Scope 2: Purchased electricity.
> - Scope 3.1: Raw materials (e.g. iron ore, scrap, alloying elements and refractories).
> - Scope 3.3: Fuels and energy carriers (e.g. coal, natural gas, biogas and hydrogen).
> - Scope 3.4: Transport of materials to the production site.
> - Scope 3.8: Operational emissions associated with upstream leased assets ...

with a **90% materiality threshold**: "at least 90% of total greenhouse gas emissions within the
defined system boundary should be accounted for."

And an upstream-decarbonisation hook directly relevant to green HBI:
> "Where independently certified product-specific upstream emission data are available, these
> should replace the corresponding default values. **This approach rewards upstream
> decarbonisation** while preserving the consistency and comparability of the classification
> methodology."

### (c) Thresholds

Derivation (Rulebook v1.1):
- Start from the IEA near-zero threshold at crude steel.
- Add a **surcharge** for the wider boundary. "The total surcharge is estimated by comparing the
  IEA approach for the Near-Zero class with values derived using a bottom-up approach for a
  reference plant with the latest technology using largely climate-neutral energy sources and a
  reduction of around 50% in Scope 3 upstream emissions. **The total impact calculated for QST is
  therefore 120 kg CO2e/t of rolled steel and 70 kg CO2e/t of rolled steel for BST.**"
- "Based on this target value, **multiplications are made by a factor for each of the levels A -
  D; all values above this are qualified as 'E'.**"

Six categories: **Near-Zero, A, B, C, D, E** — the "six-stage model" LESS attributes to the IEA.

Two product families, aligned to EU ETS free-allocation benchmark definitions:
- **QST — quality steel** (= EU ETS "EAF high alloy steel")
- **BST — reinforcing and structural steel** (= EU ETS "EAF carbon steel")

Near-Zero thresholds (JRC Table 4, confirming LESS):

| Family | 0% scrap | 100% scrap |
|---|---|---|
| **QST (quality steel)** | **< 0.52 t CO2e / t hot rolled steel** | **< 0.17 t CO2e / t hot rolled steel** |
| **BST (reinforcing & structural)** | **< 0.47 t CO2e / t hot rolled steel** | **< 0.12 t CO2e / t hot rolled steel** |

The joint RS-LESS briefing's Figure 4(ii) plots the QST scale; class boundaries read at 100% scrap as
**Near-Zero 170, A 340, B 510, C 680, D 850 kg CO2e/t hot rolled steel** — i.e. multiples 1× / 2× /
3× / 4× / 5× of the near-zero value — with E above D. Applying the same multipliers at 0% scrap
(near-zero 520) gives **A 1040, B 1560, C 2080, D 2600**.

**The exact table, from LESS Rulebook v1.1 Annex III "Threshold Value Calculation Aid"** (Quality
steel / Qualitätsstahl, kg CO2e per tonne of rolled steel — `Schwelle` = threshold, `Schrottquote` =
scrap share):

| Scrap share | Near Zero | A | B | C | D |
|---|---|---|---|---|---|
| 0% | **520** | **1 040** | **1 560** | **2 080** | **2 600** |
| 10% | 485 | 970 | 1 455 | 1 940 | 2 425 |
| 20% | 450 | 900 | 1 350 | 1 800 | 2 250 |
| 30% | 415 | 830 | 1 245 | 1 660 | 2 075 |
| 40% | 380 | 760 | 1 140 | 1 520 | 1 900 |
| 50% | 345 | 690 | 1 035 | 1 380 | 1 725 |
| 60% | 310 | 620 | 930 | 1 240 | 1 550 |
| 70% | 275 | 550 | 825 | 1 100 | 1 375 |
| 80% | 240 | 480 | 720 | 960 | 1 200 |
| 90% | 205 | 410 | 615 | 820 | 1 025 |
| 100% | **170** | **340** | **510** | **680** | **850** |

Anything above the D line is **E**. The annex states the line equation as
`E = m·(100 − x) + s`, where `E` = emissions in kg CO2e/t rolled steel, `m` = slope, `x` = scrap
share (%), `s` = emissions at 100% scrap; and shows the derivation as the **IEA approach (400 / 50)
plus the surcharge (+120 for QST, +70 for BST)** at both ends. So:

```
Near-Zero, QST:  E ≤ (400 + 120) − (350/100)·x   = 520 − 3.50·x     [x in %]
Near-Zero, BST:  E ≤ (400 +  70) − (350/100)·x   = 470 − 3.50·x
Class A = 2× Near-Zero,  B = 3×,  C = 4×,  D = 5×,  E = above D
```

**These multipliers are confirmed by the rulebook's own reference-plant worked cases (Table 3,
kg CO2e/t rolled steel):**

| Reference plant | Scope 1 | Scope 2 | Scope 3 | Sum | Class |
|---|---|---|---|---|---|
| Integrated steelworks | 1 904 | **0** (note 31) | 585 | **2 489** | **Stage E** |
| EAF QST (100% scrap, current German grid) | 207 | **294** | 304 | **805** | **Stage D** |
| DRI-EAF under ideal conditions | 16 | **0** | 529 | **544** | **Stage A** |
| EAF-QST under ideal conditions | 18 | **0** | 246 | **264** | **Stage A** |
| EAF BST (100% scrap) | — | 42% of total | 38% of total | **547** | **Stage D** |
| EAF BST under ideal conditions | — | — | — | **170** | **Stage A** |

Note 31: "The integrated steelworks covers its electricity requirements by converting the by-product
gases into electricity in its own power plant. Therefore, no electricity is purchased and the Scope 2
emissions are equal to zero."

**The most quotable sentence in the entire rulebook for this project:**

> "**The emission factor of the electricity used in the secondary route plays a decisive role in the
> classification system. This is also the largest influencing factor upon which producers can make a
> direct impact. If 50% renewable electricity is used instead of the current German electricity mix,
> the reference plant could already reach level C. With 100% green electricity level B is also
> achievable — but not level A. To reach this higher level, additional efforts such as the use of
> hydrogen or biogenic coal are necessary.**"

Arithmetic check: 805 − 294 = 511 ≈ the B boundary at 100% scrap (510). The whole distance from D to
B for a scrap EAF is *nothing but the electricity emission factor*. Under the LESS Scope 2 rule
(§(b) above), buying German GOs for green electricity "not physically delivered" moves the plant two
whole classes.

On the primary route:
> "Should direct reduction be based on climate-neutral hydrogen, then both DRI variants can be
> classified in classification level B without further changes."

and the ideal DRI-EAF case at 544 kg is Stage A, with the rulebook noting "The emissions in this
calculation mainly arise from the upstream chain of the pellets, hydrogen, dolomite lime and
aluminium used (→ Scope 3). Direct emissions are mainly caused by the carbon content of the pellets
and electrodes." — i.e. for a fully decarbonised H2-DRI-EAF, **Scope 3.1 is essentially the entire
remaining footprint**.

LESS also applies **alloying adjustment rules** so distinct limit values follow for every steel
grade, and (July 2026 recommendations) a **dimension-related material and energy adjustment**
normalising hot-rolling yield loss to a 5% reference. Worked example in that document: 800 kg
CO2e/t product → alloy adjustment to 630 → dimension adjustment to **621 kg CO2e/t** normalised
emissions, which is then compared to the threshold.

### (d) Scrap

```
scrap share = Amount_scrap / (Amount_scrap + Amount_pig_iron + Amount_DRI + Amount_HBI)
```

"It is irrelevant whether the materials mentioned are supplied solid or molten. Alloying agents,
including those containing Fe, are not to be used in determining the scrap share."

Contrast with ResponsibleSteel: LESS uses **raw tonnages**, ResponsibleSteel weights by **metallic
fraction** (98%/94%/94%). LESS's 2026 recommendation also wants Scrap Share to "consistently
include **pre-consumer scrap generated within the steel producer's own operations**" so vertical
integration does not distort the number — the *opposite* of ResponsibleSteel's exclusion of
internal scrap. These two definitions do not give the same answer for the same plant.

Scrap carries zero embodied burden (cut-off), but scrap collection/processing and its transport
ARE inside the boundary (JRC Table 2b).

### (e) Transport

Scope 3.4 in: "raw materials, auxiliary materials and operating materials used (Scope 3.1) and
**precursors including their transport to the plants (Scope 3.4)**; where applicable such
transports are to be included directly via the selected emission factors."

### (f) Imported HBI — LESS looks through, via a 95% coverage rule

> "**Intermediate products are only CS compliant if at least 95% of the total input quantity of
> pig iron, DRI, HBI and scrap is evaluated according to the CS** or if an evaluation is not
> required with respect to scrap. Should the proportion of intermediate products used that are
> not subject to CS be greater than 5%, the production is to be apportioned accordingly and only
> part of the production can be evaluated under CS."
>
> "Intermediate products from the plant operator's own or affiliated companies are to be
> determined and certified according to the CS."
>
> "The standard values defined here in the rulebook are to be used if operators are unable to
> provide the verifier with appropriate substantiation ... Alternatively, if at least 95% of an
> intermediate product's utilisation quantity is classified, **the highest emission value from
> these warranties may be used for the remaining utilisation quantity.**"

Also: "Hot-rolled steel as a product within the meaning of the CS only qualifies for
classification within this system if at least 95% of the crude steel quantities used for the
classified product are subject to the CS."

So an EU EAF melting imported HBI cannot get a LESS class unless the HBI itself is evaluated under
the classification system (or a conservative standard value is applied). The exporting HBI plant's
Scope 2 therefore does reach the EU product's class.

### (g) Credits

LESS grants credits (broader than ResponsibleSteel):
- **Granulated slag / comparable by-products sold as clinker substitute: 0.100 t CO2e/t credit.**
- Export of domed (process) gases for electricity/heat consumed outside the plant. "Electricity is
  credited using the latest ... emission factor for the **EU Member States national electricity mix
  officially published by the EU Member State's national authority (e.g. in Germany Federal
  Environment Agency UBA). For non-EU countries, verified data are required and must be accepted by
  EU as official data.**"
- Energy supplied to third parties; heat credited on the EU ETS heat benchmark.
- Carbon **offsetting is not permitted**.

### (h) Status and adoption

- Rulebook v1.0 April 2024; v1.1 March 2025 (among other changes: "Credits specification for
  electricity production for **EU- and non-EU-based companies**").
- LESS aisbl founded October 2024; certification system live; certification-body and operator
  guidelines published; FAQ on verification May 2025.
- **thyssenkrupp Steel is certified at LESS Label Level D** — verified **September 2025 by TÜV Nord**,
  for hot-rolled steel from the **blast-furnace route with 13% scrap**. Its published roadmap is
  D → C → B → A as H2-DRI comes on. (thyssenkrupp's own materials describe category A as the
  hydrogen-based endpoint and category C as achievable in the transitional phase — those are targets,
  not the current certificate.) Criticism from BDE (waste industry) and Bellona that the sliding scale
  lets high-scrap EAF steel be out-classed by low-scrap primary steel.
- **Certificates issued as of the certification-system page (https://lowemissionsteelstandard.org/certification-system) — eight, all German or Franco-German:**

  | Company | Certification body | Certified | Expires |
  |---|---|---|---|
  | Georgsmarienhütte GmbH | DNV Business Assurance B.V. | 12.09.2025 | 11.09.2026 |
  | thyssenkrupp Steel Europe AG | TÜV NORD CERT GmbH | 12.09.2025 | 11.09.2026 |
  | Salzgitter Flachstahl GmbH | TÜV NORD CERT GmbH | 12.09.2025 | 11.09.2026 |
  | Peiner Träger GmbH | TÜV NORD CERT GmbH | 12.09.2025 | 11.09.2026 |
  | Salzgitter Mannesmann Grobblech GmbH | TÜV NORD CERT GmbH | 12.09.2025 | 11.09.2026 |
  | Ilsenburger Grobblech GmbH | TÜV NORD CERT GmbH | 12.09.2025 | 11.09.2026 |
  | Saarstahl AG | TÜV NORD CERT GmbH | 27.10.2025 | 26.10.2026 |
  | Saarstahl Rail SAS | TÜV NORD CERT GmbH | 27.10.2025 | 26.10.2026 |

  Accredited certification bodies: DNV Business Assurance, GUT, proTerra, TÜV NORD CERT, TÜV SÜD.
  Two-stage verification: initial audit plus mandatory follow-up within 12 months. Individual class
  levels are not published per site. **Note the geography: LESS adoption is, so far, entirely German
  plus one French subsidiary** — which is worth weighing when judging how much an Australian or
  Brazilian HBI exporter could rely on it.
- Policy engagement is the main current activity: Public Statement on the Sliding Scale (Dec 2025),
  Public Statement on the ESPR Delegated Act (May 2026), **Technical Recommendations for a
  European Steel Classification System (23 July 2026)** aimed squarely at the IDAA voluntary label.
- Interoperability agreement with ResponsibleSteel signed Nov 2025.
- Secretary General is now **Dr. Carmen Ostwald** (LESS aisbl, c/o Wirtschaftsvereinigung Stahl,
  Rue Marie de Bourgogne 58, 1000 Brussels); Dr. Martin Theuringer signed the June 2025 joint
  briefing.

### (i) LESS's two public statements on EU policy — the live EU vehicles named

**Public Statement, 2 December 2025 — "The EU Label for Steel Should Build on the Sliding Scale
Approach":**

> "As announced in the Clean Industrial Deal, the European Commission is expected to put forward an
> **EU label for steel on 10 December [2025]** as part of the **Industrial Accelerator Act**."

> "By contrast, a purely footprint-based approach fails to guide the industry toward near-zero
> emissions… **a footprint-based approach would merely outsource the decarbonisation of primary
> production to third countries** and make the EU steel industry more vulnerable and less resilient."

> "**The EU label should include upstream scope 3 emissions (raw materials) from the start**, thereby
> ensuring that 'near-zero steel' is truly decarbonised across the entire value chain. This prevents
> carbon leakage and ensures that decarbonisation efforts are not undermined by shifting emissions
> from one part of the value chain to the other."

> "For secondary producers using the Electric Arc Furnace (EAF) route, it incentivises measures such
> as **increasing the use of green electricity** and decarbonising hot-forming processes. These
> improvements are less capital-intensive than decarbonising ore-based primary steelmaking, allowing
> EAF producers to achieve better performance classes faster. This is clearly demonstrated by the
> first certificates awarded under LESS, where secondary steel routes achieved substantially better
> classifications than primary steel routes, and by **the first steel certification by
> ResponsibleSteel, which was awarded to an EAF facility with 57% scrap inputs**."

Footnote 6 gives LESS's working definitions: "primary producers (<25% scrap), secondary producers
(>70% scrap), and producers which aren't currently classified as either (25-70% scrap)."

**Public Statement, 12 May 2026 — "A Credible and Operable Framework for Low-Carbon Steel in the ESPR
Delegated Act".** This reveals that the **ESPR (Ecodesign for Sustainable Products Regulation)
delegated act for iron and steel** is a *second, parallel* EU vehicle defining performance classes,
that the **JRC has a current proposal on the table**, and that it is **footprint-based with static
thresholds, not a sliding scale**:

> "While the current **JRC proposal** may serve transparency and comparability – valid goals of the
> ESPR – **it lacks a convincing logic for creating demand for verifiable low-carbon volumes**. Lead
> markets require more than excluding worst performers or clustering today's production into relative
> classes."

> "The JRC's current thresholds fail to sufficiently distinguish between conventional and low-carbon
> production routes, risking that **most of today's European production qualifies for lead market
> incentives** (e.g. public procurement), while significant decarbonisation investments do not improve
> the classification. For example:
> - **Scrap-based EAF: A switch to 100% renewable electricity yields no classification upgrade if
>   already in Class A.**
> - Wire rod: Natural gas-based DRI-EAF – despite billions in investments – is placed in the same
>   class as the best-performing BF-BOF routes.
> - **Hot-rolled coil: Natural gas-based DRI-EAF – despite significant emissions reductions and
>   investments – shares Class A with hydrogen-based DRI-EAF, removing incentives to transition from
>   gas to hydrogen.**"

That first bullet is directly on point for this project: **under the JRC's current ESPR proposal, a
switch to 100% renewable electricity can produce no classification change at all.**

LESS's asks: (1) "The top performance class must reflect near-zero steel production… While initially
unoccupied, it sets a clear long-term target"; (2) "**A CO2e intensity- and recycled content-based
sliding scale is the only fair way to incentivise all decarbonisation pathways – scrap-based,
hydrogen-based or renewable-powered – equally**"; (3) "Focus on key production stages (e.g. hot-rolled
steel, including **scope 1–3 emissions up to hot rolling**), covering around 90% of GHG emissions" and
"**Build on existing certification systems (e.g. LESS, ResponsibleSteel)**"; (4) "The first large-scale
European DRI plants will **start in 2027**… finalise and implement the ESPR Delegated Act, ensuring a
functioning certification system and regulation for **lead markets by 2027**."

### (j) The joint RS-LESS "Steel Decarbonisation Scale" (June 2025)

An argument for the **sliding scale** and against a pure carbon-footprint label:

> "there is currently enough scrap to satisfy just **32%** of the world's demand for new steel
> (IEA, 2021)" while the recycling rate is already ~85%; "The maximum possible post-consumer
> recycled content for steel produced in 2024 is **41%**"; IEA projects scrap meets only **46%**
> of 2050 demand even in its sustainable development scenario.
>
> "When the supply of any material is limited and the available supply is fully utilised — as is
> the case for ferrous scrap — use of that material cannot be increased, **it can only be
> redistributed** ... **Specifying steel based only on its carbon footprint is, therefore, a highly
> ineffective mechanism for driving reductions in overall GHG emissions.**"
>
> "The same flaw applies at the level of projects, customers, steelmaking sites, companies or
> national or regional boundaries: **redistributing scrap redistributes rather than reduces GHG
> emissions.**"

Signed by Dr. Martin Theuringer (Secretary General, LESS aisbl) and Annie Heaton (CEO,
ResponsibleSteel). Recommends adopting the scale for (i) the IDAA voluntary carbon-intensity
label, (ii) EU Lead Markets, (iii) other instruments.

The two scales are **not numerically identical** — RS is 4 Progress Levels on t CO2e/t *crude
steel* with full CO2e; LESS is Near-Zero + A-E on kg CO2e/t *hot rolled steel* with CO2-only
Scope 1/2. Same *shape*, hence the Nov 2025 conversion-tool agreement.

---

## 3. The IEA's near-zero emission steel definition

Source hierarchy:

| Document | Date | Role |
|---|---|---|
| IEA, *Achieving Net Zero Heavy Industry Sectors in G7 Members* — https://iea.blob.core.windows.net/assets/c4d96342-f626-4aea-8dac-df1d1e567135/AchievingNetZeroHeavyIndustrySectorsinG7Members.pdf | May 2022 | **The definition of record** (Ch. 3 + Technical Annex) |
| IEA, *Emissions Measurement and Data Collection for a Net Zero Steel Industry* — https://iea.blob.core.windows.net/assets/8f6568aa-1dd8-4578-bc61-24ceba4a07dd/EmissionsMeasurementandDataCollectionforaNetZeroSteelIndustry.pdf | Apr 2023 | Compares worldsteel / ISO / ResponsibleSteel boundaries against the definition |
| IEA, *Definitions for Near-Zero and Low-Emissions Steel and Cement…* — https://www.iea.org/reports/definitions-for-near-zero-and-low-emissions-steel-and-cement-and-underlying-emissions-measurement-methodologies | **8 Nov 2024** | Latest substantive statement; adoption status; prepared for IEA WPID and the **Climate Club** |
| Breakthrough Agenda Report 2025 / 2026 | 2025/26 | Progress tracking only; **no numeric change** |

**Nothing numeric has moved since 2022.** The 2024 paper restates 400/50 verbatim.

### (a) System boundary

Downstream cut at **crude steel, including casting, excluding rolling**:
> "The boundary of the downstream end of the supply chain is set at crude steel production —
> including casting but excluding any further semi-finishing and finishing processes because [of]
> the heterogeneity of processes at facilities producing different products." (2022, p.105)

Functional unit: **1 tonne of crude steel**.

Upstream: iron ore mining and agglomeration are IN; scrap is OUT:
> "The upstream end of the supply chain boundary encompasses the supply and processing of the main
> raw material input to steelmaking: iron ore. Mining (including extraction, transportation and
> beneficiation) and agglomeration processes are both included within the scope. **The sorting and
> transportation of steel scrap is not included, due to data constraints**, nor are the production
> processes for other material inputs to the steelmaking process (e.g. producing refractory linings
> for furnaces, ferroalloy production). The supply of limestone (to produce lime fluxes) is included
> within the boundary…"

Explicit source list (2022, pp.106-107):
- *Direct*: fossil fuel use in (i) iron ore agglomeration, (ii) ironmaking (incl. blast furnaces,
  **DRI furnaces**, innovative units, CCUS equipment), (iii) steelmaking; (iv) producing reduction
  agents (coke ovens, **on-site H2 production**); (v) lime fluxes and electrodes (process CO2);
  (vi) off-gases.
- *Indirect*: (vii) **imported electricity, heat and hydrogen**; (viii) fossil fuel supply
  (extraction, processing, transport — incl. CH4 and flaring); (ix) raw materials supply
  (extraction, beneficiation, **transportation** of iron ore and limestone).
- *Excluded*: direct CH4 and N2O; alloying elements; refractory linings; electrode *production*;
  de-oxidisers; and the transport of all of those.

**The IEA deliberately refuses to express this in GHG-Protocol scopes** — the single most important
structural point for this project:

> "While helpful for assigning varying levels of responsibility for emissions at the site or company
> level — and in widespread usage — we do not use the terminology of the [scopes]. We have chosen
> instead to be explicit about the sources themselves. **To take emissions associated with pellet
> production as an example, these could be Scope 1 emissions at one site, and Scope 3 at another
> site with a different process arrangement (e.g. using purchased pellets, produced off-site), with
> both sites having the same overall emissions intensity of steel production.**" (2022, p.108)

So **purchased pellets, sinter, coke, DRI/HBI and pig iron all carry their upstream emissions into
the number.** The boundary is drawn around *processes*, not around the fence line.

**Iron/DRI is NOT separately defined.** There is no IEA near-zero *iron* or *DRI/HBI* threshold. The
2023 report flags this as a live gap:
> "As the lines between technology arrangements may change in the future (e.g. DRI-melter-BOF,
> smelting reduction-BOF) and **trade in iron could become more common**, guidance that is either
> more detailed or more generalised will be required." (2023, p.35)

### (b) Scope 2 — counted, but the location/market question is EXPLICITLY UNRESOLVED at IEA level

Purchased electricity is definitely in — "Imported electricity, heat and hydrogen (indirect
energy-related CO2 emissions)" is a named boundary category, and 400/50 are set on a **direct +
indirect** basis.

But the IEA does not specify location-based vs market-based. The 2024 paper lists it as open:
> "Another key accounting issue that requires further analysis and discussion is **suitable
> methodologies for accounting for electricity emissions from the grid, including the level of
> regional granularity needed and use of renewable electricity credits.** The World Steel Association
> is currently undertaking analysis on this issue as part of a steel emissions methodology mapping
> exercise in support of the Steel Standards Principles discussions. **Further discussion on this
> topic is likely to be needed.**" (2024, p.30)

No grid factor is prescribed. Global averages appear only for the illustrative reference routes:
> "The global average CO2 intensity of electricity generation declines to around **140 gCO2/kWh in
> 2030** and drops below zero around 2040 in the Net Zero Emissions by 2050 Scenario, relative to a
> value of around **440 gCO2/kWh for 2020**." … "**These indirect emissions will depend on the actual
> source of the electricity used on a given site, in a given region, over a given time period.**"

And: **negative emissions are not passed through** — "The minimum emissions intensity of electricity,
hydrogen or any other vector for allocating indirect emissions, is zero."

**The direct vs direct+indirect sub-threshold trap.** The IEA splits 400/50 into "Near zero **direct**
emissions" and "Near zero **direct + indirect** emissions" at the *same numeric value*:
> "We split the threshold for steel production into two further sub-thresholds… **These two
> sub-thresholds have the same value**, as the categories of indirect emissions included must trend
> toward near zero in the long term." (2022, p.109)

For the direct-only sub-threshold, imported electricity, fossil fuel supply and raw material supply
are marked **N/A**. **So a plant can claim "near zero direct" at 400 kg while its grid electricity is
coal.** An H2-DRI-EAF plant on a dirty grid can pass "near zero direct" and fail "near zero direct +
indirect". The IEA notes the H2-DRI-EAF route "take[s] longer to reach the near zero emission
threshold with the global average parameters" precisely because its reductions are indirect.

### (c) Thresholds

Near-zero (2022, pp.108-110; restated 2024 p.17):
> "For steel, progressive according to the scrap share of metallic inputs, falling between the
> following:
> - **100% iron: 400 kg carbon dioxide equivalent (CO2-eq) per tonne of crude steel**
> - **100% scrap: 50 kg CO2-eq per tonne of crude steel.**"

Interpolation formula (Technical Annex, p.133) — strictly linear:

```
E_nz,s(s) = 400 − 350·s        [kg CO2e / t crude steel]
```

with `s` = scrap share of metallics input, 0 to 1. "In the absence of information on the share of
scrap use in a plant, country or region, **the default value of zero scrap is to be used.**"

IEA reference values (BAT energy performance, kg CO2e/t crude steel, 0% scrap unless noted):

| Emissions source | PCI BF-BOF | NG DRI-EAF | Scrap EAF (100% scrap) |
|---|---|---|---|
| Fossil fuel use in iron ore agglomeration | 235 | 40 | 0 |
| Producing reduction agents | 110 | 700 (combined with ironmaking) | 0 |
| Fossil fuel use in ironmaking | 590 | ↑ | 0 |
| Fossil fuel use in steelmaking | 0 | 25 | 30 |
| Lime fluxes and electrodes | 70 | 50 | 25 |
| Off-gases | 1 320 | 0 | 0 |
| **Imported electricity, heat, hydrogen** | **105** | **375** | **220** |
| Fossil fuel supply | 435 | 210 | 5 |
| Raw material supply | 80 | 80 | <5 |
| **Total** | **2 945** | **1 485** | **285** |

**Low-emissions (intermediate) thresholds** are set at **6× the near-zero line**, subdivided into
five bands A-E (2022 p.127; Technical Annex Table A.1):

| Band range | value at 0% scrap (kg CO2e/t) |
|---|---|
| A to E | **2 400** |
| A to D | 2 000 |
| A to C | 1 600 |
| A to B | 1 200 |
| A to A | 800 |

Full A-E line: `2400 − 2100·s`.
> "The maximum emissions intensity allowed to qualify as low emission production is set at **six
> times the near zero emission threshold**… places the low emission production threshold at around
> 10-20% below the emissions intensities of the dominant conventional process routes."

**Low-emissions is fractional, not binary.** Output is credited pro rata between the two lines:

```
P_l,s(b,s) = [E_l,s(b) − E_a,s] / [E_l,s(b) − E_nz,s(s)]   if E_a,s < E_l,s(b), else 0
```

> "Thus, the low emission production is progressively recognised, whereas the near zero emission
> threshold is binary."

### (d) Scrap

`s` = **scrap share of total metallic inputs** (not of output tonnage). The IEA offers **30% scrap**
as an optional primary/secondary dividing line for labelling only: "We propose 30% scrap use as the
cut-off below which primary near zero emission production could be explicitly recognised. The
threshold values we propose would remain the same, regardless of where this cut-off is made."
Typical primary production is 25-30% scrap; the 100%-iron endpoint is a construct — "the value for
hypothetical production with 100% iron is given here to enable calculation of the threshold for any
share of iron and scrap."

**The IEA does not define scrap precisely** (home vs manufacturing vs EOL, metallic fraction);
implementers had to (see ResponsibleSteel §1(d)). **DRI and HBI count as PRIMARY metallics, not
scrap** — so importing HBI does not raise your scrap share and does not relax your threshold.

### (e) Transport / freight — partially in

IN: transportation of **iron ore** and **limestone** (part of "raw materials supply"), and
transportation of **fossil fuels** (part of "fossil fuel supply").
OUT: "The sorting and transportation of steel scrap is not included." Also out: transport of
alloying elements, refractories, electrodes, de-oxidisers.
**Not addressed at all**: outbound transport of the crude steel product, and **transport of purchased
DRI/HBI or pig iron** — the IEA boundary is silent on this. (ResponsibleSteel and LESS fill that gap;
the IEA definition itself does not.)

FMC tightened the ore/lime side: "Transport emissions of iron ore and lime products include **all
emissions regardless of intermediary stops** between mining and steel plant."

### (f) Adoption

Applying the threshold **directly at the crude steel boundary** (2024 paper, pp.18-19):
- **ResponsibleSteel** IPS V2.1 Principle 10 — DPL 4 "aligned with the IEA's proposed threshold"
- **IDDI** Secretariat, Green Public Procurement Pledge ("noting that discussions are ongoing as to
  application by individual IDDI member countries")
- **CISA** Low Carbon Emission Steel Evaluation Method, led by Baowu, launched **October 2024**
- **First Movers Coalition** Steel Commitment
- **SteelZero** (Climate Group)

Wider boundary with adjusted numbers: **LESS** ("The threshold is adjusted slightly upwards to
account for the wider scope"); **GSCC** ("by 2045, the emission intensity required is 380 kg to 400 kg
CO2-eq per tonne of hot-rolled steel (for long and flat steel respectively), and by 2050, 120 kg
CO2-eq per tonne").

Governmental status — **recognised, not enacted**:
- **G7 2022 Communiqué**: the IEA definitions are "a robust starting point for a common understanding
  of ambitious general definitions."
- **G7 2023** Industrial Decarbonisation Agenda annex recognised the IEA **Net Zero Measurement
  Principles**.
- **COP28**: ~40 organisations endorsed them via the **Steel Standards Principles** (WTO Secretariat
  + worldsteel, 60+ endorsers).
- **Climate Club** Members Statement at **COP29** "affirmed emerging common understandings on
  definitions and convergence… on thresholds."
- **Breakthrough Agenda Report 2026 verdict**: "this progress **has not yet translated into
  formalised approaches within international agreements**, such as trade frameworks, nor been taken
  up through use of internationally interoperable approaches in domestic policy measures."

Example claim: SSAB announced in **September 2025** the "world's first near-zero CO2e steel to meet
the IEA near-zero steel and FMC thresholds" at its Montpelier, Iowa scrap EAF using "recycled scrap
metal, **fossil-free electricity**, biocoal and renewable natural gas" — i.e. a contractual-electricity
claim. https://www.ssab.com/en-us/news/2025/09/ssab-achieves-iea-threshold-for-nearzero-co2e-emissions-steel

---

## 4. worldsteel CO2 data collection, ISO 14404 series, ISO 20915

### (a) worldsteel CO2 Data Collection methodology

Sources: worldsteel *CO2 Emissions Data Collection User Guide* v5.1 (publicly retrievable copy:
https://www.jisf.or.jp/news/topics/docs/IISIDataCollectionUserGuideA5_VERSION5.1.pdf — the current
V11, May 2024, has been taken off worldsteel's site); IEA 2022 and IEA 2023 characterisations.

**Boundary — wider downstream than IEA, narrower upstream:**
> "The supply-chain boundary extends **from iron ore agglomeration through to finished steel
> products**. It does not include the emissions associated with raw materials extraction, sorting and
> transportation, nor does it incorporate upstream emissions from fossil fuel supply. The emissions
> boundary covers all direct CO2 emissions from the iron and steel sector, along with indirect
> emissions from electricity generation and the production and use of lime fluxes. **Emissions
> credits are applied when electricity and other energy carriers (like off-gases) are exported for use
> off-site.**" (IEA 2022, p.93)

**CO2 only** — no CH4, no N2O. **"Indirect emissions from transport of raw materials are not
included."**

Scopes are explicit (v5.1 Appendix 2/3):
> "Calculations incorporate **Scope 1, Scope 2 and Scope 3 emissions**, according to the GHG protocol.
> **CO2 emissions = Direct + Indirect − Credit**;
> **CO2 intensity = CO2 emissions (tonne) / crude steel (tonne)**"
> "**Scope 3** emissions… **A Scope 3 charge is applied to exported BF and BOF gas to correct their
> direct emissions credit**, taking their value in use into account…"

**Purchased electricity: A FIXED WORLD-AVERAGE FACTOR — neither location- nor market-based.** This is
the single most distinctive treatment in the whole survey:
> "Electricity — **World average value based on IEA 2006**"
> "Both energy equiv. value and emission factor of electricity are **world average values**…
> e.g. CO2 Emission factor of Electricity = Energy equiv. value × conversion factor =
> **9.8 GJ/MWh × 0.0514 tCO2/GJ = 0.504 tCO2/MWh**"

The rationale is stated in the JISF ISO 14404 User Guide
(https://www.jisf.or.jp/en/activity/climate/iso14404/documents/ISO14404UserGuide.pdf):
> "ISO 14404 applies a conversion factor that is equivalent to **world average electricity** since CO2
> emission factors of electricity depend on power supply composition of the area, **which is not
> directly related to energy saving activities of the steel plant**."

**A PPA or GO does nothing under worldsteel / ISO 14404 as written** — and neither does a genuinely
clean grid. It is a **benchmarking tool for process efficiency**, not a carbon-footprint tool. (The
guide does add: "Users are allowed to apply their own conversion factors if they are credible.")

**worldsteel v5.1 upstream factors** (t CO2/unit) — directly usable in a DRI/HBI export model:

| Item | Energy equiv. (GJ/unit) | Upstream EF (t CO2/unit) |
|---|---|---|
| Pellets (t) | 2.100 | **0.137** |
| **Gas-based DRI (t)** | 14.100 | **0.780** |
| **Coal-based DRI (t)** | 17.900 | **1.210** |
| Pig iron (t) | 20.900 | **1.855** |
| Coke (dry t) | 4.000 | 0.224 |
| Burnt lime (t) | 4.500 | 0.950 |
| **Electricity (MWh)** | 9.800 | **0.504** |
| Steam (t) | 3.800 | 0.195 |
| Oxygen (k.m3N) | 6.900 | 0.355 |
| EAF/BOF electrodes (t) | – | 0.650 |

By-product gas direct/upstream EFs (t CO2/k.m3N): COG 0.836 / 0.977; BFG 0.891 / 0.170;
BOFG 1.512 / 0.432. Applied as `Scope 1 = Direct × (Purchased − Sold)`,
`Scope 3 = (Upstream − Direct) × (Purchased − Sold)`. Slag credits: BF slag to cement 0.550,
BOF slag to cement 0.300, CO2 to external 1.000 (v5.1 notes these "are not finalised yet").

**Scrap**: counted as an input with **zero embodied burden**; "Purchased carbon steel scrap: total
external procurement of scrap (pre- and post-consumer scrap, **excluding home scrap**)". Scrap
sorting/transport out of boundary. **No sliding-scale threshold** — worldsteel publishes averages
**by process route** instead.

**Coverage**: IEA 2023 — "worldsteel's current database (2022) includes data from more than 220 sites.
These represent approximately 485 Mt of steel production, or 25% of global production." Data is
confidential; only global route averages are published.

**Sustainability Charter**: IEA 2022 p.95 — "the worldsteel Sustainability Charter, signed by 39 of
its members as of 2022, embodies nine principles… The first of these principles is climate action,
with one of the criteria being that signatories **must submit CO2 or energy consumption data to
worldsteel or national governments**." That is the whole climate obligation — a *reporting* duty, not
a performance threshold.

### (b) ISO 14404 series

| Part | Title | Edition |
|---|---|---|
| **14404-1** | …Part 1: Steel plant with blast furnace | 1st ed. 2013 → **2nd ed. 2024-09-16** |
| **14404-2** | …Part 2: Steel plant with electric arc furnace (EAF) | 1st ed. 2013 → **2nd ed. 2024-09-16** (https://www.iso.org/standard/88111.html) |
| **14404-3** | …Part 3: Steel plant with EAF and coal-based or gas-based direct reduction iron (DRI) facility | 1st ed. 2017 → **2nd ed. 2024-09-16** (https://www.iso.org/standard/85990.html) |
| **14404-4** | …Part 4: **Guidance for using the ISO 14404 series** | 2020-12-21 (https://www.iso.org/standard/77622.html) |

**On 14404-4**: it is *not* a standalone electricity-carbon-intensity standard — it is the umbrella /
guidance part. **But one of its four cross-cutting additions is exactly the electricity factor
question.** From the ISO 14404-4:2020 preview
(https://cdn.standards.iteh.ai/samples/77622/8771f743065640c7976465afd639c323/ISO-14404-4-2020.pdf):

> **Scope:** "This document provides guidance for calculating the CO2 intensity at steel plants with
> all types of process routes, by defining the boundary, CO2 emission factors and the intermediate
> products for which upstream emissions are considered… includes the **Universal Calculation Sheet**…
> i. Steel plants with different process routes from ISO 14404-1, -2 and -3 (7.2.1)
> ii. Steel plants with more than one process route (7.2.2)
> iii. **Steel plants purchasing pig iron from the outside** (7.2.3)
> iv. **Steel plants and rerollers purchasing part or all of crude steel from outside** (7.2.4)
> Moreover… additional guidance… for the following topics:
> a) Evaluation of exported slags; b) Evaluation of by-product gas; c) Evaluation of stock;
> **d) Selection of calorific values and emission factors for electricity and fuel**"

Clause structure confirms it: "9.2.2 Explanation of emission factors based on **world average
electricity equivalent**"; "9.3 Selection of calorific values and emission factors for electricity and
fuel"; "Annex D (informative) Example of template for using **different emission factors**".

**Verdict on market-based**: 14404-4 provides *guidance on selection* and an informative Annex D
template for different factors; the **default remains the world-average electricity-equivalent**. The
body text of clause 9.3 could not be obtained (only the free preview, which stops after clause 3), so
whether it explicitly authorises contractual/market-based factors is **unconfirmed** — and the JISF
rationale above (grid mix deliberately excluded because it is "not directly related to energy saving
activities of the steel plant") argues strongly that it does not.

**Relationship to worldsteel:**
> "The ISO 14404 series is based on 'CO2 Emissions Data Collection User Guide' established by the
> World Steel Association… while worldsteel method applies common boundary and CO2 emission factors
> to all steelworks regardless of their process routes, **the ISO 14404 series defines the boundary,
> CO2 emission factors and intermediate products for which upstream emissions are considered for each
> of the process routes**, such as BF-BOF (14404-1), Scrap-EAF (14404-2) and DRI-EAF (14404-3)."

**Purchased DRI/HBI is explicitly handled** via "outsourced steel production activities":
> "Intermediate products with possibilities of considering upstream emissions include the following:
> — **Electricity / steam**;
> — Substances produced in the basic activities existing in the target process route (e.g. purchased
>   coke used in the BF-BOF route);
> — Substances that substitute the iron source of the process route even if they do not exist in the
>   target process route (e.g. **purchased DRI used in the BF-BOF route**)."

Emission categories (3.1.2-3.1.4): *direct* — inside the boundary; *upstream* — "from imported
material related to outsourced steel production activities **outside the site boundary** and from
**imported electricity and steam** into the site boundary"; *credit* — "corresponds to exported
material and electricity or steam".

Boundary vs IEA (IEA 2023, p.19): direct = worldsteel list plus rolling; indirect = "electricity, heat
and hydrogen; indirect raw materials manufacture. **Indirect fossil fuel and raw materials supply are
not included.**" **CO2 only.**

**IEA's verdict:**
> "The worldsteel CO2 methodology and the ISO 14404 series **were not initially conceived with the net
> zero transition in mind**. They are designed to be route-specific… They also exclude from their
> boundary of consideration several significant categories of emissions (e.g. from fossil fuel
> extraction, mining and transportation of raw materials)."

And on credits: both use "credit… based on **global average reference values**" — versus
ResponsibleSteel, whose "credits decrease over time with progress towards net zero."

### (c) ISO 20915:2018 — LCI calculation methodology for steel products

https://www.iso.org/standard/69297.html (ISO/TC 17), developed from the worldsteel LCI methodology
(2017): https://worldsteel.org/wp-content/uploads/Life-cycle-inventory-methodology-report.pdf

**Boundary — cradle-to-gate, with optional cradle-to-gate-with-recycling / cradle-to-grave.** Four
permitted scopes. Per IEA 2023, emissions scope = direct (as worldsteel) plus off-gases; indirect =
electricity, heat, hydrogen; **indirect fossil fuel and raw material supply**; **waste treatment**.
**All GHGs** (CO2, CO, CH4, N2O). **The widest boundary of the five compared.**

**Co-product allocation = system expansion**, applied to COG, coke, benzene/tar/toluene/xylene/sulphur,
BF gas, BF slag, BOF gas, BOF slag, EAF slag:
> "System expansion is cited in section 4.3.4 of ISO 14044:2006 as one of the preferred methods to use
> since it 'avoids' allocation… **credits are given for the production (net output) of process gases
> and slags** (that are used outside the product boundary) because their production replaces the
> alternative production of similar functional products."

**Scrap credit — the closed material loop / "value of scrap" method.** The central equation:

```
LCI for 1 kg of steel product including recycling = X − (RR − S)·Y·(X_pr − X_re)
```

> "X is the cradle-to-gate LCI of the steel product.
> **RR** is the end-of-life recycling rate of the steel product.
> **S** is the scrap input to the steelmaking process — this is the net scrap consumed in the
> steelmaking process and **does not include internally generated scrap**.
> **Y** is the process yield of the EAF…
> **X_pr** is the LCI for 100% primary metal production. **This is a theoretical value for steel slab
> made in the BOF route, assuming 0% scrap input.**
> **X_re** is the LCI for 100% secondary metal production from scrap in the EAF…"

So `ScrapLCI = (X_pr − X_re)·Y` is simultaneously the **burden** on scrap consumed and the **credit**
for scrap generated — "the worldsteel methodology assumes the burdens of scrap input and the credits
for recycling the steel at the end of the life of a product are equal, per kg, and that all scrap is
treated equally." Aligned with **ISO 14044:2006 §4.3.4.3** (closed-loop recycling allocation) and
**EN 15804**. "Recycling credits should be reported separately to maximise transparency."

**This is the one framework in the survey that does NOT use the cut-off approach for scrap** — and it
is exactly the reason the JRC notes that "none of the initiatives under investigation consider credits
for the expected recyclability of finished steel products": the six *threshold* frameworks all use
cut-off, while ISO 20915 (an LCI method, not a threshold) uses the closed loop.

**Electricity in worldsteel LCI / ISO 20915 is LOCATION-BASED, no certificate provision:**
> "Wherever possible, secondary data used for inputs to the processes shall be representative of the
> region the materials are sourced from. In particular, **electricity generation must be chosen to be
> the most representative of either a specific supplier (if comparable to the actual composition of
> electricity used at a given location is known) or the most appropriate regional or national grid
> mix.**"
> "The grid electricity production associated with individual sites can have a significant effect on
> the LCI, particularly regarding CO2 emissions…"

Supplier-specific mix is allowed (equivalent to the GHG Protocol Scope 2 supplier-specific method),
but **there is no provision for RECs, GOs or PPAs.**

Two differences from the worldsteel LCI parent (IEA 2023, p.20): "First, the impact of **scrap
processing is outside the scope and boundary** of this methodology. Secondly, **ISO 20915 uses credits
for process gases, based on the savings from replacing the marginal energy for the country where the
facility is located** (including both heat production and the electricity generation mix)." That second
point matters: unlike ISO 14404's fixed world average, **ISO 20915's process-gas credit is
country-specific**.

Also: "For the worldsteel LCI methodology and ISO 20915, credits for scrap 'net exports' are included
for **cradle-to-grave assessments but not for cradle-to-gate assessments**."

Coverage: ISO/TC 17 has 27 participating + 40 observing members; full adoption would cover ~1 750 Mt
(~90% of global production). worldsteel LCI covers 160+ sites / ~400 Mt across 17 finished products.

### (d) Is there an ISO standard specifically for DRI/HBI as an intermediate product?

**No — not for emissions.**

1. **ISO/TC 17 (Steel)** publishes 14404-1/-2/-3/-4 and 20915. **ISO 14404-3 covers the *plant***
   ("Steel plant with EAF and coal-based or gas-based DRI facility"), not DRI as a traded commodity.
2. **ISO 14404-4 handles merchant iron only obliquely** — clause 7.2.3 is "Steel plants which purchase
   pig iron from outside", and the foreword names "purchased DRI used in the BF-BOF route" as an
   intermediate product with upstream emissions. That is DRI *as an input*, never a standalone DRI
   emissions-intensity standard with its own denominator.
3. **worldsteel / ISO 14404 default upstream factors exist** (gas-based DRI 0.780, coal-based DRI
   1.210 t CO2/t) but they are generic defaults inside a crude-steel calculation.
4. **ISO/TC 102** (Iron ore and direct reduced iron) publishes DRI standards — sampling, physical
   testing, HBI shipping safety (e.g. ISO 11256, ISO 15967) — but **nothing on CO2**.
5. The IEA flags the absence: "trade in iron could become more common, [so] guidance that is either
   more detailed or more generalised will be required" (2023); the IDDI proposal for "a common
   reporting point or boundary" is still a proposal.

**Practical consequence for merchant DRI/HBI**: the only operational frameworks putting a number on it
today are (i) ResponsibleSteel's Scope 3 requirement — the buyer must obtain a cradle-to-point-of-sale
declaration from the supplier, plus transport to site (10.4.5.a, .c, .d); (ii) the generic
worldsteel/ISO 14404 defaults; (iii) an EPD under ISO 20915 / ISO 14025. **There is no near-zero
*iron* threshold anywhere.**

---

## 5. GHG Protocol Scope 2 Guidance and the live revision

*(numbered 5 to follow the brief's ordering; it underpins every other entry here)*

Primary sources:
- GHG Protocol Scope 2 Guidance (2015):
  https://ghgprotocol.org/sites/default/files/2023-03/Scope%202%20Guidance.pdf
- Scope 2 Public Consultation document (Oct 2025):
  https://ghgprotocol.org/sites/default/files/2025-10/GHG-Protocol-Scope2-Public-Consultation.pdf
- Scope 2 Public Consultation Summary of Feedback (29 July 2026):
  https://ghgprotocol.org/sites/default/files/2026-07/S2-PublicConsultationSummaryofFeedback-2026.07.29.pdf
- https://ghgprotocol.org/blog/scope-2-standard-advances-isb-approves-consultation-market-and-location-based-revisions
- https://ghgprotocol.org/blog/release-ghg-protocol-opens-public-consultations-scope-2-and-electricity-sector-consequential
- https://ghgprotocol.org/blog/ghg-protocol-announces-key-standard-development-updates

### (a) The 2015 rules — still binding as of Sept 2026

**Dual reporting** (§1.5.1):
> "Companies with any operations in markets providing product or supplier-specific data in the
> form of contractual instruments **shall report scope 2 emissions in two ways** and label each
> result according to the method: one based on the location-based method, and one based on the
> market-based method. This is also termed 'dual reporting.'"

Also: "Not having contractual data for every site will not cause noncompliance."

**Contractual instruments** (Box 1.1):
> "Any type of contract between two parties for the sale and purchase of energy bundled with
> attributes about the energy generation, or for unbundled attribute claims… they can include
> **energy attribute certificates (RECs, GOs, etc.), direct contracts** (for both low-carbon,
> renewable, or fossil fuel generation), **supplier-specific emission rates**, and other default
> emission factors representing the untracked or unclaimed energy and emissions (termed the
> **residual mix**)."

**The five Scope 2 Quality Criteria** (Table 7.1). All contractual instruments shall:
1. "Convey the direct GHG emission rate attribute associated with the unit of electricity produced."
2. "Be the only instruments that carry the GHG emission rate attribute claim associated with that
   quantity of electricity generation."
3. "Be tracked and redeemed, retired, or canceled by or on behalf of the reporting entity."
4. "Be issued and redeemed **as close as possible to** the period of energy consumption to which
   the instrument is applied."
5. "Be sourced from **the same market** in which the reporting entity's electricity-consuming
   operations are located and to which the instrument is applied."

Plus criterion 6 for utility-specific factors (delivered electricity; sold-off attributes
characterised as residual mix).

The two load-bearing weaknesses: **Criterion 4 says "as close as possible", not "same hour"** — in
practice read as annual matching; **Criterion 5 says "same market", undefined** — in practice read
as national or wider.

### (b) Revision timeline and status

| Date | Event |
|---|---|
| Nov 2022 | Initial public consultation across the corporate suite; **400+ submissions** on Scope 2 alone |
| 2024–2025 | Scope 2 Technical Working Group works "for over a year" on Phase 1 topics |
| **14 July 2025** | Independent Standards Board votes: **10–1 approves** LBM + MBM revisions for consultation; **7–4 declines** the Marginal Impact Method for avoided emissions in its current form, redirecting it to the Actions and Market Instruments (AMI) workstream |
| **20 Oct 2025 – 19 Dec 2025**, extended to **31 Jan 2026** | Phase 1 public consultation: *Scope 2 inventory topics* and *Consequential Electricity-Sector Emissions Impacts* |
| **29 July 2026** | GHG Protocol publishes the **Scope 2 Public Consultation Summary of Feedback** (~1,100 responses, 56 countries) |
| **July 2026** | GHGP announces **consolidation with ISO** — a single co-branded corporate standard merging GHGP Scopes 1/2/3 + AMI with **ISO 14064-1** |
| **Q2 2027** | Integrated public consultation on the consolidated corporate standard |
| **Late 2027** | Anticipated publication of the revised Scope 2 Standard, then a multi-year phased transition |

**A revised standard has NOT been published.** It is post-consultation draft. ISB chair:
Alexander Bassen. GHGP CEO Tim Mohin framed the ISO merger as something that "will simplify
reporting, reduce duplication, and provide greater consistency across markets and jurisdictions."

**The market-based method is being retained, not abolished.** What is proposed is a large
tightening of its qualifying criteria.

### (c) What exactly is proposed — six MBM changes

1. **Quality Criteria 4 → hourly matching.** Proposed: "All contractual instruments used in the
   market-based method for scope 2 accounting shall be **issued and redeemed for the same hour**
   as the energy consumption to which the instrument is applied, except in certain cases…".
   Hierarchy where hourly instruments are unavailable: hourly instruments → monthly/annual +
   hourly production meter data → facility-specific profile → regional public profile. Demand
   side: meter data → facility load profile → market-boundary profile → time-of-use average →
   **flat average** (annual total / 8,760).
2. **Quality Criteria 5 → deliverability.** Proposed: "All contractual instruments used in the
   market-based method shall be sourced from generation that is deemed **deliverable to the
   consuming load**." Facilities are "considered as located at [their] first point of
   interconnection to a transmission network."
3. **New Criterion 6, Standard Supply Service (SSS)** — a reporter "shall not claim more than its
   **pro-rata share**" of publicly funded/mandated/shared clean resources (RPS, CES, FIT, publicly
   owned generation, regulated default service). Unclaimed pro-rata shares are **not transferable**.
4. **New Criterion 9, residual mix redefined** — "the GHG intensity of electricity, within the
   relevant market boundary and time interval, that is not claimed through contractual
   instruments, **including voluntary purchases or Standard Supply Service allocations**."
5. **Fossil-only fallback.** "Generic location-based average grid emission factors shall not be
   used." Absent a residual mix, use a fossil-only grid-average or a default fossil factor, and
   "where no clear information is available, reporters should use the most conservative applicable
   default (e.g., coal or oil) rather than the lowest-emitting option."
6. **Feasibility measures** — load profiles, exemption thresholds, legacy clause, phased
   implementation.

**Exemptions:** four options, the GWh figure still `[X]` (undecided) — (1) consumption up to [X]
GWh/yr in a deliverable market boundary; (2) SME categorisation (drawing on draft SBTi
categorisation, itself based on the EU CSDDD); (3) either; (4) both. Critically: "The exemption
would apply **only to hourly matching**; the deliverability requirement in Quality Criterion 5
still applies."

**Legacy clause:** under consideration; would let pre-existing contracts continue to count even if
they fail hourly/deliverability. Alternative floated: a single uniform effective date with lead
time plus disaggregated disclosure.

**Consequential / impact accounting is deliberately kept OUT of the inventory.** "Both have value
but address different questions and should remain distinct in reporting to preserve like-for-like
inventories." The AMI **"multi-statement" reporting approach** would report three separate
components: (i) operational emissions, (ii) market-based emissions from instruments, (iii) GHG
impact from actions using consequential methods.

**Notably absent: no incrementality/additionality requirement inside scope 2** — the single
most-repeated criticism in the feedback.

### (d) Consultation feedback (July 2026) — ~1,100 responses, 56 countries

**Q71 — support for hourly matching (QC4):**

| Respondent type | n | Support (4–5) | Neutral (3) | Low/no (1–2) |
|---|---|---|---|---|
| **Grand total** | **909** | **22%** | **7%** | **70%** |
| Company | 429 | 12% | 7% | **82%** |
| Industry group | 76 | 12% | 6% | **82%** |
| Energy supplier/utility | 81 | 25% | 12% | 62% |
| Consultant | 81 | 33% | 7% | 59% |
| Academia/research | 31 | 47% | 9% | 45% |
| NGO/civil society | 71 | 48% | 5% | 45% |
| Data/analytics or software provider | 21 | **81%** | 0% | 19% |
| Eastern Asia (region) | 156 | 15% | 3% | 81% |

**Q83 — support for deliverability (QC5):**

| Respondent type | n | Support (4–5) | Neutral (3) | Low/no (1–2) |
|---|---|---|---|---|
| **Grand total** | **875** | **30%** | **11%** | **59%** |
| Company | 413 | 19% | 10% | **71%** |
| Industry group | 72 | 13% | 17% | 71% |
| Consultant | 82 | 48% | 11% | 41% |
| NGO/civil society | 69 | 55% | 9% | 36% |
| Data/analytics provider | 20 | 80% | 5% | 15% |
| Government institution | 11 | 73% | 27% | 0% |
| Eastern Asia (region) | 144 | 9% | 10% | 81% |

Corporates and industry groups strongly opposed; NGOs, academia and software vendors in favour.
Deliverability is materially less contested than hourly matching. Some softening in the next draft
is likely, but the ISB, not the respondent majority, decides.

### (e) THE DIRECT ANSWER FOR AUSTRALIA AND BRAZIL

**Under the current (2015) rules: yes, a PPA/GO can zero out scope 2.** A plant in Australia or
Brazil can buy LGCs / I-RECs / GOs, or sign a physical or virtual PPA, and report a **market-based
scope 2 of zero** — subject only to "as close as possible" timing (read as annual) and "same
market" (read as national). The existing catch: **dual reporting means the location-based figure
still shows the coal-heavy grid**, and that number does not move.

**Under the proposed rules: substantially harder — and Australia and Brazil are named explicitly.**
The deliverability table names them in the first category, markets with zonal pricing where the
**zonal boundaries become the market boundary**:

> "1. Australia's National Electricity Market
> 2. The electricity market operated by Brazil's Chamber of Electric Energy Commercialization [CCEE]
> 3. …ENTSO-E
> 4. …Russia"

Four compounding effects for an iron/steel plant:

1. **No exemption.** A DRI/EAF plant is a several-hundred-GWh-to-TWh/yr load, far above any
   plausible `[X]` GWh threshold and failing SME categorisation under all four exemption options.
   It **must** hourly match. The design intent is stated: thresholds that exempt "a majority of CDP
   reporting companies" while keeping "the vast majority of electricity load on the grid" subject
   to hourly matching.
2. **Hourly matching against a near-baseload industrial load** requires genuine round-the-clock
   procurement — firmed supply, storage, or overbuild — not an annual certificate volume. A real
   cost, not a paperwork change.
3. **Zonal, not national, sourcing.** In Australia certificates must come from generation in the
   same NEM zone, not anywhere in the NEM, unless Alternate Methodology 1 (adjacent market, and
   average price at consumption < **1.05×** average price at generation in that hour) or Alternate
   Methodology 2 (exclusive firm transmission rights recognised by all transmission operators along
   the path, demonstrated hourly "with no direct counterbalancing reverse transactions") applies.
   Brazil is the same structure via CCEE zones.
4. **The residual/fossil fallback punishes partial matching.** Every unmatched hour can no longer be
   valued at the location-based grid average; it must use a residual mix excluding all voluntarily
   claimed and SSS generation (necessarily dirtier than grid average) or, absent a residual mix
   (likely in both Australia and Brazil today), a **fossil-only or default coal factor**. A plant
   matching 70% of hours values the other 30% at roughly coal. **The proposed rules make partially
   matched scope 2 look worse than it does today, even as the plant procures more clean energy.**

Timing: nothing binds before **late 2027 publication plus a multi-year phase-in**. A plant
commissioning in the late 2020s should be modelled against the new rules; anything reported now
runs on the 2015 criteria.

### (f) How the steel frameworks each land on the Scope 2 question

| Framework | Scope 2 rule |
|---|---|
| **ResponsibleSteel** | Default location-based per ISO 14064-1:2018 Annex E.2 ("emission factor that best characterises the pertinent grid"); **RECs / PPAs / vPPAs / green tariffs permitted** under ISO 14064-1 E.2.2, with mandatory disclosure of source and quantity. Offsets banned. |
| **LESS** | "the **real emission factor from the actual electricity mix** used by the plant operator"; **GOs from officially recognised registers permitted**, explicitly including green electricity "purchased on an accounting-related basis and **not physically delivered**". Registers named are German/EU. |
| **LESS 2026 recommendation to the EC** | "**Market-based instruments, such as Power Purchase Agreements (PPAs) and Guarantees of Origin (GoOs), should be recognised**… provided that they are substantiated by officially recognised registries and robust verification systems." |
| **IEA** | Counted; **method deliberately unresolved** — "suitable methodologies for accounting for electricity emissions from the grid, including… use of renewable electricity credits… Further discussion on this topic is likely to be needed." Negative factors floored at zero. |
| **First Movers Coalition** | Defers to GHG Protocol Scope 2; **book & claim expressly carved in for electricity** even though banned for the steel itself. The 2022 vPPA-additionality condition was **dropped** in the 2025 revision. |
| **worldsteel CO2 methodology / ISO 14404** | **Fixed world-average factor, 0.504 t CO2/MWh.** Neither location- nor market-based. A PPA does nothing; a clean grid does nothing either. |
| **ISO 20915 / worldsteel LCI** | **Location-based**: "the most appropriate regional or national grid mix", or a supplier-specific mix where known. **No REC / GO / PPA provision.** |
| **GSCC Steel Climate Standard** | **Market-based is the DEFAULT.** §6.3: "If a company is not using contractual instruments and market-based information is not available, **residual grid mix factors** shall be used… **A location-based approach using average grid emission factors shall only be used when market-based and residual mix data is unavailable.**" Instruments must meet ISO 14064-1 E.2.2, be reported as a separate line item, and **cannot** be applied to electricity from process off-gases. |
| **SBTi Steel Guidance** | Scope 1 + scope 2 for EAF power inside the iron & steel SDA boundary; method not specified. Renewable-sourcing targets are "an acceptable **alternative** to scope 2 emission reduction targets". |
| **RMI Steel Emissions Reporting Guidance** (https://rmi.org/wp-content/uploads/2022/09/steel_emissions_reporting_guidance.pdf) | Both methods allowed, but requires a **residual mix factor** to prevent double counting, and RECs/GOs may **not** be applied to off-gas-derived electricity. |
| **JRC ESPR classes (draft, Apr 2026)** | **Location-based, country-average only.** EU-average factor (EEA) for EU installations; CBAM default values reflecting "the average electricity mix of the country of origin" for imports. **No mention of PPAs, GOs or RECs anywhere in the document.** |
| **German Klimaschutzverträge (CCfD)** | **Scope 2 excluded outright** — "Die Treibhausgasemissionen des Vorhabens ergeben sich aus den Treibhausgasemissionen der geförderten Anlagen (Scope-1-Emissionen)", and the ETS benchmarks were "reduced by the indirect emissions for electricity not to be taken into account in this funding programme". Electricity enters only as a *cost* term. |
| **US GSA / Buy Clean California / Colorado / Minnesota / Canada NRC** | **Silent** — all defer to the underlying steel PCR (UL or SCS) and ISO 14025/21930; the PCR text does not resolve contractual instruments. **Unresolved.** |

The **UK government's technical annex** to its low-carbon products consultation (June 2025) is the best
neutral side-by-side, and states the LESS position verbatim:
> "**scope 2 emissions data is determined using the real emission factors from the actual electricity
> mix used by the site. However, LESS permits the use of accounting-related allocation for energy use,
> given that guarantees of origin from officially recognised registers are provided. These producers
> can claim lower energy emissions, even if such energy was not physically used at the site.**"
https://assets.publishing.service.gov.uk/media/68592299b46781eacfd71dfc/policy-framework-to-grow-the-market-for-low-carbon-industrial-products-technical-annex.pdf

The same annex records that **ResponsibleSteel applies a "burden of the doubt" multiplier of 1.2 to
default upstream data, or 1.6 for coal, coke and natural gas.**
| **EUROFER stainless proposal** | Benchmark hard-wired to a **fixed EU grid mix of 376 kg CO2/MWh**. |
| **UNESID (Spain) proposal** | EU average electricity mix for Scope 2 in the reference plants; an "**arbitrary electricity footprint**" is listed as a **knock-out (KO)** criterion. |

---

## 6. The EU's voluntary carbon-intensity label — and the parallel ESPR delegated act

Two EU vehicles are running at once and they must not be conflated:

1. **The voluntary carbon-intensity label** announced in the **European Steel and Metals Action Plan**
   (COM(2025) 125, March 2025) and to be delivered through the **Industrial Decarbonisation
   Accelerator Act (IDAA)**. Per the JRC: "the **Industrial Decarbonisation Accelerator Act will
   introduce a voluntary label on the carbon intensity of steel** (European Commission, 2025a). This
   was confirmed in the Commission's Steel and Metal Action Plan (European Commission, 2025b)…"
   The Action Plan says the methodology will be "based on a simple methodology with ETS data and
   building on the CBAM methodology", starting with steel in 2025. LESS's December 2025 statement
   expected the Commission to put it forward on **10 December 2025** as part of the Act.
   EUR-Lex: https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=celex%3A52025DC0125
   Communication PDF: https://single-market-economy.ec.europa.eu/document/download/7807ca8b-10ce-4ee2-9c11-357afe163190_en?filename=Communication+-+Steel+and+Metals+Action+Plan.pdf

2. **The ESPR delegated act for iron and steel** — Regulation (EU) 2024/1781 (Ecodesign for
   Sustainable Products Regulation). Iron and steel is "**the first intermediate product addressed
   within the ESPR framework**", with a **cradle-to-gate** boundary "from raw material extraction… to
   the factory gate". The JRC Product Bureau runs the preparatory study.
   Project page: https://susproc.jrc.ec.europa.eu/product-bureau/product-groups/642/home
   Documents: https://susproc.jrc.ec.europa.eu/product-bureau/product-groups/642/documents

**The ESPR track is where the actual numbers are.** The JRC's draft report *Classes of environmental
performance — ESPR five representative iron and steel intermediate products* (JRC144483, Seville,
**1 April 2026**) is the proposal LESS attacked in its 12 May 2026 statement.
PDF: https://susproc.jrc.ec.europa.eu/product-bureau/sites/default/files/2026-04/DraftStudyESPRSteelClassesPerformance1April_0.pdf

Other documents in the same consultation round:
- *Preparatory study iron & steel*, March 2026 (11.8 MB) — https://susproc.jrc.ec.europa.eu/product-bureau/sites/default/files/2026-03/Preparatory%20study%20iron%20%26%20steel_March2026_0.pdf
- *Task 4 LCA analysis*, March 2026 — https://susproc.jrc.ec.europa.eu/product-bureau/sites/default/files/2026-03/Task%204%20LCA%20analysis_0.pdf
- *Study on Recycled Content Steel*, March 2026 — https://susproc.jrc.ec.europa.eu/product-bureau/sites/default/files/2026-03/StudyonReCoSteel_March%202026_0.pdf
- *ESPR Steel DPP Content Proposal* (digital product passport), March 2026
- *Second Stakeholder Consultation presentation*, 13/15 April 2026
- Initial preparatory study (Tasks 1-3), 15 May 2025; first stakeholder meeting, 25 June 2024

### (a) Scope — five representative intermediate products, and NO merchant iron

> "the five representative steel intermediate products [are] **Hot Rolled Coil (HRC), Wire Rod (WR),
> Cold rolled coil Galvanized (CRCG), Electrical sheet (ES) and Stainless Steel (SS)**… These five
> representative steel intermediate products represent **50% of EU apparent consumption**."

**DRI, HBI and merchant iron do not appear anywhere in the document.** The classes attach to finished
intermediate *steel*, not to iron. So the EU labelling track — like every other framework surveyed —
gives merchant HBI no class of its own; it enters only as an input to a classified steel product.

### (b) Scope 2 — COUNTED, and STRICTLY LOCATION-BASED / COUNTRY-AVERAGE

This is the most important single finding for this project, and it runs directly against what LESS is
lobbying for. Verbatim from the JRC draft:

> **"Electricity:** The electricity modelling in the RIS follows the specifications set out in both the
> EU ETS and the CBAM. **In the CBAM regulation the emissions for the use of electricity are calculated
> using the average electricity emission factor of the exporting country for the reporting year,
> expressed in kg CO2e/kWh** [footnote 15]. **When providing results for extra-EU country starting from
> CBAM default values, the emissions of electricity are taken directly from the default values, which
> reflect the average electricity mix of the country of origin, thus the analysis presented here is
> fully aligned with the CBAM.** In the EU ETS, emissions of electricity are calculated using the **EU
> average emission factor** [footnote 16 — EEA GHG emission intensity of electricity generation in
> Europe]. Accordingly, to be in line with what [is] prescribed in the EU ETS, in the RIS the EU average
> emission factor has been used. Furthermore, **a sensitivity analysis has been carried out**. The
> electricity impact has been calculated using the emission factor of each EU Member State considered.
> Consequently, the analysis includes separate scenarios for countries such as **France, Germany, and
> Poland**. These nations have different electricity generation mixes, ranging from a high share of
> low-carbon nuclear and hydro power in France to a coal-heavy mix in Poland, resulting in distinct
> emission factors."

Footnote 15 cites **Commission Implementing Regulation (EU) 2025/2547 of 10 December 2025**, laying
down rules for the application of Regulation (EU) 2023/956 as regards the methods for calculating
embedded emissions, **Annex II, Section D.2**.

**There is no mention of market-based instruments, PPAs, Guarantees of Origin or RECs anywhere in the
document.** Electricity is modelled at country-average intensity — EU average for EU installations,
country-of-origin average for imports (via CBAM defaults), with a member-state sensitivity.

Three consequences that matter directly:
- **A Brazilian or Australian HBI/steel producer is charged its country's average grid factor.** Under
  the ESPR class methodology as drafted, buying Brazilian or Australian certificates does not help.
- **The France / Germany / Poland sensitivity is exactly the comparison in this project**, run by the
  JRC itself, and it is the reason the EU-average default matters: it flattens precisely the
  distinction between a French-nuclear EAF and a Polish-coal one.
- **This is how Scope 2 gets bolted back on to a CBAM-derived Scope 1 number.** CBAM excludes indirect
  emissions for iron and steel; the ESPR class methodology reintroduces them by applying the CBAM
  *electricity* rule (which exists for the goods where indirect emissions do count) to the steel case.

Also, on the boundary construction for imports:
> "For the primary production (BF-BOF) of products imported from outside the EU, the representativeness
> and CO2 emissions of extra-EU countries are accounted for by adopting the **CBAM default values as a
> starting point. These default values encompass direct and indirect emissions, including those
> associated with electricity consumption.** To ensure a comprehensive assessment, **Scope 3 emissions,
> particularly those linked to upstream processes and coke production, are further integrated using the
> ongoing method developed by the JRC.** For the secondary production (EAF) of products imported from
> outside the EU, the inventory defined in Task 4 is used, alongside **country-specific emission factors
> for both electricity and natural gas**…"

and for EU production:
> "Since the **EU ETS only accounts for direct emissions, these are supplemented with indirect (Scope 2)
> emissions and other relevant Scope 3 emissions** to ensure a complete life-cycle perspective."

### (c) Scrap

> "**Scrap:** The use of scrap in steelmaking is considered differently in the EU ETS, CBAM and in the
> EF. In particular, **ETS and CBAM do not account for any emission related to scrap used** [footnote 17].
> In the definition of RIS, **the EF [Environmental Footprint] approach for intermediate products has
> been used, which accounts for the impacts of scrap collection and treatment.** When declared emissions
> under CBAM are used, allocating an impact to scrap will be straightforward. A declarant in the CBAM
> registry must report **the amount of total scrap and pre-consumer scrap used** in the production
> process. However, for ETS data, such calculation is not directly possible, as EU installations are not
> required to provide scrap usage figures… **At the current stage, all the scrap has been considered
> equally, without a differentiation between pre- and post-consumer** [scrap]."

Footnote 17: "The proposal for the revision of the CBAM regulation, proposes the **inclusion of
aluminium and steel pre-consumer scrap as precursor**, thereby allowing to attribute emissions to scrap
as a precursor."

**Critically: this is NOT a sliding scale.** Scrap affects the calculated footprint, but the class
thresholds are fixed absolute values with no scrap adjustment — which is precisely the design LESS and
ResponsibleSteel are lobbying against.

### (d) Thresholds — the actual proposed classes

Two methods are used. **HRC and WR use "Method 3 (fixed percentage)"** — thresholds tuned to the
observed distribution of production volume:

> "Under the selected Method 3 (fixed percentage), class boundaries are determined by predefining the
> share of products that should fall into each class and iteratively adjusting threshold values until
> these shares… are met… **the distribution was calibrated so that approximately the top two classes
> together cover ≥30% of global production volume.** The 30% share should be interpreted as indicative
> rather than fixed."

> "The rationale for the 30% criterion is linked to the potential application of the **B/C threshold in
> Green Public Procurement (GPP)**… If applied as a procurement criterion, it should be designed in a
> manner that **does not unduly restrict competition**… Ensuring that a critical mass of products can
> comply is essential to avoid situations where contracting authorities receive few or no valid bids,
> or where prices increase disproportionately."

**Hot Rolled Coil (HRC), Table 6** (t CO2eq / t HRC; left-inclusive, right-exclusive; E is ≥):

| Class | Carbon footprint range | Description | Plants | Share of plants | Volume (Mt) | Share of volume |
|---|---|---|---|---|---|---|
| **A** | **0.00 – 1.79** | Best-in-class | 23 | 18.0% | 35.75 | 7.6% |
| **B** | **1.79 – 2.66** | High performing | 36 | 28.1% | 134.23 | 28.5% |
| **C** | **2.66 – 3.10** | Mainstream | 38 | 29.7% | 173.12 | 36.8% |
| **D** | **3.10 – 3.75** | High emitting | 22 | 17.2% | 91.89 | 19.5% |
| **E** | **≥ 3.75** | Worst-in-class | 9 | 7.0% | 35.94 | 7.6% |
| Total | | | 128 | 100% | 470.94 | 100% |

> "Class A (≤ 1.79 tCO2eq/t): Technologies falling in Class A are predominantly **EAF and DRI-based**
> technology routes. EAF plants… concentrated around ~0.4–1.2 tCO2eq/t, clearly within Class A."
> "Class A contains several plants at very low intensities (around ~0.69–0.91 tCO2eq/t)…"

**Wire Rod (WR), Table 9** (t CO2eq / t WR). Note: "For the best-in-class band (Class A), it has been
ensured that **BF-BOF plants are not included** in this [band]."

| Class | Carbon footprint range | Description | Plants | Share of plants | Volume (Mt) | Share of volume |
|---|---|---|---|---|---|---|
| **A** | **0.00 – 0.87** | Best-in-class | 21 | 20.8% | 7.23 | 9.6% |
| **B** | **0.87 – 2.43** | High performing | 17 | 16.8% | 14.94 | 19.8% |
| **C** | **2.43 – 3.10** | Mainstream | 30 | 29.7% | 22.32 | 29.6% |
| **D** | **3.10 – 3.54** | High emitting | 28 | 27.7% | 26.15 | 34.7% |
| **E** | **≥ 3.54** | Worst-in-class | 5 | 5.0% | 4.65 | 6.2% |

> "The EAF plants… span roughly ~0.2 to ~0.8 tCO2eq/t, with most points well below 0.87 tCO2eq/t…
> **H2-DRI-EAF also sit within Class A**."

**CRCG, SS and ES use "Method 2 (fixed-interval min–max)"** — five equal-width bands between a
best-available-technology minimum and a conventional-route maximum (after Senatore et al., 2025):

| Class | **CRCG** (t CO2eq/t) | **Stainless steel (SS)** | **Electrical sheet (ES)** |
|---|---|---|---|
| A | 0.77 – 2.05 | 2.41 – 3.82 | 1.28 – 2.21 |
| B | 2.05 – 3.34 | 3.82 – 5.23 | 2.21 – 3.13 |
| C | 3.34 – 4.62 | 5.23 – 6.64 | 3.13 – 4.06 |
| D | 4.62 – 5.91 | 6.64 – 8.05 | 4.06 – 4.99 |
| E | 5.91 – 7.19 | 8.05 – 9.46 | 4.99 – 5.92 |

CRCG min/max scenarios: minimum PCF (best-performing) 0.53, maximum PCF (worst-performing) 7.19
t CO2eq/t product; range computed as 7.19 − 0.77 = 6.42.

### (e) The JRC's own rationale, and LESS's objection

JRC:
> "This approach recognizes that flat products are today predominantly produced through the BF-BOF
> route, resulting in higher carbon intensity compared to the scrap-based EAF route, which is more
> commonly used to produce long products. For instance… the methodology ensures that **the
> best-performing BF-BOF installations are recognised and rewarded** relative to those with a poorer
> environmental profile (the former being classified in **class B**, while the latter fall into lower
> classes). At the same time, the scale sends a clear signal to the market that producers capable of
> manufacturing flat products through **DRI-EAF or scrap-based EAF routes, which achieve significantly
> lower carbon intensities, can reach the highest performance class (class A)**."

LESS (12 May 2026), on exactly this:
> "**Scrap-based EAF: A switch to 100% renewable electricity yields no classification upgrade if
> already in Class A.**"
> "Hot-rolled coil: Natural gas-based DRI-EAF – despite significant emissions reductions and
> investments – **shares Class A with hydrogen-based DRI-EAF**, removing incentives to transition from
> gas to hydrogen."

**Both objections are borne out by the tables above.** HRC Class A is a 0-to-1.79 t CO2eq/t band that
holds everything from a 0.4 t/t scrap EAF to a 1.7 t/t NG-DRI-EAF. A scrap EAF already inside it gains
nothing from decarbonising its power, and an H2-DRI plant is indistinguishable from an NG-DRI plant.
**Under the EU's own emerging label, an H2-DRI-EAF and an NG-DRI-EAF making the same product would
carry the same class.**

### (f) Status

- ESPR delegated act for iron and steel: **preparatory study ongoing**, second stakeholder consultation
  April 2026, classes still a **draft report** with thresholds described as "indicative and may be
  revised". No delegated act adopted.
- LESS's ask: "finalise and implement the ESPR Delegated Act, ensuring a functioning certification
  system and **regulation for lead markets by 2027**", noting "The first large-scale European DRI plants
  will start in **2027**."

---

## 7. JRC report — *Defining low-carbon emissions steel* (JRC141817, EUR 40291, 2025)

Blanco Pérez, S., Arcipowska, A., Fiorese, G., Maury, T., Napolano, L., *Defining low-carbon
emissions steel: A comparative analysis of international initiatives and standards*, Publications
Office of the EU, Luxembourg, 2025. ISBN 978-92-68-26730-1, doi:10.2760/4271464.
PDF: https://publications.jrc.ec.europa.eu/repository/bitstream/JRC141817/JRC141817_01.pdf
DOI landing: https://data.europa.eu/doi/10.2760/4271464

### (a) What it is — and what it is NOT

**It does not propose a definition.** It is a *comparative analysis* of six frameworks, written to
inform the Clean Industrial Deal, the Industrial Decarbonisation Accelerator Act's voluntary label,
and the Steel and Metals Action Plan. Its headline finding:

> "**Absence of a Universal Definition**: Despite increasing efforts to decarbonize the steel
> sector, there is currently **no universally accepted definition of 'low-carbon emissions steel.'**
> Existing initiatives vary significantly in methodologies, scopes, and emissions thresholds…"

> "The findings reveal **substantial inconsistencies in scope, system boundaries, and emissions
> accounting methodologies, which undermine comparability across frameworks**, affecting market
> transparency and fair competition. Nonetheless, despite these divergences, the study identifies a
> trend towards long-term alignment, with most initiatives targeting similar emissions intensities
> by 2050."

The six frameworks compared: **IEA**, **ResponsibleSteel**, **LESS**, **Climate Bonds Initiative**,
**Global Steel Climate Council (GSCC) Steel Climate Standard**, **Chinese Method C2F Steel (CISA)**.
Annex 1 covers, without full comparison: **First Movers Coalition**, **SBTi**, **Italian Steel
Federation**, **Indian Ministry of Steel Green Steel Taxonomy**, **UNESID (Spain)**, **Sandbag**,
**CRU**.

### (b) System boundaries — Table 2a and Annex 2

| Framework | Level | Declared unit | Processing stage | Approach | System boundary | GHGs |
|---|---|---|---|---|---|---|
| **IEA** | Site | t crude steel | Crude steel | Sliding scale, fixed | Cradle-to-gate (Scope 1 + 2 + Scope 3.1, 3.3, 3.4 partially) | CO2, indirect CH4 only from fuel supply |
| **ResponsibleSteel** | Site, Product | t crude steel | Crude steel | Sliding scale, fixed | Cradle-to-gate (Scope 1 + 2 + 3.1 + 3.3 + 3.4 + 3.5) | CO2, N2O, CH4, HFCs, PFCs, SF6, NF3 |
| **LESS** | Product | t hot-rolled steel | Hot-rolled steel | Sliding scale, fixed | Cradle-to-gate (Scope 1 + 2 + 3.1 + 3.3 + 3.4 + 3.8) | CO2 for Scope 1&2; CO2, CH4, N2O, HFCs (partial), SF6 for Scope 3 |
| **Climate Bonds Initiative** | Company | t finished steel | Finished steel | Weighted pathway, progressive | Scope 1 + 2 + 3.3 only | CO2 |
| **GSCC Steel Climate Standard** | Product, Site | t hot-rolled steel | Hot-rolled (flat/long) | Company-specific trajectory, progressive | Cradle-to-gate (Scope 1 + 2 + 3.1 + 3.3 + 3.4) | CO2, N2O, CH4, HFCs, PFCs, SF6, NF3 |
| **Chinese C2F (CISA)** | Site, Line, Product | t crude or hot-rolled steel | Crude & hot-rolled | Sliding scale, fixed | Cradle-to-gate (Scope 1 + 2 + 3.1 partially + 3.3) | CO2 |

Key JRC observations:
> "As steel products are considered intermediates, most methodologies extend their boundaries from
> cradle to gate… **Scope 1 (direct) and Scope 2 (indirect) emissions are generally included across
> the frameworks analysed.**"

> "While the IEA and ResponsibleSteel define their system boundaries at the crude steel stage (iron
> and steel making), GSCC, LESS and the Chinese Method C2F Steel extend coverage to hot-rolled
> steel, and the Climate Bonds Initiative also includes cold rolling and coating processes."

> "Regarding Scope 3 emissions upstream, **LESS, ResponsibleSteel and the Steel Climate Standard are
> the most comprehensive**… In contrast, the Climate Bonds Initiative does not include any Scope 3
> emissions upstream… the IEA includes some upstream emissions, specifically those related to iron
> ore, limestone, and coal mining but omitting non-ferrous ore mining and emissions related to iron
> and steel recycling and sorting, and to ferro-alloys."

On process gases (relevant to BF-BOF comparators, not to the H2-DRI/MOE cases):
> "While the IEA, ResponsibleSteel, the Climate Bonds Initiative, the CMC2FS and the GSCC standard
> account for these emissions, **the LESS does not include them in its system boundaries**. This
> discrepancy could lead to significant variations… accounting for **10-40% of total emissions**."

On slag and EoL:
> "According to IDDI (2023), **slag allocation can lead to a potential variation of up to 10%** in
> reported emissions… the inclusion of EoL recycling credits in product emissions footprints could
> result in a **variation of up to 14%**."

**CRITICAL FOR THIS PROJECT — Annex 2 lists "Briquetting" as an explicit process block, and ALL SIX
frameworks include it:**

| Process block | Description | IEA | GSCC | LESS | ResponsibleSteel | CBI | CMC2FS |
|---|---|---|---|---|---|---|---|
| Iron Ore Mining | Upstream extraction of iron ore | Yes | Yes | Yes | Yes | No | No |
| Limestone Quarry | Upstream extraction of limestone | Yes | Yes | Yes | Yes | No | No |
| Coal Mining | Upstream extraction of coal | Yes | Yes | Yes | Yes | No | No |
| Metal Mining (other than iron) | Upstream extraction of metals | No | Yes | Yes | Yes | No | No |
| External Scrap Collection/Processing | "Upstream collection, processing, and transport of scrap metal" | **No** | Yes | Yes | Yes | No | No |
| Other Carbon Inputs | charge/injection carbon, pet coke, charcoal | Yes | Yes | Yes | Yes | No | No |
| Pellet Plant | Iron ore processing into pellets | Yes | Yes | Yes | Yes | Yes | Yes |
| Sinter Plant / Lime Kiln / Coke Oven | — | Yes | Yes | Yes | Yes | Yes | Yes |
| Alloy/Additive Metal Processing | — | **No** | Yes | Yes | Yes | Yes | **No** |
| Refractory Production | — | No | Yes | Yes | Yes | No | No |
| **Direct Reduction Reactor** | "Ironmaking process via direct reduction of iron using CO and H2 coming from reformed natural gas, syngas or coal" | Yes | Yes | Yes | Yes | Yes | Yes |
| **Briquetting** | "Ironmaking process in which hot direct reduced iron from reduction kiln if formed into briquettes" | **Yes** | **Yes** | **Yes** | **Yes** | **Yes** | **Yes** |
| **Syngas/Hydrogen Production** | "Upstream production of hydrogen and syngas" | Yes | Yes | Yes | Yes | Yes | Yes [g] |
| **Power/Electricity Production** | "**Upstream production of electricity**" | **Yes [b]** | **Yes** | **Yes** | **Yes** | **Yes** | **Yes [g]** |
| **Upstream Materials Transport** | "Transportation of raw materials included in boundary to iron- and steel-making sites" | Partial [c] | Yes | Yes | Yes | No | No |
| Electrode/Graphite Production | "Upstream production of electrodes and graphite used in EAFs" | No | Yes | Yes | Yes | No | No |
| Electric Arc Furnace | — | Yes | Yes | Yes | Yes | Yes | Yes |
| Hot Rolling | — | No | Yes | Yes | No | Yes | Yes |
| Emissions from process gases | — | Yes | Yes | **No** | Yes | Yes | Yes |
| Export of power (credit) | — | No | Yes | Yes [e] | Yes | Yes | Yes |

Footnote [b] is directly on point for the scope 2 question:
> "**The iron and steel sector carbon budget provided in IEA Net Zero by 2050 Study (2021) accounts
> for industrial process emissions only, they do not include indirect emissions from purchased
> electricity. However, the analytical boundary used in the definition for near-zero and
> low-emissions materials in IEA Achieving Net Zero Heavy Industry Sectors in G7 Members Study
> (2022) include imported electricity.**"

Footnote [c]: IEA upstream transport is "Only for iron and limestone."
Footnote [e]: "LESS allows the export of gases to the extent permitted by the EU ETS."

### (c) Thresholds — JRC Table 4

| Organisation | Stage | Approach | Iron-ore-based (0% scrap) | Scrap-based (100% scrap) |
|---|---|---|---|---|
| **IEA** | Crude steel | Sliding scale | **0.4 t CO2e/t crude steel** | **0.05 t CO2e/t crude steel** |
| **ResponsibleSteel** | Crude steel | Sliding scale | **0.4 t CO2e/t crude steel** | **0.05 t CO2e/t crude steel** |
| **LESS — QST** | Hot rolled | Sliding scale | **0.52 t CO2/t hot rolled** | **0.17 t CO2/t hot rolled** |
| **LESS — reinforcing & structural (BST)** | Hot rolled | Sliding scale | **0.47 t CO2/t hot rolled** | **0.12 t CO2/t hot rolled** |
| **Climate Bonds Initiative** | Crude steel | Pathway by technology route | 1.81 t CO2e/t by 2030; **0.12 t CO2e/t by 2050** | 0.32 t CO2e/t by 2030; **0.12 t CO2e/t by 2050** |
| **GSCC Steel Climate Standard** | Hot rolled | Company-specific trajectory | Flat: 1.31 t CO2e/t by 2030, **0.12 by 2050**; Long: 1.11 by 2030, **0.12 by 2050** | (same — route-independent) |
| **Chinese C2F (CISA)** | Crude or hot rolled | Sliding scale | **0.4 t CO2/t** | **0.05 t CO2/t** |

CISA's five classes at 0% / 100% scrap (Figure 16): E 2.19, D 1.78, C 1.21, B 0.81, A 0.40 at 0%
scrap; and 0.41 / 0.35 / 0.24 / 0.15 / **0.05** at 100% scrap.

> "Despite the difference in their system boundaries, the IEA, ResponsibleSteel and the Chinese
> Method C2F Steel standards align in their proposed quantitative threshold… This alignment stems
> from the fact that all frameworks are based on the IEA's approach."

> "Frameworks that cover broader system boundaries, such as the GSCC and LESS, include more upstream
> processes, **leading to higher emissions thresholds**. Conversely, initiatives with narrower
> boundaries, such as the IEA and the CMC2FS, exclude specific emissions sources, leading to lower
> thresholds."

### (d) Scrap handling

> "**All the frameworks included in the analysis adopt a cut-off approach, assuming that scrap
> carries no emissions burden from its previous life.** Nonetheless, some frameworks, such as LESS,
> ResponsibleSteel and the Steel Climate Standard, account for the emissions related to collection
> and processing of scrap and incorporate these emissions within their system boundaries."

> "Regarding the recyclability of products at the End of Life (EoL) stage, **none of the initiatives
> under investigation consider credits for the expected recyclability of finished steel products**."

Methodological taxonomy (Table 3):

| Approach | Initiatives | Temporal | Scrap correlation |
|---|---|---|---|
| Sliding scale | IEA, ResponsibleSteel, LESS, CMC2FS | Fixed | Yes |
| Weighted pathway | Climate Bonds Initiative | Progressive | Yes |
| Product-based pathways | GSCC Steel Climate Standard | Progressive | **No** |

### (e) Treatment of electricity / scope 2 — what the JRC does and does NOT say

The JRC report **does not anywhere discuss location-based vs market-based Scope 2 accounting**.
This is a genuine gap in the report. It establishes only that Scope 2 is *inside the boundary* of
all six frameworks; it says nothing about whether a PPA or GO can zero it.

The **one place where the market-based question surfaces at all** is in Annex 1, in the **UNESID
(Spanish Steel Federation)** proposal, which treats a freely-chosen electricity footprint as a
knock-out criterion:

> "UNESID suggests modelling the reference plants as the average of the EU steel mills and
> **considering the EU average electricity mix when calculating Scope 2 emissions** (average
> European emission used to define the sectoral references within the EU ETS)."

and, in Figure 25's annotations, listing what is acceptable vs unacceptable:

> "Actual production route. **OK** / Inputs allocation per Steel family. **OK** / Arbitrary insumes
> reallocation, credits for by-products allocation, **arbitrary electricity footprint. KO**"

("KO" = knock-out.) UNESID also "permits carbon offsetting by purchasing emissions rights but does
not allow the application of credits for substances and energy that are dispensed beyond the
accounting limits" — a mirror-image position to LESS.

UNESID's other structural proposal: **two product families**. Type 1 (products where BF-BOF and EAF
compete: carbon hot rolled coil, seamless tubes, drawing steel, railway material) gets a sliding
scale; Type 2 (stainless, high-alloy, specialty — EAF-only in the EU) gets a **"column scale"**
based on carbon footprint alone with no recycled-content adjustment, benchmarked at a 70% recycled
/ 30% virgin reference plant. UNESID also uses **recycled metal content** rather than scrap share,
"because the production of one tonne of steel needs more than one tonne of scrap input", and adds an
**A+ class** reachable via offsetting.

### (f) Imported iron / DRI / HBI

The JRC does not treat imported iron as a distinct policy question. What it establishes indirectly
is decisive: since **briquetting and the direct reduction reactor are inside all six boundaries**,
and since ResponsibleSteel/LESS/GSCC count Scope 3.1 raw materials and 3.4 inbound transport,
imported HBI is inside the boundary of the EU melter's number under those three. The IEA and CBI
and CMC2FS treat upstream transport partially or not at all.

### (g) Other frameworks worth knowing from Annex 1

- **SBTi Steel Guidance v1.0 (2023)**: corporate-level, two separate pathways (ore-based and
  scrap-based) treated as "separate industries"; targets cover company-wide Scope 1 and Scope 2;
  Scope 3 must be in near-term targets "if relevant Scope 3 emissions are 40% or more of total
  Scopes 1, 2 and 3 emissions". JRC criticism: corporate not product/plant level.
- **Indian Ministry of Steel Green Steel Taxonomy (December 2024)**: "green steel" = **< 2.2 t
  CO2e/t finished steel**; star ratings — **5-star < 1.6**, **4-star 1.6–2.0**, **3-star 2.0–2.2**.
  Covers "Scope 1, Scope 2, and limited Scope 3 emissions up to finished steel production".
  Reviewed every three years; NISST is the MRV nodal agency.
  https://pib.gov.in/PressReleasePage.aspx?PRID=2083839
- **Sandbag** ("From niche to mainstream: Shaping demand for green steel", Zaccaro 2024): rejects
  the sliding scale, proposes separate tiered A-D scales for flat and long products on embedded
  emissions alone. "The proposal argues that the sliding scale discriminates against scrap."
  Notes that "in the EU, long steel producers already employ 100% scrap and electricity, leaving
  limited avenues for further emission reduction. **However, there are opportunities for additional
  abatement within the EAF route, depending on the electricity procurement strategy of individual
  steel plants.**"
- **EUROFER Stainless Steelmakers, "Low CO2 Stainless Steel"** — the one proposal in the JRC set
  that hard-wires a *fixed location-based grid factor*. "The calculation is based on the **average
  grid mix (376 kg CO2/MWh)** and a fixed scrap input (70% stainless scrap)." A grade-specific model
  footprint `MF_n` is calculated by linear regression over 20 grades; the threshold is
  `T_n = (100% − R)·MF_n`. Labels: A: `MF_n < T_n`; B: `T_n < MF_n < 2T_n`; C: `2T_n < MF_n < 3T_n`;
  D: `3T_n < MF_n < 4T_n`; E: `4T_n < MF_n` **or unknown**. EUROFER stainless rejects the sliding
  scale because "the production of ferroalloys (upstream Scope 3) still contributes more than 70% of
  the Product Carbon Footprint".
- **Federacciai (Italian Steel Federation)**: a dual index — a *carbon footprint index* (kg CO2e/t,
  route-agnostic, for green procurement) plus a *decarbonisation effort index* (% of the way to net
  zero, calibrated separately for the integrated route incl. DRI/EAF and for EAFs, for access to
  incentives and sustainable funds). "The absolute threshold of emissions intensity to achieve class
  A+ (net zero) is equivalent for both technology routes." Flagged by JRC as "a conceptual exercise
  rather than a fully detailed plan".
- **CRU Consulting**: defines green steel as any steel below **0.4 t CO2e**.
- **UK Construction Leadership Council, "Five Client Carbon Commitments" (early 2024)** — a green
  *public procurement* proposal. "Objective 4 specifically targets the gradual reduction and eventual
  elimination of the most carbon-intensive steel products in construction, such as reinforcing and
  structural steel." A performance scale with **seven bands**. "The scope of the carbon calculation
  includes the **quarrying, mining and transport of raw materials, and the manufacture of hot rolled
  steel. Coatings, fabrication, transport to site and assembly are excluded** from the system boundary
  considered. It is worth noting that **offsetting is not permitted**."

### (h) Producer self-definitions the JRC catalogues (relevant to "green" claims)

- **SSAB Fossil-free**: "based on iron-ore produced without fossil fuels, using HYBRIT® technology,
  which uses hydrogen instead of coal in the ore reduction process".
  **SSAB Zero**: "made of recycled steel and produced with fossil-free electricity and biogas".
- **H2 Green Steel / Stegra**: "green steel must be produced from a combination of a significant
  amount of green virgin iron and scrap in a production process that **uses electricity from
  renewable energy sources**."
- **Kloeckner Metals**: "only steel with the lowest possible carbon emissions should be called
  'green'… meaning steel produced using direct reduced iron + green hydrogen or an electric arc
  furnace + 100% scrap and renewable electricity… Other steel with lower emissions is 'carbon
  reduced.'"
- **Wang et al. (2023), on the Australian iron ore industry**: "green steel (i.e., steel produced
  using hydrogen from renewable sources as the reducing agent)."

The JRC's own comment: "A common thread among all these definitions is their focus on **carbon
intensity**… While some definitions emphasise the technology routes employed for steelmaking,
others question whether a product can be classified as 'green' if it utilises fossil feedstocks…"

### (i) JRC conclusions verbatim (the parts that matter)

> "**Uncertain role of circularity**: While the incorporation of scrap content is clearly identified
> as a key lever… the 'neutralisation' of the use of scrap in a sliding scale or similar approach
> could present a challenge for the EU steel industry as the same score would cover different levels
> of circularity and the associated CO2 footprints. In addition, with regard to the labelling
> schemes, clear methodological guidelines on scrap accounting and categorisation might be needed to
> **prevent false green claims and misleading reporting**…"

> "**Lack of consistency in terminology**: The use of terms such as 'green steel,' 'low-carbon
> emissions steel,' 'low-emissions steel,' and 'near-zero emissions steel' is inconsistent across
> initiatives, with no universally agreed definitions."

> "**Long-term convergence of thresholds**: Despite methodological differences, there is a clear
> convergence in long-term emissions intensity thresholds/targets across initiatives…"

---

## 8. Industrial Decarbonisation Accelerator Act, lead markets and public procurement criteria

### (a) EU level

The IDAA is the vehicle for the voluntary carbon-intensity label (Section 6). Its lead-market and
procurement content is still framework-level; **no numeric steel criterion has been adopted at EU
level.** The two places numbers will land are (i) the ESPR delegated act's classes (Section 6d), whose
B/C boundary the JRC explicitly designed as a potential **Green Public Procurement** criterion, and
(ii) the revision of the Public Procurement Directives.

The JRC's own framing of that link:
> "The rationale for the 30% criterion is linked to the potential application of the **B/C threshold in
> Green Public Procurement (GPP)**, which defines the boundary between the two highest populated
> classes of performance. If applied as a procurement criterion, it should be designed in a manner that
> **does not unduly restrict competition**."

So the operative EU procurement number, if adopted as drafted, would be **HRC ≥ Class B, i.e. ≤ 2.66
t CO2eq/t HRC**, and **WR ≥ Class B, i.e. ≤ 2.43 t CO2eq/t WR** — thresholds deliberately calibrated so
that ~30% of global production volume already complies.

### (b) IDDI Green Public Procurement Pledge — four Levels, no numbers

Signatories: **Austria, Canada, Germany, Japan, UAE, UK, USA.** (The domain
`industrialenergyaccelerator.org` no longer resolves; text below is from the archived pledge:
https://web.archive.org/web/20221101143302/https://www.industrialenergyaccelerator.org/the-gpp-pledge/)

> **Level One:** "Starting no later than **2025**, require **disclosure** of the embodied carbon in
> cement/concrete and steel procured for public construction projects."
> **Level Two:** "Starting no later than **2030**, conduct **whole project life cycle assessments** …
> and, by **2050**, achieve **net zero emissions** in all public construction projects."
> **Level Three:** "Starting no later than **2030**, require procurement of **low emission**
> cement/concrete and steel … applying the **highest ambition possible under national circumstances**."
> **Level Four:** "Starting in **2030**, require procurement of a share of cement and/or **crude steel
> from near zero emission material production for signature projects**."

**No numeric thresholds** — "low emission" and "near zero emission" are undefined in the pledge and
defer in practice to the IEA definition.

### (c) SteelZero (Climate Group)

- **2030**: "procure, specify or stock **lower emission steel for 50% of their steel requirement by
  2030**", qualifying via either (1) a steelmaker with a science-based target, or (2) "'lower emission
  steel' (**aligning with ResponsibleSteel Decarbonisation Progress Level 2**)" — i.e.
  **y ≤ 2.00 − 1.65·x t CO2e/t crude steel**.
- **2050**: 100% net zero steel.
- 40+ businesses, ~10 Mt committed. There is **no separate SteelZero "near zero 2030" intensity
  number**; the milestone is the 50% / PL2 formulation.

### (d) Global Steel Climate Council (GSCC) Steel Climate Standard

Included here because it is the main *rival* design to the sliding scale and the JRC compared it.

Product thresholds, **t CO2e / t hot-rolled steel**: 2026 long **1.29** / flat **1.57**; 2030
**1.11 / 1.31**; 2040 0.63 / 0.69; 2050 **0.12 / 0.12**. Company glidepath: 2026 **1.42**, 2030
**1.20**, 2050 **0.12**.

Boundary "from mining to hot rolling", explicitly **including scrap collection and processing**. All
GHGs. Scope 2 is **market-based by default** (see the table in §5(f)). It rejects the sliding scale
outright:
> "Creating a standard based on a ferrous scrap sliding scale could result in **higher emissions steel
> products being labeled as 'green'** … We must avoid a standard that enables greenwashing."

### (e) SBTi Steel Science-Based Target-Setting Guidance (v1.0, July 2023)

Corporate-level, two separate pathways, **t CO2 / t hot-rolled steel**:

| Route | 2020 | 2030 | 2040 | 2050 |
|---|---|---|---|---|
| 100% ore-based | 2.42 | **1.71** | 0.77 | 0.11 |
| 100% scrap-based | 0.50 | **0.37** | 0.24 | 0.11 |

A company's pathway is the **scrap-weighted average** of the two endpoints. Scope 1 + scope 2 for EAF
power inside the iron & steel SDA boundary; **purchased HBI counts inside**; purchased power for cold
rolling and purchased ferroalloys fall outside. Mandatory scope 3 category 3 target. SBTi deliberately
avoids the phrase "sliding scale" — it is "used for setting **product-level** standards, which is not
the goal of the SBTi". Stakeholder survey open to **11 September 2026**.

---

## 9. Member-state and national procurement / funding criteria with explicit numbers

### (a) Germany — Klimaschutzverträge (carbon contracts for difference)

**Förderrichtlinie (FRL KSV), BAnz AT 10.04.2024 B1** —
https://www.bundesanzeiger.de/pub/publication/z26rd4pzhmkLYKeETr6/content/z26rd4pzhmkLYKeETr6/BAnz%20AT%2010.04.2024%20B1.pdf

Round 1 criteria (Nr. 4.15, verbatim):
> "Sie beträgt **mindestens 10 kt CO2-Äquivalente pro Kalenderjahr**"
> "**Spätestens ab dem dritten vollständigen Kalenderjahr** … muss die relative
> Treibhausgasemissionsminderung gegenüber dem Referenzsystem **mindestens 60 %** betragen"
> "Eine relative Treibhausgasemissionsminderung von **mindestens 90 %** … **in den letzten zwölf Monaten
> der Laufzeit** … (**Zugangskriterium Klimaneutralität**)"

Projects under €15 m total funding ineligible; term 15 years.

**Round 2** — Vorverfahren 6 Oct – 1 Dec 2025; amended FRL cleared by the Commission 24 March 2025;
round approved **7 May 2026 (SA.122980)**; Förderaufruf published 5 May 2026, revised 3 Sept 2026,
deadline **7 September 2026**. Terms loosened:
- Volume **€5 bn** (€3 bn basic + €2 bn additional); max basic project €700 m; projects > €2.5 bn ineligible
- Minimum project size **10 kt → 5 kt CO2e/yr**
- Max bid price **€600/t → €550/t CO2e**
- Milestones relaxed to **≥50% from the 4th full calendar year and ≥85% in the final year** (was 60%/90%)
- Latest start of operation 1 Jan 2031; CCU/CCS now eligible

Förderaufruf (EN):
https://www.co2-differenzvertraege.info/lw_resource/datapool/systemfiles/agent/ewbpublications/e9bcbbc2-7ac7-11f1-9305-fa163ed847d2/live/document/CCfD_GV2026_F%C3%B6rderaufruf_260903_EN.pdf

**Steel reference systems:**

| Reference system | Specific GHG | Product |
|---|---|---|
| **2-4 Primary steel** | **1.321 t CO2-eq/t** | "tonnes of **liquid pig iron**"; the funded product may be **pig iron (DRI alone) or crude steel** |
| 5 EAF carbon steel | 0.050 t CO2-eq/t | secondary crude steel |
| 6 EAF high-alloy steel | 0.103 t CO2-eq/t | secondary crude steel |

The primary-steel reference assumes a **20% scrap share**, with per-tonne energy carriers: electricity
0.10 MWh, natural gas 0.67 MWh, coking coal 2.83 MWh, steam coal 0.86 MWh.

**Implied absolute ceilings** for primary steel: **≤0.661 t CO2e/t** (50% cut) and **≤0.198 t CO2e/t**
(85% cut) under the 2026 round; ≤0.528 and ≤0.132 under 2024 terms.

**Scope 2: EXCLUDED OUTRIGHT.** FRL 7.1(e): "Die Treibhausgasemissionen des Vorhabens ergeben sich aus
den Treibhausgasemissionen der geförderten Anlagen (**Scope-1-Emissionen**)". The Förderaufruf states
the ETS benchmark values from Implementing Reg. (EU) 2021/447 "were **reduced by the indirect emissions
for electricity not to be taken into account in this funding programme**". Electricity enters only as a
*cost* term (the Dynamisierung is priced off SMARD day-ahead). **So a hydrogen-DRI plant on ordinary
grid power books zero charge emissions for its electricity under CCfD.** The location/market-based
question does not arise.

**Note the boundary**: this is the one instrument in the whole survey that names **liquid pig iron /
DRI alone** as a fundable product with its own reference value. It is the closest thing to a numeric
"green iron" threshold in force anywhere — but it is a *funding* criterion, not a label, and it counts
Scope 1 only.

15 first-round contracts were concluded (paper, sugar, chemicals); **no steel project appears among
them** — the large steel plants went through separate individually notified aid. *(Indicated, not
fully verified.)*

### (b) Germany — individual steel funding decisions, 2026 status

- **thyssenkrupp tkH2Steel**: ~€2 bn public (federal + NRW) of ~€3 bn; 2.5 Mt DRI/yr; "up to 3.5 Mt CO2
  per year" avoided; DR tower assembly began Feb 2026; startup 2027. **The hydrogen condition is being
  unwound.** 18 Aug 2026: thyssenkrupp "in advanced discussions to modify funding conditions"; CFO Axel
  Hamann says the **European Commission has approved the proposed modification as state-aid
  compatible**. 1 Sept 2026: reported as a concluded "DRI-Kompromiss mit Berlin und Brüssel" allowing
  the plant to **run initially on natural gas**. *(No formal amended Zuwendungsbescheid found.)*
  https://www.marketsteel.de/news-details/thyssenkrupp-koennte-dri-anlage-zunaechst-ohne-wasserstoff-betreiben.html
- **Salzgitter SALCOS**: proceeding. 9 July 2026 Salzgitter took 100% of HKM; 29 July 2026 "startet
  Transformation der HKM mit Deutschlands größtem Elektrolichtbogenofen". **Peiner Träger's published
  EAF product number: 299 kg CO2e/t** (100% scrap + renewable electricity).
- **ArcelorMittal Bremen + Eisenhüttenstadt — WITHDRAWN.** Announced 19/20 June 2025; the **~€1.3 bn
  federal grant (plus ~€251 m from Land Bremen) was returned**. Reiner Blaschek: "**Wir haben uns
  entschieden, das gemeinsame Dekarbonisierungsprojekt der ArcelorMittal-Standorte Bremen und
  Eisenhüttenstadt nicht umzusetzen.**" Reasons cited: uncompetitive electricity/energy prices,
  insufficient competitively priced green hydrogen, **and the absence of green lead markets**. Original
  terms: ~€2.5 bn total investment, two blast furnaces replaced by EAF+DRI by 2030, up to 5.8 Mt CO2/yr
  avoided, 3.8 Mt/yr of CO2-reduced steel.
- **Stahl-Holding-Saar "Power4Steel"**: one DRI plant plus **two EAFs at Dillingen and Völklingen**;
  per BT-Drs. 21/3563, "**2,6 Mrd. Euro an öffentlicher Förderung**" plus a **€1.7 bn project financing
  finalised October 2025**. Per-site capacities, emission targets and hydrogen contracts are classified
  **VS-Vertraulich** and deposited in the Bundestag's Geheimschutzstelle — **no public numeric
  CO2-intensity table exists** for these projects.

### (c) Germany — the "Grüner Stahl" definition (BMWK Leitmärkte concept, 22 May 2024)

https://www.bundeswirtschaftsministerium.de/Redaktion/DE/Publikationen/Klimaschutz/leitmaerkte-fuer-klimafreundliche-grundstoffe.pdf

> "Die Definition bezieht sich auf **warmgewalzten Stahl**. Sie berücksichtigt alle relevanten
> (**über 90 %**) Quellen von Treibhausgasemissionen"
> "Der '**cradle-to-gate**' Bilanzrahmen umfasst neben den direkten THG-Emissionen der produzierenden
> Anlage (**Scope 1**), die THG-Emissionen aus der **Energiebereitstellung (Scope 2)** sowie indirekte
> THG-Emissionen der Vorkette … Eisenerz sowie weiteren eingesetzten Stoffen wie Kalk oder
> Legierungsmittel (**Scope 3**)"
> "Neben der Kategorie '**near zero**' sind vier weitere Kategorien von **A bis D**"
> "**Für Baustahl zusätzlich 70 kg CO2äq/t und für Qualitätsstahl 120 kg CO2äq/t**" [the surcharge over IEA]

Values: **Qualitätsstahl** near-zero **0–520 kg CO2äq/t**, emissionsarm in four steps up to
**2 600 kg CO2äq/t**; **Baustahl** near-zero **0–120 kg CO2äq/t**, emissionsarm up to **600 kg
CO2äq/t**. *(Caution: 520 is the QST 0%-scrap extrapolation and 120 is the BST 100%-scrap value —
these are not like-for-like. Use the full LESS table in §2(c).)*

Illustrative EU ladder proposed: **class D or better from 2030; C from 2035; A from 2040; Near Zero
from 2045.** **Scope 2 method is not stated** — it defers to the WV Stahl rulebook, i.e. LESS, which
permits GOs.

### (d) Germany — is there a Leitmarktgesetz? An empowerment with no number

No standalone act. The hook is a clause inserted into the GWB by the **Gesetz zur Beschleunigung der
Vergabe öffentlicher Aufträge** (passed 12 May 2026; **BGBl. 2026 I Nr. 137, 18 May 2026; in force
1 July 2026**) — https://www.recht.bund.de/bgbl/1/2026/137/regelungstext.pdf

**§ 113 Abs. 1 Nr. 9 GWB (new), verbatim:**
> "**verpflichtender Anforderungen an die Klimafreundlichkeit bei der Beschaffung von Leistungen,
> insbesondere hinsichtlich der Verwendung von emissionsarmen Grundstoffen wie Stahl und Zement.**"

**No number.** The implementing Rechtsverordnung must be laid before the Bundestag before the
Bundesrat.

**Government's own position, BT-Drs. 21/3563 (8 January 2026)** —
https://dserver.bundestag.de/btd/21/035/2103563.pdf
> "Für Stahl und Zement liegen bereits mit den privaten Label-Initiativen **Low Emission Steel Standard
> (LESS)** und Cement Carbon Class (CCC) solche Grundlagen vor. … Das **BMWE wird die Umsetzung zeitnah
> im Jahr 2026 anstoßen.**"
> "**Es sind keine Quotenvorgaben für klimafreundlichen Stahl und Zement in der öffentlichen Beschaffung
> vorgesehen.**"
> "**Aktuell ist kein Förderprogramm seitens BMWE für die Kompensation von Mehrkosten klimafreundlicher
> Grundstoffe geplant.**"

Green premium at Deutsche Bahn "etwa **10 bis 15 Prozent**"; whole-product cost increase from H2/NG-DRI
primary steel "zwischen **0,04 Prozent und 0,31 Prozent bei Automobilen** sowie **0,5 Prozent bis 4,0
Prozent bei Nichtwohngebäuden**". *(The "25% emission-reduced steel from 2026" figure sometimes quoted
is a Germanwatch demand, not policy.)*

Other German instruments, none with a steel number: **AVV Klima (2021)** — lifecycle GHG estimate plus a
**CO2-Schattenpreis** in the award decision; **§ 15(1) KSG** — federal administration climate-neutral by
2030; **QNG / BNB** — building-level lifecycle reporting. **Deutsche Bahn**, 5 Nov 2025: DB InfraGO +
Saarstahl Rail, **~1 000 t rails ≈ 22 track-km**, made at **Ascoval (France)** EAF, "**bis zu 70 Prozent
weniger CO2**" — relative only, no threshold. DB ramp-up targets (LV InfraGO, July 2025): 2026 50 km
(~5 800 t); 2027 100 km (~11 600 t); 2028-30 500 km/yr (~58 000 t/yr) — volume targets, not binding, no
CO2 threshold. DB baseline demand for rails and switches is ~270 000 t/yr (2020-25).

**Structural point.** Germany's two instruments use *incompatible* boundaries and reach the same
practical outcome by opposite routes: **CCfD pays on Scope 1 only, with electricity emissions
deliberately subtracted out** — grid power is free; **LESS/Leitmärkte classify on cradle-to-gate Scope
1+2+3-upstream, but Scope 2 is market-based via registry GOs** — so contractual green power zeroes the
electricity anyway. **In both German regimes, electricity carbon can be discharged contractually rather
than physically.**

### (e) Other member states — no numeric steel thresholds anywhere

**France.** ArcelorMittal Dunkerque **€850 m** (approved 20 July 2023). Conditions are **absolute, not
intensity-based**: "un potentiel de réduction des émissions de CO2 d'**au moins 4,4 millions de tonnes
par an**" and "environ **70 millions de tonnes** … sur le cycle de vie du projet, d'une durée de **15
ans**". **No tCO2/t steel condition.** ArcelorMittal has since paused/relaunched consultation on the
Dunkerque DRI project. **RE2020** sets *building-level* carbon limits (Ic_construction in kgCO2e/m²),
not per tonne of steel. No "acier bas-carbone" décret with a tCO2/t number.

**Netherlands.** Tata Steel IJmuiden — **Joint Letter of Intent signed 29 September 2025**: "**circa 5,4
Mton CO2 per jaar**", rising "**tot circa 7,2 megaton**"; total investment **€4-6.5 bn**, state
contribution "**maximaal € 2 miljard**", payable only after a final binding agreement; particulates
~−35%, lead −68%, nitrogen −44%. **No intensity condition.** The definitive maatwerkafspraak is not
concluded (Lodewijk Asscher appointed government envoy). **MKI/DuboCalc** is a *monetised* environmental
cost indicator used as an award criterion, not a tCO2/t threshold; the **"Sturende MKI" bill submitted
29 July 2026** mandates minimum environmental performance for **concrete and asphalt only — steel gets
no numeric limits**.

**Sweden.** Upphandlingsmyndigheten confirms it has **no steel criteria yet** — concrete criteria are in
development and "**efter det står stålprodukter på tur**". Trafikverket's klimatkrav require
project-specific emission factors for steel but **no per-tonne limit**. Its updated 2026 building
criteria carry numbers, but at **building level** ("siffersatta värden för maximal klimatpåverkan i CO2e
per m² BTA") in three tiers (Bas / Avancerad / Spjutspets). No national green-steel procurement rule
attached to HYBRIT or Stegra.

**Spain.** €460 m approved for ArcelorMittal (Gijón/Sestao), IP/23/849 — **conditions could not be
read** (PDF would not decode). ArcelorMittal has **suspended the Gijón green-hydrogen DRI plan** and
will start on natural gas while retaining the PERTE claim. **PERTE numeric conditions: gap.**

**Italy.** Acciaierie d'Italia in extraordinary administration; a ~€1 bn DRI d'Italia agreement for
Taranto, with those funds **redirected** by the government (July 2026). "Green plants from 2033"; total
decarbonisation cost ~€10 bn. **No tCO2/t condition found.**

### (f) United Kingdom

**UK Steel Strategy, published March 2026** —
https://assets.publishing.service.gov.uk/media/69bbd096f7b1c24d8e23ce06/uk-steel-strategy.pdf
**No numeric carbon-intensity procurement standard.** Procurement policy is **origin**-focused: PPN 022
(June 2025) requires in-scope buyers to "consult UK Steel's digital catalogue"; Wales has WPPN 008.
Green procurement is tied to IDDI: "The proposals will help the UK to fulfil its COP28 commitment to
the Industrial Deep Decarbonisation Initiative's (IDDI) Green Public Procurement pledge."

**UK CBAM — indirect emissions delayed, verbatim:**
> "The government has legislated in Finance Bill 2025-26 to introduce the CBAM from **1 January 2027**.
> The **inclusion of indirect emissions within scope of the UK CBAM will be delayed until 2029 at the
> earliest.** This is to reflect continued support for the Energy Intensive Industries (EII)
> Compensation Scheme."

"The CBAM charge will be based on the **direct emissions** embodied within the imported CBAM goods."
Registration threshold £50 000. Draft secondary legislation first published 10 February 2026 for a
6-week technical consultation. **So the UK, like the EU, is blind to the grid the iron was made on —
and has now put a date on when it might stop being so (2029 at the earliest).**

**UK low-carbon products framework — decision deferred.** Consultation 23 June – 29 Sept 2025 (109 valid
responses); **government response 28 May 2026** —
https://assets.publishing.service.gov.uk/media/6a18151a916cd732dcdaac4a/low-carbon-products-government-response.pdf
- Concept was A-G classes with procurement commitments "(e.g., **Class D steel by 2030**)" — **no class
  boundaries set**
- The **EPA model was dropped**: "Given the minimal support on the suitability of the U.S. Environmental
  Protection Agency's (EPA) approach to setting limits for low embodied carbon steel, the government
  will no longer consider this model as a potential option." Will explore **CARES Sustainable
  Construction Steel (SCS)** and the EU **ESPR** classifications instead.
- The pivotal question, verbatim: "The government notes the **lack of clear consensus on whether to use
  a steel product classification that applies a scrap sliding scale**. Determining whether to use a
  sliding scale approach is a **priority** for government … The government will therefore assess the
  suitability of a sliding scale approach as a first step."

**Port Talbot EAF**: £500 m Grant Funding Agreement within a £1.25 bn package; planning approved Feb
2025. **No CO2-intensity condition found.**

### (g) North America — where the hard numbers actually are

**US GSA, IRA Low Embodied Carbon Steel Requirements (Dec 2023)**, kg CO2e per metric tonne. (The PDF
has been pulled from gsa.gov; archived at
https://web.archive.org/web/20250807100014/https://www.gsa.gov/system/files/Steel%20-%20GSA%20IRA%20Low%20Embodied%20Carbon%20Requirements%20%28Dec.%202023%29_508.pdf)

| Category | Top 20% | Top 40% | Better Than Average |
|---|---|---|---|
| Rebar (fabricated) | 728 | 794 | 850 |
| Rebar (unfabricated) | 611 | 716 | 760 |
| Hollow Structural Sections (fabricated) | 1 778 | 1 854 | 1 898 |
| HSS from **Electric Arc Furnaces** (unfabricated) | 1 580 | 1 620 | 1 652 |
| HSS from **Integrated Mills** (unfabricated) | TBD | TBD | TBD |
| Hot-Rolled Sections (fabricated) | 1 022 | 1 128 | 1 163 |
| Hot-Rolled Sections (unfabricated) | 686 | 713 | 869 |
| Cold-Formed and Galvanized | 2 228 | 2 324 | 2 408 |
| Structural Plate from **EAF** (unfabricated) | 987 | 1 152 | 1 190 |
| Structural Plate from **Integrated Mills** (unfabricated) | TBD | TBD | TBD |

**No scrap sliding scale — GSA segregates by ROUTE instead**, setting EAF limits and leaving
integrated-mill limits "TBD", and additionally requiring integrated-mill steel to come from a mill with
an **ENERGY STAR Energy Performance Score**. Assemblies qualify if ≥80% of cost or weight is compliant
steel. Scope 2 sits inside A1-A3 by construction but **location vs market is unresolved** — GSA defers
to the UL or SCS steel PCRs and ISO 14025/21930.
**Status: defunded going forward.** §60021 of the One Big Beautiful Bill Act 2025 (P.L. 119-21)
rescinded unobligated IRA §60503 funds; already-obligated money (~$1.6 bn as of Jan 2025) runs to
30 Sept 2026.

**US FHWA §60506**: GWP thresholds for concrete, asphalt, glass and steel were published Dec 2024, but
the steel page now returns 403/404 and the only Wayback capture archived a server error. **The FHWA
steel numbers could not be obtained — genuine gap.** Unobligated funds also rescinded by the OBBBA;
$1.2 bn had been awarded to 39 state DOTs in Nov 2024.

**US EPA IRA label programme**: §60112 is the EPD *assistance* programme; the **label** sits under
**§60116**. EPA's Label Program Approach and PCR Criteria (both Aug 2024) are still online; **2026
operational status not verified.** The UK has dropped this model (above).

**Buy Clean California**: effective **1 January 2025**. Concrete reinforcing steel limits were **890
(unfabricated) / 920 (fabricated)** effective 1 Jan 2022 and were **revised down to 755 / 778**. A
compliance-report update was issued 23 April 2026 for PCR changes. No scrap sliding scale, no EAF/BF
split, market vs location not specified.

**Colorado — binding caps as of January 2026.** Document EE-5.1, 01/2026:
https://osa.colorado.gov/sites/osa/files/documents/OSA%20Max.%20GWP%20Limits%202026-01.pdf
kg CO2e/t: fabricated rebar **1 030**; hot-rolled **1 220**; plate **1 730**; HSS **1 990**; cold-formed
framing **2 843**; roof/floor deck **2 350**; open-web joist/joist girders **1 450**; post-tension steel
"**No sufficient data to set a valid threshold at this time**". A1-A3; projects over **$500 000** from
1 Jan 2024; separate A4 transport reporting for products sourced >100 miles. Statutory ratchet: "By
January 1, 2026, and every four years thereafter, the OSA shall review … **The OSA shall not adjust the
number upward** for any eligible material."

**Minnesota**: https://mn.gov/admin/assets/2026-MNBuyCleanGWPLimits_2026-01-15_tcm36-724796.pdf
(15 Jan 2026, updated 12 Feb 2026). **Rebar in buildings: 0.755 t CO2e/t** cradle-to-mill-gate;
fabricated-equivalent **0.854**. Applies to projects advertised on or after **15 July 2026**. All other
steel products are **EPD-disclosure only; limits to be implemented in 2028**. The fabricated conversion
assumes "a scrap rate of 3.0% (A1 multiplier), transportation (A2) of 0.0490 tons/ton, and fabrication
(A3) of 0.0270 tons/ton". The statute defines "integrated" and "secondary" steel production but the
definitions do not yet drive differentiated limits.

Other states: Washington — EPD/reporting only. Oregon — "working toward benchmarks", no numbers.
New York — **concrete only; steel is not covered**. Maryland — cement/concrete only. New Jersey — no
Buy Clean steel GWP limits found (gap).

**Canada — the most interesting design in the whole survey, because the scrap adjustment runs
backwards.** NRC, *National tiered greenhouse gas emissions limits for steel construction products*,
published **26 Aug 2025**, DOI 10.4224/34dr-wm27 —
https://nrc-publications.canada.ca/eng/view/author/version/?id=e33c6800-6f6a-4fc8-b0e4-a0e10113c4ba

kg CO2e/t:

| Product | Best 20% | Best 40% | Best 50% |
|---|---|---|---|
| Rebar | 703 | 773 | 793 |
| Rebar (unfabricated) | 622 | 684 | 702 |
| Hot-rolled sections | 1 008 | 1 105 | 1 140 |
| Hot-rolled sections (unfabricated) | 820 | 898 | 927 |
| Hollow structural sections | 1 542 | 1 810 | 1 908 |
| HSS (unfabricated) | 1 329 | 1 560 | 1 645 |
| Structural plate | 1 476 | 1 646 | 1 709 |
| Structural plate (unfabricated) | 1 262 | 1 407 | 1 461 |
| Hot-rolled sheet, strip, plate | 1 158 | 1 401 | 1 505 |
| Cold-rolled sheet, strip, plate | 1 848 | 2 021 | 2 102 |

Carbon and low-alloy steel only. Procurement is a cascade: Best-20% unless not reasonably sourceable →
Best-40% → Best-50% → exemption.

**The primary steel adjustment factor — an inverted sliding scale, verbatim:**
> "**2.1. Primary steel adjustment factor** … an adjustment factor that reflects the difference between
> best-in-class primary and secondary steel production of **up to 585 kgCO2e/t** can be deducted from
> the embodied carbon of that product **in proportion to its primary steel content**. This deduction
> acts as a **primary steel credit**, enabling products that contain cleaner primary steel to be
> compliant with EC limits."

`EC_with_credit = EC_whole_product − (F_primary × 585)`, with `F_primary = 1 − F_recycled`, EC taken from
A1-A3 in the EPD. The 585 figure derives from Argonne National Laboratory DRI modelling (Zang et al.,
2023). And it is ring-fenced: "The credit … **shall not be used in any embodied carbon disclosure** of
the project's structural materials or as part of whole-building life cycle assessments."

This is worth dwelling on for the export case. Rather than *loosening the threshold* as scrap falls
(IEA / ResponsibleSteel / LESS), Canada *subtracts a credit* proportional to primary content — same
intent, opposite mechanics, and quarantined so it cannot leak into reported footprints. **A
green-HBI-fed EAF product would take the full 585 kg/t credit on its primary fraction.**

---

## 10. First Movers Coalition — near-zero steel commitment

Primary sources:
- FMC Steel commitment, 2025 revision:
  https://reports.weforum.org/docs/WEF_First_Movers_Coalition_Steel_Commitment.pdf
- FMC Steel commitment, 2022 original (superseded):
  https://www3.weforum.org/docs/WEF_FMC_Steel_2022.pdf
- FMC Impact Brief 2026:
  https://reports.weforum.org/docs/WEF_First_Movers_Coalition_Impact_Brief_2026.pdf
- https://initiatives.weforum.org/first-movers-coalition/steel (403 to automated fetch)

### (a) The commitment and threshold

> **"At least 10% (by volume) of all our steel purchased per year will be near-zero emissions
> (as per FMC definition) by 2030"**

> "Crude steel from breakthrough technology production facilities. Per IEA guidance, the steel
> should emit **<0.4 (0% scrap inputs) to <0.05 t (100% scrap inputs) of CO2e per tonne of crude
> steel produced**"

Scrap-adjusted sliding scale: threshold = **0.4 − 0.35·s** tCO2e/t crude steel, `s` = scrap fraction.

**The threshold was tightened in the 2025 review.** The 2022 text read "<0.4 t (with 0% scrap
inputs) to **<0.1 t** (with 100% scrap inputs) of **CO2**". The 100%-scrap end moved from 0.1 to
**0.05**, and the metric from CO2 to **CO2e**. Many secondary sources still quote 0.1 — use 0.05.
Document labelled "Commitment reviewed in 2025", revised "through a biennial Commitment Review
process".

### (b) Boundary — scope 2 IS included

Cradle-to-crude-steel, aligned to IEA (2022). Footnote 1:
> "FMC sets ambitious standards, including a supply chain boundary inclusive of **raw material
> preparation (iron ore and limestone) and fossil fuel supply (including extraction,
> transportation, and beneficiation) through steelmaking and casting** (including all iron ore and
> limestone processing transportation emissions; **does not include sorting and transportation of
> steel scrap**). Transport emissions of iron ore and lime products include all emissions
> regardless of intermediary stops between mining and steel plant."

Also: "Lime-based slag formers are included within the IEA boundary and therefore should be
included in the calculation of emissions." Covers **both flat and long steel**.

Three points for the model:
- **Scrap is boundary-free upstream** — sorting and transport of scrap explicitly excluded. A real
  thumb on the scale for scrap-heavy routes beyond the sliding scale itself.
- **Upstream fossil fuel supply chain is in** — methane from gas extraction counts, which matters
  for NG-DRI.
- Under the IEA definition, indirect (electricity/heat) emissions are approximated at **regional
  average grid intensity** — IEA's own default is essentially *location-based*.

### (c) Chain of custody — and the electricity carve-out

**For the steel itself**: identity preserved, segregated, and (with registry tracking and
third-party verification) controlled blending and site-level mass balance. Explicitly excluded:
"Group-level mass balance and **book & claim** chain of custody models are not applicable to the
FMC steel commitment at this time."

**But for electricity, book & claim survives.** Footnote 8:
> "Exceptions apply to the use of **book and claim market-based approach to carbon accounting for
> energy attribute certificates required to prove the origin of purchased energy per Greenhouse Gas
> (GHG) Protocol Scope 2 Guidance**"

So **FMC defers to the GHG Protocol for electricity** — whatever the Scope 2 Standard permits, FMC
permits. Section 5 above flows straight through to FMC compliance.

Note the 2022 version had a *stricter* condition the 2025 version drops: "FMC permits the use of
virtual PPA to satisfy Scope 2 emission thresholds, **if additionality is confirmed by independent
expert third party**." That sentence is gone; the 2025 text points generically at GHG Protocol
Scope 2, which contains no additionality requirement. On its face a **loosening on the electricity
side at the same time as the headline threshold tightened**.

### (d) Status (2025–26)

WEF still runs it. FMC Impact Brief 2026, end-2025 figures:
- **101 members** (from 35 founding members at COP26 Glasgow, 2021)
- **14 government partners** (Spain most recent)
- **>$19 billion** in committed green product demand by 2030
- **~26 Mt CO2** estimated annual emissions reductions in 2030
- **>130 offtake agreements and investments** signed
- 7 sectors: aviation, shipping, trucking, aluminium, cement/concrete, steel, CDR

Steel-relevant 2025–26 activity: SSAB Montpelier, Iowa now meets both IEA and FMC near-zero
thresholds by adding HYBRIT hydrogen-reduced iron (material contracted to GE Vernova wind towers);
Ecolab joined a **€60 million** round for **GravitHy** (H2-DRI) via the Japan Hydrogen Fund with Rio
Tinto, Siemens Financial Services and Engie New Ventures; FMC convened **workshops on green iron and
near-zero materials in Australia** and met in China for the first time. The Impact Brief flags the
problem directly: "Clear, harmonized, industry definitions of near-zero steel are also needed to
accelerate action and build scale", and lists "limited availability of high-quality scrap" as a
binding constraint on reaching the low end of the sliding scale.

---

## 11. The asymmetry nobody in the steel standards has fixed: RFNBO hydrogen rules vs Scope 2 rules

Not on the original list, but it is the single sharpest place where a **hydrogen-DRI plant is
treated differently from a molten-oxide-electrolysis plant** — and it comes from EU energy law,
not from any steel standard.

**Commission Delegated Regulation (EU) 2023/1184** sets the conditions under which hydrogen counts
as a renewable fuel of non-biological origin (RFNBO):
https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R1184

- **Additionality (Art. 5(a))** — the renewable electricity installation must have "come into
  operation not earlier than **36 months** before the installation producing the [RFNBO]".
  **Art. 11** transitional: installations "that come into operation before **1 January 2028**" are
  exempt from Art. 5(a)-(b) until **1 January 2038** — not applicable to capacity added after
  1 January 2028.
- **Temporal correlation (Art. 6)** — until **31 December 2029**, production must occur "during the
  same **calendar month** as the renewable electricity produced"; from **1 January 2030**, "during
  the same **one-hour period** as the renewable electricity produced". Waived where "the clearing
  price of electricity … is lower or equal to **EUR 20 per MWh**".
- **Geographic correlation (Art. 7(1))** — the generator must be "located … in the **same bidding
  zone** as the electrolyser" (a), or in an interconnected bidding zone where prices are "equal or
  higher than in the bidding zone" of the producer (b), or in "an offshore bidding zone that is
  interconnected" (c).

**Why this matters here.** Additionality + hourly + bidding-zone is exactly the package the GHG
Protocol Scope 2 revision is *proposing* (minus additionality) and does not yet have. So today:

- A **hydrogen-DRI** plant whose hydrogen must qualify as RFNBO (for EU regulatory purposes, EU
  funding, or an EU customer's RED III industry target) already faces hourly matching, bidding-zone
  correlation and a 36-month additionality test on its electricity.
- A **molten-oxide-electrolysis** or **aqueous-electrowinning** plant, whose entire decarbonisation
  claim is a Scope 2 claim, faces **only** the 2015 GHG Protocol criteria — "as close as possible"
  (annual) and "same market" (national). No additionality, no hourly, no zone.

For a *non-EU* H2-DRI plant exporting HBI, RFNBO rules do not bind at all unless the hydrogen (or
the resulting fuel) is imported into the EU; the HBI's embodied hydrogen is accounted under the
steel standard's Scope 3.1/3.3 rules, not under 2023/1184.

Also relevant: **ResponsibleSteel's default embodied GHG value for hydrogen is 13.916 kg CO2e/kg H2**
(Annex 5, Table A1; 2024.1 GaBi/Sphera) — a grey-hydrogen figure. At roughly 50-60 kg H2 per tonne
of DRI, using the default would load **~0.7-0.8 t CO2e per tonne of DRI** onto the account. An
H2-DRI producer that does not supply verified primary data for its hydrogen is therefore scored as
if it were running on SMR hydrogen. (Natural gas default: 1.064 t CO2e / kNm3, worldsteel CO2
methodology.)

---

## 12. Synthesis — direct answers to the three questions

### Q1. Is there ANY standard under which exported HBI made with coal-grid electricity would be disqualified?

**Yes — three, and they work by different mechanisms.**

**(1) ResponsibleSteel** is the strongest, because it forces the exporting HBI plant's own Scope 2
into the EU melter's number. Criterion 10.4.5.a requires the EU EAF to count the HBI's cradle-to-gate
**Scope 1 + Scope 2 + Scope 3**; 10.4.5.b makes supplier primary data mandatory where it exists and
forbids reverting to the default even if the default is lower; and 10.4.5.d makes the purchaser
responsible for the freight. So a coal-grid HBI plant cannot hide behind the fence line: either it
supplies verified data showing its high Scope 2 (and drags the EU EAF's number up), or no data is
supplied and the default 1.219 t CO2e/t HBI applies. Either way, the EU melter's Progress Level is
degraded.

*But note the loophole*: the coal-grid HBI plant can buy Australian LGCs or Brazilian I-RECs and
report a market-based Scope 2 of zero in the primary data it supplies. ResponsibleSteel expressly
permits "renewable energy certificates, power purchase agreements, virtual power purchase agreements,
or green tariffs" under ISO 14064-1 E.2.2, with no hourly, deliverability or additionality test. The
only defence is **disclosure**: 10.7.1.b(iii) requires the site to publish whether such instruments
were used and "a description of the source and quantity". So the standard makes the claim *visible*
rather than *impossible*.

**(2) worldsteel CO2 methodology / ISO 14404** disqualify coal-grid HBI in a completely different and
in some ways cleaner way: they use a **fixed world-average electricity factor of 0.504 t CO2/MWh** and
neither reward a clean grid nor accept a PPA. A coal-grid HBI plant and a hydro-powered one get the
*same* electricity charge. That makes the framework immune to greenwashing on electricity, but also
useless for rewarding genuine renewables. These are process-efficiency benchmarking tools, not
carbon-footprint tools, and the IEA says so: they "were not initially conceived with the net zero
transition in mind."

**(3) LESS** disqualifies via its **95% coverage rule** — an EU EAF cannot obtain a LESS class unless
at least 95% of its pig iron / DRI / HBI / scrap input is itself evaluated under the classification
system, and an unevaluated remainder must be valued at "the highest emission value" among the
warranties held. A non-EU HBI producer with no LESS evaluation therefore blocks the classification.
But LESS also permits Guarantees of Origin, and its Annex VI carries an explicit "Electricity
(renewable) — 0.0000" Scope 2 factor, so an *evaluated* coal-grid plant holding accepted GOs would
still report zero Scope 2. The friction for Australia/Brazil is that the registers LESS names are
**German** ("German Environment Agency Register", "German Biogas Register") and its declared
geographical area is the **EU** — a non-EU producer currently has no named register to point at.

**(4) The EU's own emerging ESPR class methodology is, unexpectedly, the strictest of all on this exact
point.** The JRC draft models electricity **only at country-average intensity** — the EU average (EEA)
for EU installations, and CBAM default values "which reflect the average electricity mix of the country
of origin" for imports — with a France/Germany/Poland member-state sensitivity, and **no mention of
PPAs, GOs or RECs anywhere in the document**. Under this methodology a Brazilian or Australian producer
is charged its country's grid factor and **certificates do not help**. This is the one place in the
survey where the answer to "can a PPA make grid electricity count as zero?" is currently *no*. It is
also, notably, the opposite of what LESS is lobbying the Commission for.

*Caveat*: the JRC text describes the electricity rule for the **Reference Installation Scenarios** used
to *calibrate* the thresholds. Whether the declared PCF of a real product would be allowed to use a
market-based factor is not stated. The document's silence on contractual instruments should be read as
unresolved rather than as an explicit prohibition — but the RIS calibration is country-average, which
sets the frame.

**Frameworks under which coal-grid HBI is essentially invisible:**
- **CBAM** — excludes indirect emissions for iron, steel and hydrogen (the premise of this work).
- **IEA near-zero, "direct emissions" sub-threshold** — imported electricity is marked N/A. A plant can
  claim "near zero **direct**" at 400 kg while running on coal. Only the "direct + indirect"
  sub-threshold sees the grid, and the two carry the same number, so which one is being claimed must
  be stated. The IEA itself has *not resolved* location vs market: "suitable methodologies for
  accounting for electricity emissions from the grid, including… use of renewable electricity
  credits… **Further discussion on this topic is likely to be needed.**"
- **First Movers Coalition** — bans book-and-claim for the steel but **expressly carves it back in for
  electricity** (footnote 8), and dropped the 2022 vPPA-additionality condition in the 2025 revision.
- **Climate Bonds Initiative** — no upstream Scope 3 at all, so purchased HBI's embodied emissions are
  outside the boundary entirely.
- **Chinese C2F (CMC2FS)** — Scope 3.1 only "partially" (energy preparation), no external scrap
  collection, no upstream transport.
- **GSCC Steel Climate Standard** — market-based is the *default* and location-based the last resort;
  the framework most permissive of a contractual zero.
- **Germany's Klimaschutzverträge** — Scope 2 excluded outright, with the ETS benchmarks explicitly
  "reduced by the indirect emissions for electricity not to be taken into account in this funding
  programme". A hydrogen-DRI plant on coal grid power books nothing for its electricity.
- **UK CBAM** — same blind spot as the EU, now with a date on it: "**the inclusion of indirect
  emissions within scope of the UK CBAM will be delayed until 2029 at the earliest**."
- **US GSA, Buy Clean California, Colorado, Minnesota, Canada NRC** — all defer to the underlying steel
  PCR, which does not resolve contractual instruments. Unresolved rather than permissive, but in
  practice an EPD prepared on a market-based basis would likely pass.

### Q2. Which standard best rewards genuinely renewable-powered iron?

**ResponsibleSteel, on balance — but only if the producer supplies verified primary data, and only
because its disclosure rule makes the difference between a real renewable supply and a certificate
purchase legible to the buyer.** Its combination of (i) full CO2e rather than CO2-only, (ii) mandatory
supplier primary data that cannot be overridden by a lower default, (iii) freight inside the boundary,
(iv) an outright ban on offsets, and (v) mandatory disclosure of REC/PPA source and quantity, is the
most complete package. A genuinely renewable-powered Australian or Brazilian HBI plant supplying
verified primary data can beat the 1.219 t CO2e/t default substantially and hand its EU customer a
better Progress Level; a certificate-only claim is legal but must be labelled as such.

**LESS rewards it most *legibly*, because its scale is calibrated so that electricity moves you two
whole classes.** The rulebook's own worked case: a 100%-scrap EAF QST reference plant on the German
grid sits at 805 kg CO2e/t rolled steel = **Stage D**; "If 50% renewable electricity is used instead of
the current German electricity mix, the reference plant could already reach level **C**. With 100%
green electricity level **B** is also achievable — but not level A." And its ideal H2-DRI-EAF case
reports **Scope 2 = 0** and 544 kg total = **Stage A**. So LESS is the standard that most visibly pays
for renewable power. Its weakness is that it accepts GOs for electricity "not physically delivered",
so the reward attaches to the certificate rather than to the physical supply.

**In the coming regime, the GHG Protocol Scope 2 revision is the thing that will actually separate
real renewables from certificates** — and it names Australia's NEM and Brazil's CCEE explicitly in its
zonal-deliverability list. Under the proposed rules a DRI/HBI plant (far too large for any exemption)
would have to hourly-match within its own NEM or CCEE zone, and every unmatched hour would be valued
at a **residual mix or a default coal factor** rather than the grid average — which makes partially
matched procurement look *worse* than today. Publication is anticipated **late 2027** with a
multi-year phase-in, so this binds a plant commissioning in the late 2020s but nothing reported now.

**And the standard with the *purest* physical signal is the one nobody uses for claims**: ISO 20915 /
the worldsteel LCI methodology, which requires "the most appropriate regional or national grid mix" or
a known supplier-specific mix and **has no REC/GO/PPA provision at all**. Likewise EUROFER's stainless
proposal (fixed 376 kg CO2/MWh EU mix) and UNESID's, which lists an "arbitrary electricity footprint"
as a **knock-out criterion** — and, most consequentially, **the JRC's ESPR class methodology, which is
country-average throughout**.

**Two designs actively reward *primary* iron rather than merely tolerating it, and both are relevant to
an HBI exporter:**
- **Canada's NRC "primary steel adjustment factor"** subtracts a credit of **up to 585 kg CO2e/t in
  proportion to primary steel content** from the product's embodied carbon before comparing to the
  limit. A green-HBI-fed EAF product takes the full credit on its primary fraction. Unlike a sliding
  scale, this is a *deduction from the number*, ring-fenced so it "shall not be used in any embodied
  carbon disclosure".
- **Germany's CCfD** is the only instrument anywhere that names **liquid pig iron / DRI alone** as a
  fundable product with its own reference value (**1.321 t CO2-eq/t liquid pig iron at 20% scrap**) —
  the closest thing to a numeric "green iron" threshold in force. But it counts **Scope 1 only**, so it
  rewards the *reductant* switch and is completely indifferent to the grid.

### Q3. Which standards treat a hydrogen-DRI plant differently from a molten-oxide-electrolysis plant?

**None of the steel standards do — by design, and in at least two cases the drafters have been told so
and declined to change it. The sliding-scale family is explicitly technology-neutral**: LESS says the classification is set "in a technology-open manner"; ResponsibleSteel,
IEA, FMC and CISA all score t CO2e per t crude steel regardless of route. Whether the electrons went
into an electrolyser making H2 or directly into an MOE cell does not change the number *if* both are
counted at the same electricity factor.

**But four asymmetries fall out of the accounting mechanics, and they all favour or disfavour one
route:**

1. **Where the burden lands — Scope 2 vs Scope 3.** For **MOE and aqueous electrowinning**, the entire
   decarbonisation claim is a **Scope 2** claim: reduction is done by electrons at the cell. For
   **H2-DRI**, the burden splits — on-site electrolysis puts it in Scope 2, but *purchased* hydrogen
   puts it in **Scope 3.1 / 3.3**, where a different set of rules and default factors apply.
   ResponsibleSteel's hydrogen default is **13.916 kg CO2e/kg H2** (grey), so an H2-DRI plant without
   verified hydrogen data is scored as if running on SMR — roughly 0.7-0.8 t CO2e per tonne of DRI.
   An MOE plant has no such trap; it also has no such shelter.

2. **The RFNBO regime (Section 11) binds hydrogen and not electricity.** An H2-DRI plant whose hydrogen
   must qualify as RFNBO already faces **36-month additionality**, **hourly temporal correlation from
   1 January 2030** (monthly until then) and **same-bidding-zone geographic correlation** under
   Delegated Regulation (EU) 2023/1184. An MOE plant making the equivalent claim through Scope 2 faces
   only the 2015 GHG Protocol criteria — annual, national, no additionality. **The hydrogen route is
   held to a materially stricter electricity-provenance test than the direct-electrification route**,
   for the same underlying physical question. This will partially converge if the Scope 2 revision
   lands as drafted, but additionality is proposed to stay *outside* Scope 2 (in the AMI
   multi-statement), so it will not converge fully.

3. **Residual Scope 1.** H2-DRI-EAF retains genuine Scope 1 from carburisation carbon, graphite
   electrodes and pellet carbon — LESS's ideal case attributes its remaining direct emissions to
   "the carbon content of the pellets and electrodes". MOE with an inert anode has essentially no
   process CO2 at all. Under the **IEA "near zero direct emissions" sub-threshold** (where imported
   electricity is N/A), MOE therefore looks near-perfect and H2-DRI does not — an artefact of the
   sub-threshold, not of real climate performance.

4. **Scrap-share classification.** Both routes produce **primary** metallics: ResponsibleSteel counts
   "pig iron, DRI, HBI, ferro-alloys" as `p`, LESS counts "scrap, pig iron, DRI and HBI" in its
   denominator. MOE/electrowinning iron is not named in either list, which is a live drafting gap —
   it would presumably be treated as primary iron, but neither standard says so. **ResponsibleSteel's
   default metallic fraction table (98% scrap / 94% DRI-HBI / 94% pig iron) has no MOE entry either.**
   Anyone modelling MOE against these standards has to make an assumption the standards do not yet
   support.

5. **The EU's own draft ESPR classes actively collapse the distinction.** LESS's objection is that
   "Hot-rolled coil: Natural gas-based DRI-EAF – despite significant emissions reductions and
   investments – **shares Class A with hydrogen-based DRI-EAF, removing incentives to transition from
   gas to hydrogen**", and the JRC's Class A band for HRC (0.00-1.79 t CO2eq/t) confirms it. For wire
   rod the JRC observes that EAF plants span ~0.2-0.8 t CO2eq/t and "**H2-DRI-EAF also sit within Class
   A**". So under the EU's emerging label, **an H2-DRI-EAF and an NG-DRI-EAF making the same product
   would carry the same class** — while under Germany's CCfD (Scope 1 only, ≥85% reduction against
   1.321 t CO2-eq/t liquid pig iron in the final year) they would not be remotely equivalent. These are
   the two EU-adjacent instruments a European H2-DRI project actually faces, and they disagree.

6. **Route-segregation is the alternative some regimes chose instead.** US GSA sets separate limits for
   **EAF** and **integrated mills** (the latter left "TBD" and additionally gated on an ENERGY STAR
   score); Canada uses the primary-steel credit; SBTi runs two entirely separate corporate pathways and
   calls the routes "separate industries". None of these separates H2-DRI from MOE either — they all
   split on scrap/ore, not on how the ore was reduced.

**One further point worth stating plainly:** the IEA's refusal to use the scope vocabulary is what
makes all of this comparable in the first place — "these could be Scope 1 emissions at one site, and
Scope 3 at another site with a different process arrangement… with both sites having the same overall
emissions intensity of steel production." That is exactly the property you need for an
Australia/Brazil-HBI-into-EU-EAF comparison against Germany/France/Spain domestic production. **A
model reporting scope-partitioned numbers needs a mapping layer before it can be compared to an
IEA-boundary threshold.**

---

## 13. Master comparison table

| Standard | Boundary (start → end) | Iron/DRI separate? | Scope 2 counted? | Location or market-based? | Headline threshold | Scrap | Transport in? | Status |
|---|---|---|---|---|---|---|---|---|
| **ResponsibleSteel v2.1.1** | cradle → crude steel (site); PCF cradle→gate A1-A3 | No graded scale for iron; stand-alone iron plants can hold Core Site Cert. | **Yes**, full CO2e | **Location default (ISO 14064-1 E.2 grid factor); MARKET-BASED PERMITTED** (REC/PPA/vPPA/green tariff, ISO 14064-1 E.2.2) with mandatory disclosure. Offsets banned. | DPL4 near-zero `y ≤ 0.40 − 0.35x` t CO2e/t crude steel; DPL1 `2.80 − 2.30x` | Sliding scale on **metallic input**; DRI/HBI = primary; cut-off (zero burden) | **Yes** — Scope 3.4 incl. purchaser's own carriage | 47 certificates, 45 Core / **2 Certified Steel**; Brazil 10, Australia 1 |
| **LESS v1.1** | cradle → **hot-rolled steel** (incl. 1st heat) | No | **Yes**, but **CO2 only** for Scope 1&2 | "real emission factor from the actual electricity mix" BUT **GOs from officially recognised registers expressly permitted**, incl. green power "not physically delivered". Annex VI: renewable electricity = **0.0000** Scope 2. Registers named are **German**. | Near-Zero QST `520 − 3.5x` kg/t HR; BST `470 − 3.5x`; A/B/C/D = 2×/3×/4×/5× | Sliding scale on **raw tonnage** (scrap ÷ scrap+pig iron+DRI+HBI); cut-off | **Yes** — Scope 3.4 | 8 certificates, all German/Franco-German; tk Steel at **Level D** |
| **IEA (2022, restated 2024)** | cradle → crude steel incl. casting, excl. rolling | **No** — flagged as a gap | **Yes** ("imported electricity, heat and hydrogen") — but a **direct-only sub-threshold** exists at the *same value* | **UNRESOLVED** — "use of renewable electricity credits… Further discussion on this topic is likely to be needed" | Near-zero `400 − 350s` kg CO2e/t crude steel; low-emissions line `2400 − 2100s`, A-E, credited fractionally | Sliding scale on metallic input; DRI/HBI = primary; scrap sorting/transport **excluded** | **Partial** — ore + limestone + fossil fuels only; **not** purchased DRI/HBI | G7/Climate Club "robust starting point"; not enacted anywhere |
| **worldsteel CO2 / ISO 14404** | ore agglomeration → finished steel | 14404-3 covers the DRI-EAF *plant*; purchased DRI handled as "outsourced activity" | **Yes**, Scope 1+2+3, CO2 only | **NEITHER** — fixed **world-average 0.504 t CO2/MWh**. A PPA does nothing; a clean grid does nothing. | No threshold — route averages only | Cut-off, zero burden; sorting/transport out | **No** | 220+ sites, ~25% of global production; ISO 14404-1/2/3 2nd ed. 2024 |
| **ISO 20915 / worldsteel LCI** | cradle → gate (optionally cradle→grave) | No | Yes, all GHGs | **LOCATION-BASED** — "most appropriate regional or national grid mix" or known supplier mix. **No REC/GO/PPA provision.** | No threshold — an LCI method | **Closed-loop / "value of scrap"**: `X − (RR − S)·Y·(X_pr − X_re)` | Yes | ISO/TC 17; ~90% of production if fully adopted |
| **GHG Protocol Scope 2** | n/a (accounting rule) | n/a | n/a | **Dual reporting mandatory.** 2015: annual, "same market". **Proposed**: hourly matching, zonal deliverability (Australia NEM and Brazil CCEE named), SSS pro-rata, fossil-only fallback, **no additionality** | n/a | n/a | n/a | Consultation closed 31 Jan 2026; feedback published 29 July 2026; ISO consolidation; **publication anticipated late 2027** + multi-year phase-in |
| **EU voluntary label / IDAA** | "building on the CBAM methodology" + ETS data | Not addressed | Expected yes, via the ESPR track | Not yet stated at IDAA level | None published | Not stated | Not stated | Announced ESMAP March 2025; expected 10 Dec 2025 |
| **JRC ESPR classes (draft 1 Apr 2026)** | cradle → gate, five intermediate products (HRC, WR, CRCG, ES, SS) | **No — merchant iron absent entirely** | **Yes** — ETS/CBAM Scope 1 + Scope 2 + Scope 3 | **LOCATION-BASED, COUNTRY-AVERAGE.** EU average (EEA) for EU; CBAM defaults = "average electricity mix of the country of origin" for imports. **No PPA/GO/REC provision.** | **HRC**: A ≤1.79, B <2.66, C <3.10, D <3.75, E ≥3.75 t CO2eq/t. **WR**: A ≤0.87, B <2.43, C <3.10, D <3.54, E ≥3.54 | **NOT a sliding scale** — fixed absolute bands; scrap affects the footprint only, EF approach counts collection/treatment | Yes (cradle-to-gate) | Draft; 2nd stakeholder consultation Apr 2026; thresholds "indicative and may be revised" |
| **IDAA lead markets / GPP** | — | — | — | — | Implied **HRC ≥ Class B (≤2.66)** / **WR ≥ Class B (≤2.43)** if the B/C boundary is used, calibrated so ~30% of global volume complies | — | — | Not adopted |
| **Germany CCfD (KSV)** | the funded installations only | **YES — liquid pig iron / DRI is a fundable product in its own right** | **NO — excluded outright** ("Scope-1-Emissionen"; benchmarks "reduced by the indirect emissions for electricity") | n/a | Primary steel reference **1.321 t CO2-eq/t liquid pig iron** at 20% scrap → **≥50% from yr 4 and ≥85% in final yr** (2026 round) ⇒ ≤0.661 / ≤0.198 t CO2e/t | 20% scrap assumed in the reference; no sliding scale | No | Round 2 closed **7 Sept 2026**; €5 bn; max bid €550/t CO2e |
| **Germany Leitmärkte / § 113(1) Nr. 9 GWB** | cradle → hot-rolled (>90% of emissions) | No | **Yes**, Scope 1+2+3 upstream | Not stated — defers to the WV Stahl rulebook (LESS ⇒ GOs permitted) | Near-zero + A-D; QST 0-520, BST 0-120 kg CO2äq/t; ladder D 2030 / C 2035 / A 2040 / NZ 2045 | LESS sliding scale | Yes | Enabling clause in force 1 July 2026, **no number**; "keine Quotenvorgaben"; BMWE to start implementation in 2026 |
| **First Movers Coalition (2025 rev.)** | cradle → crude steel incl. casting (IEA boundary) | No | **Yes** | **Defers to GHG Protocol; book & claim expressly carved in for electricity** (banned for the steel itself). 2022 vPPA-additionality condition **dropped**. | `<0.4 − 0.35s` t CO2e/t crude steel (was <0.1 at 100% scrap, tightened to **0.05**, CO2→CO2e) | Sliding scale; **scrap sorting/transport excluded** | Ore + lime "regardless of intermediary stops"; scrap excluded | 101 members, >$19 bn, 10% of purchases by 2030 |
| **GSCC Steel Climate Standard** | mining → hot rolling, incl. scrap collection | No | Yes, all GHGs | **MARKET-BASED BY DEFAULT**; residual mix next; "location-based… **only when market-based and residual mix data is unavailable**". Not applicable to off-gas electricity. | t CO2e/t HR: 2026 long 1.29 / flat 1.57; 2030 1.11/1.31; 2050 **0.12/0.12** | **No sliding scale by design** — "could result in higher emissions steel products being labeled as 'green'" | Yes | Product certification live |
| **Climate Bonds Initiative** | crude/finished steel, **no upstream Scope 3** | No | Yes (Scope 1+2+3.3 only), CO2 only | Not specified | Primary 1.81 (2030) → **0.12 (2050)**; secondary 0.32 (2030) → 0.12 (2050) t CO2e/t | Weighted pathway on ore:scrap ratio | **No** | Tier 1 / Tier 2 certification |
| **CISA C2F (China)** | cradle → crude or hot-rolled | No | Yes, CO2 only, Scope 3.1 partial | Not specified | A **0.4 → 0.05** t CO2/t; E 2.19 → 0.41 | Sliding scale, 5 classes | No | Launched Oct 2024; RS interoperability agreement Nov 2025 |
| **SBTi Steel v1.0** | company-level SDA | Purchased HBI counts **inside** | Yes (scope 1+2) | Not specified; renewable-sourcing targets an "acceptable alternative" | t CO2/t HR: ore-based 1.71 (2030) → 0.11 (2050); scrap-based 0.37 → 0.11 | Two pathways, scrap-weighted average | — | Stakeholder survey to 11 Sept 2026 |
| **IDDI GPP Pledge** | — | — | — | — | **Four Levels, no numbers** (defers to IEA) | — | — | AT, CA, DE, JP, UAE, UK, US; pledge site offline |
| **SteelZero** | — | — | — | — | 50% of requirement at **ResponsibleSteel DPL2** (`≤2.00 − 1.65x`) by 2030; 100% net zero 2050 | via RS | — | 40+ businesses, ~10 Mt |
| **US GSA (IRA)** | A1-A3 | No | Inside A1-A3 | **Unresolved** — defers to UL/SCS PCR | Rebar unfab. **611/716/760**; hot-rolled sections unfab. **686/713/869**; EAF plate **987/1152/1190** kg CO2e/t | **No sliding scale — route segregation** (EAF set, integrated "TBD" + ENERGY STAR gate) | A1-A3 only | Defunded going forward (OBBBA 2025) |
| **Buy Clean California** | A1-A3 | No | Inside A1-A3 | Unresolved | Rebar **755 (unfab.) / 778 (fab.)** kg CO2e/t, tightened from 890/920 | No | A1-A3 | Effective 1 Jan 2025 |
| **Colorado (EE-5.1, 01/2026)** | A1-A3 | No | Inside A1-A3 | Unresolved | Fab. rebar **1030**; hot-rolled **1220**; plate **1730**; HSS **1990**; cold-formed **2843** kg CO2e/t | No | A1-A3 (+ separate A4 report >100 miles) | **Binding caps**, projects >$500k; ratchet-only-downward |
| **Minnesota** | cradle → mill gate | No | Inside | Unresolved | **Rebar 0.755 t CO2e/t** (fab. 0.854) | No | — | Projects advertised from 15 July 2026; other products EPD-only until 2028 |
| **Canada NRC (Aug 2025)** | A1-A3 | No | Inside A1-A3 | Unresolved | Rebar **703/773/793**; hot-rolled sections **1008/1105/1140**; cold-rolled sheet **1848/2021/2102** kg CO2e/t | **INVERTED sliding scale** — a **primary steel credit of up to 585 kg CO2e/t × primary fraction** deducted before comparison; ring-fenced from disclosure | A1-A3 | Published 26 Aug 2025 |
| **India Green Steel Taxonomy (Dec 2024)** | up to finished steel | No | Yes (Scope 1+2 + limited 3) | Not specified | "green steel" **< 2.2**; 5-star **<1.6**; 4-star **1.6-2.0**; 3-star **2.0-2.2** t CO2e/t finished steel | No sliding scale | — | NISST as MRV agency; reviewed every 3 years |
| **UK** | — | — | — | — | **None.** Low-carbon products framework deferred 28 May 2026; sliding-scale question "a priority"; **UK CBAM indirect emissions delayed to 2029 at the earliest** | — | — | Steel Strategy March 2026; procurement is origin-based (PPN 022) |

---

## 14. Gaps and caveats

**Genuinely unresolved after this research:**
- **FHWA §60506 steel GWP numbers** — pages removed (403/404), no usable Wayback capture.
- **EPA IRA §60116 label programme** — 2026 operational status not verified.
- **ISO 14404-4 clause 9.3** — only the free preview (through clause 3) was obtainable; whether it
  explicitly authorises contractual/market-based electricity factors is unconfirmed. The JISF rationale
  argues strongly that it does not.
- **Spanish PERTE and Italian numeric conditions** — Spain's IP/23/849 PDF would not decode; Italy's
  conditions not found.
- **Canada Treasury Board Standard on Embodied Carbon in Construction** and Contracting Policy Notice
  2025-6 — fetch blocked.
- **Whether the US steel PCRs (UL / SCS) permit RECs** — not resolved; this determines the scope 2
  treatment for GSA, BCCA, Colorado, Minnesota and Canada NRC alike.
- **LESS's DRI standard emission factor** — the public annex publishes the source (Duarte et al. 2021)
  but not the number.
- **The `[X]` GWh exemption threshold** in the GHG Protocol Scope 2 draft is genuinely undecided.
- **FMC per-company steel signatory list** — weforum.org returns 403 to automated fetches.
- IDDI 2026 signatory count; OECD steel standards work; Salzgitter's funding split; New Jersey Buy Clean
  steel limits.

**Things to treat with care:**
- The **JRC ESPR classes are a draft** and the thresholds are described in the document itself as
  "indicative and may be revised". LESS is actively lobbying to replace them with a sliding scale, and
  the UK has made the same question "a priority". The class boundaries above should be treated as the
  current proposal, not settled law.
- The **LESS 450/170 and 520/170 pairs are both correct** — 450 is the 20%-scrap reference point in the
  Rulebook, 520 the 0%-scrap extrapolation the JRC quotes. `E = m(100−x) + s` with `m = (p−s)/80`
  reconciles them exactly.
- The **BMWK Leitmärkte "0-520 (QST) / 0-120 (BST)"** pair is *not* like-for-like — 520 is the QST
  0%-scrap value and 120 the BST 100%-scrap value. Use the full LESS table.
- **FMC's 100%-scrap threshold is 0.05, not 0.1** — the 2022 figure was tightened in the 2025 review and
  many secondary sources still quote the old number.
- Column 3 of the LESS Annex VI factor table is *inferred* to be the 50%-reduced Scope 3 variant used in
  the near-zero reference plant, on the basis of the exact 2:1 ratio across every row.
- The **thyssenkrupp hydrogen-condition unwinding** (running the DRI plant initially on natural gas) is
  reported and Commission-approved as state-aid compatible, but no formal amended Zuwendungsbescheid was
  located.

