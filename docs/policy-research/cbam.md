# EU CBAM — which emissions count, for iron/steel and hydrogen

Research note prepared 2026-09-07. Focus: green iron/steel produced in Australia and
Brazil, exported as HBI (hot briquetted iron, CN 7203) into an EU electric arc furnace.

**Sourcing note.** `eur-lex.europa.eu` is behind an AWS WAF challenge and cannot be
fetched directly (both `WebFetch` and `curl` receive a `202` challenge). All EU legal
texts below were retrieved from the EU Publications Office cellar endpoint, which serves
the identical authenticated OJ text:

```
curl -sL -H "Accept: application/xhtml+xml" -H "Accept-Language: eng" \
     "https://publications.europa.eu/resource/celex/<CELEX>"
```

Citations give both the ELI (`http://data.europa.eu/eli/...`) and the CELEX id used.

---

## 0. The legal stack as it actually stands (September 2026)

| Instrument | What it is | ELI / CELEX |
|---|---|---|
| Regulation (EU) 2023/956 | The CBAM Regulation (base act) | http://data.europa.eu/eli/reg/2023/956/oj — CELEX `32023R0956` |
| Regulation (EU) 2025/2083 | Amends 2023/956, "simplifying and strengthening the CBAM" (the *Omnibus*), of 8 Oct 2025, OJ 17.10.2025 | http://data.europa.eu/eli/reg/2025/2083/oj — CELEX `32025R2083` |
| Consolidated 2023/956 | Post-Omnibus consolidation | http://data.europa.eu/eli/reg/2023/956/2025-10-20 |
| Implementing Reg. (EU) **2025/2547** | **The Methodology Act** — "methods for the calculation of emissions embedded in goods", of 10 Dec 2025, OJ 22.12.2025 | http://data.europa.eu/eli/reg_impl/2025/2547/oj — CELEX `32025R2547` |
| Implementing Reg. (EU) 2025/2546 | Verification principles, of 10 Dec 2025 | http://data.europa.eu/eli/reg_impl/2025/2546/oj |
| Delegated Reg. (EU) **2025/2551** | **Accreditation of verifiers** — *not* the methodology act, of 20 Nov 2025 | http://data.europa.eu/eli/reg_del/2025/2551/oj |
| Implementing Reg. (EU) **2025/2620** | Free allocation adjustment + **CBAM benchmarks**, of 16 Dec 2025 | http://data.europa.eu/eli/reg_impl/2025/2620/oj — CELEX `32025R2620` |
| Implementing Reg. (EU) **2025/2621** | **Default values**, of 16 Dec 2025, OJ 31.12.2025 | http://data.europa.eu/eli/reg_impl/2025/2621/oj — CELEX `32025R2621` |
| Implementing Reg. (EU) **2026/1740** | **Corrects 2025/2621, replacing Annexes I and IV in full**, of 20 July 2026, OJ 31.7.2026, applies retroactively from 1 Jan 2026 | CELEX `32026R1740` |
| Guidance No. 3 | CBAM methods for the calculation of emissions embedded in goods | https://taxation-customs.ec.europa.eu/document/download/29b9eec7-1a4b-4eb6-ab85-96a0c9e35fd0_en |
| Guidance No. 4 | CBAM calculation of the free allocation adjustment | https://taxation-customs.ec.europa.eu/document/download/3aa2c730-4f1b-4524-8f71-d9cfe4880ac7_en |
| Guidance No. 5b | Sector-specific guidance: **hydrogen** | https://taxation-customs.ec.europa.eu/document/download/04ceca3b-466e-4ed9-a1a5-a47c2e5e5be0_en |
| Guidance No. 5d | Sector-specific guidance: **iron and steel** | https://taxation-customs.ec.europa.eu/document/download/3eb64513-8255-4c71-875a-bd0a8c7ff772_en |

Landing page for all of the above:
https://taxation-customs.ec.europa.eu/carbon-border-adjustment-mechanism/cbam-legislation-and-guidance_en

**Correction to the premise in the brief.** Delegated Regulation (EU) 2025/2551 is *not*
the Methodology Act. Its official title is: *"Commission Delegated Regulation (EU)
2025/2551 of 20 November 2025 supplementing Regulation (EU) 2023/956 … by specifying the
conditions for granting accreditation to verifiers, for the control and oversight of
accredited verifiers, for the withdrawal of accreditation and for mutual recognition and
peer evaluation of accreditation bodies."* The Methodology Act is **Implementing
Regulation (EU) 2025/2547**.

**Second correction.** **Implementing Regulation (EU) 2026/1740** (July 2026) replaced
Annexes I and IV of 2025/2621 wholesale and applies retroactively from 1 January 2026. Any
default-value table quoted from the December 2025 version of 2025/2621 is superseded. In
particular the per-year "including mark-up" columns were **deleted**; the mark-up is now
applied inside the CBAM Registry to the "total emissions" column (see §4).

---

## 1. ANNEX II to Regulation (EU) 2023/956 — verbatim

### 1.1 Article 7(1), the operative provision

> **Article 7 — Calculation of embedded emissions**
>
> 1. Embedded emissions in goods shall be calculated pursuant to the methods set out in
> Annex IV. **For goods listed in Annex II only direct emissions shall be calculated and
> taken into account.**
>
> 2. Embedded emissions in goods other than electricity shall be determined based on the
> actual emissions in accordance with the methods set out in points 2 and 3 of Annex IV.
> Where the actual emissions cannot be adequately determined, as well as in the case of
> indirect emissions, the embedded emissions shall be determined by reference to default
> values in accordance with the methods set out in point 4.1 of Annex IV.

(CELEX `32023R0956`, Article 7.)

### 1.2 Annex II verbatim (original 2023 text)

> **ANNEX II**
>
> **List of goods for which only direct emissions are to be taken into account, pursuant
> to Article 7(1)**
>
> **Iron and steel**
>
> | CN code | Greenhouse gas |
> |---|---|
> | **72 – Iron and steel**<br>Except:<br>7202 2 – Ferro-silicon<br>7202 30 00 – Ferro-silico-manganese<br>7202 50 00 – Ferro-silico-chromium<br>7202 70 00 – Ferro-molybdenum<br>7202 80 00 – Ferro-tungsten and ferro-silico-tungsten<br>7202 91 00 – Ferro-titanium and ferro-silico-titanium<br>7202 92 00 – Ferro-vanadium<br>7202 93 00 – Ferro-niobium<br>7202 99 – Other:<br>  7202 99 10 – Ferro-phosphorus<br>  7202 99 30 – Ferro-silico-magnesium<br>  7202 99 80 – Other<br>7204 – Ferrous waste and scrap; remelting scrap ingots and steel | Carbon dioxide |
> | 7301 – Sheet piling of iron or steel, whether or not drilled, punched or made from assembled elements; welded angles, shapes and sections, of iron or steel | Carbon dioxide |
> | 7302 – Railway or tramway track construction material of iron or steel … | Carbon dioxide |
> | 7303 00 – Tubes, pipes and hollow profiles, of cast iron | Carbon dioxide |
> | 7304 – Tubes, pipes and hollow profiles, seamless, of iron (other than cast iron) or steel | Carbon dioxide |
> | 7305 – Other tubes and pipes … the external diameter of which exceeds 406,4 mm, of iron or steel | Carbon dioxide |
> | 7306 – Other tubes, pipes and hollow profiles … of iron or steel | Carbon dioxide |
> | 7307 – Tube or pipe fittings (for example, couplings, elbows, sleeves), of iron or steel | Carbon dioxide |
> | 7308 – Structures (excluding prefabricated buildings of heading 9406) and parts of structures … of iron or steel … | Carbon dioxide |
> | 7309 00 – Reservoirs, tanks, vats and similar containers … of a capacity exceeding 300 l … | Carbon dioxide |
> | 7310 – Tanks, casks, drums, cans, boxes and similar containers … of a capacity not exceeding 300 l … | Carbon dioxide |
> | 7311 00 – Containers for compressed or liquefied gas, of iron or steel | Carbon dioxide |
> | 7318 – Screws, bolts, nuts, coach screws, screw hooks, rivets, cotters, cotter pins, washers (including spring washers) and similar articles, of iron or steel | Carbon dioxide |
> | 7326 – Other articles of iron or steel | Carbon dioxide |
>
> **Aluminium**
>
> | CN code | Greenhouse gas |
> |---|---|
> | 7601 – Unwrought aluminium | Carbon dioxide and perfluorocarbons |
> | 7603 – Aluminium powders and flakes | " |
> | 7604 – Aluminium bars, rods and profiles | " |
> | 7605 – Aluminium wire | " |
> | 7606 – Aluminium plates, sheets and strip, of a thickness exceeding 0,2 mm | " |
> | 7607 – Aluminium foil … of a thickness (excluding any backing) not exceeding 0,2 mm | " |
> | 7608 – Aluminium tubes and pipes | " |
> | 7609 00 00 – Aluminium tube or pipe fittings | " |
> | 7610 – Aluminium structures … and parts of structures … | " |
> | 7611 00 00 – Aluminium reservoirs, tanks, vats and similar containers … exceeding 300 litres … | " |
> | 7612 – Aluminium casks, drums, cans, boxes and similar containers … not exceeding 300 litres … | " |
> | 7613 00 00 – Aluminium containers for compressed or liquefied gas | " |
> | 7614 – Stranded wire, cables, plaited bands and the like, of aluminium, not electrically insulated | " |
> | 7616 – Other articles of aluminium | " |
>
> **Chemicals**
>
> | CN code | Greenhouse gas |
> |---|---|
> | **2804 10 00 – Hydrogen** | Carbon dioxide |

### 1.3 The one amendment: electricity added

Regulation (EU) 2025/2083, Article 1, point (29):

> (29) in Annex II, the following table is added:
>
> ‘**Electricity**
>
> | CN code | Greenhouse gas |
> |---|---|
> | 2716 00 00 – Electrical energy | Carbon dioxide’ |

Rationale, recital (87) of 2025/2083:

> Annex II to Regulation (EU) 2023/956 lists the goods for which only direct emissions
> should be taken into account in the calculation of embedded emissions. For goods not
> listed in that Annex, both direct and indirect emissions should be taken into account.
> Since indirect emissions are not relevant in the case of electricity generation,
> electricity should be added to the list of goods in that Annex.

2025/2083 made **no other change to Annex II**. It amended Annex I only to narrow the
kaolinic-clay entry (`2507 00 80` → `ex 2507 00 80 – Other kaolinic clays except
non-calcined kaolinic clays`, point (28)), amended Annex IV (point (30) + its own Annex I),
and added a new Annex VII (point (33) + its own Annex II).

### 1.4 Complete current Annex II scope, by sector

| Sector | Goods covered by Annex II |
|---|---|
| Iron and steel | All of CN Chapter **72**, except ferro-silicon/-silico-manganese/-silico-chromium/-molybdenum/-tungsten/-titanium/-vanadium/-niobium/-phosphorus/-silico-magnesium and "other" ferro-alloys of 7202, and except **7204** (ferrous waste and scrap); plus **7301, 7302, 7303 00, 7304, 7305, 7306, 7307, 7308, 7309 00, 7310, 7311 00, 7318, 7326** |
| Aluminium | 7601, 7603, 7604, 7605, 7606, 7607, 7608, 7609 00 00, 7610, 7611 00 00, 7612, 7613 00 00, 7614, 7616 |
| Chemicals | **2804 10 00 — Hydrogen** |
| Electricity (added 2025) | 2716 00 00 |

Everything else in Annex I — i.e. cement (2507 00 80 ex, 2523 xx), fertilisers (2808, 2814,
2834 21 00, 3102, 3105) and **sintered ore (2601 12 00)** — carries **direct + indirect**.

### 1.5 CONFIRMED: hydrogen (2804 10 00) is on Annex II

Yes — under the heading "Chemicals". It has been there since the original 2023 text and
was not touched by the Omnibus. **Hydrogen's embedded emissions for CBAM purposes are
direct emissions only; its indirect (electricity) emissions are excluded.**

### 1.6 CONFIRMED: the sintered-ore exception

Annex **I** (scope) lists, under "Iron and steel":

> **2601 12 00 – Agglomerated iron ores and concentrates, other than roasted iron
> pyrites** | Carbon dioxide

Annex **II** does **not** list 2601 12 00 — its iron-and-steel entry begins at Chapter 72.
Sintered ore / iron ore pellets are therefore the one iron-and-steel-sector good for which
**both direct and indirect emissions count**.

This is visible in the default-value tables: every CN 72xx row carries `N/A` in the
indirect column, whereas 2601 12 00 carries a real number (see §4).

Guidance 5d states it explicitly (p. 8):

> For sintered ore (CN 2601 12 00) – not listed in Annex II to the CBAM Regulation,
> indirect emissions are to be included.

**Practical consequence for the HBI case.** DR-grade pellets fall under CN 2601 12 00. When
they are bought in as a precursor, their *indirect* emissions flow into the HBI's embedded
emissions even though the HBI itself is Annex II. If pelletising happens inside the same
installation and the pellets are never sold or transferred out, Article 4(9) of the
Methodology Act allows a **joint production process**, in which case there is no separate
precursor and the pelletising electricity never enters the calculation (see §3.3).

---

## 2. THE METHODOLOGY ACT — Implementing Regulation (EU) 2025/2547

Full title:

> **COMMISSION IMPLEMENTING REGULATION (EU) 2025/2547 of 10 December 2025 laying down
> rules for the application of Regulation (EU) 2023/956 of the European Parliament and the
> Council as regards the methods for the calculation of emissions embedded in goods**
> (OJ L, 2025/2547, 22.12.2025). Adopted under Article 7(7) of 2023/956. Replaces
> Implementing Regulation (EU) 2023/1773, which governed the transitional period.

Structure: 5 chapters, 16 articles, Annexes I–V.

- **Annex I** — Definitions, functional unit and **system boundaries**. Point 2 = mapping
  of CN codes to aggregated goods categories. Point 3 = the boundaries: 3.1 cross-sectoral,
  3.6 hydrogen, 3.12 sintered ore, 3.13 pig iron, 3.14 DRI, 3.15 crude steel, 3.16 iron or
  steel products, 3.17–3.19 aluminium.
- **Annex II** — Monitoring rules (A: principles and monitoring plan; B: direct emissions,
  incl. **B.3.2 mass balance** and B.3.3 biomass zero-rating; C: heat flows; D: electricity;
  **E: precursors**; F: activity levels).
- **Annex III** — Rules for attributing emissions to goods (A: attribution to production
  processes; **B: calculation of specific embedded emissions of complex goods**).
- **Annex IV** — Templates for the operator's emissions report.
- **Annex V** — Alternative default values.

### 2.1 The controlling article: Article 3

> **Article 3 — System boundaries**
>
> 1. In order to quantify and calculate specific embedded emissions of goods, the processes
> within an installation that occur within the system boundaries, defined per aggregated
> goods category in accordance with Annex I, shall be taken into account.
>
> 2. **The system boundaries shall cover direct emissions, indirect emissions for goods not
> listed in Annex II to Regulation (EU) 2023/956, and the embedded emissions of any
> precursor.**

### 2.2 Cross-sectoral rule — Annex I, point 3.1

> Specific embedded emissions shall be calculated as the emissions of the production
> process and, for complex goods, the embedded emissions of the precursors to produce the
> functional unit of the good during the reporting period.
>
> The system boundaries are defined per aggregated goods categories and cover the direct
> emissions, the indirect emissions from electricity consumption where relevant under
> Regulation (EU) 2023/956, emitted by all processes directly or indirectly linked to the
> production processes, and the embedded emissions of precursors, **independently of
> whether these precursors are produced in the installation or acquired from a different
> installation**. In addition to these general rules, the specific details of each
> aggregated goods category are set out in points 3.2 to 3.19. Any CBAM goods produced by
> means of a production route not listed in points 3.2 to 3.19 is subject to the
> cross-sectoral rules described in this point, and to the sector-specific rules if the
> production route is a combination of the production routes listed in points 3.2 to 3.19.
>
> **The purchase and maintenance of infrastructure and equipment are excluded from the
> system boundaries.**
>
> **When the production process of complex goods listed in Annex II to Regulation (EU)
> 2023/956 includes one or more precursors not listed in that Annex, the indirect emissions
> of those precursors will be included in the calculation of the embedded emissions of the
> complex goods. When the production process of complex goods not listed in that Annex
> includes one or more precursors listed in that Annex, the indirect emissions of these
> precursors will not be included in the calculation of the embedded emissions of the
> complex goods.**

That last paragraph is the single most load-bearing sentence for this project. Applied to
HBI:

- HBI/DRI (CN 7203) **is** an Annex II good.
- Sintered ore / pellets (2601 12 00) **is not** → **its indirect emissions ARE included in
  HBI's embedded emissions.**
- Hydrogen (2804 10 00) **is** → **its indirect emissions are NOT included.**

The sentence also means the fallback clause matters: a production route not enumerated in
points 3.2–3.19 (e.g. water electrolysis, or a novel H2 shaft) is not out of scope — it
simply falls under the cross-sectoral rules of point 3.1.

### 2.3 Point 3.6 — Hydrogen

> **3.6. Hydrogen**
>
> **3.6.1. Special provisions**
>
> Only the production of pure hydrogen or mixtures of hydrogen with nitrogen usable in
> ammonia production shall be considered. Not covered are the consumption of synthesis gas
> or of hydrogen as precursor within refineries or organic chemical installations, where
> hydrogen is exclusively used within those plants and not used for the production of goods
> listed in Annex I to Regulation (EU) 2023/956.
>
> **3.6.2. System boundary**
>
> **3.6.2.1. Steam reforming and partial oxidation**
>
> For those production routes, direct emissions monitoring shall take into account:
> - all processes directly or indirectly linked to hydrogen production and the separation
>   of hydrogen and carbon monoxide, and flue gas cleaning;
> - all fuels used in the hydrogen production process irrespective of their energetic or
>   non-energetic use, and fuels used for other combustion processes including for the
>   purpose of producing hot water or steam.
>
> **3.6.2.2. Steam cracking**
>
> For that production route, direct emissions monitoring shall take into account:
> - all processes directly linked to hydrogen production;
> - all processes directly or indirectly linked to the production processes, and from flue
>   gas cleaning.

A gasification route follows: *"That production route applies where hydrogen is produced by
gasification of coal, heavy refinery fuels or other fossil feedstock. Input materials may
include biomass, for which the provisions of point B.3.3 of Annex II shall be taken into
account."*

**No electrolysis route is defined.** Guidance 5b (p. 10) confirms this is deliberate:

> Note that chlor-alkali electrolysis and water electrolysis are not mentioned in the
> Methodology Act, as these production routes include minimal direct emissions.

Electrolytic hydrogen therefore falls back on the cross-sectoral rule in Annex I point 3.1,
and its CBAM embedded emissions are essentially **zero**: no fuel combustion, no process
emissions from carbonates, and — because hydrogen is Annex II — no indirect emissions from
the electrolyser's electricity, regardless of how that electricity was generated.

Guidance 5b, §2.3.1 and footnote 2:

> In the definitive period, for CBAM purposes the hydrogen sector has to account for direct
> emissions only, since hydrogen is listed in Annex II to the CBAM Regulation as a good for
> which only direct emissions are taken into account. Indirect emissions were reported
> separately during the transitional period but are not included in the embedded emissions
> of hydrogen in the definitive period. Emissions must be reported in metric tonnes CO2
> equivalent (t CO2e) emissions per tonne of output.

> Hydrogen is defined as a simple good in accordance with Article 1 and Annex I of the
> Methodology Act, as there are no possible precursors listed for the production of
> hydrogen and the embedded emissions are determined solely from the direct emissions of
> the hydrogen production process.
>
> There are no precursors relevant for hydrogen. However, hydrogen may itself be a
> precursor for other processes, where it is separately produced for use as a chemical
> feedstock to produce ammonia, or to produce pig iron or direct reduced iron (DRI).

Note the "not covered" carve-out in 3.6.1 is about *captive* hydrogen inside refineries and
organic-chemical plants. Captive hydrogen produced and consumed inside a DRI installation
is *not* carved out — it is used to produce an Annex I good. In practice, if the electrolyser
sits inside the same installation as the shaft and the H2 is never sold or transferred, the
operator can use Article 4(9) and treat the whole thing as one joint production process.

### 2.4 Point 3.13 — Pig iron

> **3.13. Pig Iron**
>
> **3.13.1. Special provisions**
>
> This aggregated goods category includes non-alloyed pig iron from blast furnaces as well
> as alloy-containing pig irons (e.g., spiegeleisen), irrespective of the physical form
> (e.g. ingots, granules). NPI (nickel pig iron) is included if the nickel content is lower
> than 10 %. In integrated steel plants, liquid pig iron ('hot metal') directly charged to
> the oxygen converter is the product which separates the production process for pig iron
> from the production process of crude steel. Where the installation does not sell or
> transfer pig iron to other installations, a joint production process including crude
> steel can be established making subject to the rules of Article 4.
>
> **3.13.2. System boundary**
>
> **3.13.2.1. Blast furnace route** — direct emissions monitoring shall encompass:
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from **fuels and reducing agents** such as coke, coke dust, coal, fuel oils, plastic
>   wastes, natural gas, wood wastes, charcoal, as well as from **waste gases** such as coke
>   oven gas, blast furnace gas or converter gas;
> - where biomass is used, the provisions of point B.3.3 of Annex II shall be taken into
>   account;
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from **process materials** such as limestone, magnesite, and other carbonates, carbonate
>   ores; materials for flue gas cleaning;
> - **carbon remaining in the product or in slags or wastes is taken into account by using a
>   mass balance method in accordance with point B.3.2 of Annex II.**
>
> **3.13.2.2. Smelting reduction** — same four bullets, with "waste gases from the process
> or converter gas" in place of the coke-oven/blast-furnace gas list.

### 2.5 Point 3.14 — DRI (Direct Reduced Iron) — the one that matters here

> **3.14. DRI (Direct Reduced Iron)**
>
> **3.14.1. Special provisions**
>
> **There is only one production route defined**, although different technologies may use
> different qualities of ores, which may require pelletisation or sintering, and different
> reducing agents (natural gas, diverse fossil fuels or biomass, **hydrogen**). **Therefore,
> precursors sintered ore or hydrogen may be relevant.** As products, iron sponge, **hot
> briquetted iron (HBI)** or other forms of direct reduced iron may be relevant, including
> DRI which is immediately fed to electric arc furnaces or other downstream processes.
>
> Where the installation does not sell or transfer DRI to other installations, a joint
> production process including steel can be established making subject to the rules of
> Article 4.
>
> **3.14.2. System boundary**
>
> For that production route, direct emissions monitoring shall encompass:
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from **fuels and reducing agents** such as coal, natural gas, fuel oils, waste gases
>   from the process or converter gas, etc.;
> - where biogas or other forms of biomass are used, the provisions of point B.3.3 of Annex
>   II shall be taken into account;
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from **process materials** such as limestone, magnesite, and other carbonates, carbonate
>   ores, materials for flue gas cleaning;
> - **carbon remaining in the product or in slags or wastes is taken into account by using a
>   mass balance method in accordance with point B.3.2 of Annex III.**

**Drafting error flagged.** The last bullet cites "point B.3.2 of **Annex III**". B.3.2
("Mass balance method") is in **Annex II**, not Annex III. Points 3.6, 3.13.2.1 and 3.13.2.2
all correctly cite Annex II, whereas 3.14.2, 3.15.2.1 and 3.15.2.2 cite Annex III. Annex III
of 2025/2547 is "Rules for attributing emissions to goods" and contains no point B.3.2 —
its point B is "Calculation of specific embedded emissions of complex goods". Guidance 5d
reproduces both variants inconsistently (quotes 3.14 as "Annex II" on p. 24, quotes 3.15 as
"Annex III" on p. 30). Treat as a typographical slip: the mass balance method is Annex II
point B.3.2 throughout.

**Briquetting is inside the boundary.** Guidance 5d §2.2.3.5 lists the production steps
within DRI installation system boundaries:

> - Raw material handling and pre-treatment.
> - Fuel storage and preparation – coal, natural gas or hydrogen etc.
> - **Direct reduction process for iron production – all steps for the DRI process, forming
>   into hot briquetted iron (HBI) if applicable.**
> - Emissions control – in particular flue gas cleaning.
>
> […] Precursors (if used in the process) are: **sintered ore; hydrogen**; pig iron or DRI
> from other installations or production processes; and ferro-alloys FeMn, FeCr, FeNi, if
> used. During the transitional period, indirect emissions that result from electricity
> consumed by the production process also had to be monitored and reported. **In the
> definitive period, indirect emissions are not taken into account for the calculation of
> embedded emissions of iron and steel CBAM goods listed in Annex II to the CBAM
> Regulation.**

And on the wider context:

> Direct reduction involves the production of solid primary iron from high grade iron ores
> (pellets, sinter or concentrates). There are different technologies that may use
> different qualities of ores (which may require pelletisation or sintering) and different
> fuels and reducing agents (natural gas, diverse fossil fuels or biomass, hydrogen). The
> solid product is called direct reduced iron (DRI). Different types of DRI are produced,
> for example 'iron sponge' and hot briquetted iron (HBI). Some DRI is used directly as a
> feedstock in EAFs or for other downstream processes. **It is expected that production
> routes using hydrogen will play a major role in decarbonising the steel industry in
> coming years.**

So HBI briquetting is **not** a separate CBAM production process — it is part of the DRI
aggregated goods category. Its energy use sits inside the DRI boundary, direct emissions
only; the electricity used for briquetting is indirect and therefore excluded.

Guidance 5d Table (§2.2.2, precursor mapping) confirms the precursor set:

| Aggregated goods category | Possible precursors |
|---|---|
| Pig iron | Hydrogen, sintered ore, ferro alloys, pig iron/DRI |
| **Direct Reduced Iron (DRI)** | **Hydrogen, sintered ore, ferro alloys, pig iron/DRI** |
| Crude steel | Ferro alloys, pig iron, DRI, crude steel |
| Iron or steel products | Ferro alloys, pig iron, DRI, crude steel, iron or steel products |

### 2.6 Point 3.15 — Crude steel

> **3.15.1. Special provisions**
>
> The system boundary shall cover all necessary activities and units for obtaining crude
> steel:
> - if the process starts from hot metal (liquid pig iron), the system boundary shall
>   include the basic oxygen converter, vacuum degassing, secondary metallurgy, argon oxygen
>   decarburisation / vacuum oxygen decarburisation, continuous casting or ingot casting,
>   where relevant hot-rolling or forging, and all necessary auxiliary activities such as
>   transfers, re-heating, and flue gas cleaning;
> - **if the process uses an electric arc furnace, the system boundary shall include all
>   relevant activities and units such as the electric arc furnace itself, secondary
>   metallurgy, vacuum degassing, argon oxygen decarburisation / vacuum oxygen
>   decarburisation, continuous casting or ingot casting, where relevant hot-rolling or
>   forging, and all necessary auxiliary activities such as transfers, heating of raw
>   materials and equipment, re-heating, and flue gas cleaning;**
> - only primary hot-rolling and rough shaping by forging to obtain the semi-finished
>   products under CN codes 7207, 7218 and 7224 are included in this aggregated goods
>   category. All other rolling and forging processes are included in the aggregated goods
>   category 'iron or steel products'.
>
> **3.15.2.1. Basic oxygen steelmaking** — direct emissions monitoring shall encompass:
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from fuels such as coal, natural gas, fuel oils, waste gases such as blast furnace gas,
>   coke oven gas or converter gas;
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from process materials such as limestone, magnesite, and other carbonates, carbonate
>   ores; materials for flue gas cleaning;
> - **carbon entering the process in scrap, alloys, graphite etc. and carbon remaining in
>   the product or in slags or wastes is taken into account by using a mass balance method
>   in accordance with point B.3.2 of Annex III [read: Annex II].**
>
> **3.15.2.2. Electric arc furnace** — direct emissions monitoring shall take into account:
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from fuels such as coal, natural gas, fuel oils, as well as from waste gases such as
>   blast furnace gas, coke oven gas or converter gas;
> - **all processes directly or indirectly linked to the production processes emitting CO2
>   from the consumption of electrodes and electrode pastes;**
> - all processes directly or indirectly linked to the production processes emitting CO2
>   from process materials such as limestone, magnesite, and other carbonates, carbonate
>   ores; materials for flue gas cleaning;
> - **carbon entering the process, e.g. in the form of scrap, alloys and graphite, and
>   carbon remaining in the product or in slags or wastes is taken into account by using a
>   mass balance method in accordance with point B.3.2 of Annex III [read: Annex II].**

Note: an EAF's electricity is *indirect* and therefore excluded for crude steel. What *is*
charged: carbon in electrodes and electrode paste, carbon in the metallic charge that
oxidises to CO2, carbonate fluxes, and any fossil fuel burned (burners, ladle heating).

For the modelled configuration this matters twice over. The **EU** EAF is inside the EU ETS,
not CBAM, so the carbon carried in the imported HBI oxidises on the EU side and is paid for
by the EU steelmaker's ETS allowances. Only the HBI's *own* embedded emissions cross the
border.

### 2.7 What counts as a "source stream" — Article 1 definitions

> (3) **'system boundary'** means the group of chemical or physical processes included in
> the calculation of embedded emissions of goods under the same aggregated goods category;
>
> (6) **'production route'** means a specific technology used in a production process to
> produce goods;
>
> (7) **'precursor'** means any input material into a production process included in the
> list of goods set out in Annex I to Regulation (EU) 2023/956;
>
> (8) **'source stream'** means either of the following:
> (a) a specific fuel type, raw material or product giving rise to emissions of relevant
> greenhouse gases at one or more emission sources as a result of its consumption or
> production;
> (b) **a specific fuel type, raw material or product containing carbon and included in the
> calculation of greenhouse gas emissions using a mass balance method;**
>
> (9) **'emission source'** means a separately identifiable part of an installation or a
> process within an installation, from which relevant greenhouse gases are emitted;
>
> (12) **'activity data'** means the amount of fuels or materials consumed or produced by a
> process relevant for the calculation-based methodology, expressed in terajoules, mass in
> tonnes or, for gases, volume in normal cubic metres, as appropriate.

For a hydrogen-DRI plant the relevant source streams are: any natural gas or other fossil
fuel burned (start-up heaters, reformer trim, top-gas fuel), carbonate fluxes, carbon
contained in the ore, carbon in the product (negative), carbon in slag/dust (negative), plus
any carburising carbon added to bring HBI to EAF spec. Hydrogen as reducing agent yields
water, not CO2, so it contributes nothing to direct emissions of the DRI process; it enters
only as a precursor carrying its own SEE (≈0 if electrolytic).

### 2.8 The mass balance method — Annex II point B.3.2

> **B.3.2. Mass balance method**
>
> The CO2 quantities relevant for each source stream shall be calculated based on the carbon
> content in each material, **without distinguishing fuels and process materials**. **Carbon
> leaving the installation in products instead of being emitted is taken into account by
> output source streams, which have therefore negative activity data.**
>
> The emissions corresponding to each source stream shall be calculated as follows:
> (Equation 12)
>
> Where:
> - AD_k is the activity data [t] of material k; **for outputs, AD_k is negative**;
> - f is the ratio of the molar masses of CO2 and C: **f = 3,664 t CO2/t C**; and
> - CC_k is the carbon content of material k (dimensionless and positive).
>
> If the carbon content of a fuel k is calculated from an emission factor expressed in
> t CO2/TJ, the following equation shall be used: (Equation 13)
> If the carbon content of a material or fuel k is calculated from an emission factor
> expressed in t CO2/t, the following equation shall be used: (Equation 14)
> For mixed fuels, the zero-rated biomass fraction may be taken into account, provided that
> the criteria provided in point B.3.3 are met as follows: (Equation 15)
>
> Where CC_pre,k is the preliminary carbon content of fuel k (i.e. emission factor assuming
> the total fuel is fossil) and BF_k is the zero rated biomass fraction of fuel k
> (dimensionless).
>
> **For fossil fuels or materials and where the biomass fraction is not known, BF shall be
> set to the conservative value zero. Where biomass is used as input material or fuel, and
> output materials contain carbon, the overall mass balance shall treat the biomass fraction
> conservatively**, meaning that the total mass of carbon corresponding to the zero-rated
> carbon fractions of the carbon contained in all relevant output materials is not lower
> than the total mass of zero-rated fractions of the carbon contained in input materials and
> fuels, except if the operator provides evidence of a lower biomass fraction in the output
> materials by a 'trace the atom' (stoichiometric) method or by carbon-14 analyses.

In words: emissions = 3.664 × Σ_k (mass_k × carbon_content_k), with **outputs entering
negatively**.

Two modelling consequences:

1. **Carbon in the HBI is subtracted.** HBI is commonly carburised to 1–4 % C for EAF
   melting. Under B.3.2 that carbon leaves as a negative output stream and therefore
   *reduces* the DRI installation's reported direct emissions. It reappears downstream in
   the EU EAF's own mass balance, where it oxidises — but that is EU ETS territory, not
   CBAM.
2. **Biochar carburisation is conservatively treated.** If a bio-based carbon source is
   used, the operator must show the biogenic fraction in the *outputs* is not overstated,
   via stoichiometric tracing or C-14 analysis. Without such evidence, the biogenic fraction
   assigned to outputs must be at least as large as that in the inputs — which caps the
   negative-output credit you can claim for biogenic carbon.

The general B.3 hierarchy: Annex II point B.2 offers a calculation-based methodology
(standard method or mass balance) and a measurement-based methodology (CEMS), and requires
*"the monitoring methodology that gives the most accurate and reliable results"* unless the
sector-specific rules in B.9 mandate one. For DRI, points 3.14.2 mandate the mass balance
for product/slag/waste carbon.

### 2.9 Is transport inside or outside the boundary? — **OUTSIDE**

Three independent confirmations:

1. **Methodology Act, Annex II, point B.1 (completeness of source streams):**
   > The boundaries of the installation and its production processes shall be clearly known
   > to the operator and defined in the monitoring plan […] The following principles shall
   > apply:
   > (a) as a minimum, all relevant greenhouse gas emissions emission sources and source
   > streams associated directly or indirectly with the production of goods listed in point 2
   > of Annex I shall be covered;
   > (b) all emissions from regular operations shall be included, as well as from abnormal
   > events, including start-up, shut-down and emergency situations, over the reporting
   > period;
   > **(c) emissions from mobile machinery for transportation purposes shall be excluded.**

2. **Methodology Act, Annex I, point 3.1:** *"The purchase and maintenance of infrastructure
   and equipment are excluded from the system boundaries."*

3. **Regulation 2023/956, Article 30(2)(a)(ii)** lists *"embedded emissions in the transport
   of the goods listed in Annex I and transportation services"* as something the Commission
   must *assess the possibility of extending the scope to* — i.e. it is presently outside.

**Consequence for the model:** shipping HBI from Australia or Brazil to an EU port carries
**zero CBAM liability**. Neither bunker fuel nor port handling nor inland haulage enters
embedded emissions, and on-site mobile machinery (haul trucks, loaders, stackers) is
explicitly excluded even inside the fence line. EU-leg maritime emissions are separately in
the EU ETS since 2024 under the maritime extension, but that is the shipping company's
obligation, priced into the freight rate, not part of the good's embedded emissions and not
part of the CBAM declarant's surrender obligation. Any freight-emission term in an LCOI/LCOS
model is therefore a *voluntary* accounting choice, not a CBAM cost — though it does show up
indirectly through ETS-inflated EU-leg freight rates.

---

## 3. PRECURSORS — how embedded emissions flow into complex goods

### 3.1 The recursion — Annex III, point B

> **B. CALCULATION OF SPECIFIC EMBEDDED EMISSIONS OF COMPLEX GOODS**
>
> In accordance with Annex IV to Regulation (EU) 2023/956, the specific embedded emissions
> SEE_g of complex goods g shall be calculated as follows: (Equations 59, 60)
>
> Where:
> - SEE_g are the specific **direct or indirect** embedded emissions of (complex) goods g
>   expressed in t CO2e per functional unit;
> - AttrEm_g are the attributed direct or indirect emissions of the production process
>   yielding goods g determined in accordance with point A.3 of this Annex for the reporting
>   period, expressed in t CO2e;
> - AL_g is the activity level of the production process yielding goods g …;
> - EE_InpMat are the embedded direct or indirect emissions of all precursors consumed
>   during the reporting period, expressed in t CO2e;
> - M_i is the mass of precursor i used in the production process yielding g during the
>   reporting period, expressed in functional units of precursor i, and
> - SEE_i are the specific direct or indirect embedded emissions of precursor i expressed in
>   t CO2e per functional unit of precursor i.
>
> In this calculation, only precursors not covered by the same production process as goods g
> are taken into account. Where the same precursor is obtained from different production
> processes, the precursor from each installation shall be treated separately.
>
> **If a precursor i originates in the Union or in one of the countries or territories
> exempted pursuant to point 1 of Annex III to Regulation (EU) 2023/956 the specific direct
> or indirect embedded emissions of that precursor shall be counted as zero.**
>
> Where a precursor i itself has precursors, those precursors are first taken into account
> using the same calculation method in order to calculate the embedded emissions of the
> precursor i before they are used for calculating the embedded emissions of goods g. **This
> method is used recursively to all precursors which are complex goods.**
>
> **The parameter M_i refers to the total mass of precursor required to produce the amount
> AL_g. It also includes quantities of the precursor which do not end up in the complex goods
> but may be spilt, cut off, combusted, chemically modified, etc. in the production process
> and leave the process as by-products, scrap, residues, wastes, or emissions.**
>
> In order to provide data which can be used independently of activity levels, the specific
> mass consumption m_i for each precursor i shall be determined … (Equation 61)
> Thereby the specific embedded emissions of complex goods g may be expressed as: (Equation 62)
>
> Where ae_g are the specific attributed direct or indirect emissions of the production
> process yielding goods g, expressed in t CO2e per tonne of g, being equivalent to specific
> embedded emissions without precursors' embedded emissions … and m_i is the specific mass
> consumption of precursor i used in the production process yielding one functional unit of
> goods g …

In compact form: **SEE_g = ae_g + Σ_i (m_i × SEE_i)**, where `ae_g` is the process's own
specific attributed emissions and `m_i` the specific mass consumption of precursor i per
functional unit of g. The whole thing is evaluated twice in parallel — once for the direct
channel, once for the indirect channel — and for Annex II goods the indirect channel is
simply not evaluated for the good itself, and is dropped for any Annex II precursor.

The `M_i` rule is worth flagging for hydrogen: **all** hydrogen fed to the shaft counts,
including recycled/purged excess that does not chemically end up in the iron. For a shaft
operating at, say, 1.8–2.2× stoichiometric H2, `M_i` is the full quantity fed unless
unreacted hydrogen is recovered *within the same production process*, in which case it never
crosses a process boundary and is not double counted. This only matters if `SEE_hydrogen` is
non-zero (i.e. fossil hydrogen); for electrolytic hydrogen the whole term vanishes.

### 3.2 What value does hydrogen carry as a DRI reductant?

`SEE_hydrogen` = its **direct emissions only**, per:
- Article 7(1) of Regulation (EU) 2023/956 + Annex II (hydrogen 2804 10 00 listed);
- Article 3(2) of the Methodology Act;
- and the third paragraph of Annex I point 3.1 of the Methodology Act:

> When the production process of complex goods listed in Annex II … includes one or more
> precursors listed in that Annex, the indirect emissions of these precursors will not be
> included …

Since both DRI (CN 7203, inside Chapter 72) and hydrogen (2804 10 00) are Annex II goods, the
hydrogen precursor contributes **direct emissions only**.

For **electrolytic hydrogen** those direct emissions are effectively zero: water electrolysis
has no fuel combustion and no carbonate process emissions, and — per Guidance 5b — the
Methodology Act does not even define an electrolysis route because it involves "minimal
direct emissions". The electricity driving the electrolyser is indirect and is excluded twice
over: once because hydrogen is Annex II, and again because DRI is Annex II.

For **SMR/ATR/coal-gasification hydrogen** the direct emissions are large and flow straight
through into the HBI. The default (see §4) is 10.82 t CO2e per tonne of H2 for both Australia
and Brazil — roughly SMR-with-no-CCS territory.

**Order-of-magnitude sanity check for the model.** At roughly 55–60 kg H2 per tonne of DRI
(≈ 0.055–0.060 t/t), fossil hydrogen at the country default of 10.82 t CO2e/t H2 would add
about **0.60–0.65 t CO2e per tonne of HBI** on the precursor term alone. Electrolytic
hydrogen adds ≈ 0. That difference is the single largest lever in the CBAM cost of the
imported HBI.

**Confirmed: indirect emissions of the hydrogen precursor are NOT included.**

### 3.3 Joint production process — the pelletising escape hatch

> **Article 4(9)** — Where precursors relevant for complex goods are produced in the same
> installation as the complex goods, and where the respective precursors are not transferred
> out for sale or use in other production processes, the production of precursors and complex
> goods may be covered by a **joint production process**. In that case, monitoring and
> calculation of embedded emissions of the precursors and complex goods shall be carried out
> jointly.

Annex I point 3.14.1 repeats this for DRI ("Where the installation does not sell or transfer
DRI to other installations, a joint production process including steel can be established").
Guidance 5d §2.2.3.5 repeats it again.

This is directly relevant. An integrated pellet plant + H2-DRI + briquetting site that does
not sell pellets can define one joint production process, so that CN 2601 12 00 never appears
as a precursor and its indirect (electricity) emissions never enter the calculation.
Conversely, a project that **buys** DR-grade pellets from a third party inherits those
pellets' indirect emissions — and pelletising is electricity-hungry (induration fans, grinding),
so on a grid-average emission factor this is not negligible.

Order of magnitude from the Commission's own defaults: for Australia, 2601 12 00 carries
0.070 direct + **0.020 indirect** = 0.090 t CO2e/t pellets; for Brazil, 0.190 direct +
**0.010 indirect** = 0.200 t CO2e/t. At ≈1.45 t pellets per tonne of DRI the *indirect*
component alone is ≈0.029 t CO2e/t HBI (AU) or ≈0.015 (BR) — small, but non-zero and it does
not disappear just because the shaft runs on green hydrogen.

Also note **Article 4(7)** — the anti-splitting rule:

> Splitting an installation into different installations, with the result that production
> routes otherwise pertaining to a single production process are carried out in separate
> installations, shall only be allowed where the operators demonstrate valid commercial
> reasons for this split that are related to their economic activity. Commercial reasons shall
> be considered as valid where circumventing Regulation (EU) 2023/956 is not their main
> purpose or one of their main purposes.

And **Article 4(6)**: *"Where goods to which the same functional unit applies are produced
using different production routes within an installation, a single production process shall be
used encompassing all production routes."* — i.e. you cannot run a green shaft and a grey
shaft side by side on the same site and export only the green tonnes as a separate product.
They collapse into one production process with one blended SEE. (See §9.)

### 3.4 Actual vs default for precursors — Annex II, point D and point E

> for those goods whose production processes include precursors, making those goods 'complex
> goods', the embedded emissions of the precursor shall be determined in accordance with
> point E of this Annex, and shall be added to the embedded emissions of the complex goods
> produced, by applying the rules provided in point B of Annex III. **Where precursors are
> themselves complex goods, that process shall be repeated recursively until no more
> precursors are at stake.**

> For precursors produced outside the installation and originating in third countries and
> territories that are not exempted pursuant to point 1 of Annex III to Regulation (EU)
> 2023/956, **actual data obtained from the operator of the installation producing the
> precursor shall be used only if the following conditions are met**:
> (a) the data must be taken from a verification report that has been issued by a verifier
> having an accreditation in accordance with Article 18 of Delegated Regulation (EU)
> 2025/2551 valid at the time of issuing the verification report and for the sectoral scope
> required for the aggregated goods category of the precursor under consideration; and
> (b) the verification report must cover the reporting period during which the precursor was
> produced.
>
> **Where the operator does not have a verification report meeting conditions (a) and (b), the
> relevant default values, made available in accordance with Annex IV of Regulation (EU)
> 2023/956, for the precursor shall be used.**

Annex II point E ("Monitoring of precursors") requires, for a precursor obtained from another
installation:
> - specific embedded direct and indirect emissions of the precursor as average over the
>   reporting period, expressed in tonnes CO2e per tonne of precursor;
> - quantity of the precursor used in each production process of the installation.

Article 15 allows mixing: *"The specific embedded emissions for complex goods may be
calculated by determining actual emissions for the production processes within the
installation producing the complex goods, and default values for one or more precursors of
the complex goods."*

Articles 13 and 14 fix the reporting period and averaging:
> **Article 13** — The default reporting period of a precursor shall be the year of production
> of the complex good. However, where operators provide the verifier with sufficient evidence
> to identify the actual time of production, the reporting period shall be the period during
> which the precursor was produced.

> **Article 14(1)** — Where an installation producing complex goods receives, from another
> installation, precursors under a given CN code produced during different reporting periods,
> the embedded emissions of the complex goods shall … be determined as the **weighted average**
> …
> **(2)** Where an installation producing complex goods receives precursors under a given CN
> code from multiple installations, the embedded emissions … shall … **by default** be
> determined as the **weighted average** …
> **(3)** Where operators provide the verifier with sufficient evidence demonstrating that,
> out of the precursors under a given CN code received from multiple installations, the
> installation producing the complex goods used, for a given production process, only
> precursors from a single installation, or from a subset of installations, the embedded
> emissions of those precursors … shall be determined … based on the embedded emissions of the
> precursors obtained from that single installation …

Article 14(3) is the precursor analogue of the electricity rule in Article 9(2), and is the
main structural opening for physical-allocation choices (see §9).

**Verification is the gate.** The whole "green HBI" case rests on getting verified actual
values. Without an accredited verifier's report for the reporting period, the declarant is
pushed onto the default of 1.325 t CO2e/t (Australia and Brazil both fall back to "Other
Countries and Territories"), plus mark-up — which would destroy the economics of a green
project. Accreditation must be by a National Accreditation Body in the European Accreditation
network under EN ISO/IEC 14065, per Delegated Regulation (EU) 2025/2551.

### 3.5 THE GUIDANCE 5d ERROR — confirmed

Guidance No. 5d, "Sector-specific guidance document on iron and steel", p. 8, states:

> Where iron and steel CBAM goods are complex goods produced using precursors (such as
> sintered ore, pig iron, DRI or ferro-alloys), the embedded emissions of the iron and steel
> goods must include the full embedded emissions of those precursors, including any indirect
> emissions of electricity associated with the precursors where relevant (i.e. where the
> precursor is **not** listed in Annex II to the CBAM Regulation), **for example, hydrogen
> where used as a precursor**.

**This example is wrong.** Hydrogen, CN 2804 10 00, *is* listed in Annex II to Regulation (EU)
2023/956, under "Chemicals" — it has been since the original 2023 text, and Regulation (EU)
2025/2083 did not remove it (it only *added* electricity). The correct example of a
*non*-Annex-II iron-and-steel precursor is **sintered ore, CN 2601 12 00**, which the very same
guidance document correctly identifies two paragraphs earlier:

> For sintered ore (CN 2601 12 00) – not listed in Annex II to the CBAM Regulation, indirect
> emissions are to be included.

and again, correctly, in the general statement immediately preceding the faulty example:

> In the definitive period, in line with Article 7(1) and Annex II ("List of goods for which
> only direct emissions are to be taken into account") to the CBAM Regulation and Article 3(2)
> of the Methodology Act, only direct emissions are taken into account when determining the
> embedded emissions of iron and steel CBAM goods (**except for sintered ore**). Where iron and
> steel goods listed in Annex II are used as precursors in the production of other CBAM goods,
> their indirect emissions are likewise not included in the embedded emissions of those complex
> goods, in accordance with the cross-sectoral rules in the Methodology Act.

**Which is right:** Annex II of the Regulation, and Guidance 5b. Guidance 5b is unambiguous and
repeats the point at least six times, e.g.:

> In the definitive period, for CBAM purposes the hydrogen sector has to account for direct
> emissions only, since hydrogen is listed in Annex II to the CBAM Regulation as a good for
> which only direct emissions are taken into account.

> **Notes to figure:** For hydrogen (Annex II goods), only direct emissions are taken into
> account for embedded emissions; indirect emissions from electricity are not included.

> Where hydrogen is used as a precursor in the production of ammonia and the conditions laid
> down in the Methodology Act for disregarding indirect emissions of precursors are met (e.g.
> where the precursor itself is listed in Annex II), the indirect emissions of electricity
> associated with the hydrogen precursor may be disregarded when determining the embedded
> emissions of the ammonia CBAM goods.

Note that 5b's drafting is itself hedged ("where the conditions … are met", "may be
disregarded"), whereas the Methodology Act's Annex I point 3.1 is mandatory in form ("will not
be included") and Article 3(2) is the binding rule. Guidance documents have no legal force —
the CBAM guidance PDFs carry the standard disclaimer that only the Court of Justice of the
European Union can give an authoritative interpretation of Union law. **On any conflict,
Article 7(1) of 2023/956 + Annex II + Article 3(2) of 2025/2547 govern, and hydrogen carries
direct emissions only.**

The independent cross-check is the Commission's own numbers: in Annex I of IR 2025/2621 (as
replaced by IR 2026/1740), the row for `2804 10 00 – Hydrogen` shows **`N/A` in the
indirect-emissions column** for every country listed, exactly as for CN 72xx and unlike CN
2601 12 00. The Commission's default-value tables treat hydrogen as an Annex II good.

---

## 4. DEFAULT VALUES AND BENCHMARKS

Two different numbers, often confused, both needed to compute a CBAM bill:

- **Default value** (Implementing Regulation (EU) 2025/2621, as corrected by (EU) 2026/1740)
  = the fallback *embedded emissions* of the imported good, per tonne, when verified actual
  data is unavailable. This is the **gross liability**.
- **CBAM benchmark** (Implementing Regulation (EU) 2025/2620) = the EU-ETS-derived
  free-allocation equivalent, per tonne, which is **subtracted** (scaled by the CBAM factor
  and the cross-sectoral correction factor). This is the **relief**.

Certificates due ≈ (embedded emissions − free allocation adjustment) × mass, where the free
allocation adjustment = CBAM_factor(y) × CSCF(y) × CBAM_benchmark.

### 4.1 The default values regulation

> **COMMISSION IMPLEMENTING REGULATION (EU) 2025/2621 of 16 December 2025 laying down rules
> for the application of Regulation (EU) 2023/956 of the European Parliament and the Council
> as regards the establishment of default values** (OJ L, 2025/2621, 31.12.2025).
> http://data.europa.eu/eli/reg_impl/2025/2621/oj

> **Article 1 — Default values**
>
> 1. Where the embedded emissions in imported goods are determined on the basis of default
> values in accordance with Article 7(2), point (b), of Regulation (EU) 2023/956, the default
> values laid down in **Annex I** to this Regulation shall be used.
>
> 2. Where the embedded emissions of complex goods are determined on the basis of actual
> values, and the embedded emissions of precursors used in the production of those complex
> goods are determined on the basis of default values in accordance with Article 15 of
> Commission Implementing Regulation (EU) 2025/2547, the default values laid down in Annex I
> shall be used for those precursors.
>
> 3. Where the specific indirect embedded emissions are determined on the basis of default
> values in accordance with Article 7(4) of Regulation (EU) 2023/956, the default values laid
> down in **Annex II** to this Regulation shall be used.
>
> 4. Where the embedded direct emissions in electricity imported into the customs territory of
> the Union are determined on the basis of default values in accordance with Article 7(3) of
> Regulation (EU) 2023/956, the default values laid down in **Annex III** to this Regulation
> shall be used.
>
> 5. By way of derogation from paragraph 2, **where a country of production cannot be
> identified for a precursor, the default values laid down in Annex IV shall be used.**

> **Article 2** — … It shall apply from 1 January 2026. **It shall be revised in 2027 at the
> latest.**

The Methodology Act says the same at Article 11(5): *"The Commission shall conduct a review of
the default values by December 2027 at the latest."*

**Superseded by the correcting act.** Implementing Regulation (EU) 2026/1740 of 20 July 2026
*"correcting Implementing Regulation (EU) 2025/2621 as regards Annexes I and IV thereto"*
(OJ L, 2026/1740, 31.7.2026, CELEX `32026R1740`) replaces both annexes in full and applies
retroactively from 1 January 2026:

> **Article 1** — 1. Annex I to Implementing Regulation (EU) 2025/2621 is replaced by the text
> in Annex I to this Regulation. 2. Annex IV to Implementing Regulation (EU) 2025/2621 is
> replaced by the text in Annex II to this Regulation.

Its recital (10) explains the most visible change:

> Rounding rules have been applied when determining the default values. In order to avoid
> potential inconsistencies as a result of rounding, **the columns with default values that
> include mark-ups in all the tables for 2026, 2027 and 2028 in Annexes I and IV … should be
> deleted**. The methodology for calculating the default values with the mark-ups should be
> stated separately in opening paragraphs in Annexes I and IV. **The final default values that
> include the mark-ups should then be calculated within the CBAM Registry** and this
> calculation should be based on the default values for total emissions as provided for in
> Annex I. The default values for direct emissions and indirect emissions in Annex I have only
> been provided for information.

Also worth knowing, recital (8), which is a direct confirmation of the sintered-ore point:

> The default value for indirect emissions for CN code 2601 12 00 for Angola in one of the
> tables in Annex I … is stated as being not applicable as a result of a transcription mistake.
> However, **it is clear that indirect emissions for this CN code are within the CBAM's scope,
> in accordance with Article 7(1) of Regulation (EU) 2023/956**.

### 4.2 The mark-up

Opening paragraph of Annex I as replaced by 2026/1740:

> For the calculation of the number of CBAM certificates, the default values of the column
> 'total emissions' shall be selected and increased as follows:
>
> For goods in the **cement, iron and steel, aluminium and hydrogen** sectors, the mark-up shall
> be **10 % for the year 2026, 20 % for the year 2027 and 30 % for the year 2028 and onwards**.
>
> For goods in the **fertiliser** sector, the mark-up shall be **1 % for the year 2026 and
> onwards**.

Rationale, recitals (23) and (25) of 2025/2621:

> To ensure the environmental integrity of the CBAM, default values for embedded emissions of
> goods other than electricity should include a mark-up to account for the deviations of an
> individual installation with emission levels higher than the relevant average emission
> intensity of the producer country. Given the difficulties to verify that installation-specific
> data from third countries is of a sufficiently high quality, a suitable proxy should be
> applied to estimate the variations of individual installations compared to the average.
> Therefore, the proposed mark-up is based on existing deviations of Union installations with
> respect to Union average altogether.

> To avoid immediate disproportionate impacts on prices of goods, and to give economic operators
> time to adapt, the mark-up should be gradually phased-in.

### 4.3 The fallback chain

Annex I opening paragraph:

> Where a country or territory is not explicitly listed, the default value for the respective
> good from the table "Other countries and territories" needs to be selected. **Where a country
> or territory is explicitly listed but no value is provided or the relevant field shows "–",
> the default value for the respective good from the table "Other countries and territories"
> needs to be selected.**

And the "Other countries and territories" values themselves are set by Annex IV, point 4.1 of
Regulation 2023/956 (as replaced by the Omnibus):

> Default values shall be set at the average emission intensity of each exporting country and
> for each of the goods listed in Annex I other than electricity, increased by a proportionately
> designed mark-up. This mark-up shall be determined in the implementing acts adopted pursuant
> to Article 7(7) … **When reliable data for the exporting country cannot be applied for a type
> of goods, the default values shall be based on the average emission intensity of the 10
> exporting countries with the highest emission intensities for which reliable data can be
> applied for that type of goods.**

That is why the fallback numbers are punitive.

### 4.4 THE NUMBERS — DRI (CN 7203), Australia and Brazil

**Both Australia and Brazil show "–" for CN 7203.** Neither country has a country-specific DRI
default. Both therefore take the "Other Countries and Territories" value.

| Table | CN 7203 direct | indirect | total | 2026 (×1.10) | 2027 (×1.20) | 2028+ (×1.30) |
|---|---|---|---|---|---|---|
| **Australia** | – | N/A | – | (fallback) | (fallback) | (fallback) |
| **Brazil** | – | N/A | – | (fallback) | (fallback) | (fallback) |
| **Other Countries and Territories** (applies to AU and BR) | **1,325** | N/A | **1,325** | **1,458** | **1,590** | **1,723** |
| Annex IV — precursor, country of production unknown | **4,200** | — | 4,200 | 4,620 | 5,040 | 5,460 |

All in t CO2e per tonne of good. The 2026/2027/2028 columns are the mark-up-inclusive values;
in the current (2026/1740) text these columns are gone from the annex and are computed by the
CBAM Registry from the "total emissions" column using the mark-up schedule above — the
arithmetic is unchanged.

**Comparison anchors from the same tables (direct = total, since Annex II):**

| CN code | Description | Australia | Brazil | Other countries |
|---|---|---|---|---|
| 2601 12 00 | Agglomerated iron ores (sintered ore / pellets) | 0,070 direct + **0,020 indirect** = 0,090 | 0,190 direct + **0,010 indirect** = 0,200 | (see table) |
| 7201 | Pig iron and spiegeleisen | 2,790 | 1,478 | — |
| 7203 | DRI / spongy ferrous products | – (→1,325) | – (→1,325) | **1,325** |
| 7205 | Granules and powders of pig iron etc. | 2,810 | 1,430 | 3,715 |
| **2804 10 00** | **Hydrogen** | **10,820** | **10,820** | **17,740** |

Note that 2601 12 00 is the **only** row in the iron-and-steel block with a number in the
indirect column — every 72xx row and the hydrogen row show `N/A`. This is the cleanest possible
empirical confirmation of the Annex II boundary.

### 4.5 THE NUMBERS — Hydrogen (CN 2804 10 00)

| Table | direct | indirect | total | 2026 | 2027 | 2028+ |
|---|---|---|---|---|---|---|
| **Australia** | **10,820** | **N/A** | **10,820** | 11,902 | 12,984 | 14,066 |
| **Brazil** | **10,820** | **N/A** | **10,820** | 11,902 | 12,984 | 14,066 |
| Other Countries and Territories | 17,740 | N/A | 17,740 | 19,514 | 21,288 | 23,062 |
| Annex IV — precursor, country unknown | 26,640 | — | 26,640 | 29,304 | 31,968 | 34,632 |

t CO2e per tonne of hydrogen. `N/A` in the indirect column for every country — hydrogen is an
Annex II good, so there is no indirect default to publish.

10.82 t CO2e/t H2 is roughly unabated steam-methane-reforming plus some upstream slack. It is
*not* a plausible number for electrolytic hydrogen, so a green project must produce verified
actual values for its hydrogen precursor (or, better, define a joint production process under
Article 4(9) so the hydrogen never appears as a separate precursor at all).

### 4.6 Indirect-emission electricity factors (Annex II of 2025/2621)

Relevant only for non-Annex-II goods, i.e. for this project **only for bought-in pellets**:

| Country | Emission factor (t CO2eq/MWh) |
|---|---|
| Australia | **0,645** |
| Brazil | **0,096** |
| Other countries and territories | 0,465 |
| European Union (Annex III, for imported electricity) | 0,612 |

Source note in the regulation: *"The data provided are based on data sourced from the
International Energy Agency (IEA) and is subject to a Creative Commons Non-Commercial
Share-Alike 4.0 CC BY NC SA licence."*

### 4.7 CBAM BENCHMARKS — Implementing Regulation (EU) 2025/2620 — verified

> **COMMISSION IMPLEMENTING REGULATION (EU) 2025/2620 of 16 December 2025 laying down rules for
> the application of Regulation (EU) 2023/956 … as regards the calculation of the free
> allocation adjustment to the number of CBAM certificates to be surrendered** (OJ L, 2025/2620,
> 22.12.2025). http://data.europa.eu/eli/reg_impl/2025/2620/oj

The Annex, point 5.3, is the benchmark table. Two columns:
- **Column A — BMg\*** = the *process-related* CBAM benchmark (the production process yielding
  the good, excluding precursors). Used when computing SEFA from **actual** data.
- **Column B — BMg** = the *default* CBAM benchmark (whole good including precursors). Used when
  the declaration uses **default** values.

Indicator legend (point 5.3):

> (C) Carbon Steel based on BF/BOF
> (D) Carbon Steel based on DRI/EAF
> (E) Carbon Steel based on Scrap/EAF
> (F) Low alloy Steel based on BF/BOF
> (G) Low alloy Steel based on DRI/EAF
> (H) Low alloy Steel based on scrap/EAF
> (J) High alloy Steel (based on EAF)
> (K) primary Aluminium
> (L) secondary Aluminium
> (A) grey clinker / cement, (B) white clinker / cement
> (1) Value is to be used for production years 2026-27
> (2) Value is to be used for production years 2028-30

**Verification of the 1,370 / 0,481 / 0,072 triple.** These numbers are real and they are
Column B (default CBAM benchmark, whole good) — but they are **not** the crude steel row. They
are the values for **hot-rolled flat products of CN 7208** (and the other 72xx flat/long
product codes). The crude-steel and semi-finished rows carry different numbers:

| CN code | Column A (BMg\*, process-related) | Column B (BMg, default, whole good) |
|---|---|---|
| **7206 10 00 / 7206 90 00** (ingots, crude steel) | 0,150 (C) / 0,027 (D) / 0,027 (E) | **1,288 (C) / 0,424 (D) / 0,027 (E)** |
| **7207 11 11 / 11 14 / 11 16 / 12 10 …** (semi-finished, rolled or continuous cast) | 0,188 (C) / 0,065 (D) / 0,065 (E) | **1,364 (C) / 0,475 (D) / 0,066 (E)** |
| 7207 11 90 (semi-finished, forged) | 0,453 (C) / 0,330 (D) / 0,330 (E) | 1,629 (C) / 0,740 (D) / 0,331 (E) |
| **7208 10 00, 7208 25 00, 7208 26 00, 7208 27 00 …** (hot-rolled flat) | 0,044 | **1,370 (C) / 0,481 (D) / 0,072 (E)** |
| **7203 10 00 / 7203 90 00 (DRI)** | **0,295** | **0,397** |
| **2804 10 00 (Hydrogen)** | **5,089** | **5,089** |

t CO2e per tonne. So: **the "1,370 / 0,481 / 0,072" attribution to Implementing Regulation
2025/2621 in the brief is wrong on two counts.** The numbers come from Implementing Regulation
**2025/2620** (the free-allocation-adjustment act), not 2025/2621 (default values); and within
2025/2620 they are the Column B benchmarks for hot-rolled flat steel, not for crude steel. The
"BF-BOF / DRI-EAF / scrap-EAF" reading of the (C)/(D)/(E) indicators is correct.

**The DRI benchmark, 0,295 / 0,397, is the number that matters for HBI imports.** Recital (47)
of 2025/2620 explains why a dedicated DRI benchmark had to be invented:

> Currently, direct reduced iron (DRI) is covered by the hot metal benchmark of the EU ETS and,
> without further differentiation, imports of steel based on natural-gas DRI would receive a free
> allocation adjustment that exceeds their embedded emissions for the first years in which a CBAM
> obligation is due, which means that no CBAM certificates would be due for DRI-based goods.
> Compared to this, secondary steel imports would face a CBAM obligation, despite having lower
> actual embedded emissions than DRI. In addition, the potential free allocation adjustment
> stemming from the hot metal benchmark would create a situation in which more carbon-intensive
> natural gas based DRI imports would receive more free allocation adjustment than secondary
> steel producers, which receive comparably less free allocation based on the electric arc
> furnace (EAF) benchmarks, therefore increasing the carbon leakage risk for secondary steel
> producers in the Union. In line with the principles to be applied for the allocation of free
> allowances under the EU ETS, to ensure the environmental integrity of the CBAM and to address
> the potential risk of carbon leakage of the production of secondary steel in the Union, **a
> dedicated CBAM benchmark for natural gas-based DRI should be created. Taking into account the
> relative level of embedded emissions, the level of the DRI benchmark should be chosen to ensure
> that the CBAM obligation for primary natural gas-based DRI imports is lower than for primary
> blast furnace steel, but higher than for secondary steel.**

And recital (49) — the reason CBAM benchmarks for steel are *lower* than the corresponding ETS
benchmarks:

> The rules for the exchangeability of fuel and electricity have been removed for the
> determination of free allocation under the EU ETS starting in 2026. This means that free
> allocation granted under some ETS product benchmarks in the steel sector will cover indirect
> emissions to a certain extent. **As the CBAM scope currently only covers direct emissions in
> the steel sector, only the direct emission share of the respective ETS benchmarks should be
> considered when determining the corresponding CBAM benchmarks.** For these benchmarks,
> improvement rates in accordance with [points (c) and (d) of the third subparagraph of] Article
> 10a(2) of Directive 2003/87/EC do not appropriately reflect the direct emissions to be covered
> by the CBAM benchmarks. Therefore, the average direct emissions of the 10 % best installations
> under these ETS benchmarks in the new baseline years 2021 and 2022 should serve as proxy
> instead of the EU ETS benchmarks for the purpose of calculating the CBAM benchmarks.

**Benchmarks are provisional for 2026.** Recital (35):

> EU ETS benchmarks for the period from 2026 to 2030 will only become available in early 2026.
> However, based on the already collected data, estimates with a high degree of certainty can
> already be made … **the CBAM benchmarks in 2026 should be based on the estimated ETS benchmarks
> to be applied in the period from 2026 to 2030 and apply from 1 January 2026. These CBAM
> benchmarks should be reviewed at the latest one month after the final EU ETS benchmarks for the
> period from 2026 to 2030 are published. The updated CBAM benchmarks based on the final ETS
> benchmarks for the period 2026 to 2030 should apply to goods imported from 1 January 2027.**

The Commission publishes a working spreadsheet of the benchmarks:
https://taxation-customs.ec.europa.eu/document/download/9877523c-2a02-4926-a211-aefae7cf6d0d_en
(`CBAM Benchmarks_20260206.xlsx`)

### 4.8 The free allocation adjustment formula

> **Article 1** — 1. The adjustment to the number of CBAM certificates, as referred to in
> Article 31 of Regulation (EU) 2023/956 ('free allocation adjustment'), shall be calculated in
> accordance with point 2 of the Annex to this Regulation. 2. **The free allocation adjustment
> for electrical energy (CN code 2716 00 00) shall be zero.**

Annex, point 2:

> FAA_g = SEFA_g,y × M_g
>
> where SEFA_g,y is the specific embedded free allocation of good g in year y (t CO2e / tonne),
> M_g the mass imported, y the reporting period per Article 7 of IR 2025/2547.

Annex, point 3.1 (actual data):

> SFA_Proc_g,y = CBAM_y × CSCF_y × BM_g\*
>
> where CBAM_y is **the CBAM factor referred to in Article 10a(1a) of Directive 2003/87/EC** for
> year y (dimensionless); CSCF_y is the cross-sectoral correction factor for year y determined
> pursuant to Article 14(6) of Delegated Regulation (EU) 2019/331 and published under Article
> 10a(5) of Directive 2003/87/EC; BM_g\* is the process-related CBAM benchmark from Column A.

Point 4 does the same with Column B when defaults are used. Point 5.1:

> Where default values are used to determine SEFA of a final good or of a precursor, the same
> production route shall be used as indicated in Annex I to Commission Implementing Regulation
> (EU) 2025/2621 for the country of origin of that good or precursor. Where different alloy
> grades for steel are given in the table for the same CN code, **the highest benchmark value
> given for the relevant production year is used.**

Article 4 mirrors the Methodology Act's Article 14 for precursors from multiple installations
(weighted average by default; single-installation attribution on sufficient evidence).

### 4.9 A structural asymmetry worth modelling: precursor benchmarks

Point 3.3 of the Annex gives the free allocation adjustment for a **complex** good under actual
data:

> For a complex good, the calculation of the SEFA shall take into account the production process as
> well as the SEFA of each precursor and shall be calculated using the following equation:
> (Equation 4)
>
> Where … SFA_Proc_g,y is the process-related specific free allocation … m_i is the specific mass of
> precursor i consumed for the production of one tonne of good g, as determined in accordance with
> the rules set out in Implementing Regulation (EU) 2025/2547 … SEFA_i,y' is the specific embedded
> free allocation of precursor i …
>
> Where precursors are themselves complex goods, the calculation of SEFA_i,y' shall be repeated
> recursively using Equations 2, 3 and 4, as appropriate, until no more precursors are relevant.

and, where the producer has not supplied a verified value for the precursor, SEFA_i is derived from
**Column B** of the benchmark table for that precursor's CN code and country of origin.

So: **SEFA_HBI = CBAM_y × CSCF_y × 0,295  +  Σ_i m_i × SEFA_i.**

The benchmark values for the two candidate precursors are:

| Precursor | CN code | Column A | Column B |
|---|---|---|---|
| Sintered ore / pellets | 2601 12 00 | 0,086 | 0,086 |
| Hydrogen | 2804 10 00 | 5,089 | 5,089 |

**The hydrogen benchmark is very large relative to green hydrogen's actual embedded emissions
(≈0).** At a typical ≈0,057 t H2 per tonne of DRI, treating hydrogen as a *bought-in precursor*
adds roughly `0,057 × 5,089 ≈ 0,29 t CO2e/t HBI` of free-allocation relief on top of the DRI
process benchmark of 0,295 — i.e. it roughly **doubles** the relief, to ≈0,585 t CO2e/t.

By contrast, folding the electrolyser into a **joint production process** under Article 4(9) of the
Methodology Act means hydrogen is not a precursor at all, and the whole good sits under the single
DRI benchmark of 0,295.

**Why this matters, and why it usually does not.** For a genuinely near-zero-emissions green HBI
plant the surrender obligation is floored at zero either way, so the choice is economically
irrelevant. It becomes material in three cases the model should be able to represent:

1. **Blended or partially fossil reductant.** A shaft on a blended H2/natural-gas reductant has
   real direct emissions; whether hydrogen is a separate precursor then determines whether the
   0,29 t/t of extra relief is available. Buying the hydrogen from a separate installation is worth
   money.
2. **After 2030**, as the CBAM factor decays, the *absolute* size of the relief shrinks
   proportionally, so the gap between the two structures narrows in euro terms even as the gross
   liability rises.
3. **Sintered ore cuts the other way.** Its benchmark (0,086) is small, while buying it in drags in
   its indirect emissions (0,020 t CO2e/t for Australia, 0,010 for Brazil, per §4.4). At ≈1,45 t
   pellets/t DRI the relief is ≈0,125 t CO2e/t and the added liability ≈0,029 (AU) — so on these
   numbers buying pellets is *net favourable* on paper, though that conclusion is sensitive to the
   actual verified pellet emissions rather than the defaults.

**Caveat.** The extra precursor relief is only available where the precursor's SEFA can be
established — either a producer-supplied verified value or the Column B benchmark for that
precursor's CN code and country of origin. And the CBAM benchmarks for 2026 are provisional
(recital (35) of 2025/2620); they are to be revised once the final 2026-2030 EU ETS benchmarks are
published, applying to imports from 1 January 2027. Do not treat 0,295 / 0,397 / 5,089 / 0,086 as
fixed through 2030.

---

## 5. THE 2025–2026 CHANGES

Two separate things are often run together. **Regulation (EU) 2025/2083** (October 2025) is
enacted law — the "Omnibus"/simplification package. **COM(2025) 989** (December 2025) is a
*proposal* on scope extension and anti-circumvention; it is not law and is subject to the
ordinary legislative procedure.

### 5.1 The Omnibus — Regulation (EU) 2025/2083

> **Regulation (EU) 2025/2083 of the European Parliament and of the Council of 8 October 2025
> amending Regulation (EU) 2023/956 as regards simplifying and strengthening the carbon border
> adjustment mechanism** — OJ L, 2025/2083, 17.10.2025.
> http://data.europa.eu/eli/reg/2025/2083/oj — CELEX `32025R2083`. In force **20 October 2025**.

#### 5.1.1 The 50-tonne de minimis — new Article 2a (not Article 2(3))

> **‘Article 2a — De minimis exemption**
>
> 1. An importer, including any importer with the status of an authorised CBAM declarant, shall
> be exempted from the obligations under this Regulation, where **the net mass of the imported
> goods in a given calendar year does not cumulatively exceed the single mass-based threshold
> laid down in point 1 of Annex VII** (the “single mass-based threshold”). That threshold shall
> apply to the total net mass of goods under all CN codes aggregated per importer and per
> calendar year. In such a case, the importer, including an importer with the status of an
> authorised CBAM declarant, shall declare that exemption in the relevant customs declaration.
>
> 2. Where, within the relevant calendar year, an importer … **exceeds** the single mass-based
> threshold, the importer or the authorised CBAM declarant **shall be subject to all obligations
> under this Regulation in respect of all emissions embedded in all goods imported in that
> calendar year.**
>
> 3. By 30 April of each calendar year, the Commission shall assess, on the basis of the import
> data for the preceding 12 calendar months, whether the single mass-based threshold ensures that
> paragraph 1 of this Article applies to no more than 1 % of the emissions embedded in the
> imported goods and processed products. The Commission shall adopt delegated acts … to amend the
> single mass-based threshold … **where the value of the resulting threshold deviates from the
> applicable threshold by more than 15 tonnes**. The amended single mass-based threshold shall
> apply from 1 January of the following calendar year.
>
> **4. This Article shall not apply to imports of electricity or hydrogen.’**

The 50 t figure is in the new Annex VII:

> **‘ANNEX VII — The single mass-based threshold**
> 1. The single mass-based threshold referred to in Article 2a shall be set at **50 tonnes of net
> mass.**

Annex VII point 2 sets out the recalculation methodology (choose the threshold so ≥99 % of
embedded emissions stay in scope). Its footnote is directly relevant:

> The emission intensities E_j are based on default values (without mark-up) for emissions
> published for the transitional period. For cement and fertiliser products, direct emissions and
> indirect emissions are considered; **for aluminium and iron and steel products, only direct
> emissions are considered.**

**Hydrogen and electricity are carved out of the de minimis entirely** (Article 2a(4)). Recital
(4) of 2025/2083 explains:

> In the electricity and hydrogen sectors, key features such as quantity of imports, trade
> patterns, customs information and emission intensities differ substantially from those in the
> iron and steel, aluminium, fertilisers and cement sectors. Those differences imply that making
> electricity and hydrogen imports subject to a single mass-based threshold would require
> introducing complex adjustments that would not allow for the substantial reduction of
> administrative costs for importers in those sectors. **Imports of electricity or hydrogen should
> therefore not be included under the de minimis exemption.**

The threshold applies **cumulatively across the four goods sectors** (iron and steel, aluminium,
fertilisers, cement) — one aggregate figure per importer per calendar year, not per CN code and
not per consignment.

**Monitoring and enforcement** — new Article 25a. Highlights:

> 1. The Commission shall monitor the imports of goods for the purpose of monitoring the
> compliance with the single mass-based threshold. … The Commission shall periodically and
> automatically exchange with competent authorities the information necessary for the monitoring
> of importers via the CBAM registry. **Such information shall include a list of importers that
> exceed 90 % of the single mass-based threshold.**
>
> 3. Where the competent authority concludes that an importer that is not an authorised CBAM
> declarant has exceeded the single mass-based threshold, it shall without undue delay adopt a
> decision to that effect. … **The submission of an appeal against a decision determining that
> the importer has exceeded the single mass-based threshold shall not have suspensive effect.**
>
> 4. For the purpose of determining whether an importer has exceeded the single mass-based
> threshold, a competent authority shall **disregard a practice, arrangement or a series thereof
> which has been put into place for the main purpose or one of the main purposes of falling below
> the single mass-based threshold and which is non-genuine.**

Article 27(2)(b) was correspondingly rewritten to read:

> (b) **artificially splitting imports, including via non-genuine arrangements, to avoid
> exceeding the single mass-based threshold.**

On exceedance: retroactive full-year liability (Art. 2a(2)); customs must stop further imports
until authorisation is obtained; a new penalty in Article 26(2a), reducible where the threshold
was exceeded by no more than 10 %; a grace period for the quarterly holding obligation
(Art. 22(2a): "by the end of the quarter following that in which the … threshold is exceeded");
and a repurchase safety valve (Art. 23(2)) for anyone who bought certificates in anticipation
and then stayed under.

**Relevance to this project:** an HBI cargo is measured in tens of thousands of tonnes. The
50 t threshold is irrelevant to a bulk HBI trade — it exempts small parcel importers, not
commodity flows.

#### 5.1.2 Authorised declarant rules and dates

Article 5 rewritten:

> ‘1. Any importer established in a Member State shall, **prior to importing goods into the
> customs territory of the Union**, apply for the status of authorised CBAM declarant …
>
> 1a. An **indirect customs representative shall obtain the status of authorised CBAM declarant
> prior to importing goods** … irrespective of whether the importer is exempted from the
> obligations under this Regulation pursuant to Article 2a …
>
> 1b. Where Article 2a applies, the importer shall submit the application for an authorisation
> **in cases where that importer expects to exceed the single mass-based threshold**.’

> ‘7a. An authorised CBAM declarant **may delegate the submission of CBAM declarations** … to a
> person acting on behalf and in the name of that authorised CBAM declarant. The authorised CBAM
> declarant shall remain responsible for compliance …’

The transitional arrangement — new Article 17(7a). (Note: the phrase "occasional importer"
appears nowhere in 2025/2083; the actual mechanism is a one-off 2026 grace period.)

> ‘7a. By way of derogation from Article 4, where an importer or an indirect customs
> representative has submitted an application in accordance with Article 5 **by 31 March 2026**,
> such an importer or indirect customs representative **may provisionally continue to import goods
> until the competent authority takes a decision** under this Article.
>
> Where the competent authority refuses to grant the authorisation …, the competent authority
> shall establish, within one month of the date of the decision, the emissions embedded in the
> goods imported between 1 January 2026 and the date of that decision … **by reference to default
> values** … Those established emissions shall be used for the calculation of penalties in
> accordance with Article 26(2a).’

**Date shifts — both confirmed.**

Article 6(1) as replaced:
> ‘1. **By 30 September of each year, and for the first time in 2027 for the year 2026**, each
> authorised CBAM declarant shall use the CBAM registry … to submit a CBAM declaration for the
> preceding calendar year.’

(Original date was 31 May. So 31 May 2027 → **30 September 2027**.) Article 22(1) moves
surrender to the same date.

Article 20(1) as replaced:
> ‘1. **From 1 February 2027**, a Member State shall sell CBAM certificates on a common central
> platform to authorised CBAM declarants established in that Member State.’

Backed by Article 36(2), to which two points are added:
> ‘(c) Article 22(2) shall apply from 1 January 2027;
> (d) Article 20(1), (3), (4) and (5) shall apply from 1 February 2027.’

**So there is no certificate purchasing at all during 2026**; the whole 2026 liability is
settled in 2027. Recital (132):

> To provide authorised CBAM declarants sufficient time to prepare for compliance with the
> amended obligations …, Member States should start selling CBAM certificates in 2027 for
> emissions embedded in goods imported during the year 2026. **The price of CBAM certificates
> purchased in 2027 and corresponding to emissions embedded in goods imported into the Union in
> 2026 should reflect the prices of EU ETS allowances in 2026.**

Operationalised by new Article 21(1a):
> ‘1a. By way of derogation from paragraph 1, the Commission shall calculate the price of CBAM
> certificates that corresponds to the embedded emissions declared in respect of the year 2026 …
> as the **quarterly average** of the closing prices of EU ETS allowances on the auction platform
> … **of the quarter of importation of the goods** in which those emissions are embedded.’

Companion rules: certificates bought in 2027 for 2026 emissions may only be repurchased in 2027
(Art. 23(2a)), and are cancelled without compensation on 1 November 2027 (Art. 24(2)).

#### 5.1.3 The quarterly holding requirement: 80 % → 50 %

Article 22(2) as replaced:

> ‘2. **From 2027**, the authorised CBAM declarant shall ensure that the number of CBAM
> certificates on its account in the CBAM registry at the end of each quarter corresponds to **at
> least 50 %** of the embedded emissions in all goods it has imported since the beginning of the
> calendar year determined by reference to either of the following:
>
> (a) default values in accordance with the methods set out in Annex IV **without the mark-up** as
> referred to in point 4.1 of that Annex; or
>
> (b) the number of CBAM certificates surrendered in accordance with paragraph 1 for the calendar
> year preceding the year of the surrender, provided that the customs declaration for the import
> of goods refers to the same goods by CN code and countries of origin as the CBAM declaration
> submitted in the calendar year preceding the current year.
>
> **For the purpose of this paragraph, the adjustment for free allocation referred to in Article
> 31 shall be taken into account.**’

Four distinct reliefs are bundled here: the drop from 80 % to 50 %; the start pushed to 2027
(no quarterly obligation at all in 2026); the default-value basis now **without** the mark-up;
and option (b) letting a stable importer anchor on last year's actual surrender.

#### 5.1.4 Changes touching embedded emissions, precursors and default values

**Important correction to the Article 7 and Article 9 texts quoted earlier.** Both were amended
by 2025/2083; the versions in §1.1 and §8 below are updated accordingly.

Article 7(2) as replaced:
> ‘2. Embedded emissions in goods other than electricity shall be determined:
> (a) based on the actual emissions in accordance with the methods set out in points 2 and 3 of
> Annex IV; **or**
> (b) by reference to default values in accordance with the methods set out in point 4.1 of Annex
> IV.’

i.e. the declarant now has a genuine *choice*, rather than defaults being a fallback only where
actual emissions "cannot be adequately determined". This matters: it is a legal option, and a
green producer will always choose (a).

Article 8(1) as replaced — **verification is required only where actual emissions are used**:
> ‘1. **Where the embedded emissions are determined on the basis of actual emissions**, the
> authorised CBAM declarant shall ensure that the total embedded emissions declared in the CBAM
> declaration … are verified by a verifier accredited pursuant to Article 18 …’

Article 7(7)(a) now requires the implementing acts to determine
> **“system boundaries of production processes, which shall be aligned with those covered by the
> EU ETS”**, and relevant input materials (precursors) …

with recital (16) explaining the practical effect — finishing steps drop out of the boundary:

> The embedded emissions of some aluminium and steel goods currently included in the scope of
> Regulation (EU) 2023/956 are primarily determined by the embedded emissions of input materials
> (precursors), while the emissions arising during the production steps of those goods are
> typically relatively low. Those production steps consist of finishing processes that are carried
> out by separate installations not covered by the EU emissions trading system … **the embedded
> emissions of those production processes should be excluded from the system boundaries** for the
> calculation of emissions, by aligning the system boundaries of production processes with those
> covered by the EU ETS.

Recital (18) + the replaced Annex IV point 3 — **EU/exempted-origin precursors count as zero**:

> Only input materials (precursors) listed in Annex I and originating in third countries and
> territories that are not exempted pursuant to point 1 of Annex III are to be considered.

Annex IV point 4.1 as replaced (the fallback rule already quoted at §4.3), and Annex IV point 7
second paragraph, which makes region-specific adaptations **downward only**:

> Where declarants for goods produced in a third country, a group of third countries or a region
> within a third country can demonstrate, on the basis of reliable data, that alternative
> region-specific adaptations of default values are **lower** than the default values determined
> by the Commission, such region-specific adaptations can be used.

Other scope changes:
- New **Article 2(3a)** carves out offshore electricity and hydrogen: *"This Regulation shall not
  apply to: (a) electricity generated on the continental shelf or in the exclusive economic zone
  of a Member State or of a country or territory listed in points 1 and 2 of Annex III;
  (b) hydrogen originating on the continental shelf or in the exclusive economic zone of a Member
  State or of a country or territory listed in point 1 of Annex III."*
- Annex I: `2507 00 80` narrowed to `ex 2507 00 80 – Other kaolinic clays except non-calcined
  kaolinic clays` — the only Annex I goods change.
- Annex II: electricity added (see §1.3).
- New **Article 10a**: registration of accredited verifiers in the CBAM registry, requests
  *"within two months of the date on which the accreditation was granted, but not before
  1 September 2026"*.
- Article 3(31) "operator" redefined to include *"a parent company that controls an installation
  in a third country"*.
- Article 30(6) review must now also cover *"the application of the single mass-based threshold,
  including the possibility of increasing that threshold and of introducing a supplementary
  consignment-based threshold."*

### 5.2 The 17 December 2025 package — three documents, not one

| Doc | Date | Title | CELEX |
|---|---|---|---|
| **COM(2025) 783 final** | 16.12.2025 | Report from the Commission … on the application of the CBAM Regulation — **the Article 30(2) review** | `52025DC0783` |
| **COM(2025) 989 final** / 2025/0419(COD) | 17.12.2025 | Proposal … amending Regulation (EU) 2023/956 **as regards the extension of its scope to downstream goods and anti-circumvention measures** | `52025PC0989` |
| **COM(2025) 990 final** / 2025/0418(COD) | 17.12.2025 | Proposal … **establishing the Temporary Decarbonisation Fund** | `52025PC0990` |
| SWD(2025) 988 / 989 / 987 | 17.12.2025 | Impact assessment (part 1, 65 pp.), its executive summary, subsidiarity grid | `52025SC0988` etc. |

**Naming correction.** Several law-firm notes call the proposal "COM 2025/0419". That is the
interinstitutional **procedure** number 2025/0419(COD). The document number is **COM(2025) 989
final**.

#### 5.2.1 Downstream extension — what and when

**Date: 1 January 2028.** Article 2 of the proposal:

> This Regulation shall enter into force on the third day following that of its publication …
> Points 1 and 6 of Annex II, shall apply from 1 January 2026.
> However, Article 1(6), point (a), Article 1(8), points (a), (b) and (c), **Article 1(21), (23),
> and (24), and point 2 of Annex II shall apply from 1 January 2028.**

Article 1(21) is the Annex I goods-list amendment. Explanatory memorandum:

> … The changes requiring implementation in the CBAM registry or a launch at the start of the
> calendar year, **including the extension of scope to downstream products, will apply on
> 1 January 2028.**

**How many products.** Counted directly from the proposal's annexes:

- The **Iron and steel** table grows from **15 to 22** entries. New: `7312 10` (stranded wire,
  ropes, cables), `7314 39 00` (welded grill/netting/fencing), `7320 20 89` and `7320 90 90`
  (springs), `7323 94 00` and `7323 99 00` (table/kitchen/household articles), `7325` (other cast
  articles of iron or steel).
- A new table, **"Combined metal products"**, adds **107 CN-code entries** — machinery (84xx),
  electrical equipment (85xx), motor vehicles and parts (8704, 8706, 8707, 8708, 8716), plus
  `7314 31/41/49`, `7317 00`, `ex 7415 10 00`, `ex 8302 42/49`, `ex 8309 90 90`, `9018 32 10`,
  `9401 79 00`, `9403 10`. 102 of these carry "Carbon dioxide and perfluorocarbons" (aluminium-
  bearing); 5 carry "Carbon dioxide" only.
- The Aluminium table in Annex I is **not** amended.

Total: **114 new CN-code entries.** **The widely reported "around 180 products" figure appears in
no annex count and is stated nowhere in COM(2025) 989 or SWD(2025) 989.** It is presumably an
expansion of headings such as `7325`, `8704 21`, `8708 40` down to 8-digit level. Do not treat
"180" as a legal figure.

Named examples matching press coverage: `8418 10` combined refrigerator-freezers, `8450 11/12/19`
washing machines, `8451 21 00` drying machines, `8504 31 80`/`8504 33 00` transformers,
`8428 70 00` industrial robots, `8432 80 00`/`8432 90 00` agricultural machinery, `9403 10` metal
office furniture.

Selection method (explanatory memorandum):
> First, the trade intensity of goods was taken as a proxy for their tradability. … Second, a cost
> push indicator captures how much the carbon cost of CBAM inputs drives a downstream good's
> overall costs compared to its overall value added. In addition, to ensure that only products
> with the highest climate relevance are included, goods below a specified floor of total embedded
> emissions at sectoral level were excluded from the selection.

Scale: *"the estimated reduction in yearly GHG emissions is approximately 0,7 Mt of CO2 equivalent
emissions (CO2e) by 2030"*; *"around EUR 0.58 billion of annual revenues by 2030 … reaching an
estimated EUR 0.69 billion by 2035"*; *"around 3,800 – 3,900 additional SMEs facing CBAM
obligations"*.

Two structural novelties for downstream goods:
- **A mark-up waiver.** New subparagraph in Article 7(7): *"The implementing acts referred to in
  the first subparagraph may provide a list of downstream goods for which, due to the complexity
  of the supply chain and without prejudice to the environmental integrity of the CBAM, **no
  mark-up is to apply**."*
- **Composition-based precursor mass.** Annex IV amendment: *"However, for goods listed in
  sections 'Iron and Steel', 'Aluminium' and 'Combined Metal Goods' of Annex I, **M_i is a
  function of the content of goods used as input materials (precursors) in the manufacturing of
  the good**."*

#### 5.2.2 Indirect emissions — CONFIRMED: deliberately kept OUT for metals

**The December 2025 proposal does not extend CBAM to indirect emissions for iron/steel,
aluminium or hydrogen.** Three independent confirmations:

**(a) The proposal text.** "Indirect emissions" appears in COM(2025) 989 only in a narrow
methodological clarification for *electricity* default values (recital 47):

> To ensure a consistent methodological approach with respect to the default values applied for
> indirect emissions, it should be clarified that the alternative default value for indirect
> emissions that a third country, or a group of third countries, may demonstrate to be lower than
> the one established by the Commission, should be based on the same calculation method as the
> default values for indirect emissions determined by the Commission.

**There is no amendment to Annex II** — the annex that decides which sectors' indirect emissions
count.

**(b) SWD(2025) 988** — the impact assessment — does **not** assess an indirect-emissions
extension at all. Its three problem strands are downstream leakage, avoidance/circumvention, and
electricity import rules. Secondary sources pointing to SWD(2025) 988 as the location of the
Commission's indirect-emissions reasoning are wrong.

**(c) The reasoning is in COM(2025) 783 final (16.12.2025), Chapter 5.1** — the Article 30(2)
review report. This is the authoritative statement:

> Carbon leakage protection measures in the EU ETS cover both direct emissions and, more
> selectively, indirect emissions. For indirect emissions, the main instrument is **indirect cost
> compensation (ICC)**, where Member States can compensate a share of electricity-related carbon
> costs for electro-intensive industries listed in the State aid guidelines. …
>
> Some Member States have opted to grant ICC, while others have not. ICC covers up to a maximum
> allowed aid intensity of 75 % of eligible costs. … **Cement, fertiliser and agglomerated iron
> ore imports were included in the CBAM scope in 2023 for both direct and indirect emissions, as
> EU producers in these sectors were not eligible for ICC. By contrast, aluminium, steel and
> hydrogen are exempt from CBAM charges for indirect emissions, since EU producers in these
> sectors were eligible for ICC. This approach was designed to avoid double carbon leakage
> protection, meaning a situation where EU producers would receive compensation for their indirect
> carbon costs through ICC while importers were simultaneously charged under the CBAM for indirect
> emissions in those sectors.**

**This is the answer to "why".** The Annex I / Annex II split is not a judgement about whether
electricity emissions matter. It is a mirror of the EU's own domestic instrument set: sectors
whose EU producers get state-aid indirect cost compensation are excluded from CBAM's indirect
scope, to avoid double protection. It follows that the indirect-emissions question is legally
tied to the future of ICC and of the EU ETS state aid guidelines, not to CBAM design in
isolation.

The report canvasses five technical solutions for a future extension (immediate full coverage
with ICC kept or removed; covering only the ICC-uncompensated share; gradual CBAM-in/ICC-out;
ICC phase-out then CBAM; ICC retained with adjusted amounts), and concludes:

> **In conclusion, the analysis conducted so far confirms that technical solutions for expanding
> the indirect emissions scope of the CBAM can be conceived, beyond the approach enshrined in the
> current text of the CBAM Regulation, according to which ICC and CBAM are mutually exclusive.
> However, further analysis is needed, on one hand, to fully explore the viability in practice of
> several of these solutions and, on the other hand, to perform a more developed assessment of the
> impacts that such solutions would have.**

The explicit deferral, stated twice:

> **Step 2 provides for a report to follow (in 2027) with an evaluation of ways to extend the
> scope further (i) to indirect emissions from further CBAM goods (iron and steel, aluminium and
> hydrogen) and (ii) to other sectors.**

> Step 2: in 2027 the Commission will consider the possibility to provide for further extensions
> to other downstream products and other ETS sectors such as chemicals, indirect emissions.

The Commission's standing Q&A (updated 27 May 2026) states the present position flatly:

> Indirect emissions are taken into account only for CBAM goods for which indirect emissions fall
> within the scope of the CBAM, namely **cement and fertilisers (and agglomerated iron ore)**.

Also noted in the review: *"Indirect emissions typically do not exceed 10 % of direct emissions"*
for the sectors where they *are* counted (cement, fertilisers) — which is why the Commission
regards the current asymmetry as low-stakes for those sectors, though obviously not for
electricity-intensive green iron.

And, on why electricity downstream products are not addressed (COM(2025) 989 footnote 9):

> Downstream products of electricity are not considered given that electricity is used in the
> production process of virtually all goods, thus rendering the determination of the input share
> and embedded emissions of electricity in all possible imported goods unfeasible.

**Bottom line for this model: through at least 2028, imported iron/steel, aluminium and hydrogen
carry a CBAM charge on direct emissions only. The earliest realistic legislative vehicle for a
change is a proposal following the 2027 report, which would then need to clear the ordinary
legislative procedure — so 2030 at the earliest in practice.**

#### 5.2.3 Anti-circumvention measures proposed

- **New definition, Article 3(35):** *"‘abusive practices’ are practices pursued by an actor for
  the purpose of gaining a benefit by unduly avoiding, wholly or partially, the CBAM financial
  liability and thereby undermining the effectiveness of the CBAM to address the risk of carbon
  leakage in the EU."*
- **New Article 6(7)** — a delegated power to attach extra evidence conditions to the use of
  actual emissions for high-risk good/origin combinations:
  > Where the Commission … finds that there is sufficient evidence pointing towards a high risk
  > of abusive practices **for a combination of goods and origins**, it may inform importers and
  > authorised CBAM declarants about these risks, … and **it is empowered to adopt delegated acts**
  > … by laying down the methods for the identification of the combination of goods and origins,
  > the information to be declared for the use of actual emissions for those combinations …
  > **The Commission shall adopt the delegated acts … within three months of finding that there is
  > sufficient evidence** …

  with a matching new declaration item: *"evidence demonstrating that the high risk of abusive
  practices has not materialised."*
- **Traceability — new Article 6(2)(e):** *"where applicable for the purpose of addressing the
  risk of misdeclaration resulting from the lack of supply chain traceability, **evidence that the
  goods imported during the preceding calendar year were produced at the declared installation and
  at the actual time of production** referred to in the CBAM declaration"*. Reviewers may demand
  that evidence.
- **New circumvention practice, Article 27(2)(c):** *"**artificially adjusting the supply chains to
  make the goods benefit from lower default values.**"*
- **Composition granularity:** implementing powers to *"further detail CN codes to better capture
  the specific composition of the different products falling within any given CN code under the
  CBAM scope"*.
- **Operator registration required** in order to use actual values.
- **New Article 27a — "Serious and unforeseen circumstances"**, an escape hatch allowing the
  Commission to remove a good from Annex I by delegated act where its inclusion *"causes severe
  harm to the Union internal market due to serious and unforeseen circumstances related to the
  impact on the prices of goods"*. DG TAXUD's Q&A of 8 January 2026 confirms removal could operate
  **retroactively**.

#### 5.2.4 Pre-consumer scrap — added as a precursor, not as a product

**Correction to the common framing.** Pre-consumer scrap is not added as a separate CBAM
*product*. It is added as an **input material (precursor)** in a brand-new **Annex VIII**, "List
of non-CBAM goods and greenhouse gases considered as input materials (precursors)":

> **Iron and steel** — **ex 7204** Ferrous waste and scrap; remelting scrap ingots and steel
> **except post-consumer scrap** — Carbon dioxide
> **Aluminium** — **ex 7602** Aluminium waste and scrap **except post-consumer scrap** — Carbon
> dioxide

7204 and 7602 remain **excluded** from Annex I, so importing scrap as such still triggers no CBAM
liability. The change bites only when scrap is an input to a CBAM good.

Recitals (19)–(20):

> (19) Emissions from the production of pre-consumer scrap in the Union are subject to a carbon
> price since, under the EU ETS, emissions are measured at installation level. Since pre-consumer
> aluminium and pre-consumer steel scrap under Regulation (EU) 2023/956 are assigned
> zero-emissions, imported goods using pre-consumer aluminium and pre-consumer steel scrap as
> input material are subject to a lower carbon price compared to goods produced in the Union, thus
> weakening the effectiveness of the CBAM …
>
> (20) … Since pre-consumer scrap is a co-product generated unintentionally in the production
> process of metal goods and immediately reusable in a production process, **it is not considered
> at risk of carbon leakage in its own right. Therefore, the emissions of pre-consumer aluminium
> scrap and pre-consumer steel scrap should only be taken into account when used as a precursor
> for goods listed in Annex I** … The Commission should ensure that the monitoring, reporting and
> verification of emissions embedded in pre-consumer scrap used as input material (precursor) is
> not circumvented, including by misreporting pre-consumer scrap as post-consumer scrap to lower
> the determination of embedded emissions.

**Post-consumer scrap keeps its zero rating**, deliberately:

> In particular, it was considered that the inclusion of post-consumer scrap as CBAM precursor, as
> proposed under option 2, **could disincentivise the circular economy and would not be consistent
> with several EU policies in this area.**

#### 5.2.5 Hydrogen and DRI/HBI — no change proposed

**Hydrogen.** COM(2025) 989 makes no hydrogen-specific amendment. The explanatory memorandum
defers:

> A potential extension to downstream products in other CBAM sectors, namely those related to
> cement, fertilisers and hydrogen, is discussed in the Commission's review report set out under
> Article 30(2) of the CBAM Regulation. **An extension to these goods will be considered in a
> future legislative revision.**

**DRI/HBI (CN 7203).** The replaced Iron and steel table in COM(2025) 989's Annex I still opens
with "72 – Iron and steel / Except: [ferro-alloys] / 7204 – Ferrous waste and scrap …". **7203 is
not in the exception list**, so sponge iron, DRI and HBI remain fully in scope exactly as today.
The only ferro-alloy edit is a cosmetic split of `7202 2` into `7202 21 00, 7202 29`.

One DRI-adjacent change: DRI/HBI used as a *precursor* would fall under the new Annex IV rule
that, for Iron and Steel / Aluminium / Combined Metal Goods, `M_i` is a function of the precursor
*content* of the good.

#### 5.2.6 Exports: the Temporary Decarbonisation Fund (COM(2025) 990) — and it is not a rebate

There **is** a separate act addressing the export side of the EU steel problem, but it is not an
export rebate and the word "export" does not appear anywhere in COM(2025) 990. It is framed
entirely as addressing the "remaining risk of carbon leakage" during the free-allocation
phase-out. Recital (4):

> Energy-intensive industries covered by Directive 2003/87/EC progressively internalise the cost
> of their greenhouse gas emissions. The reduced Union-wide emissions cap, combined with the
> gradual phase-out of free allocation …, requires cost-intensive and rapid adaptations …
> **That remaining risk of carbon leakage is not fully prevented by Regulation (EU) 2023/956** …
> and should therefore be addressed through additional measures supporting the transition and
> promoting the decarbonisation of industrial sectors.

Key parameters:

- **Duration:** *"The Fund shall provide financial support in the period 2028-2029 to address the
  remaining risk of carbon leakage associated with carbon intensive goods produced by eligible
  operators of installations in the period 2026-2027."* (Art. 2(2))
- **Financing — 25 % of CBAM revenues.** Art. 3(2): *"Those contributions shall correspond to
  **25 % of the revenues that each Member State has collected from the sale of CBAM certificates**
  … in relation to embedded emissions declared for 2026 and 2027."*
- **Support formula** (Art. 9): based on the amount of free allocation phased out, computed per
  Article 16(8) of Delegated Regulation (EU) 2019/331, adjusted to the production share of listed
  goods, multiplied by the 2026–2027 average EUA closing price.
- **Conditionality:** support is conditional on decarbonisation investments (Art. 7).
- **Eligibility test:** production of Annex-listed goods, or (Member State opt-in) goods with
  *"a low ratio of value to weight"* at heightened remaining carbon-leakage risk — **not** exports.
- **No-precedent clause**, recital (9): *"**The transitory character of the Fund precludes any
  interpretation that it may constitute a precedent, a model or a reference point for the EU ETS
  review.**"*
- **The Annex lists 142 CN codes** — 7 aluminium, 8 fertilisers, 127 iron and steel, beginning
  with `26011200` (agglomerated iron ore) and **including `72039000` (DRI/sponge iron)**.

Secondary sources describing the Fund as "primarily targeting export-related carbon leakage" (and
EUROFER's framing about the share of steel exports covered) are reading the economic intent, not
the legal text. The economics are real, but the Commission has not stated exports as the purpose.

**Relevance to this project:** the Fund subsidises *EU* installations, including EU producers of
DRI/sponge iron under 7203 90 00. It is a competitive factor for an importer of Australian or
Brazilian HBI in 2028–2029, not a benefit available to them.

---

## 6. THE REVIEW CLAUSE — Article 30 of Regulation (EU) 2023/956

Full text of the parts that matter:

> **Article 30 — Review and reporting by the Commission**
>
> 1. The Commission, in consultation with relevant stakeholders, shall collect the information
> necessary **with a view to extending the scope of this Regulation as indicated in and pursuant
> to paragraph 2, point (a)**, and to developing methods of calculating embedded emissions based
> on environmental footprint methods.
>
> 2. **Before the end of the transitional period referred to in Article 32**, the Commission shall
> present a report to the European Parliament and to the Council on the application of this
> Regulation.
>
> The report shall contain an assessment of:
>
> (a) **the possibility to extend the scope to:**
> **(i) embedded indirect emissions in the goods listed in Annex II;**
> (ii) embedded emissions in the transport of the goods listed in Annex I and transportation
> services;
> (iii) goods at risk of carbon leakage other than those listed in Annex I, and specifically
> organic chemicals and polymers;
> (iv) other input materials (precursors) for the goods listed in Annex I;
>
> (b) the criteria to be used to identify goods to be included in the list in Annex I … based on
> the sectors at risk of carbon leakage identified pursuant to Article 10b of Directive
> 2003/87/EC; that assessment shall be accompanied by **a timetable ending in 2030 for the
> gradual inclusion of the goods within the scope of this Regulation** …;
>
> (c) the technical requirements for calculating embedded emissions for other goods to be
> included in the list in Annex I;
>
> (d) the progress made in international discussions regarding climate action;
>
> (e) the governance system, including the administrative costs;
>
> (f) the impact of this Regulation on goods listed in Annex I imported from developing countries
> with special interest to the least developed countries as identified by the United Nations
> (LDCs) and on the effects of the technical assistance given;
>
> (g) **the methodology for the calculation of indirect emissions pursuant to Article 7(7) and
> point 4.3 of Annex IV.**
>
> 3. **At least one year before the end of the transitional period**, the Commission shall present
> a report to the European Parliament and to the Council that **identifies products further down
> the value chain of the goods listed in Annex I** that it recommends to be considered for
> inclusion within the scope of this Regulation. To that end, the Commission shall develop, in a
> timely manner, a methodology that should be based on relevance in terms of cumulated greenhouse
> gas emissions and risk of carbon leakage.
>
> 4. The reports referred to in paragraphs 2 and 3 shall, where appropriate, be accompanied by a
> legislative proposal by the end of the transitional period, including a detailed impact
> assessment, in particular with a view to extending the scope of this Regulation on the basis of
> the conclusions drawn in those reports.
>
> 5. Every two years from the end of the transitional period … the Commission shall assess the
> effectiveness of the CBAM in addressing the carbon leakage risk of goods produced in the Union
> for export to third countries which do not apply the EU ETS or a similar carbon pricing
> mechanism. …
>
> 6. The Commission shall monitor the functioning of the CBAM with a view to evaluating the
> impacts and possible adjustments in its application.
>
> **Before 1 January 2028, as well as every two years thereafter**, the Commission shall present a
> report to the European Parliament and to the Council on the application of this Regulation and
> functioning of the CBAM. The report shall contain at least the following:
>
> (a) an assessment of the impact of the CBAM on:
> (i) carbon leakage, including in relation to exports;
> (ii) the sectors covered;
> (iii) internal market, economic and territorial impact throughout the Union;
> (iv) inflation and the price of commodities;
> (v) the effect on industries using goods listed in Annex I;
> **(vi) international trade, including resource shuffling; and**
> (vii) LDCs;
>
> (b) an assessment of:
> (i) the governance system, including an assessment of the implementation and administration of
> the authorisation of CBAM declarants by Member States;
> (ii) the scope of this Regulation;
> **(iii) practices of circumvention;**
> (iv) the application of penalties in Member States;
>
> (c) results of investigations and penalties imposed;
>
> (d) **aggregated information on the emission intensity for each country of origin for the
> different goods listed in Annex I.**

### 6.1 The exact mandate on indirect emissions

**The article reference is Article 30(1) + Article 30(2)(a)(i)** of Regulation (EU) 2023/956.
The mandate is:

> the possibility to extend the scope to … **embedded indirect emissions in the goods listed in
> Annex II**

and the Commission "shall collect the information necessary with a view to extending the scope"
to exactly that.

Recital (148) of the Regulation reinforces it:

> The Commission should, as part of that reporting, collect the information necessary **with a
> view to the further extension of the scope of this Regulation to embedded indirect emissions in
> the goods listed in Annex II as soon as possible**, as well as to other goods and services that
> could be at risk of carbon leakage, such as downstream products, and to developing methods of
> calculating embedded emissions based on the environmental footprint methods …

### 6.2 On the "2027 report" claim

**Careful with dates.** Two distinct reporting duties exist and neither is a "2027 report" in
terms:

- **Article 30(2)**: due *"before the end of the transitional period"*, i.e. **before 31 December
  2025**. This is the report that had to assess extension to indirect emissions in Annex II
  goods. This obligation was discharged as part of the December 2025 package (see §5).
- **Article 30(6)**: due *"before 1 January 2028, as well as every two years thereafter"*. This
  is the recurring application report, and it is the one that must cover **resource shuffling**
  (para 6(a)(vi)) and **practices of circumvention** (6(b)(iii)) and publish per-country emission
  intensities (6(d)).

So the *report due in 2027* is the Article 30(6) report (deadline 1 January 2028), and its
mandate is the impact/functioning review — including resource shuffling — **not** specifically a
fresh decision on indirect emissions. The indirect-emissions extension mandate sits in Article
30(2)(a)(i), whose deadline has already passed.

There is also a separate, dated review commitment in the secondary legislation, which is the one
most likely to move numbers relevant to this project: **default values and mark-ups must be
revised by December 2027 at the latest** (Article 2 of IR 2025/2621, and Article 11(5) of IR
2025/2547), with the Commission stating in recital (37) of 2025/2621 that it *"should make all
necessary efforts … to ensure that a revision of the default values can already be carried out in
2026."*

---

## 7. PHASE-IN — the CBAM factor schedule

### 7.1 The legal source

The schedule is **not** in the CBAM Regulation. Article 31 of Regulation (EU) 2023/956 only
establishes the principle:

> **Article 31 — Free allocation of allowances under the EU ETS and obligation to surrender CBAM
> certificates**
>
> 1. The CBAM certificates to be surrendered in accordance with Article 22 of this Regulation
> shall be adjusted to reflect the extent to which EU ETS allowances are allocated free of charge
> in accordance with Article 10a of Directive 2003/87/EC to installations producing, within the
> Union, the goods listed in Annex I to this Regulation.
>
> 2. The Commission is empowered to adopt implementing acts laying down detailed rules for the
> calculation of the adjustment … Such detailed rules shall be elaborated by reference to the
> principles applied in the EU ETS for the free allocation of allowances … taking account of the
> different benchmarks used in the EU ETS for free allocation with a view to combining those
> benchmarks into corresponding values for the goods concerned, and **taking into account relevant
> input materials (precursors)**.

The numbers live in **Directive 2003/87/EC, Article 10a(1a)**, as inserted by Directive (EU)
2023/959 (the ETS revision). Consolidated text (CELEX `02003L0087-20240301`), Article 10a(1a),
second subparagraph:

> By way of derogation from the first subparagraph of this paragraph, for the first years of
> application of Regulation (EU) 2023/956, the production of goods listed in Annex I to that
> Regulation shall benefit from free allocation in reduced amounts. **A factor reducing the free
> allocation for the production of those goods shall be applied (CBAM factor). The CBAM factor
> shall be equal to 100 % for the period between the entry into force of that Regulation and the
> end of 2025 and, subject to the application of provisions referred to in Article 36(2), point
> (b), of that Regulation, shall be equal to 97,5 % in 2026, 95 % in 2027, 90 % in 2028, 77,5 % in
> 2029, 51,5 % in 2030, 39 % in 2031, 26,5 % in 2032 and 14 % in 2033. From 2034, no CBAM factor
> shall apply.**
>
> The reduction of free allocation shall be calculated annually as the average share of the demand
> for free allocation for the production of goods listed in Annex I to Regulation (EU) 2023/956
> compared to the calculated total free allocation demand for all installations, for the relevant
> period referred to in Article 11(1) of this Directive. The CBAM factor shall be applied in this
> calculation.

### 7.2 VERIFIED — but mind the terminology

The user's schedule ("2.5 % in 2026 rising to 48.5 % in 2030 and 100 % by 2034") is **numerically
correct** as the *CBAM phase-in*, i.e. the share of benchmark-level emissions that is **not**
shielded by free allocation. But it is the **complement** of the legally defined "CBAM factor",
which is the share of free allocation **retained**:

| Year | "CBAM factor" (Art. 10a(1a)) = free allocation retained | CBAM phase-in = 1 − CBAM factor |
|---|---|---|
| to end 2025 | 100 % | 0 % |
| **2026** | **97,5 %** | **2,5 %** |
| **2027** | **95 %** | **5 %** |
| **2028** | **90 %** | **10 %** |
| **2029** | **77,5 %** | **22,5 %** |
| **2030** | **51,5 %** | **48,5 %** |
| **2031** | **39 %** | **61 %** |
| **2032** | **26,5 %** | **73,5 %** |
| **2033** | **14 %** | **86 %** |
| **2034 onwards** | none (0 %) | **100 %** |

This distinction is load-bearing for anyone reading IR 2025/2620: in the free-allocation-adjustment
equations, `CBAM_y` is the *statutory* CBAM factor (97.5 %, 95 %, …), and it multiplies the
benchmark to give the **relief**, not the charge. Writing 2.5 % into that slot would invert the
result.

### 7.3 How it tracks the ETS free allocation phase-out for steel

The two are the same instrument seen from two sides:

1. EU steel installations lose free allocation on exactly this schedule — from 2026 an EU BF-BOF
   plant retains only 97.5 % of its benchmark-based free allowances, and by 2034 none.
2. Imports get a *matching* deduction — the free allocation adjustment — so that in any year the
   effective carbon charge on an importer and on an EU producer at the same emissions intensity is
   the same.

Concretely, for a good g in year y:

```
certificates due per tonne  ≈  SEE_g  −  CBAM_factor(y) × CSCF(y) × BM_g
```

For **HBI at CN 7203** with a default declaration in 2026:
```
1.325 × 1.10 (mark-up)  −  0.975 × CSCF × 0.397   ≈  1.458 − ~0.38  ≈  1.08 t CO2e/t
```
and in 2030:
```
1.325 × 1.30            −  0.515 × CSCF × 0.397   ≈  1.723 − ~0.20  ≈  1.52 t CO2e/t
```
(CSCF ≈ 1 for illustration; it is published annually and is typically slightly below 1.)

For **verified green HBI at, say, 0.05 t CO2e/t direct** the bracket goes negative, i.e. the free
allocation adjustment exceeds the embedded emissions, and **no certificates are due**. Note the
floor at zero is *implicit* rather than stated: Article 6(2)(c) frames the surrender obligation as
the total embedded emissions "after the reduction that is due on the account of the carbon price
paid in a third country in accordance with Article 9 and the adjustment necessary to reflect the
extent to which EU ETS allowances are allocated free of charge in accordance with Article 31", and
neither the Regulation nor IR 2025/2620 provides for a negative surrender or a transferable credit
— the only express zero-setting is Article 1(2) of IR 2025/2620 for electricity. Treat "adjustment
capped at the embedded emissions, no credit generated" as the correct modelling assumption, and
flag it as an inference rather than a quoted rule. This holds comfortably through 2030 on these
numbers; as the CBAM factor decays the shield shrinks, but a genuinely near-zero-direct-emission
HBI still owes essentially nothing even once the CBAM factor reaches zero in 2034.

The strategic point for the model: **CBAM does not reward green HBI relative to EU-produced green
iron; it removes the penalty relative to fossil imports.** The competitive delta a green Australian
or Brazilian HBI exporter captures against, say, a natural-gas DRI competitor is roughly
(1.325 − 0.05) × mark-up × EUA price × phase-in-independent-of-benchmark — because both receive the
same benchmark deduction. At an EUA price of €80/t and the 2026 mark-up, that is on the order of
€110/t HBI in 2026, rising as the mark-up climbs; against a BF-BOF-based competitor the gap is
larger still. Against **scrap-EAF** it is negative — scrap is not even a CBAM good (7204 is excluded
from Annex I's Chapter 72 sweep).

---

## 8. ARTICLE 9 — CARBON PRICE PAID ABROAD

### 8.1 The definition — Article 3(29) of Regulation (EU) 2023/956

> ‘**carbon price**’ means the monetary amount paid in a third country, under a carbon emissions
> reduction scheme, in the form of a tax, levy or fee or in the form of emission allowances under a
> greenhouse gas emissions trading system, calculated on greenhouse gases covered by such a
> measure, and released during the production of goods;

Unamended by Regulation (EU) 2025/2083, which touched only Article 3 points (15) and (31).

### 8.2 Article 9 as replaced by Regulation (EU) 2025/2083

The heading changed from "carbon price paid in the **country of origin**" to "carbon price paid in
**a third country**" — a deliberate widening, per recital (63) of 2025/2083: *"Since the carbon
price can be paid in a third country other than the country of origin of the imported goods, such a
carbon price should also be eligible for deduction."*

> **‘Article 9 — Carbon price paid in a third country**
>
> 1. **Where the embedded emissions are determined on the basis of actual emissions**, an
> authorised CBAM declarant may claim in the CBAM declaration a reduction in the number of CBAM
> certificates to be surrendered in order to take into account the carbon price paid in a third
> country for the declared embedded emissions. **The reduction may be claimed only if the carbon
> price has been effectively paid in a third country. In such a case, any rebate or other form of
> compensation available in that country that would have resulted in a reduction of that carbon
> price shall be taken into account.**
>
> 2. The authorised CBAM declarant shall keep records of the documentation required to demonstrate
> that the declared embedded emissions were subject to a carbon price in a third country that has
> been effectively paid as referred to in paragraph 1. The authorised CBAM declarant shall in
> particular keep evidence related to **any rebate or other form of compensation available**, in
> particular the references to the relevant legislation of that country. The information contained
> in that documentation shall be **certified by a person that is independent from the authorised
> CBAM declarant and from the authorities of the third country**. The name and contact information
> of that independent person shall appear on the documentation. The authorised CBAM declarant shall
> also keep evidence of the actual payment of the carbon price.
>
> 3. The authorised CBAM declarant shall keep the records referred to in paragraph 2 until the end
> of the fourth year after the year during which the CBAM declaration has been or should have been
> submitted.
>
> 4. By way of derogation from paragraphs 1, 2 and 3, an authorised CBAM declarant may claim, in
> the CBAM declaration, a reduction … **by reference to yearly default carbon prices**. In such a
> case, any rebate or other form of compensation available in that country that would have resulted
> in a reduction of that default carbon price shall be taken into account. **The reduction may be
> claimed only where a carbon price was set by the rules applicable in the third country and a
> yearly default carbon price can be determined, including on a conservative basis, for that third
> country. Where the embedded emissions are determined on the basis of default values, a reduction
> may be claimed only by a reference to yearly default carbon prices.**
>
> **As from 2027**, the Commission **may**, for third countries where carbon pricing rules are in
> place, determine and make available, in the CBAM registry referred to in Article 14, **the default
> carbon prices for those third countries** and publish the methodology for their calculation. …
>
> **5. The Commission is empowered to adopt implementing acts** concerning the conversion of the
> yearly average carbon price effectively paid in accordance with paragraph 1 … and of the yearly
> default carbon prices determined in accordance with paragraph 4 … into a corresponding reduction
> of the number of CBAM certificates to be surrendered. Those acts shall also govern the conversion
> … into euro at the yearly average exchange rate, **the evidence required of the actual payment of
> the carbon price, examples of any relevant rebate or other form of compensation** …, the
> qualifications of the independent person … and the conditions to ascertain that person’s
> independence. …’

**Note the renumbering.** The implementing-act empowerment is now **Article 9(5)**, not 9(4); 2025/2083
inserted the new default-carbon-price paragraph as 9(4) and pushed the old empowerment down. The
Commission's own Q&A (updated 27 May 2026) still cites "Article 9(4)" for the implementing act —
**the Q&A is stale on this point**; the draft act itself cites 9(5).

**Two structural changes matter for modelling.** First, **actual-value declarations are now a
precondition** for claiming a real, certified carbon price (Art. 9(1) opening words). A declarant on
default emission values is confined to Commission-published default carbon prices, which do not yet
exist and may never for a given country. Second, the separate treaty route survives:

> **Article 2(12)** — The Union may conclude agreements with third countries or territories with a
> view to taking into account carbon pricing mechanisms in such countries or territories for the
> purposes of the application of Article 9.

### 8.3 The Article 9(5) implementing act — DRAFT PUBLISHED, NOT YET ADOPTED

**Full title** (Ref. Ares(2026)4841230 – 13/05/2026):

> COMMISSION IMPLEMENTING REGULATION (EU) …/… of XXX **laying down rules for the application of
> Regulation (EU) 2023/956 as regards the conversion of the carbon price paid in a third country
> into a corresponding reduction in the number of CBAM certificates to be surrendered, the evidence
> of payment of that carbon price, the qualifications of the independent person and conditions to
> ascertain its independance and qualifications** (Text with EEA relevance)

*(the misspelling "independance" is in the original)*

**Status.** No CELEX or ELI exists — it is unadopted. Better Regulation initiative **14830**: call for
evidence 28 Aug – 25 Sep 2025; **draft published 13 May 2026**; feedback closed 10 June 2026 (177
submissions); current stage `ISC_WORKFLOW`, `adoptionDate: None`, `ADOPTION_WORKFLOW` marked
`UPCOMING`. It carries the standard banner *"This draft has not been adopted or endorsed by the
European Commission."* The Commission's CBAM legislation page (updated August 2026) lists no such
regulation among the adopted acts.

- Initiative page: https://ec.europa.eu/info/law/better-regulation/have-your-say/initiatives/14830
- Draft text: https://ec.europa.eu/info/law/better-regulation/api/download/090166e52d80a42e
- Annexes: https://ec.europa.eu/info/law/better-regulation/api/download/090166e52d80a42f
- News item: https://taxation-customs.ec.europa.eu/news/carbon-price-paid-third-countries-2026-05-13_en
- Call-for-evidence synopsis: https://taxation-customs.ec.europa.eu/document/download/c4049ce0-48c2-43ff-9414-1dd11b3941f4_en

**Article 32** of the draft: enters into force on the third day after OJ publication, and *"shall
apply from **1 January 2026**"* — retroactive to the first definitive-regime year, whose declaration
falls due 30 September 2027. So the timetable is workable even though the act is late.

**None of the December 2025 acts covers carbon price paid abroad.** Verified titles: 2025/2546
(verification principles), 2025/2547 (methodology), 2025/2548 (certificate price calculation),
2025/2619 (customs information), 2025/2620 (free allocation adjustment), 2025/2621 (default values),
DR 2025/2551 (accreditation). 2025/2620 is the Article 31 free-allocation adjustment — a different
mechanism from the Article 9 deduction and easily confused with it.

### 8.4 What the draft act says — the parts that decide the Australian question

**Definitions — Article 2:**

> (3) ‘**carbon price mechanism**’ means carbon tax, carbon fee or carbon levy, or an emissions
> trading system;
>
> (5) ‘**baseline-and-credit emission trading system**’ means the form of emission trading system
> where a baseline is established as an emission limit under which no carbon price is due and
> tradable emission credits are issued to entities that emit less than the baseline, and where
> emissions credits must be purchased by entities emitting more than the baseline;

**Recital (10):**
> For the purpose of ensuring equivalence with the carbon price paid under the EU ETS, **credits or
> other emission units purchased under a baseline-and-credit emissions trading system should be
> considered equivalent to allowances paid under an emissions trading system.**

**Recital (8) — the qualifying filter** (recital only; no corresponding operative article):
> …only carbon prices paid on specific embedded emissions under a carbon price mechanism in a third
> country should give rise to a reduction … where that scheme takes the form of a tax, levy or fee
> or of emission allowances under a greenhouse gas emissions trading system that is **binding in
> nature and imposes compliance obligations on all operators active in the relevant sectors covered
> by that mechanism without discrimination**.

**Recital (9)** extends eligibility to fuel-based carbon taxes not paid directly by the operator,
provided the rate is consistent with the fuel's emission factor.

> **Legal tension worth flagging.** Article 3(29) of the parent Regulation confines the carbon price
> to "a tax, levy or fee or … **emission allowances** under a greenhouse gas emissions trading
> system". ACCUs and SMCs are *credits*, not allowances. The draft bridges this by defining
> baseline-and-credit as a *form of* ETS (Art. 2(5)) and declaring its credits *equivalent to*
> allowances (recital 10). Defensible, but it is the implementing act doing interpretive work the
> parent act does not obviously authorise. If the act is challenged, this is where.

**How rebates and free allocation are netted off — Article 8(1):**

> When preparing the operator's carbon price report, the operator shall identify and take into
> account the following rebates or other forms of compensation:
> **(a)** a reduced tax rate under a carbon tax, levy or fee;
> **(b)** any exemption of the emissions coverage in the carbon price mechanism, including:
>   (i) emissions associated with **free allowances** received by the operator;
>   (ii) **emissions that are below an emission intensity baseline under a baseline-and-credit
>   emissions trading system**; and
>   (iii) emissions exempted from the application of the carbon tax, levy or fee;
> **(c)** a **refund in monetary value** that partially or totally compensates the carbon price paid,
> including forms of **indirect cost compensation** due to a carbon price mechanism;
> **(d)** any other rebate or form of compensation that is based on any relevant parameters
> establishing the effective carbon price to be paid on the emissions covered by a carbon price
> mechanism.

**Recital (14):** *"**Any modification of a parameter that lowers the obligation to pay the carbon
price should be regarded as a compensation** for that purpose."*

**Article 8(2) / recital (15) — the one carve-out.** Revenue recycling is *not* a rebate, provided
cumulatively: all covered installations are eligible irrespective of the price each paid; all must
apply; the granting decision is public; and the subsidy's stated objective is reducing the
beneficiary installation's emissions.

**Article 8(3) — a hard trap:**
> Where an operator … is entitled to a rebate or other form of compensation but the operator is not
> able to provide evidence of the amount to be deducted, the operator shall **not be entitled to any
> deduction** from the carbon price effectively paid.

(Softened in the following subparagraphs: an officially established or ascertainable maximum may be
used where the rebate has not yet been received; a rebate the operator can prove it never requested,
or was refused, is disregarded.)

**The arithmetic — Annex I.** Section 5 defines the rebated quantum to include *"emissions
associated with free allowances … for which no carbon price has been paid"* and *"**emissions that
are below an emission intensity baseline that are exempted from payment of a carbon price under a
baseline-and-credit emission trading system**"*. Section 6.1, equations 12a/12b:

> **EFF_CP_DIR_g = (EM_DIR_g − Rebated_EM_DIR_g) × EFF_CP_DIR**
> **EFF_CP_IND_g = (EM_IND_g − Rebated_EM_IND_g) × EFF_CP_IND**

and, for monetary refunds, §3.4 equations 6a/6b: `EFF_CP_DIR = CP_DIR − RC_DIR`. **Rebates bite
twice — on the tonnage base and on the price rate.**

**Conversion to certificates — Article 6(1):**
> Reduction_ActualCarbonPrice,g = (€EFF_CP_g / Ref Price CBAM) × Q_g

and Article 6(2) for the default route, `(€DCP × SEE_g / Ref Price CBAM) × Q_g`, where `Ref Price
CBAM` is the average published CBAM certificate price for the year of import under IR 2025/2548.

**Currency — Article 5:** yearly average exchange rate for the reporting year, published by the
Commission using ECB rates or, where appropriate, Eurostat.

**Evidence and the independent person.** Article 7: an electronic **carbon price report**, in
English, via the CBAM registry. Recital (7): emissions subject to the carbon price must be
attributed to goods per the system boundaries and production processes of **IR 2025/2547**, with a
**5 % tolerance** between the mechanism's emission boundary and CBAM's. Recital (23) / Annex I
§3.3.1: for an ETS, normally the **weighted average auctioning price**; failing that the yearly
average secondary-market price published by the market's managing authority; alternatively records
of individual purchases matched to allowances surrendered. Article 9(1): the independent person must
be **accredited** for the CBAM activity group in Annex III and comply with **EN ISO/IEC 17029:2019**.
Article 9(2): independence from the operator, the Member State competent authority, the Commission,
**and any third-country authority regulating or supervising the carbon price mechanism**. Article 10:
**reasonable assurance**. Article 13: reliance on the emissions verification report is permitted only
where accreditation was valid, the opinion satisfactory and the periods match; for bought-in
precursors, on the precursor installation's certification report on equivalent conditions — otherwise
the operator falls back to the default carbon price.

**Paris Agreement Article 6 credits.** Recital (12): domestic carbon credits used for compliance are
recognised without extra criteria. Recital (13) and Annex I §3.3.4/§3.5.4: international credits must
be Article 6.2 ITMOs registered on the UNFCCC CARP, or Article 6.4 units, and are capped:

> …international carbon credits pursuant to Article 6(2) and 6(4) of the Paris Agreement … may only
> be claimed to a maximum of **10 % of the reported and confirmed CPM emissions** covered by the
> third-country carbon price mechanism … Where more than 10 % … are covered by such international
> carbon credits, **a price of zero shall be assigned** to CPM emissions covered by international
> carbon credits in excess of this 10 % threshold.

**The Commission Q&A (27 May 2026), Q1.5**, is narrower and reads as hostile to intensity schemes:
> The deduction concerns compliance schemes, **essentially carbon taxes explicitly levied on the
> embedded emissions and allowances under an emission trading system**. The carbon price should be
> effectively paid, in the sense that it should not be rebated or compensated. This would for example
> be the case if free allowances under an ETS are granted, or there are exemptions under a carbon tax.

The Q&A never mentions baseline-and-credit. The draft implementing regulation, published two weeks
earlier, is later, more specific and unambiguous — **where they diverge, the draft governs, but it is
a draft**, and the Q&A's "Article 9(4)" citation is out of date.

### 8.5 Australia's Safeguard Mechanism — the instrument qualifies; almost none of the tonnes do

**Eligible in principle.** Article 2(5) of the draft is a near-verbatim description of the reformed
Safeguard Mechanism, and recital (10) makes its credits equivalent to allowances. The drafting looks
deliberate, with schemes like Australia's and Canada's OBPS in view.

**But the baseline is a rebate.** Article 8(1)(b)(ii) and Annex I §5(b) classify "emissions that are
below an emission intensity baseline under a baseline-and-credit emissions trading system" as a
rebate. Feed that into equation 12a: for a facility at or below baseline, `Rebated_EM = EM`, so
**EFF_CP = 0**. A facility above baseline gets a deduction only on the excess slice, priced at the
ACCU/SMC price actually paid.

**The scheme.** NGER Act 2007; NGER (Safeguard Mechanism) Rule 2015; *Safeguard Mechanism (Crediting)
Amendment Act 2023* (No. 14 of 2023, assented 11 April 2023). Facilities above 100,000 t CO2-e scope 1;
production-adjusted baselines declining **4.9 %/yr to 30 June 2030**; excess managed by 1 April with
ACCUs, SMCs or a flexibility measure. **Emissions up to the baseline attract no compliance cost at
all.** SMCs are issued to facilities *below* baseline — the Clean Energy Regulator: *"SMCs are an
accounting tool within the Safeguard Mechanism's regulated emissions limit. **They are not
offsets.**"* and, because baselines are intensity-based, *"SMCs may even be earned where a facility's
total emissions increase."* **TEBA** determinations cut the decline rate to 1–2 % for three years
where compliance cost exceeds 3 % of revenue/EBIT; the CER states their purpose is to *"mitigate risk
of cross-border carbon leakage by supporting the competitiveness of trade-exposed businesses"* — which
under Article 8(1)(d) is itself a compensating parameter, cutting against Australia a second time.

**The Australian Government agrees with this reading.** DCCEEW's submission to the EU consultation,
Ref. Ares(2026)5894923 – 10/06/2026:

> **Australia supports Article 2(3) …; and Article 2(5) …. These definitions clarify how the Safeguard
> Mechanism would be recognised by the EU CBAM.** Australia also supports **Article 8(1)(b)(ii) and
> Recital 10, which clarifies how the Implementing Regulation would treat Safeguard Mechanism
> baselines.**

On recital (8): *"Our understanding is that this would apply to the Safeguard Mechanism, because 'all
operators… without discrimination' means equal treatment **within the covered scope** … Australia
encourages the European Commission to clarify this in the Implementing Regulation or in explanatory
material."* Since neither ACCUs nor SMCs are auctioned, Australia supports the secondary-market route
and offers its published **Default Prescribed Unit Price**.

**Australia's actual ask is a timing fix, not a recognition fix.** The Safeguard runs July–June, NGER
reports fall due 31 October, and final FY2026-27 compliance is 31 March 2028 — after the CBAM
declaration deadline:

> …this could result in a deduction for the carbon price paid not applying for the second half of each
> calendar year for affected goods – **potentially around half of the annual emissions for those
> goods.** … This would be inconsistent with the objective of Article 9 …, which is to ensure that a
> carbon price is not paid twice on the same emissions.

**Australia's own Carbon Leakage Review defines the effective price identically.** Prof. Frank Jotzo's
final report (Feb 2025, released Feb 2026) does not address EU recognition at all — its three
recommendations concern an *Australian* BCA — but footnote 114 is directly on point:

> 'Effective price paid' is defined for jurisdictions with an ETS or carbon tax as the explicit carbon
> price … **after taking account of free allocations**. **In the context of the Safeguard Mechanism,
> the equivalent 'effective price paid' is defined as the carbon compliance costs (per tonne) which
> arise on emissions above Safeguard regulated baselines, i.e. from ACCU or SMC surrender.**

and at p. 24: *"Facilities can emit up to their baseline level. **This is equivalent to the role that
free allowances play under emissions trading schemes.**"*

**The UK has already made the call, the same way.** HMRC's list of qualifying carbon pricing schemes
(published 27 August 2026) names the **Australia Safeguard Mechanism** among 16 recognised schemes —
but with the same netting: *"**Emissions covered by free allowances do not qualify for relief because
no effective carbon price has been paid on them.** … **In some cases, no carbon price relief will be
available.**"* SI 2026/809, reg. 7(2)(a)(iii) lists *"the threshold above which a carbon price is
charged"* among the elements netted off — the direct UK analogue of Article 8(1)(b)(ii). **So the UK
listing is not evidence that Australian exporters get meaningful relief, and should not be cited that
way.**

> **Secondary-source conflict.** Clayton Utz (August 2026) writes that whether Safeguard costs "will
> be creditable against CBAM liabilities remains to be determined" and that "its recognition will
> depend on how the European Commission classifies it." That was overtaken three months earlier by the
> draft's Articles 2(5) and 8(1)(b)(ii). The residual uncertainty is narrower: adoption of the act, the
> recital (8) clarification, and the timing mismatch.

**Upshot for the model.** Do **not** credit Australian production with an Article 9 deduction of ACCU
price × embedded emissions. For a baseline-compliant facility — and a green HBI plant would be far
below baseline, *earning* SMCs — the deduction is approximately **zero**. The real Australian carbon
cost in the model is the Safeguard's own above-baseline liability, which for a green plant is also
approximately zero (and is arguably a small revenue from SMC sales).

### 8.6 Brazil — SBCE (Lei 15.042/2024): no price is paid, and none will be for years

**The law.** *Lei nº 15.042, de 11 de dezembro de 2024*, DOU 12.12.2024, in force on publication
(Art. 58), sanctioned without vetoes. It creates the **CBE** (*Cota Brasileira de Emissões*, granted
*"de forma gratuita ou onerosa"*, Art. 2 VI) and the **CRVE** (verified reduction certificate,
Art. 2 III).

**Phasing — Article 50 verbatim:**
> I – **fase I**: período de 12 (doze) meses, **prorrogável por mais 12 (doze) meses**, para a edição
> da regulamentação desta Lei …;
> II – **fase II**: período de 1 (um) ano para operacionalização … dos instrumentos para relato de
> emissões;
> III – **fase III**: período de 2 (dois) anos, no qual os operadores estarão sujeitos **somente** ao
> dever de submissão de plano de monitoramento e de apresentação de relato …;
> IV – **fase IV**: vigência do primeiro Plano Nacional de Alocação, com **distribuição não onerosa de
> CBEs** …;
> V – **fase V**: implementação plena do SBCE …

Two provisions settle when anyone pays. **Art. 50(IV):** the first National Allocation Plan
distributes CBEs **free of charge only**. **Art. 11 §3º:** *"O início da cobrança pela outorga onerosa
das CBEs seguirá as fases de implementação do SBCE"* — paid allocation begins only at Phase V.

**Status, September 2026.** The permanent *Órgão Gestor* **does not exist**. Decree 12.677 of
15 October 2025 created a *Secretaria Extraordinária do Mercado de Carbono* (SEMC) inside the Finance
Ministry, exercising the managing body's functions *"de forma temporária até que se crie e entre em
funcionamento o órgão gestor do SBCE"*. Critically, SEMC was **not** given the Art. 8 competences
VIII–IX (allocation plan), **XI (issue CBEs)**, **XII (auctions)** or XIII–XVI (approve monitoring
plans, receive reports, run reconciliation, price stabilisation). Every function capable of generating
a payable carbon price is reserved to a body that does not yet exist. ICAP records the ETS status as
**"under development"**.

**Scope covers steel, from ~2031.** SEMC's sectoral proposal (19 May 2026) puts **iron and steel in
the first wave**, submitting monitoring plans in **2027**; SEMC's 2026 calendar carries the line item
*"Workshop MRV **ferro e aço**"*. On the draft MRV ordinance's four-year run-in (consulted 28 July –
28 August 2026), **surrender obligations start in the fifth year — 2031 for the first wave.**
Thresholds (Art. 30 read with Art. 29): >10,000 tCO2e/yr = monitoring plan and reporting only;
**>25,000 tCO2e/yr** adds the *conciliação periódica de obrigações*, i.e. surrender. Primary
agriculture is excluded outright (Art. 1 §2º).

**No other Brazilian carbon price qualifies.** There is no federal carbon tax, and state-level pricing
on SBCE sectors is legally foreclosed — Art. 22: *"vedadas a dupla regulação institucional e **qualquer
tributação sobre emissões de GEE** por atividades, por instalações ou por fontes reguladas pelo
SBCE."* Fuel excises (CIDE, ICMS) fall outside both Art. 3(29) of the CBAM Regulation and Art. 2(3) of
the draft implementing act.

**No EU–Brazil CBAM recognition instrument exists.** The EU–Mercosur interim Trade Agreement (signed
17 January 2026, provisionally applied 1 May 2026) contains no CBAM chapter, annex or joint
declaration.

**Upshot for the model: zero Article 9 deduction for Brazilian HBI throughout any horizon ending
before ~2031**, and even then only for above-free-allocation emissions.

*Caveat: no decree formally invoking the Art. 50(I) 12-month extension was located. The circumstantial
evidence (SEMC's 2026 calendar, the July 2026 consultation) is conclusive that Phase I is still
running, which is only lawful if the extension was taken — treat as a strong inference, not a verified
fact.*

### 8.7 State of play on recognition generally

**There is no Commission "list of recognised carbon prices", and no legal concept of "recognition" in
the Regulation.** Two distinct routes exist and are routinely conflated:

**(a) Full exemption — Article 2(4) + Annex III point 1.** Entry conditions, Article 2(6):
> (a) **the EU ETS applies to that third country or territory or an agreement has been concluded …
> fully linking the EU ETS and the emission trading system of that third country or territory**;
> (b) the carbon price paid in the country in which the goods originate is effectively charged on the
> greenhouse gas emissions embedded in those goods **without any rebates beyond those also applied in
> accordance with the EU ETS**.

Annex III point 1, unchanged since adoption: **Iceland, Liechtenstein, Norway, Switzerland**, plus
Büsingen, Heligoland, Livigno, Ceuta, Melilla. Iceland/Liechtenstein/Norway qualify because the EU ETS
applies via the EEA; **Switzerland is the only country ever exempted via the linking limb** (agreement
signed 7 December 2017, in force 1 January 2020). Annex III point 2 (electricity) remains an empty
placeholder. Per Q&A 7.19, Norway and Iceland "aim to" incorporate CBAM into the EEA Agreement **as of
1 January 2027**, after which they apply CBAM at their own external borders.

**(b) Partial deduction — Article 9.** Declarant-driven, per-consignment, independently certified.
**No list, no country approval.** The only forward-looking list-like object is the set of **default
carbon prices** the Commission "**may**" publish in the CBAM registry **as from 2027** (Art. 9(4)).
None has been published; the enabling act is the unadopted draft above; and **the draft names no
countries anywhere in its text or annexes.**

**EU–UK ETS linkage: negotiating, nothing signed.** UK ETS policy overview, last updated 7 September
2026:
> In May 2025, the UK Government and European Union agreed … to work towards linking the UK ETS and
> EU ETS. … The May 2025 Summit **Common Understanding sets out the parameters for a potential linking
> agreement**. **Following approval of the EU negotiating mandate in November 2025, the UK and EU have
> now begun negotiations** … In December 2025, the UK and EU issued a joint statement noting an aim to
> conclude negotiations by the time of the next UK-EU summit.

**The UK is not in Annex III; no CBAM exemption applies in either direction.** Article 2(6)(a) makes
full linkage a legal precondition, not a policy choice. Commentary describing May 2025 as a concluded
agreement to link is wrong.

**UK CBAM starts 1 January 2027** (Finance Act 2026), covering aluminium, cement, fertilisers,
hydrogen, iron and steel; first accounting period calendar 2027, return due 31 May 2028; registration
threshold £50,000. Its policy summary: *"**There are currently no international arrangements or
agreements in place.** … **There are currently no jurisdictions exempt from CBAM.**"*

**Other jurisdictions.** **Türkiye**: Climate Law No. 7552 (Resmî Gazete 09.07.2025); implementing ETS
regulation published 27 August 2026; notably includes a **voluntary complementary carbon price** so
exporters can demonstrate a higher effective carbon cost for CBAM. **Ukraine**: draft ETS law released
15 May 2026, Phase 1 from 2028 with no cap. **Western Balkans / Energy Community**: no exemption; the
only route is the conditioned electricity carve-out in Art. 2(7), requiring an ETS for electricity
"with a price equivalent to the EU ETS … to be finalised by **1 January 2030**" — and Annex III point 2
is still empty. **China/India/Brazil**: no bilateral instrument. The live multilateral channel is the
**Open Coalition on Compliance Carbon Markets** (endorsed COP30 November 2025, launched Florence
7 May 2026; Brazil chair, China co-chair, EU founding participant) — but its remit is MRV, accounting
methodology and credit integrity, and **CBAM is not mentioned in its launch materials**. Groundwork,
not recognition.

**Modelling recommendation:** treat Article 9 as a **zero-value line** in the base case for both
Australia and Brazil. If a sensitivity is wanted, parameterise it as a fraction of the CBAM liability
rather than as a carbon price — the rebate-netting rule makes the effective recognised price highly
uncertain, and for a low-emissions facility it is structurally near zero by construction.

---

## 9. RESOURCE SHUFFLING

**Short answer: there is no operative anti-resource-shuffling rule in EU CBAM law today.** The
Regulation names the practice in a recital and in a biennial reporting mandate, and nowhere else. The
Methodology Act closes it *within a single installation* and expressly *permits* the precursor-level
version of it. The December 2025 proposal does not close it. The live development is Parliament's
July 2026 ENVI report, which would insert a statutory definition.

### 9.1 Article 27 — verbatim, in full

> **Article 27 — Circumvention**
>
> 1. The Commission shall take action in accordance with this Article, based on relevant and objective
> data, to address practices of circumvention of this Regulation.
>
> 2. **Practices of circumvention shall be defined as a change in the pattern of trade in goods, which
> stems from a practice, process or work, for which there is insufficient due cause or economic
> justification other than to avoid, wholly or partially, any of the obligations laid down in this
> Regulation.** Such practice, process or work may consist of, but is not limited to:
> **(a)** slightly modifying the goods concerned to make those goods fall under CN codes which are not
> listed in Annex I, except where the modification alters their essential characteristics;
> **(b)** artificially splitting shipments into consignments the intrinsic value of which does not
> exceed the threshold referred to in Article 2(3).
>
> 3. The Commission shall continuously monitor the situation at Union level with a view to identifying
> practices of circumvention, including by way of market surveillance or on the basis of any relevant
> source of information, such as submissions by, and reporting from, civil society organisations.
>
> 4. A Member State or any party that has been affected by, or has benefited from, any of the
> situations referred to in paragraph 2 may notify the Commission … Interested parties other than
> directly affected or benefited parties, such as environmental organisations and non-governmental
> organisations, which find concrete evidence of practices of circumvention may also notify the
> Commission.
>
> 5. … The Commission shall conclude the investigation within nine months from the date of
> notification. …
>
> 6. Where the Commission … has sufficient reasons to believe that the circumstances referred to in
> **paragraph 2, point (a)** … are occurring in one or more Member States by way of an established
> pattern, it is empowered to adopt delegated acts in accordance with Article 28 **to amend the list of
> goods in Annex I by adding the relevant slightly modified products** …, for anti-circumvention
> purposes.

Two observations. First, the illustrative list in paragraph 2 is about product modification and
shipment splitting, not emissions allocation — though the chapeau is open-ended. Second, and
decisively, **the only remedial power (paragraph 6) is confined to point (a)**. Even if resource
shuffling fell within the *definition* in paragraph 2, Article 27 gives the Commission **no operative
remedy** for it.

Regulation (EU) 2025/2083 amended only point (b), to read: *"artificially splitting imports, including
via non-genuine arrangements, to avoid exceeding the single mass-based threshold."*

### 9.2 Recital 150 — the legislator saw it coming and chose "monitored"

> Practices of circumvention of this Regulation should be **monitored** and addressed by the
> Commission, including where operators could slightly modify their goods without altering their
> essential characteristics, or artificially split shipments … Situations where goods would be sent to
> a third country or region prior to their importation to the Union market …, **or where operators in
> third countries would export their less greenhouse gas emissions intensive products to the Union and
> keep their more greenhouse gas emissions intensive products for other markets, or reorganisation by
> exporters or producers of their patterns and channels of sale and production, or any other kinds of
> dual production and dual sale practices**, with the aim of avoiding the obligations under this
> Regulation, **should also be monitored**.

This describes resource shuffling exactly. It is a recital; its verb is "monitored"; and it feeds no
operative remedy. Recitals are interpretive aids, not obligations.

### 9.3 Article 30(6)(a)(vi) — filed under trade effects, not circumvention

> 6. … Before 1 January 2028, as well as every two years thereafter, the Commission shall present a
> report … The report shall contain at least the following:
> **(a)** an assessment of the impact of the CBAM on: (i) carbon leakage, including in relation to
> exports; (ii) the sectors covered; (iii) internal market, economic and territorial impact throughout
> the Union; (iv) inflation and the price of commodities; (v) the effect on industries using goods
> listed in Annex I; **(vi) international trade, including resource shuffling;** and (vii) LDCs;
> **(b)** an assessment of: (i) the governance system …; (ii) the scope of this Regulation;
> **(iii) practices of circumvention;** (iv) the application of penalties in Member States; …

Note the drafting: resource shuffling sits under *"impact on international trade"*, **not** under
*"practices of circumvention"* in point (b)(iii). The legislator placed it as a trade effect to be
observed, not a practice to be policed. COM(2025) 989 Article 1(20) reproduces (vi) unchanged.

**COM(2025) 783, the Article 30(2) review of December 2025, never uses the term.** Its §4.3 "CBAM
avoidance: circumvention & other practices to unduly lower the CBAM liability" identifies only:
misclassification of goods, under-declaration of quantities, missing declarations, de minimis
misreporting, **misdeclaration of emission intensities**, "abusive practices", and the **scrap
loophole**.

### 9.4 The attribution rules — where the answer actually lies

**The anti-shuffling rule — Article 4(6) of IR 2025/2547:**
> Where goods to which the same functional unit applies are produced using **different production
> routes within an installation, a single production process shall be used encompassing all production
> routes.**

Recital (7): *"the emissions attributable to goods to which the same functional unit applies should be
the **weighted average** of the emissions of all the production routes used within the installation."*
Guidance No. 3 confirms: *"If more than one production route co-exists at your installation for
producing goods which share the same functional unit, then a single production process shall be
defined, so that the different production routes are jointly monitored."*

**For steel this is directly on point.** Annex I §3.15 defines two crude-steel routes — basic oxygen
steelmaking and electric arc furnace. **An integrated mill cannot designate its EAF output for the EU
while averaging away its BF-BOF output — within one installation.** COM(2025) 783 frames Article 4(6)
explicitly as anti-circumvention: *"This measure will reduce the risk of circumvention practices."*

**The escape hatch — Article 4(7):**
> **Splitting an installation into different installations**, with the result that production routes
> otherwise pertaining to a single production process are carried out in separate installations,
> **shall only be allowed where the operators demonstrate valid commercial reasons for this split that
> are related to their economic activity. Commercial reasons shall be considered as valid where
> circumventing Regulation (EU) 2023/956 is not their main purpose or one of their main purposes.**

This is the only explicit anti-circumvention rule in the Methodology Act, and it is a
**subjective-purpose test with the burden on the operator**. It does not reach the far commoner case: a
group that *already* operates a BF-BOF mill and a separate EAF/DRI mill as distinct installations, and
simply routes the EAF/DRI output to the EU. No splitting occurs, so Article 4(7) is never engaged.

**Article 9 — electricity source attribution — DOES NOT APPLY TO STEEL:**
> 1. Where an installation producing goods listed in Annex I … **and not listed in Annex II** …,
> receives, during a reporting period, electricity from multiple sources, and where actual emissions
> are reported …, the embedded indirect emissions … shall be determined **by default** [as the weighted
> average of source emission factors].
> 2. However, where operators provide the verifier with sufficient evidence demonstrating that the
> installation … used, **for a given production process, only electricity from one single source, or
> from a subset of sources**, the embedded indirect emissions … shall be determined … based on the
> emission factor of that single source …

Chapter 72 is in Annex II, so Article 9(1) applies only to cement, fertilisers and sintered ore. **For
iron and steel the electricity-attribution channel is closed by scope, not by anti-abuse design.**

**Article 14(3) — precursor source attribution — DOES apply, and is the live channel:**
> Where operators provide the verifier with sufficient evidence demonstrating that, out of the
> precursors under a given CN code received from multiple installations, the installation producing the
> complex goods used, **for a given production process, only precursors from a single installation, or
> from a subset of installations**, the embedded emissions of those precursors … shall be determined …
> based on the embedded emissions of the precursors obtained from **that single installation** …

Recital (28) frames this purely as proportionality — *"To ensure proportionality with respect to this
default method"* — with **no anti-abuse qualifier at all**, in contrast to Article 4(7). A steelmaker
buying DRI/HBI or pig iron from several suppliers may, on evidence, attribute the cleanest supplier's
precursors to the production process feeding EU-bound goods. Article 15 compounds it: actual values may
be combined with default values precursor-by-precursor.

Note also the granularity: Articles 9(2) and 14(3) both say *"for a given production process"* — not
"for a given consignment" or "for a given customer". You cannot attribute inputs to *export tonnes*;
only to a *production process*.

**The net position.** CBAM's accounting boundary has moved *downward* over time — installation →
production process → single electricity source / subset of precursor installations. Each step narrows
the averaging pool and therefore widens the shuffling headroom. Article 4(6) is the one countervailing
rule, and it operates only within a single installation.

### 9.5 COM(2025) 989 does not close it

**The entire new circumvention limb is one sentence.** Article 1(16):
> in Article 27 (2), the following point (c) is added:
> **‘(c) artificially adjusting the supply chains to make the goods benefit from lower default
> values.’**

This targets **default-value arbitrage** — routing through a jurisdiction with a favourable default —
not actual-value resource shuffling. It is the mirror image of the problem: it addresses producers
*fleeing* actual values, not producers cherry-picking among them.

**The broader hook is "abusive practices"** — new Article 3(35), plus Article 6(2)(f) and Article 6(7)
letting the Commission, on finding "sufficient evidence pointing towards a high risk of abusive
practices **for a combination of goods and origins**", adopt delegated acts within three months
requiring extra evidence "demonstrating that the high risk of abusive practices has not materialised".
Also new: Article 6(2)(e), evidence "that the goods imported during the preceding calendar year **were
produced at the declared installation** and at the actual time of production". Broad enough to *reach*
shuffling in principle — but a power to act later, per goods/origin pair, with no criteria and no test
case.

**Three things cut the other way:**
- It **adds pre-consumer scrap as a CBAM precursor**, closing the scrap channel — which the impact
  assessment says is *"around 40 % of the total scrap intake for both the Aluminium and the Steel
  making processes based on industry intelligence."*
- It **loosens the electricity actual-value conditions**: recital (48) removes the physical-congestion
  criterion and the direct-connection alternative and permits indirect PPAs. Recital 122 of the original
  Regulation had justified those very conditions as being *"To avoid the risk of circumvention"*.
- **SWD(2025) 988 contains zero occurrences of "shuffl".** It frames the issue only as *"abusive
  practices that could occur when actors exploit the possibility of using actual emissions"*, and
  concedes it cannot quantify: *"avoidance practices can typically not be directly observed … given
  that CBAM adjustment only applies from 2026 onwards, avoidance strategies are probably not yet being
  employed by operators."*

**Parliament is going further than the Commission.** The ENVI report adopted **9 July 2026**
(rapporteur Mohammed Chahim, S&D/NL) expands scope to 457 products and, per EPRS, *"would … define
terms such as 'pre-consumer aluminium scrap' and **'resource shuffling'**. It details assessment
criteria for goods and origins at high risk of abusive practices."* Council's general approach (12 June
2026) covers 200 additional goods and does not define the term. **If Parliament's definition survives
trilogue, resource shuffling moves from a reporting mandate into operative law.** Watch trilogue — this
is the single most consequential open item for an exporter relying on installation-level accounting.

### 9.6 What analysts say

**ERCST — the dedicated monograph.** Marcu, Mehling, Cosbey & Fleury, *"The EU CBAM in Practice: The
Challenge of Resource Shuffling"*, 22 October 2025:

> …the term broadly denotes **a reallocation of production or trading relationships that alters the
> attributed carbon intensity of traded goods or energy without changing actual physical emissions**.
> In its simplest form, a regulated entity substitutes low-carbon output for delivery into the
> jurisdiction subject to carbon constraints, while redirecting higher-emission output to unregulated
> markets. The net effect is a statistical reduction in measured emissions, even though global
> atmospheric concentrations remain unchanged.

They classify it as **avoidance, not circumvention** — *"circumvention involves deliberate
non-compliance or falsifying information, whereas avoidance involves legally leveraging the rules of
the mechanism, without technically breaking the letter of the law."* On balance they are **sceptical
that CBAM is badly exposed on steel**:

> …the potential for resource shuffling tends to be largest when emissions intensities vary
> significantly across producers and their installations, as is the case in the power sector … **When
> it comes to industrial process emissions, the variability in carbon intensities of production tends
> to be more limited, at least across comparable production routes.**

**But footnote 32 names this project's exact case:**

> Resource shuffling is, in special cases, also a serious scope 1 problem. **The CBAM assigns zero
> carbon emissions to scrap in aluminum and steel**, which may create a significant risk of resource
> shuffling by the use, or falsely declared use, of scrap in foreign production. … **In the case of
> steel, resource shuffling may arise also from the use of other low carbon input materials, such as
> directly reduced iron (DRI) and/or hot briquetted iron (HBI).**

They assess six options and recommend **none of the in-CBAM fixes**, warning that mandatory defaults
carry WTO exposure (*US–Superfund* 1987, *US–Reformulated Gasoline* 1996: default-based border measures
"may constitute de facto discrimination if domestic producers are assessed on actual performance").
Their response to COM(2025) 989 (17 March 2026) is the direct answer on whether the proposal closes the
gap — and they hope it does not:

> **In the worst-case scenario, this provision would be used to address resource shuffling** – that is,
> used to prevent clean foreign producers from responding to the carbon price imposed by CBAM. Such a
> scenario would disincentivize foreign producers from decarbonizing investment. … it could even result
> in **higher overall global emissions**.
> **They should be elaborated in such a way to ensure that they will not be used to address resource
> shuffling, which ultimately is avoidance—normal market behaviour—not circumvention.**

**Sandbag — the monetary quantification, steel-heavy.** *"A Scrap Game"* (Assous et al., 2024). Their
shuffling scenario is defined physically: *"Imported steel products are all made from electric
furnaces, with **50 % scrap and 50 % DRI or pig iron for flat products** vs. 100 % scrap for long
products."* CBAM fees on EU imports of Chinese goods at full phase-out:

| Sector | Business-as-usual | Resource shuffling | Defaults only |
|---|---|---|---|
| Iron and steel | €447.7m | **€259.2m** | €528.5m |
| Aluminium | €136.7m | €47.7m | €160.7m |
| **Total** | **€590.1m** | **€312.6m** | **€695.2m** |

A **42 % erosion of the steel charge**, China-only. They note *"as no enforcement mechanism is
provided, it is possible that some resource shuffling will occur"*, and that only **5 %** of
transitional-period declarants chose actual data when given the option. Their remedy is systematic
country-level default values.

**CRU Consulting for the German BMWi** — *"Assessing the drivers and scale of potential resource
shuffling under a CBAM"* (ref. ST2254-21, 2021) — is the strongest steel-specific quantification.
Plant-level, two anonymised EU trading partners. It splits the phenomenon into **"output shuffling"**
(imports switching to low-emissions suppliers) and **"input shuffling"** (designating lower-carbon
inputs). Against a 5 %-of-GVA carbon-leakage trigger at €50/tCO2: output shuffling alone gives flats
**5–13 % of GVA**, longs **35–58 %**; output plus input shuffling gives flats **18–31 %**, longs
**45–61 %**. Directly on the BF-BOF-plus-EAF case:

> Resource shuffling might also result from additional use of less carbon intensive EAF steel for
> exports of flat products, which are normally produced via the BF-BOF route. **However, company
> ownership and quality constraints due to technology and customer will limit the potential switch from
> BOF to EAF**, as not all qualities of steel can currently be produced via the EAF route.

**Agora Industry** (Sartor, Cosbey & Shawkat, 2022) reaches the opposite, reassuring conclusion — but
only under scope-1-only coverage, and with a named exception: *"assuming the CBAM only covers scope 1
and not scope 2 emissions, then the risks of resource shuffling leading to carbon leakage from the EU
will be very low for these products by 2030. **The one possible exception may be some steel
production**, although here the risks would fade with time."*

**Fastmarkets** modelling (reported in ERCST): *"In iron and steel … model results show **price
increases up to 50 % higher under national averages** compared to facility-level reporting"* — but also
that switching to national averages *"may decrease total value added in CBAM-covered sectors in the EU
and **increase overall emissions**"*.

**Bruegel** (WP 05/2026, Bahí, Fuchs & Reverdy): *"Countries with carbon pricing policies in place are
**less prone to such trade reshuffling**, as carbon prices typically apply to the production of goods
that serve both domestic and export markets, independent of their ultimate destination."*

**Carbon Market Watch dissents**, treating the framing as an industry artefact: *"**Claims of 'resource
shuffling'** … **have luckily not led to a rehauling or an abuse of default values** – the possibility
to declare actual values is what makes the CBAM a climate measure, instead of a trade protection one."*
**EUROFER** says the opposite: *"the proposed measures are highly uncertain in time and effectiveness
since they do not provide convincing deterrents but only potential ex-post fixes."*

**AEGIS Europe** (2022) proposed the group-level fix aimed precisely at this: *"**This loophole is due
to the definition in Article 7 and Annex III … that provides that the emissions to be measured are
those of the 'Installation'**"*, proposing a switch to "a 'group-responsibility basis'". **Not
adopted.** **CITP** (Winters & Zhang, April 2025) gives the sharpest argument against firm-average
accounting: *"They may be willing to pay €X to halve the carbon intensity of exports to the EU and so
halve their CBAM payment. But if, say, sales to the EU are only 10 % of their output, **averaging means
the investment will now reduce intensity and payment by just 5 % … and so not be worthwhile.**"*

**The California lineage.** CARB's 2011 definition (Cal. Code Regs. tit. 17, §95802(a)(250)): *"any
plan, scheme, or artifice to receive credit based on emissions reductions that have not occurred…"*.
Cullenward & Weiskopf (2013) argued the safe harbours "swallowed the rule". But the ex-post record is
deflationary — Fowlie, Petersen & Reguant (*AEA P&P*, 2021): *"Simulations suggest significant
potential for leakage via resource shuffling. Realized emissions outcomes indicate that **this
potential has not been fully realized**."* Mehling & Ritz (*Oxford Rev. Econ. Pol.* 39(1), 2023) state
the trade-off canonically.

**Terminology is a position marker.** Authors who call it *avoidance* (ERCST, CITP) oppose closing the
gap; authors who call it *circumvention* (Sandbag, EUROFER) want it closed. ERCST distinguishes
electricity-type shuffling (pure re-labelling, emissions-neutral) from goods-market shuffling (which
does move market shares and can move global emissions). Bruegel says "trade reshuffling"; CRU splits
"output" and "input" shuffling; AEGIS says "source shuffling".

*Coverage note: no CEPS publication on CBAM circumvention or resource shuffling was found. IISD has
published nothing under its own name using the term — its* State of Border Carbon Adjustments 2026
*(July 2026) contains zero hits for "shuffl"; the IISD link runs through Aaron Cosbey personally, who
publishes this work via ERCST.*

### 9.7 What this means for the project

For a **dedicated greenfield green HBI plant in Australia or Brazil**, resource shuffling is not the
relevant question — there is no dirty sister production to shuffle against, and the plant's verified
actual emissions are simply low. The relevant provisions are all permissive and cut in your favour:

- **Article 14(3)** lets the plant attribute a *specific* low-carbon pellet or hydrogen supplier to its
  production process, rather than being averaged with a higher-carbon supplier, on evidence to the
  verifier. **Contract for the named installation and keep the evidence trail.** Note this provision
  has no anti-abuse qualifier.
- **Article 4(9)** lets an integrated site fold pelletising, electrolysis and briquetting into one
  joint production process, removing the sintered-ore indirect-emissions exposure entirely.
- **Article 4(6) is the trap to avoid:** do **not** plan a site where a green shaft and a natural-gas
  shaft both produce CN 7203. They collapse into a single production process with one blended
  emissions figure, and the green premium is averaged away.
- **Article 4(7)** means separating them into two legally distinct installations to escape 4(6)
  requires demonstrable commercial reasons unrelated to CBAM avoidance.

**Two forward risks worth a sensitivity:**
1. **Parliament's ENVI position** (9 July 2026) would put a statutory definition of "resource
   shuffling" into the Regulation. Depending on drafting, that could reach a producer that supplies the
   EU from its cleanest of several installations. ERCST's public warning is precisely that this power
   would be used against clean exporters responding rationally to a carbon price.
2. **COM(2025) 989 Article 6(7):** if the Commission designates "green HBI from [country]" as a
   high-risk good/origin combination, extra evidence obligations follow by delegated act within three
   months. A documentation-burden risk, not an emissions-accounting risk.

Note the asymmetric irony flagged by ERCST footnote 32: **HBI is named in the literature as a shuffling
vector**, because it is a low-carbon input that can be selectively directed. A genuine green HBI
exporter is on the right side of that argument substantively, but should expect to carry the
documentary burden of proving it.

---

## 10. SUMMARY OF CORRECTIONS TO THE STARTING PREMISES

| Premise as stated | Verdict |
|---|---|
| Iron, steel and hydrogen are all in CBAM Annex II; indirect emissions excluded in the definitive period | **Correct**, and confirmed by Article 3(2) and Annex I point 3.1 of IR 2025/2547, by Guidance 5b, and by the `N/A` indirect column in the Commission's own default-value tables. Note Annex II now also includes electricity (added by Reg. (EU) 2025/2083). |
| "Delegated Regulation (EU) 2025/2551" is the Methodology Act | **Wrong.** 2025/2551 is the verifier-accreditation delegated act. The Methodology Act is **Implementing Regulation (EU) 2025/2547** of 10 December 2025. |
| Implementing Regulations 2025/2620 and 2025/2621 | **Correct.** 2025/2620 = free allocation adjustment + CBAM benchmarks; 2025/2621 = default values. But 2025/2621's Annexes I and IV were **replaced retroactively by IR (EU) 2026/1740** of 20 July 2026 — use the corrected annexes. |
| Sintered ore (2601 12 00) exception | **Correct.** It is in Annex I but not Annex II, so direct **and** indirect emissions count. Confirmed by Guidance 5d and by recital (8) of IR 2026/1740. |
| Guidance 5d errs in offering hydrogen as a non-Annex-II precursor example | **Correct — the guidance is wrong.** Hydrogen (2804 10 00) is in Annex II. The correct example is sintered ore. Guidance 5b and the Regulation govern. |
| Benchmarks BF-BOF 1.370 / DRI-EAF 0.481 / scrap-EAF 0.072, attributed to IR 2025/2621 | **Numbers real, attribution wrong on two counts.** They are from **IR 2025/2620** (not 2025/2621), and they are the Column B (default) benchmarks for **hot-rolled flat products of CN 7208**, not for crude steel. Crude steel (7206) is 1,288 / 0,424 / 0,027; semi-finished (7207 rolled) is 1,364 / 0,475 / 0,066. The DRI benchmark itself is **0,295 (Column A) / 0,397 (Column B)**. |
| CBAM factor 2.5 % in 2026 → 48.5 % in 2030 → 100 % by 2034 | **Numerically correct as the CBAM phase-in**, but the legally defined "CBAM factor" in Article 10a(1a) of Directive 2003/87/EC is the **complement** (97,5 % / 95 % / 90 % / 77,5 % / 51,5 % / 39 % / 26,5 % / 14 % / none). IR 2025/2620's formulas use the statutory factor. |
| A 2027 Commission report on extending to indirect emissions | **Half right.** The Article 30(2) mandate (which covers indirect-emissions extension) was due *before the end of 2025* and was discharged by **COM(2025) 783** of 16 December 2025, which explicitly deferred the question to a 2027 report. The recurring Article 30(6) report — the one covering **resource shuffling** — is due **before 1 January 2028**. |
| The Omnibus 50-tonne de minimis | **Correct**, but it is a new **Article 2a**, not Article 2(3); it sits at 50 t in a new **Annex VII**; it applies **cumulatively across iron & steel, aluminium, fertilisers and cement**; and it **does not apply to hydrogen or electricity** at all. |
| December 2025 proposal kept indirect emissions out for metals | **Correct.** COM(2025) 989 contains no Annex II amendment. The reasoning is in **COM(2025) 783 §5.1**, not in SWD(2025) 988: metals are excluded from CBAM's indirect scope because EU producers receive **indirect cost compensation (ICC)** state aid, and charging importers too would be double protection. |
| "~180 downstream products" | **Not a legal figure.** The proposal's annexes contain **114 new CN-code entries** (7 added to Iron and steel + a new 107-entry "Combined metal products" table). No product count appears in COM(2025) 989 or SWD(2025) 989. |
| Anti-resource-shuffling rule | **No operative rule exists.** Article 27(2) lists only slight modification and artificial import splitting, and the only remedy (Art. 27(6)) is confined to point (a). **Recital 150 describes resource shuffling almost verbatim and says it "should … be monitored"**; Article 30(6)(a)(vi) files it under *trade impacts*, not *circumvention practices*. What *is* closed is within-installation shuffling (Art. 4(6)–(7) of IR 2025/2547). **Article 14(3) expressly permits single-supplier precursor attribution with no anti-abuse qualifier** — that is the open channel for steel. Art. 9 of IR 2025/2547 (electricity source attribution) does not apply to steel at all, since Chapter 72 is in Annex II. COM(2025) 989 does not close it; **Parliament's ENVI report of 9 July 2026 would insert a statutory definition** — watch trilogue. |
| Article 9 carbon price deduction | Article 9 was **replaced** by Reg. (EU) 2025/2083 (now "carbon price paid in **a third country**"), and the implementing-act power renumbered to **9(5)** — the Commission's own Q&A still cites 9(4) and is stale. The implementing act is **a published draft of 13 May 2026** (Ref. Ares(2026)4841230), feedback closed 10 June 2026, still in inter-service consultation; it would apply retroactively from 1 January 2026. It **does** recognise baseline-and-credit systems (Art. 2(5), recital 10) — so Australia's Safeguard Mechanism qualifies as an *instrument* — but Art. 8(1)(b)(ii) and Annex I §5(b) net out below-baseline emissions as a rebate, so a baseline-compliant facility gets **zero** deduction. DCCEEW's own submission endorses that reading. Brazil's SBCE has no operative price (first-wave surrender ~2031). No Commission list of recognised carbon prices exists; the only full exemption route is Annex III, which lists only Iceland, Liechtenstein, Norway and Switzerland. |

---

## 11. THE SHORT VERSION FOR THE MODEL

For **HBI (CN 7203) produced in Australia or Brazil with green hydrogen and imported into an EU
EAF**, in the definitive period:

**What is charged**
- Direct CO2 from fuels and reducing agents burned at the DRI/briquetting installation.
- Direct CO2 from carbonate process materials (fluxes, flue-gas cleaning reagents).
- Net carbon by mass balance (Annex II B.3.2), with carbon leaving in the HBI, slag and dust
  entering **negatively** at 3.664 t CO2/t C.
- **Full** embedded emissions of any bought-in **sintered ore / pellets (CN 2601 12 00)** —
  **direct AND indirect**, because that CN code is not in Annex II.
- **Direct-only** embedded emissions of any bought-in **hydrogen (CN 2804 10 00)**.

**What is not charged**
- Any electricity consumed at the DRI/briquetting plant (Annex II good → indirect excluded).
- Any electricity consumed making the hydrogen (Annex II precursor → indirect excluded).
- Ocean freight, port handling, inland haulage, and on-site mobile machinery.
- Capex and infrastructure.
- Downstream EU EAF emissions (those are EU ETS).

**Key numbers**
- DRI default if actual data is unavailable: **1,325 t CO2e/t** (both countries fall back to
  "Other Countries and Territories"), × 1.10 in 2026, 1.20 in 2027, 1.30 from 2028.
- Hydrogen default: **10,820 t CO2e/t H2** for both Australia and Brazil (17,740 for "other
  countries", 26,640 if the country of production cannot be identified).
- Free allocation adjustment for 7203: CBAM benchmark **0,295 (actual data)** or **0,397
  (default data)** × CBAM factor × CSCF.
- CBAM factor: 97,5 % (2026) → 51,5 % (2030) → 14 % (2033) → none (2034).

**Design implications**
1. **Verified actual data is existential.** The gap between a verified ~0,05 t CO2e/t and the
   1,325 default (× mark-up) is the whole business case. Budget for an EN ISO/IEC 14065-accredited
   verifier from a European Accreditation member NAB, and for annual re-verification.
2. **Integrate the pellet plant, or contract a named low-carbon pellet supplier.** Bought-in pellets
   drag their *indirect* emissions in. Article 4(9) (joint production process) eliminates the
   exposure; Article 14(3) (single-supplier attribution on evidence) mitigates it.
3. **Do not co-locate a fossil DRI route producing the same CN code.** Article 4(6) would average
   them into one number.
4. **Freight is not a CBAM cost.** Any transport-emissions term in the model is a voluntary
   accounting choice, not a regulatory liability.
5. **Assume no Article 9 carbon-price credit** for Australia or Brazil in the base case. The draft
   implementing act recognises Australia's Safeguard Mechanism as an eligible *instrument*, but nets
   below-baseline emissions out as a rebate — for a low-emissions plant the deduction is zero by
   construction. Brazil pays nothing under the SBCE until roughly 2031.
6. **Contract for named precursor installations.** Article 14(3) lets you attribute a single named
   pellet or hydrogen supplier to the production process, rather than being averaged across suppliers,
   on evidence to the verifier — and unlike Article 4(7) it carries no anti-abuse qualifier.
7. **The Annex II shelter is stable but not permanent.** It rests on EU state-aid indirect cost
   compensation, which the 2026 EU ETS review will revisit; the Commission has explicitly parked
   the indirect-emissions extension for a 2027 report. A model running past 2030 should carry a
   sensitivity in which the EAF/electrolyser electricity emission factor starts to count.

---

## 12. SOURCE INDEX

**EU primary law** (all via `https://publications.europa.eu/resource/celex/<CELEX>`)
- Regulation (EU) 2023/956 — `32023R0956` · consolidated `02023R0956-20251020` · ELI http://data.europa.eu/eli/reg/2023/956/oj
- Regulation (EU) 2025/2083 — `32025R2083` · ELI http://data.europa.eu/eli/reg/2025/2083/oj
- Directive 2003/87/EC consolidated — `02003L0087-20240301` (Article 10a(1a), CBAM factor)
- IR (EU) 2025/2546 (verification) — `32025R2546`
- **IR (EU) 2025/2547 (Methodology Act)** — `32025R2547` · ELI http://data.europa.eu/eli/reg_impl/2025/2547/oj
- IR (EU) 2025/2548 (certificate price) — `32025R2548`
- IR (EU) 2025/2619 (customs information)
- **IR (EU) 2025/2620 (free allocation adjustment + benchmarks)** — `32025R2620`
- **IR (EU) 2025/2621 (default values)** — `32025R2621`
- **IR (EU) 2026/1740 (corrects 2025/2621 Annexes I and IV)** — `32026R1740`
- DR (EU) 2025/2551 (verifier accreditation) — `32025R2551`
- COM(2025) 783 final (Article 30(2) review) — `52025DC0783` · English text http://publications.europa.eu/resource/cellar/05f0b7f5-da86-11f0-8da2-01aa75ed71a1.0001.03/DOC_1
- COM(2025) 989 final (downstream + anti-circumvention proposal) — `52025PC0989` · http://publications.europa.eu/resource/cellar/a837cf93-db4d-11f0-8da2-01aa75ed71a1.0001.03/DOC_1
- COM(2025) 990 final (Temporary Decarbonisation Fund) — `52025PC0990`
- SWD(2025) 988 (impact assessment) — `52025SC0988` · http://publications.europa.eu/resource/cellar/b2be0e26-db4d-11f0-8da2-01aa75ed71a1.0001.03/DOC_1

**Draft Article 9(5) implementing act (unadopted)**
- Initiative 14830 — https://ec.europa.eu/info/law/better-regulation/have-your-say/initiatives/14830
- Draft text — https://ec.europa.eu/info/law/better-regulation/api/download/090166e52d80a42e
- Annexes — https://ec.europa.eu/info/law/better-regulation/api/download/090166e52d80a42f
- News item — https://taxation-customs.ec.europa.eu/news/carbon-price-paid-third-countries-2026-05-13_en
- Call-for-evidence synopsis — https://taxation-customs.ec.europa.eu/document/download/c4049ce0-48c2-43ff-9414-1dd11b3941f4_en
- DCCEEW (Australia) submission, Ares(2026)5894923 — https://ec.europa.eu/info/law/better-regulation/api/download/090166e52ec41661

**Commission guidance (DG TAXUD)**
- Legislation and guidance index — https://taxation-customs.ec.europa.eu/carbon-border-adjustment-mechanism/cbam-legislation-and-guidance_en
- Guidance No. 3 (methods) — https://taxation-customs.ec.europa.eu/document/download/29b9eec7-1a4b-4eb6-ab85-96a0c9e35fd0_en
- Guidance No. 4 (free allocation adjustment) — https://taxation-customs.ec.europa.eu/document/download/3aa2c730-4f1b-4524-8f71-d9cfe4880ac7_en
- **Guidance No. 5b (hydrogen)** — https://taxation-customs.ec.europa.eu/document/download/04ceca3b-466e-4ed9-a1a5-a47c2e5e5be0_en
- **Guidance No. 5d (iron and steel)** — https://taxation-customs.ec.europa.eu/document/download/3eb64513-8255-4c71-875a-bd0a8c7ff772_en
- CBAM Questions and Answers (27 May 2026) — https://taxation-customs.ec.europa.eu/document/download/013fa763-5dce-4726-a204-69fec04d5ce2_en
- Q&A on new Article 27a (8 January 2026) — https://taxation-customs.ec.europa.eu/document/download/412387b4-bf0c-4316-af9e-958857b3dea9_en
- CBAM benchmarks spreadsheet — https://taxation-customs.ec.europa.eu/document/download/9877523c-2a02-4926-a211-aefae7cf6d0d_en
- Corrected default values spreadsheet — https://taxation-customs.ec.europa.eu/document/download/1c05d211-80cb-4aaa-8ef0-e08005a95d7e_en

**Australia**
- Clean Energy Regulator, Safeguard Mechanism — https://cer.gov.au/schemes/safeguard-mechanism · baselines https://cer.gov.au/schemes/safeguard-mechanism/safeguard-baselines · SMCs https://cer.gov.au/schemes/safeguard-mechanism/safeguard-mechanism-credit-units
- Safeguard Mechanism (Crediting) Amendment Act 2023 — https://www.legislation.gov.au/C2023A00014/asmade/text
- Carbon Leakage Review final report (Jotzo) — https://www.dcceew.gov.au/sites/default/files/documents/carbon-leakage-review-final-report.pdf
- ICAP — https://icapcarbonaction.com/en/ets/australian-safeguard-mechanism

**Brazil**
- Lei 15.042/2024 — https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2024/lei/l15042.htm
- Decreto 12.677/2025 (SEMC) — https://www.planalto.gov.br/ccivil_03/_ato2023-2026/2025/decreto/d12677.htm
- SEMC 2026 calendar — https://www.gov.br/fazenda/pt-br/composicao/orgaos/mercado-de-carbono/calendario-previsto
- ICAP — https://icapcarbonaction.com/en/ets/brazilian-greenhouse-gas-emissions-trading-system

**United Kingdom**
- UK ETS policy overview (7 Sep 2026) — https://www.gov.uk/government/publications/uk-emissions-trading-scheme-uk-ets-policy-overview/uk-emissions-trading-scheme-uk-ets-a-policy-overview
- UK CBAM policy summary — https://www.gov.uk/government/publications/carbon-border-adjustment-mechanism-cbam-policy-summary/carbon-border-adjustment-mechanism-cbam-policy-summary
- Qualifying carbon pricing schemes list (27 Aug 2026) — https://www.gov.uk/government/publications/uk-cbam-current-qualifying-carbon-pricing-schemes/carbon-border-adjustment-mechanism-list-of-current-qualifying-carbon-pricing-schemes
- SI 2026/809 — https://www.legislation.gov.uk/uksi/2026/809/made

**Analyst and academic sources (resource shuffling)**
- ERCST, *The EU CBAM in Practice: The Challenge of Resource Shuffling* (22 Oct 2025) — https://ercst.org/the-eu-cbam-in-practice-the-danger-of-resource-shuffling/ (PDF https://ercst.org/?wpdmdl=20246)
- ERCST, response to COM(2025) 989 (17 Mar 2026) — https://ercst.org/response-to-eu-cbam-q4-2025-legislative-proposal/
- Sandbag, *A Scrap Game* (2024) — https://sandbag.be/wp-content/uploads/Sandbag-CBAM-Scrap-Game-2024.pdf
- Sandbag, position on default values — https://sandbag.be/wp-content/uploads/For-a-systematic-use-of-default-value-in-the-CBAM_Position-Paper_Sandbag-1.pdf
- **CRU Consulting for BMWi**, *Assessing the drivers and scale of potential resource shuffling under a CBAM* (2021) — https://www.bundeswirtschaftsministerium.de/Redaktion/DE/Downloads/A/assessing-drivers-and-scale-of-pot-resource-shuffling-under-CBAM.pdf
- Agora Industry, *Getting the transition to CBAM right* (2022) — https://www.agora-industry.org/publications/getting-the-transition-to-cbam-right
- Bruegel WP 05/2026 — https://www.bruegel.org/sites/default/files/2026-03/WP%2005%202026.pdf
- IISD, *State of Border Carbon Adjustments 2026* — https://www.iisd.org/system/files/2026-07/state-of-border-carbon-adjustments-2026.pdf
- Carbon Market Watch on COM(2025) 989 — https://carbonmarketwatch.org/2025/12/18/proposed-cbam-reforms-serve-industrial-lobbies/
- EUROFER response — https://www.eurofer.eu/press-releases/cbam-proposals-single-out-key-loopholes-but-fall-short-of-ensuring-comprehensive-and-structural-solutions-warns-eurofer
- AEGIS Europe on source shuffling — https://aegiseurope.squarespace.com/s/AEGIS-Europe-on-source-shuffling.pdf
- CITP (Winters & Zhang) — https://citp.ac.uk/publications/the-cbam-evolves-the-eus-omnibus-regulation-and-anti-shuffling-measures
- EPRS briefing on the ENVI report — https://www.europarl.europa.eu/RegData/etudes/ATAG/2026/791461/EPRS_ATA(2026)791461_EN.pdf
- Mehling & Ritz, *Oxford Rev. Econ. Pol.* 39(1) 2023 — https://doi.org/10.1093/oxrep/grac043
- Fowlie, Petersen & Reguant, *AEA P&P* 2021 — https://doi.org/10.1257/pandp.20211073
- Cullenward & Weiskopf (Stanford, 2013) — https://law.stanford.edu/index.php?webauth-document=publication/440262/doc/slspublic/Resource%20Shuffling%20-%20Cullenward%20and%20Weiskopf.pdf

---

## 13. RESIDUAL UNCERTAINTIES

1. **The Article 9(5) implementing act is a draft.** Its Articles 2(5), 8(1)(b)(ii) and Annex I §5(b) —
   the provisions that decide the Australian question — could change before adoption.
2. **Brazil's Phase I extension is inferred, not verified.** No decree formally invoking Art. 50(I) of
   Lei 15.042/2024 was located.
3. **Australian government sources** were reachable only via archived copies in the parallel research
   pass; the absence of a DFAT or Senate position on EU recognition is unconfirmed rather than
   established.
4. **CBAM benchmarks for 2026 are provisional** and are to be revised once the final 2026-2030 EU ETS
   benchmarks are published, applying to imports from 1 January 2027 (recital (35), IR 2025/2620).
   Default values must be revised by December 2027 at the latest.
5. **COM(2025) 989 is a proposal.** Council's general approach (12 June 2026) and Parliament's ENVI
   report (9 July 2026) differ from it and from each other on scope (200 vs 457 additional goods) and
   on whether "resource shuffling" gets a statutory definition. Trilogue outcome unknown.
6. **The zero floor on the free allocation adjustment is an inference**, not an express provision (see
   §7.3).
7. **Two drafting errors in IR 2025/2547** are noted in §2.5: points 3.14.2, 3.15.2.1 and 3.15.2.2 cite
   "point B.3.2 of Annex III" where B.3.2 sits in Annex II.
8. **Guidance No. 5d contains a substantive error** (§3.5) offering hydrogen as a non-Annex-II precursor
   example. It has not been corrected as of the 14 August 2026 version.
