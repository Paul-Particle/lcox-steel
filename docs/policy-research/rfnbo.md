# EU rules on renewable (RFNBO) and low-carbon hydrogen — which emissions count, and how electricity is valued

Research note, 7 September 2026. Prepared for a model of H2-DRI production in DE, FR, ES, AU, BR.

## 0. Sourcing note — what is primary and what is not

EUR-Lex blocks direct automated fetching (HTTP 202 / empty body) both on the
`legal-content` URLs and on the ELI redirects. **All EU legal texts quoted below were
obtained verbatim through the `r.jina.ai` text-extraction proxy in front of the EUR-Lex
HTML pages**, i.e. they are the official EUR-Lex English text, machine-transcribed. Where
the proxy dropped an image (the Part C formulas are rendered as images on EUR-Lex) this
is flagged explicitly.

Primary texts obtained in full:

| Text | URL fetched |
|---|---|
| CDR (EU) 2023/1184 — RFNBO "additionality" act | https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R1184 |
| CDR (EU) 2024/1408 — amends 1184 (term alignment) | https://eur-lex.europa.eu/eli/reg_del/2024/1408/oj/eng |
| CDR (EU) 2023/1185 — GHG methodology (RFNBO + RCF) | https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R1185 |
| CDR (EU) 2025/2359 — low-carbon fuels GHG methodology (8 Jul 2025) | https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=OJ:L_202502359 |
| Directive (EU) 2018/2001 consolidated at 20.11.2023 (= RED III) | https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:02018L2001-20231120 |

Official Commission interpretation (soft law, not binding):

| Document | URL |
|---|---|
| Commission Q&A "Implementation of hydrogen delegated acts", v. 26.07.2023 (19 pp.) | https://energy.ec.europa.eu/system/files/2023-07/2023_07_26_Document_Certification_questions.pdf — **now 404 on the live site**; retrieved from the Internet Archive snapshot of 2024-03-01: http://web.archive.org/web/20240301030803if_/https://energy.ec.europa.eu/system/files/2023-07/2023_07_26_Document_Certification_questions.pdf . A later version (14.03.2024) exists on CIRCABC (https://circabc.europa.eu/ui/group/8f5f9424-a7ef-4dbf-b914-1af1d12ff5d2/library/ca8efd4d-cb44-4aec-914d-3d95f95ea293/details) but that page is not machine-fetchable; **quotes below are from the July 2023 version.** |
| Commission Communication C/2025/2983 — *Guidance on the targets ... in Articles 22a, 22b and 25* (adopted as C(2024) 5042 final, 2.9.2024) | https://energy.ec.europa.eu/document/download/0c574279-b71d-4aa0-9403-daf9ea5a8491_en?filename=C_2024_5042_1_EN_ACT_part1_v8.pdf |

Secondary, used only where flagged:

- RSB EU Standard for RFNBOs and RCFs, RSB-STD-11-103 v1.0, Dec 2025 — https://rsb.org/wp-content/uploads/2025/11/rsb-std-11-103-v.1.0-eu-standard-for-rfnbos-and-rcfs.pdf
- CertifHy, *Assessment on Bidding-Zone Equivalence: Brazil*, v1.0, 3 Nov 2025 — https://www.certifhy.eu/content/uploads/2026/04/CertifHy-Decision-on-Brazil_V2.pdf
- PtX Hub / Ecologic Institute, *FAQ on EU requirements for renewable hydrogen and its derivatives* (2024) — https://www.ecologic.eu/sites/default/files/publication/2024/60022-FAQ-EU-requirements-green-hydrogen-and-PtX.pdf
- Commission press release on the low-carbon hydrogen act — https://ec.europa.eu/commission/presscorner/detail/en/ip_25_1743
- EPRS briefing, *Delegated act on low-carbon hydrogen* (2025) — https://www.europarl.europa.eu/RegData/etudes/BRIE/2025/777921/EPRS_BRI(2025)777921_EN.pdf

---

## 0.1 A numbering trap you must not fall into

DR 2023/1184 and DR 2023/1185 cite **"Article 27(3)" of Directive (EU) 2018/2001**.
That is the RED II numbering. RED III (Directive (EU) 2023/2413) renumbered it: the
provision is now **Article 27(6)**. Concretely, in the consolidated RED III:

- Art. 27(6) 1st subpara = the "partially renewable" default (country RES-E share of year n−2);
- Art. 27(6) 2nd subpara = direct connection (what 1184 Art. 3 implements, cited there as "Art 27(3) fifth subparagraph");
- Art. 27(6) 3rd subpara = grid electricity (what 1184 Art. 4 implements, cited there as "sixth subparagraph").

CDR (EU) 2024/1408 (14.3.2024) confirms this by re-basing 1184 on "Article 27(6), fourth
subparagraph". Later texts (2025/2359, RSB, ISCC) all cite Art. 27(6).

**2024/1408 also matters substantively**: it replaced "renewable liquid and gaseous
*transport* fuels of non-biological origin" with "renewable fuels of non-biological
origin" throughout 1184. Recital (1):

> "Those changes expanded the scope of the replaced term in Directive (EU)2018/2001,
> which previously referred only to liquid and gaseous fuels used in the transport sector
> but, following the amendment, it also refers to liquid and gaseous fuels used in the
> electricity sector, **in non-energy purposes in the industrial sector** and in the
> heating and cooling sector."

So the RFNBO electricity rules apply to hydrogen used to reduce iron ore. Everywhere
below I write "RFNBO" for the post-2024/1408 term.

---

## 1. DR 2023/1184 — every route to "fully renewable" electricity

### 1.0 Scope (Art. 1, as amended)

> "This Regulation lays down detailed rules for determining when electricity used for the
> production of renewable fuels of non-biological origin can be considered fully
> renewable. These rules shall apply to the production of renewable fuels of
> non-biological origin via electrolysis and analogously for less common production
> pathways.
> They shall apply **regardless of whether the renewable fuel of non-biological origin is
> produced inside or outside the territory of the Union**."

### 1.1 Definitions that bind the routes (Art. 2)

- (1) **bidding zone** = as in Art. 2(65) Reg. (EU) 2019/943 for Member States, "**or an
  equivalent concept for third countries**".
- (2) **direct line** = Art. 2(41) Dir. (EU) 2019/944.
- (3) **installation generating renewable electricity** = "individual units, or groups of
  units, producing electricity in one or several locations from the same or from different
  renewable sources ... **excluding units producing electricity from biomass and storage
  units**".
- (5) **come into operation** = "starting production of renewable fuels of non-biological
  origin or renewable electricity for the first time or **following a repowering** ...
  requiring investments exceeding 30 % of the investment that would be needed to build a
  similar new installation".
- (7) **imbalance settlement period** = Art. 2(15) Reg. 2019/943, "or an equivalent concept
  for third countries".

### 1.2 Route A — Direct connection (Art. 3)

All three of the following must be evidenced:

- **(a)** the RES installation is connected to the fuel installation **via a direct line**,
  *or* both take place within the same installation;
- **(b)** the RES installation "came into operation **not earlier than 36 months before**
  the installation producing [the RFNBO]"; added RFNBO capacity counts as part of the
  existing installation if at the same site and within 36 months of the initial
  installation coming into operation;
- **(c)** the RES installation is **not connected to the grid**, *or* it is grid-connected
  but "a **smart metering system** that measures all electricity flows from the grid shows
  that **no electricity has been taken from the grid** to produce [the RFNBO]".

Final subparagraph: "If the fuel producer also uses electricity from the grid, it may
count it as fully renewable if it complies with the rules set out in Article 4." So the
routes combine.

Commission Q&A confirmations (Q9–Q11):
- Q9: the **no-public-support exclusion does NOT apply** to directly connected installations.
- Q10: **no PPA required** for a direct connection.
- Q11: both the RES installation and the electrolyser **may also be grid-connected**; then
  Art. 3 governs the direct-line electricity and Art. 4 the grid electricity.

Note what is *absent* from Art. 3: no additionality-style support exclusion, no temporal
correlation, no geographical correlation, no PPA. The only "additionality-like" condition
is the 36-month age rule in (b).

### 1.3 Route B — >90 % renewable bidding zone (Art. 4(1))

> "Fuel producers may count electricity taken from the grid as fully renewable if the
> installation ... is located in a bidding zone where the average proportion of renewable
> electricity **exceeded 90 % in the previous calendar year** and the production ... does
> **not exceed a maximum number of hours** set in relation to the proportion of renewable
> electricity in the bidding zone.
>
> This maximum number of hours shall be calculated by multiplying the total number of
> hours in each calendar year by the share of renewable electricity reported for the
> bidding zone ... The average share of renewable electricity shall be determined by
> dividing the gross final consumption of electricity from renewable sources in the
> bidding zone calculated by analogy to the rules set out in Article 7(2) of Directive
> (EU)2018/2001 by the gross electricity production from all energy sources as defined in
> Annex B to Regulation (EC) No 1099/2008, except from water previously pumped uphill,
> plus imports minus exports of electricity to the bidding zone. **Once the average share
> of renewable electricity exceeds 90 % in a calendar year, it shall be continued to be
> considered to be higher than 90 % for the subsequent five calendar years.**"

Conditions, in full:
1. Location in a bidding zone whose RES-E share exceeded 90 % in the previous calendar year.
2. Annual RFNBO production capped at `8760 h x RES-E share` (full-load-hour cap).
3. Nothing else — **no PPA, no additionality, no temporal correlation, no geographical
   correlation**.
4. GO housekeeping still applies via Art. 19 RED (Q&A Q26: if GOs were issued for that
   electricity they must be cancelled).

Q&A Q12 (data sources for the RES-E share): "Imports and exports are not considered in the
numerator. Where bidding zones are identical to countries, the latest data on the RES-E
that has been published by **Eurostat** are to be used for EU Member States and the latest
data ... published by the **IEA for third countries**. When IEA data is not available, data
from the national statistical institutes may be used. Where bidding zones are not identical
to countries, data from official national statistics have to be used that have been derived
in line with the methodology applied for determining the RES-E share in the **SHARES** tool."

Q&A Q14 (breach of the hour cap): "the hydrogen produced during the maximum number of hours
... (8760 hours x RES-E share) would count as renewable (RFNBO) and hydrogen produced
outside of these hours would count as non-renewable." — i.e. a partial, not a total, loss.

Rationale, recital (5): "Adding additional installations producing renewable electricity is
not necessary given that it can be reasonably assumed that producing renewable hydrogen in
a bidding zone where the share of renewable energy exceeds 90 % allows meeting the 70 %
greenhouse gas saving criterion ... and it may create challenges for the operation of
electricity system."

### 1.4 Route C — bidding zone below 18 gCO2eq/MJ (Art. 4(2)) — the France/Sweden route

> "**Where the conditions set out under paragraph 1 are not met**, fuel producers may count
> electricity taken from the grid as fully renewable if the installation ... is located in
> a bidding zone where the emission intensity of electricity is **lower than
> 18 gCO2eq/MJ**, provided that the following criteria are met:
>
> (a) the fuel producers have **concluded directly, or via intermediaries, one or more
> renewables power purchase agreements** with economic operators producing renewable
> electricity in one or more installations generating renewable electricity **for an amount
> that is at least equivalent to the amount of electricity that is claimed as fully
> renewable** and the electricity claimed is **effectively produced** in this or these
> installations;
>
> (b) the conditions on **temporal and geographical correlation in accordance with Articles
> 6 and 7** are met.
>
> The emission intensity of electricity shall be determined following the approach for
> calculating the average carbon intensity of grid electricity in the methodology ... set
> out in the delegated act adopted pursuant to Article 28(5) of Directive (EU)2018/2001
> **based on latest available data**.
>
> Once the emission intensity of electricity is lower than 18 gCO2eq/MJ in a calendar year,
> the average emission intensity of electricity shall be continued to be considered to be
> lower than 18 gCO2eq/MJ for the subsequent five calendar years."

**Direct answer to the question asked:** yes — **additionality (Article 5) is waived; the
PPA requirement, temporal correlation (Art. 6) and geographical correlation (Art. 7) all
still bind.** The PPA under Art. 4(2)(a) is *not* subject to the Art. 5(a) 36-month
newness test nor the Art. 5(b) no-public-support test — those live only in Art. 5, which
Art. 4(2) does not cross-reference. Confirmed by recital (6):

> "Similarly, in bidding zones, where the emission intensity of electricity is below
> 18 gCO2eq/MJ, **adding further installations producing renewable electricity is not
> required** to achieve the 70 % emissions savings for renewable hydrogen. In such cases,
> it is appropriate to consider electricity taken from the grid as fully renewable
> **provided that the renewable properties of electricity are demonstrated with renewables
> power purchase agreements and by applying criteria for temporal and geographic
> correlation.** Lack of compliance with these conditions and criteria would prevent
> electricity ... from being considered as fully renewable."

Note also the chapeau: Route C is available only "where the conditions set out under
paragraph 1 are not met" — the routes are drafted as an ordered cascade, though in
practice they operate as alternatives.

Q&A Q18 sharpens the PPA condition: a PPA with a **retailer/supplier** does **not**
satisfy Art. 4(2)(a) — "Fuel producers are required to have concluded directly, or via
intermediaries, one or more renewables PPAs with economic operators **producing** renewable
electricity. While electricity suppliers could act as intermediaries (i.e. facilitators of
the contracting), the fuel producer would need to conclude renewables PPAs with economic
operators producing renewable electricity."

Q&A Q17: the GOs attached to the PPA "need to comply with the general requirements in
Article 19 of RED and furthermore carry **the same attributes as the physical installation**
producing the electricity. This includes e.g. the location of the installation, the age of
the installation, and the time of the production. The associated GOs need to be cancelled
before the expiry of the validity period and the volume cancelled shall correspond to that
claimed under the PPA."

**Why 18?** It is the electrolyser back-calculation of the 70 % threshold. 28.2 gCO2eq/MJ
of H2 divided by ~1.56 MJ_e per MJ_H2 (a ~64 % LHV-efficient electrolyser, ~52 kWh/kg)
gives ~18.1 gCO2eq/MJ_e. In other units, 18 gCO2eq/MJ = **64.8 gCO2eq/kWh**.

### 1.5 Route D — curtailment / redispatch (Art. 4(3))

> "Electricity taken from the grid ... may also be counted as fully renewable if the
> electricity ... is consumed **during an imbalance settlement period** during which the
> fuel producer can demonstrate, **based on evidence from the national transmission system
> operator**, that:
> (a) power-generating installations using renewable energy sources were **redispatched
> downwards** in accordance with Article 13 of Regulation (EU)2019/943;
> (b) the electricity consumed for the production of [the RFNBO] **reduced the need for
> redispatching by a corresponding amount**."

No PPA, no additionality, no temporal/geographical correlation. Q&A Q15: "The delegated act
does not set conditions for the **reason** of the redispatch." Q&A Annex: implementing this
in a third country "will ... only be feasible if it set out entities adopting the tasks of
national transmission system operators as well as rules for redispatching."

### 1.6 Route E — general additionality + correlation (Art. 4(4) -> Arts. 5, 6, 7)

> "Where the conditions in paragraphs 1, 2 and 3 are not met, fuel producers may count
> electricity taken from the grid as fully renewable if it complies with the conditions on
> **additionality, temporal correlation and geographic correlation** in accordance with
> Articles 5, 6 and 7."

**Article 5 — Additionality (full text):**

> "The additionality condition ... shall be considered complied with if fuel producers
> **produce an amount of renewable electricity in their own installations that is at least
> equivalent** to the amount of electricity claimed as fully renewable, **or** have
> concluded directly, or via intermediaries, one or more **renewables power purchase
> agreements** with economic operators producing renewable electricity ... for an amount of
> renewable electricity that is at least equivalent to the amount of electricity that is
> claimed as fully renewable and the electricity claimed is **effectively produced** in this
> or these installations, provided that the following criteria are met:
>
> **(a)** The installation generating renewable electricity **came into operation not
> earlier than 36 months before** the installation producing the [RFNBO].
> [2nd subpara] Where an installation ... complied with the requirements ... under a
> renewables PPA with a fuel producer **that has ended**, it shall be considered to have
> come into operation at the same time as the installation producing the [RFNBO] under a
> **new** renewables PPA.
> [3rd subpara] Where **additional production capacity is added** to an existing
> installation producing [the RFNBO], the added capacity shall be considered to have come
> into operation at the same time as the initial installation, provided that the capacity is
> added at the same site and the addition takes place **no later than 36 months** after the
> initial installation came into operation.
>
> **(b)** The installation generating renewable electricity **has not received support in
> the form of operating aid or investment aid**, excluding support received by installations
> **before their repowering**, financial support for **land** or for **grid connections**,
> support that **does not constitute net support**, such as support that is fully repaid and
> support for installations ... supplying installations producing [RFNBOs] used for
> **research, testing and demonstration**."

Q&A Q19 on what counts as aid: "any payments received from public authorities for the
construction of the installations ... and any benefits received from public authorities for
the production of renewable electricity, including **feed-in tariffs, feed-in premiums,
reductions applying for the production, contracts for difference or any direct payments**
linked to the production of renewable electricity. Operating aid or investment aid does
**not** include obligations or restrictions placed on energy consumers, producers or
suppliers such as **renewable energy obligations**." A CfD must be shown ex-ante unlikely to
result in net support and verified ex-post.

Q&A Q20 on "come into operation": step-by-step commissioning => the **first** date counts.

**Transitional exemption (Art. 11):**

> "Article 5, points (a) and (b) shall **not apply until 1 January 2038** to installations
> producing [RFNBOs] that **come into operation before 1 January 2028**. This exemption
> shall not apply to capacity added after 1 January 2028."

So: an electrolyser online before 2028 is free of both the 36-month newness test and the
no-subsidy test for its first decade — but temporal and geographical correlation still
apply throughout.

### 1.7 Common rules, certification, review (Arts. 8–10)

**Art. 8** requires hourly-resolution reporting where relevant, disaggregating electricity
into six buckets: (i) grid electricity not fully renewable + its renewable proportion,
(ii) direct-connection (Art. 3), (iii) Art. 4(1), (iv) Art. 4(2), (v) Art. 4(3),
(vi) Art. 4(4); plus (b) all renewable electricity generated by the contracted
installations regardless of use, and (c) renewable and non-renewable fuel volumes produced.

**Art. 9**: producers inside or outside the EU may use national schemes or Commission-
recognised voluntary schemes under Art. 30(4) RED; Member States must then not require
further evidence.

**Art. 10**: Commission report to EP and Council by **1 July 2028** assessing the impact,
including of temporal correlation. (Mirrored in RED III Art. 27(6), with an explicit power
to amend the methodology "in order to facilitate the ramp-up of the hydrogen industry" —
this is the legal hook for the relaxation debate currently running.)

### 1.8 Combining routes — Q&A Q24 (important for modelling)

> "Fuel producers may combine those options to source renewable electricity provided the
> way the electricity is sourced is fully documented ... **This applies also for electricity
> sourced during the same time interval.** For each way of sourcing electricity, the
> dedicated rules apply. If for instance an electrolyser is fed with 50 % electricity that
> counts as fully renewable and 50 % electricity that is only 40 % renewable, **70 % of the
> total hydrogen produced will be renewable**. The remaining 30 % cannot be made renewable
> by applying the rules of the RFNBO delegated act. The hydrogen produced from the remaining
> 30 % electricity may count as **low carbon hydrogen** under the ... Hydrogen and Gas
> Market Decarbonisation Package."

---

## 2. Temporal and geographical correlation

### 2.1 Temporal correlation (Art. 6) — full text

> "**Until 31 December 2029** the temporal correlation condition referred to in Article 4(2)
> and (4), shall be considered complied with if the [RFNBO] is produced **during the same
> calendar month** as the renewable electricity produced under the renewables power purchase
> agreement **or** from renewable electricity **from a new storage asset that is located
> behind the same network connection point as the electrolyser or the installation
> generating renewable electricity, that has been charged during the same calendar month**
> in which the electricity under the renewables PPA has been produced.
>
> **From 1 January 2030**, the temporal correlation condition shall be considered complied
> with if the [RFNBO] is produced **during the same one-hour period** as the renewable
> electricity produced under the renewables PPA or from renewable electricity from a new
> storage asset ... charged during the same one-hour period ... **Following a notification
> to the Commission, Member States may apply the rules set out in this paragraph from
> 1 July 2027** for [RFNBOs] produced in their territory.
>
> The temporal correlation condition shall **always be considered complied with** if the
> [RFNBO] is produced during a one-hour period where the clearing price of electricity
> resulting from single day-ahead market coupling in the bidding zone ... is **lower or
> equal to EUR 20 per MWh** or **lower than 0,36 times the price of an allowance to emit
> 1 tonne of carbon dioxide equivalent** during the relevant period [EU ETS]."

Confirmed: **same calendar month until 31.12.2029; same hour from 1.1.2030**, with an
optional Member-State opt-in to hourly from **1 July 2027**.

The EUR 20/MWh-or-0.36 x EUA price escape is unconditional — it does not depend on the
route, and it is why a model should treat low-price hours specially. At an EUA price of
EUR 80/t the second limb gives EUR 28.8/MWh; the applicable trigger is whichever limb the
producer can satisfy (they are alternatives, joined by "or").

**Storage.** There is no "storage exemption" in the sense of relaxing correlation. Storage
is an *alternative source of the correlated electricity*, subject to three cumulative
conditions: the asset must be (i) **new**, (ii) **located behind the same network
connection point** as the electrolyser or as the RES installation, and (iii) **charged in
the same calendar month (later: hour)** as the PPA electricity was produced. Separately,
Art. 2(3) excludes storage units from the definition of "installation generating renewable
electricity", so a battery can never itself be the PPA counterparty installation.

Recital (11) sets out the logic: the three ways to show "renewable electricity is
available" are same-period production, stored renewable electricity, and prices so low that
"fossil-based electricity generation is not economically viable".

Q&A Q25: hourly data recording is only strictly needed where hourly checks apply (hourly
correlation from 2030; interconnected-bidding-zone sourcing; curtailment). Q&A Q55: where
hourly correlation applies, the **GHG emission intensity itself must be computed hourly**,
not monthly.

### 2.2 Geographical correlation (Art. 7) — full text

> "1. The geographical correlation condition referred to in Article 4(2) and (4) shall be
> considered complied with if **at least one** of the following criteria relating to the
> location of the electrolyser is fulfilled:
> **(a)** the installation generating renewable electricity under the renewables PPA is
> located, **or was located at the time when it came into operation**, in the **same bidding
> zone** as the electrolyser;
> **(b)** the installation ... is located in an **interconnected bidding zone**, including
> in another Member State, **and electricity prices in the relevant time period on the
> day-ahead market ... in the interconnected bidding zone is equal or higher than in the
> bidding zone where the [RFNBO] is produced**;
> **(c)** the installation ... is located in an **offshore bidding zone that is
> interconnected** with the bidding zone where the electrolyser is located.
>
> 2. Without prejudice to Articles 14 and 15 of Regulation (EU)2019/943, **Member States may
> introduce additional criteria** concerning the location of electrolysers and the
> installation producing renewable electricity ... in order to ensure compatibility of
> capacity additions with the national planning of the hydrogen and electricity grid. Any
> additional criteria shall have no negative impact on the functioning of the internal
> electricity market."

Q&A Q21: interconnected bidding zones **need not be adjacent** (e.g. a sub-sea cable).
Q&A Q22: an "offshore bidding zone" is one comprising only offshore areas; **none had been
created** as of the Q&A. Q&A Q23: **no requirement to monitor physical flows or to book
interconnector capacity.**

Note the "was located at the time when it came into operation" clause in (a): it protects a
PPA against a later bidding-zone reconfiguration (relevant to Germany, where a zone split
is a live proposal).

Practically, limb (b) means an out-of-zone PPA only delivers in hours when the generating
zone is at least as expensive as the electrolyser's zone — i.e. exactly the hours when the
interconnector is *not* exporting into the electrolyser's zone under congestion.

---

## 3. Which bidding zones pass the 18 gCO2eq/MJ test

### 3.1 Is there an official list?

**No.** There is no Commission list of qualifying bidding zones. The regime is:

- Art. 4(2), 2nd subpara: the intensity "shall be determined following the approach for
  calculating the average carbon intensity of grid electricity in the methodology ... set
  out in the delegated act adopted pursuant to Article 28(5)" — i.e. **Part C of the Annex
  to DR 2023/1185** — "**based on latest available data**".
- Art. 9 + recital (14): demonstration runs through **recognised voluntary schemes** and
  their auditors. So it is scheme-verified self-assessment against a prescribed method,
  with a **default table** (Part C, Table A) that the Commission said it would refresh.
- Q&A Q53: "How often does the Commission plan to update the factors provided for Emission
  intensity of electricity in the European Union, provided in Table A ...? Reply: The
  objective is to **update them annually**. Data will be made available on the website of
  the Commission." **I found no evidence that such an annual RFNBO update has ever been
  published.** The Commission energy pages I fetched list no such table.

### 3.2 The two published tables — and they disagree

**(i) DR 2023/1185, Annex, Part C, Table A — "Emission intensity of electricity in the
European Union 2020", source JRC 2022.** This is the legally-designated default for RFNBO
purposes: "If the greenhouse gas emission intensity of electricity is determined at country
level, these values shall be used for electricity sourced in the European Union **until more
recent data becomes available**". It covers **generated** electricity only (no net imports).

| Country | gCO2eq/MJ | Country | gCO2eq/MJ |
|---|---|---|---|
| Austria | 39,7 | Latvia | 39,4 |
| Belgium | 56,7 | Lithuania | 57,7 |
| Bulgaria | 119,2 | Luxembourg | 52,0 |
| Croatia | 55,4 | Malta | 133,9 |
| Cyprus | 206,6 | Netherlands | 99,9 |
| Czechia | 132,5 | Poland | 196,5 |
| Denmark | 27,1 | Portugal | 61,6 |
| Estonia | 139,8 | Romania | 86,1 |
| Finland | 22,9 | Slovakia | 45,6 |
| **France** | **19,6** | Slovenia | 70,1 |
| **Germany** | **99,3** | **Spain** | **54,1** |
| Greece | 125,2 | **Sweden** | **4,1** |
| Hungary | 72,9 | Ireland | 89,4 |
| Italy | 92,3 | | |

On these numbers **only Sweden (4,1) is below 18. France at 19,6 fails.** The RSB EU
Standard (Dec 2025) still reproduces exactly this table for RFNBO work, so certification
practice as of late 2025 is anchored on the 2020 JRC values.

**(ii) DR (EU) 2025/2359, Annex, Part C, Table 5 — "Emission intensity of generated and net
imported electricity in Member States from 2019 to 2023", source JRC 2025 from Eurostat
data.** This is the table for **low-carbon** fuels, not formally for RFNBO — but it is
computed by the same Part C method, it is the *latest available data* on that method, and
it now includes **net imports**. Rule attached to it: "**One of the five most recent
available annual values may be selected** for electricity sourced in the respective
countries."

| Country | 2019 | 2020 | 2021 | 2022 | 2023 |
|---|---|---|---|---|---|
| Austria | 65,2 | 55,6 | 62,7 | 65,3 | 43,8 |
| Belgium | 57,0 | 58,2 | 47,9 | 53,2 | 48,2 |
| Bulgaria | 136,7 | 117,6 | 129,4 | 149,7 | 100,5 |
| Croatia | 76,1 | 63,0 | 79,9 | 87,8 | 64,3 |
| Cyprus | 203,4 | 199,3 | 194,3 | 191,7 | 184,6 |
| Czechia | 146,5 | 132,0 | 142,5 | 146,7 | 127,6 |
| Denmark | 37,1 | 22,6 | 27,5 | 26,3 | 15,9 |
| Estonia | 162,6 | 88,8 | 111,0 | 135,4 | 78,0 |
| Finland | 24,3 | 18,7 | 21,5 | 18,9 | **12,5** |
| **France** | 18,8 | **17,8** | 18,3 | 25,0 | **15,4** |
| **Germany** | 110,5 | 99,7 | 110,2 | 117,2 | 103,8 |
| Greece | 158,3 | 127,9 | 115,5 | 115,4 | 101,1 |
| Hungary | 80,2 | 73,0 | 70,8 | 71,3 | 54,6 |
| Ireland | 100,0 | 92,2 | 110,5 | 101,4 | 85,6 |
| Italy | 97,6 | 92,4 | 97,0 | 108,1 | 87,9 |
| Latvia | 84,7 | 57,5 | 68,4 | 85,9 | 44,6 |
| Lithuania | 33,8 | 31,8 | 35,6 | 32,1 | 19,1 |
| Luxembourg | 86,2 | 76,5 | 76,1 | 87,1 | 70,6 |
| Malta | 122,7 | 129,8 | 120,4 | 121,7 | 115,7 |
| Netherlands | 123,9 | 99,7 | 101,8 | 96,0 | 77,8 |
| Poland | 211,9 | 198,1 | 211,2 | 202,8 | 174,8 |
| Portugal | 81,0 | 64,4 | 53,1 | 56,9 | 39,1 |
| Romania | 108,0 | 91,3 | 88,1 | 93,9 | 73,1 |
| Slovakia | 85,8 | 79,1 | 86,6 | 93,2 | 60,9 |
| Slovenia | 72,3 | 66,4 | 68,8 | 67,9 | 54,2 |
| **Spain** | 69,4 | 54,7 | 52,6 | 60,8 | 47,3 |
| **Sweden** | 4,3 | 3,3 | 3,7 | 3,6 | **3,4** |

**Conflict to flag.** Table A (2020, generation only) gives France **19,6** — fails. Table 5
gives France **17,8** for the same year 2020 and **15,4** for 2023 — passes. Same country,
same nominal method, different numbers, because Table 5 adds net imports and uses a 2025
JRC vintage of the Eurostat data. Since 1184 Art. 4(2) says "based on **latest available
data**", the defensible reading is that a French project today uses the 2023 value of 15,4
and passes; but a conservative auditor working only from the RFNBO delegated act's own
annex would use 19,6 and refuse. **I would not model France's 18 g route as legally certain
without a scheme-level confirmation.** Sweden passes on every published vintage. Germany
(99–118) and Spain (47–69) are nowhere near.

Practically this matters less than it looks, because a French or Swedish electrolyser would
in any case need a PPA plus temporal and geographical correlation under Route C — the only
thing Route C buys you is escaping **additionality** (36-month newness and the no-subsidy
rule), which is exactly what lets a French project sign a PPA with an existing, subsidised
wind or solar farm.

### 3.3 The five-year persistence rule

Text (Art. 4(2), 3rd subpara): "Once the emission intensity of electricity is lower than
18 gCO2eq/MJ in a calendar year, the average emission intensity of electricity shall be
continued to be considered to be lower than 18 gCO2eq/MJ for the subsequent five calendar
years." The parallel wording sits in Art. 4(1) for the 90 % test.

The Commission reads this as a *loss-of-status* rule, not a five-year grace window
(Q&A Q13):

> "A bidding zone is **no longer considered** under Article 4(1) ... to have a share of
> renewable electricity higher than 90 % if the actual share **drops below 90 % for more
> than 5 consecutive years**. **The same principle applies to the calculation of the
> emission intensity of electricity in the bidding zone in the context of the application of
> Article 4(2).**"

RSB operationalises this as: verify the zone was below 18 gCO2eq/MJ "in at least one of the
past six years". Note France's 2022 value (25,0, the nuclear-outage year) is therefore
harmless as long as an adjacent year qualifies.

### 3.4 Third countries with no EU bidding zone (Australia, Brazil)

Two separate questions: (a) what is a "bidding zone" there; (b) what is its emission
intensity.

**(a) Bidding-zone equivalence.** 1184 recital (3): third-country producers "rely on
equivalent concepts provided the objective of this Regulation is maintained ... In case of
bidding zones such concept could be **similar market regulations, the physical
characteristics of the electricity grid, notably the level of interconnection or as a last
resort the country**." The Q&A Annex gives the operational cascade:

> 1. "Certifiers should assess whether at the location of the electrolyser, market
>    regulations applied are **similar to the rules set out for bidding zones in Regulation
>    (EU) 2019/943** ... 'similar' means that there are rules requiring **establishing hourly
>    prices for electricity in a geographical area**. If such rules are in place, the
>    geographical area for which the prices are established should be considered as a bidding
>    zone ...
> 2. If such rules are not in place, certifiers should assess whether the electricity network
>    ... is integrated or whether there are **several separated networks**. If there are
>    several networks, **each network should be considered as a bidding zone** ...
> 3. If the electricity network of the country is integrated and there are no geographically
>    differentiated electricity prices, **the whole country may be considered as one bidding
>    zone** ...
> 4. Where the methodology requires certain conditions ... e.g., on the average proportion of
>    renewable electricity (Article 4(1)), **the emission intensity of electricity (Article
>    4(2))** or the price of electricity (Articles 6 and 7(1)), the conditions **can only be
>    considered as fulfilled if compliance can be demonstrated based on reliable data from
>    official sources**."

*Brazil.* CertifHy has issued a formal determination (3 Nov 2025, v1.0): the **four ONS/CCEE
sub-markets — South (S), Southeast/Central-West (SE/CO), Northeast (NE), North (N) — are
treated as EU-bidding-zone equivalents** for geographic correlation under DR 2023/1184.
Reasons given: formally established market areas; own hourly spot price (*PLD horario*) per
sub-market; no explicit capacity allocation; transmission limits produce locational price
divergence; ONS/CCEE publish audit-grade data. Consequence stated: "hydrogen producers must
source electricity **within the same sub-market** as their electrolyser installation to
comply with RFNBO geographic correlation requirements."

*Australia.* I found **no** published scheme determination. On the Q&A cascade, the NEM
falls squarely in step 1: it sets **regional reference prices per region on a sub-hourly
settlement basis**, so the five NEM regions (QLD, NSW, VIC, SA, TAS) would be the
bidding-zone equivalents, and the **SWIS/WEM (Western Australia) and the NT are separate
networks** => separate zones under step 2. That is my inference from the cascade, not a
scheme decision — treat it as such.

**(b) Emission intensity for a third country.** DR 2023/1185 Part C: "Data on electricity
production and fuel consumption shall be sourced from **IEA Data and statistics** ... For EU
Member States, Eurostat data are more detailed and can be used instead. Where the greenhouse
gas emission intensity is established at the level of bidding zones, **data from official
national statistics of the same level of detail as the IEA data** shall be used." So an
Australian NEM-region or Brazilian sub-market figure must be built bottom-up from AEMO/ONS
generation-by-fuel data using the Part C formula and the Tables 1–3 factors — there is no
default table for non-EU countries. Neither Australia (any region) nor Brazil comes close to
18 gCO2eq/MJ nationally; Brazil's *national* RES-E share (roughly high-80s to low-90s %)
puts the **Art. 4(1) >90 % route** genuinely in play for some sub-markets in some years —
that is the route worth checking for a Brazilian project, not the 18 g route. (The RES-E
share must be evidenced from IEA data, per Q&A Q12.)

---

## 4. DR 2023/1185 — the GHG methodology

### 4.1 The formula (Annex, Part A, point 1)

> "**E = e_i + e_p + e_td + e_u − e_ccs**
> where:
> **E** = total emissions from the use of the fuel (gCO2eq/MJ fuel)
> **e_i** = e_i,elastic + e_i,rigid − e_ex-use : emissions from supply of inputs
> **e_i,elastic** = emissions from elastic inputs
> **e_i,rigid** = emissions from rigid inputs
> **e_ex-use** = emissions from inputs' existing use or fate
> **e_p** = emissions from processing
> **e_td** = emissions from transport and distribution
> **e_u** = emissions from combusting the fuel in its end-use
> **e_ccs** = emission savings from carbon capture and geological storage
> **Emissions from the manufacture of machinery and equipment shall not be taken into
> account.**"

**Confirmed: capital goods / manufacture of machinery and equipment are excluded.** The
identical sentence appears in DR 2025/2359 Part A point 1. Part C adds, for grid
electricity: "The emissions from the **construction and decommissioning and waste
management** of electricity producing facilities are not considered."

Term-by-term, from the Annex:

- **e_i, elastic inputs (points 4, 7, 8).** "Rigid inputs are those whose supply cannot be
  expanded to meet extra demand ... In principle, elastic inputs are those whose supply can
  be increased to meet extra demand." Elastic inputs from an *incorporated process* -> actual
  supply-chain data, "including emissions arising from the extraction of the primary energy
  required to make the input, processing of the input and transportation"; **"Combustion
  emissions related to the carbon content of fuel inputs shall not be included"** (fn 3: "If
  carbon intensities are taken from the table in part B, combustion emissions shall not be
  considered. This is because combustion emissions are counted in processing or in the
  combustion emissions of the final fuel."). Elastic inputs **not** from an incorporated
  process -> the **Part B standard values**. If not listed: "the latest version of the
  JEC-WTW report, the ECOINVENT database, official sources such as the IPCC, IEA or
  government, other reviewed sources such as the E3 and GEMIS database and peer reviewed
  publications."
  ("Incorporated process", fn 2: processes "in the same industrial complex, or that supply
  the input via a dedicated supply infrastructure, or that supply more than half of the
  energy of all inputs" to the fuel production.)
- **e_i, rigid inputs (point 9).** Diversion/opportunity-cost emissions: the lost production
  of electricity, heat or products, valued at the relevant emission factor; for lost
  electricity, "the emission factors to consider are for grid electricity generation in the
  country where the displacement occurred determined according to the appropriate methodology
  set out under points 5 or 6". Baseline for the first 20 years = average output from that
  input over the 3 years before the fuel production started; after 20 years = BAT minimum
  energy performance standards. Explicit examples: coke oven gas, blast furnace gas in a
  steelworks, refinery gas.
- **e_ex-use (point 10).** Avoided emissions of the input's existing use/fate, including
  "the CO2 equivalent of the carbon incorporated in the chemical composition of the fuel that
  would have otherwise been emitted as CO2 into the atmosphere", subject to five alternative
  CO2-source conditions (ETS-activity CO2 priced upstream and bound before **2036** — **2041**
  for non-power-sector CO2; DAC; CO2 from compliant biofuels/bioliquids/biomass fuels not
  already credited; CO2 from compliant RFNBO/RCF combustion; naturally-released geological
  CO2). Excluded: CO2 from fuel deliberately combusted to make CO2, and CO2 already credited
  elsewhere. Point 11 makes the 2036/2041 dates reviewable against the 2040 climate target.
- **e_p (point 12).** "direct atmospheric emissions from the processing itself, from waste
  treatment and from leakages." Q&A Q54: H2 leakage counts as an **energy loss** (raising
  intensity proportionally) until a GWP for H2 is added to the RED annex.
- **e_td (point 16).** "emissions from the storage and distribution of the finished fuels.
  Emissions attributed to inputs e_i shall include emissions from their associated transport
  and storage."
- **e_u (point 13).** "the total combustion emissions of the fuel in use." For pure hydrogen
  this is zero-carbon; for H2-DRI the hydrogen is a reductant, not combusted for energy, and
  the DRI is not itself an RFNBO (see section 7).
- **e_ccs (point 17).** Credit for CO2 permanently stored under Dir. 2009/31/EC; the storage
  and CO2 transport operations' own emissions go into e_p.

Other operative rules:
- **Point 14 / GWPs.** "The greenhouse gases taken into account ... shall be the same as
  specified in **Annex V, part C, point 4** of Directive (EU)2018/2001" — i.e. CO2 = 1,
  **CH4 = 25, N2O = 298** (IPCC AR4 100-year). (DR 2025/2359 instead uses the GWPs in
  Delegated Regulation (EU) 2020/1044 — a divergence between the two acts.)
- **Point 1, averaging.** "The greenhouse gas emissions intensity may be calculated as an
  average for the entire production of fuels occurring during a period of **at most one
  calendar month** ... Where electricity qualifying as fully renewable ... is used as input
  that enhances the heating value of the fuel or intermediate products, the time interval
  shall be **in line with the requirements applying for temporal correlation**" — hourly from
  2030, and each interval must individually clear 70 %.
- **Point 1, mixing.** "If a fuel is a mix of [RFNBOs], recycled carbon fuels and other
  fuels, **all (fuel) types shall be considered to have the same emission intensity**",
  except in co-processing, where a proportional "virtual process" split by energetic value of
  inputs applies.
- **Point 15, co-products.** Allocation at the end of the producing process; by **physical
  causality** where the co-product ratio is variable; by **energy content** where fixed and
  all co-products are fuels/electricity/heat; by **economic value** (3-year average
  factory-gate price) where some co-products are materials with no energy content.
- **Point 3, RFNBO fraction.** "the fraction of [RFNBOs] shall be determined by dividing the
  relevant renewable energy input into the process by the total relevant energy inputs";
  relevant energy for material inputs is the LHV of what enters the fuel's molecular
  structure; **for electricity used to enhance the heating value it is the electricity
  energy**; useful heat = total heat x Carnot efficiency (Annex V, Part C, point 1(b) RED).

### 4.2 Savings and the fossil comparator (Annex, Part A, point 2)

> "**Savings = (E_F − E)/E_F** where E = total emissions from the use of [the RFNBO or RCF];
> E_F = total emissions from the fossil fuel comparator. **For all renewable liquid and
> gaseous transport fuels of non-biological origin and recycled carbon fuels, the total
> emissions from the fossil fuel comparator shall be 94 gCO2eq/MJ.**"

### 4.3 PART B — 'STANDARD VALUES' FOR GREENHOUSE GAS EMISSION INTENSITIES OF ELASTIC INPUTS

Verbatim from the Annex. "The GHG intensities of inputs other than electricity are shown in
the table below". Note fn 3: **when you take a value from this table you use the upstream
column only** — the combustion column is already counted in e_p or in e_u of the final fuel.
Decimal commas are as in the OJ.

**B.1 — Energy carriers [gCO2eq/MJ]**

| Input | Total emissions | Upstream emissions | Combustion emissions |
|---|---|---|---|
| Natural gas | 66,0 | 9,7 | 56,2 |
| Diesel | 95,1 | 21,9 | 73,2 |
| Gasoline | 93,3 | 19,9 | 73,4 |
| Heavy fuel oil | 94,2 | 13,6 | 80,6 |
| Methanol | 97,1 | 28,2 | 68,9 |
| Hard coal | 112,3 | 16,2 | 96,1 |
| Lignite | 116,7 | 1,7 | 115,0 |

**B.2 — Material inputs [gCO2eq/kg]**

| Input | Total emissions |
|---|---|
| Ammonia | 2 351,3 |
| Calcium chloride (CaCl2) | 38,8 |
| Cyclohexane | 723,0 |
| Hydrochloric acid (HCl) | 1 061,1 |
| Lubricants | 947,0 |
| Magnesium sulphate (MgSO4) | 191,8 |
| Nitrogen | 56,4 |
| Phosphoric acid (H3PO4) | 3 124,7 |
| Potassium hydroxide (KOH) | 419,1 |
| Pure CaO for processes | 1 193,2 |
| Sodium carbonate (Na2CO3) | 1 245,1 |
| Sodium chloride (NaCl) | 13,3 |
| Sodium hydroxide (NaOH) | 529,7 |
| Sodium methoxide (Na(CH3O)) | 2 425,5 |
| SO2 | 53,3 |
| Sulphuric acid (H2SO4) | 217,5 |
| Urea | 1 846,6 |

*(Source line in the OJ for the material table in the equivalent 2025 act: "JEC-WTW report
and Renewable Energy Directive calculations".)*

**Caveat if you intend to use B.1 as a general emission-factor table for power generation:**
Part B is the table of *elastic input* intensities, not the table used to build grid
intensities. The grid-intensity build in Part C uses **its own** factors — IPCC 2006
stationary-combustion factors (Table 1/Table 2) plus JEC WTW v5 upstream factors (Table 3),
reproduced in section 4.4 below. The two sets are close but not identical (e.g. hard coal
upstream: 16,2 in Part B vs **15,9** in Part C Table 3; natural gas upstream 9,7 in Part B vs
**12,7** in Part C Table 3; lignite/brown coal upstream 1,7 in both). **Use Part C Tables
1–3 for electricity; use Part B for process inputs.**

### 4.4 PART C — GHG EMISSION INTENSITY OF ELECTRICITY (DR 2023/1185)

Verbatim, method section:

> "The greenhouse gas emission intensity of electricity shall be determined **at the level of
> countries or at the level of bidding zones**. The ... intensity may be determined at the
> level of bidding zones **only if the required data are publicly available**. The calculation
> [of] the carbon intensity of electricity, expressed as gCO2eq/kWh electricity, shall
> consider all potential primary energy sources for electricity generation, type of plant,
> conversion efficiencies and own electricity consumption in the power plant.
>
> The calculation shall consider all carbon equivalent emissions, associated with the
> **combustion and supply** of the fuels used for electricity production ...
>
> Greenhouse Gases other than CO2 shall be converted to CO2eq by multiplying their GWP ...
> over the 100-year time horizon as set out in Annex V, part C, point 4 to Directive
> (EU)2018/2001. **Because of their biogenic origin, CO2 emissions from the combustion of
> biomass fuels are not accounted for, but emissions of CH4 and N2O shall be accounted for.**
>
> For the calculation of the GHG emissions from fuels combustion, the **IPCC default emission
> factors for stationary combustion in the energy industries** shall be used (IPCC 2006). The
> upstream emissions shall include emissions from all the processes and phases required to
> make the fuel ready to supply the power production ... extraction, refining and transport ...
>
> In addition, all the upstream emissions from the cultivation, harvesting, collection,
> processing and transport of biomass shall be considered. **Peat and the components of waste
> materials that are from fossil origins shall be treated as a fossil fuel.**
>
> The fuels used for gross electricity production in electricity only plants are determined
> based on the electricity production and the efficiency of conversion to electricity. In the
> case of **Combined Heat and Power (CHP)**, the fuels used for heat produced in CHP shall be
> counted by considering **alternative heat production with average overall efficiencies of
> 85 %**, while the rest shall be attributed to electricity generation.
>
> For **nuclear** power plants, the conversion efficiency from nuclear heat shall be assumed
> to be **33 %** or data provided by Eurostat or a similar, accredited source.
>
> **No fuels are associated with electricity production from renewables that include hydro,
> solar, wind and geothermal. The emissions from the construction and decommissioning and
> waste management of electricity producing facilities are not considered.** Thus, the carbon
> equivalent emissions associated with the renewable electricity (wind, solar, hydro and
> geothermal) production are considered to be equal to zero."

Formulas (rendered as images on EUR-Lex; the surrounding legend is reproduced verbatim, and
the algebra is written out explicitly in the equivalent Part C of DR 2025/2359, from which
the symbol definitions below are taken):

- `e_gross_prod = SUM_i (c_i,ups + c_i,comb) x B_i` — CO2eq emissions from gross production,
  where `c_i,ups` = upstream CO2eq emission factor of fuel *i*, `c_i,comb` = combustion
  CO2eq emission factor of fuel *i* (including CH4 and N2O expressed as CO2eq), `B_i` =
  consumption of fuel *i* for electricity generation [MJ], `i = 1…k` over the fuels used.
- `E_net = E_gross − E_own − E_pump` — net electricity production, gross minus power-plant
  own consumption minus pumped-storage losses.
- `CI = e_gross_prod / E_net` — "The carbon intensity of net produced electricity shall be
  the total gross GHG emissions for producing or using the net electricity", in gCO2eq/MJ.

Data: "Data on electricity production and fuel consumption shall be sourced from **IEA** Data
and statistics ... For EU Member States, **Eurostat** data are more detailed and can be used
instead. Where the ... intensity is established at the level of bidding zones, data from
**official national statistics of the same level of detail as the IEA data** shall be used.
Fuel consumption data shall include available data at the highest level of detail available
from national statistics: solid fossil fuels, manufactured gases, peat and peat products, oil
shale and oil sands, oil and petroleum products, natural gas, renewables and biofuels,
non-renewable waste and nuclear."

**Note: DR 2023/1185 Part C contains no net-import step.** Net trade was added only in DR
2025/2359 ("Once the national electricity production and its carbon intensity calculated,
net yearly imports from other countries shall be taken into account ... this calculation
should be carried out **iteratively until values converge, at least three times**"). This is
the mechanical reason Table A and Table 5 disagree for France.

**Table 1 — Default emissions factors for stationary combustion [g/MJ fuel, NCV], IPCC 2006.**
Multiply by the GWPs (CO2 1, CH4 25, N2O 298).

| Fuel | CO2 | CH4 | N2O |
|---|---|---|---|
| **Solid fossil fuels** | | | |
| Anthracite | 98,3 | 0,001 | 0,0015 |
| Coking coal | 94,6 | 0,001 | 0,0015 |
| Other bituminous coal | 94,6 | 0,001 | 0,0015 |
| Sub-bituminous coal | 96,1 | 0,001 | 0,0015 |
| Lignite | 101 | 0,001 | 0,0015 |
| Patent fuel | 97,5 | 0,001 | 0,0015 |
| Coke oven coke | 107 | 0,001 | 0,0015 |
| Gas coke | 107 | 0,001 | 0,0001 |
| Coal tar | 80,7 | 0,001 | 0,0015 |
| Brown coal briquettes | 97,5 | 0,001 | 0,0015 |
| **Manufactured gases** | | | |
| Gas works gas | 44,4 | 0,001 | 0,0001 |
| Coke oven gas | 44,4 | 0,001 | 0,0001 |
| Blast furnace gas | 260 | 0,001 | 0,0001 |
| Other recovered gases | 182 | 0,001 | 0,0001 |
| Peat and peat products | 106 | 0,001 | 0,0015 |
| Oil shale and oil sands | 73,3 | 0,003 | 0,0006 |
| **Oil and petroleum products** | | | |
| Crude oil | 73,3 | 0,003 | 0,0006 |
| Natural gas liquids | 64,2 | 0,003 | 0,0006 |
| Refinery feedstocks | 73,3 | 0,003 | 0,0006 |
| Additives and oxygenates | 73,3 | 0,003 | 0,0006 |
| Other hydrocarbons | 73,3 | 0,003 | 0,0006 |
| Refinery gas | 57,6 | 0,001 | 0,0001 |
| Ethane | 61,6 | 0,001 | 0,0001 |
| Liquefied petroleum gases | 63,1 | 0,001 | 0,0001 |
| Motor gasoline | 69,3 | 0,003 | 0,0006 |
| Aviation gasoline | 70 | 0,003 | 0,0006 |
| Gasoline-type jet fuel | 70 | 0,003 | 0,0006 |
| Kerosene-type jet fuel | 71,5 | 0,003 | 0,0006 |
| Other kerosene | 71,5 | 0,003 | 0,0006 |
| Naphtha | 73,3 | 0,003 | 0,0006 |
| Gas oil and diesel oil | 74,1 | 0,003 | 0,0006 |
| Fuel oil | 77,4 | 0,003 | 0,0006 |
| White spirit and SBP | 73,3 | 0,003 | 0,0006 |
| Lubricants | 73,3 | 0,003 | 0,0006 |
| Bitumen | 80,7 | 0,003 | 0,0006 |
| Petroleum coke | 97,5 | 0,003 | 0,0006 |
| Paraffin waxes | 73,3 | 0,003 | 0,0006 |
| Other oil products | 73,3 | 0,003 | 0,0006 |
| **Natural gas** | 56,1 | 0,001 | 0,0001 |
| **Waste** | | | |
| Industrial waste (non-renewable) | 143 | 0,03 | 0,004 |
| Non-renewable municipal waste | 91,7 | 0,03 | 0,004 |

**Table 2 — Default emissions factors for stationary combustion of fuels of biomass origin
[g/MJ, NCV], IPCC 2006.**

| Fuel | CO2 | CH4 | N2O |
|---|---|---|---|
| Primary solid biofuels | 0 | 0,03 | 0,004 |
| Charcoal | 0 | 0,2 | 0,004 |
| Biogases | 0 | 0,001 | 0,0001 |
| Renewable municipal waste | 0 | 0,03 | 0,004 |
| Pure biogasoline | 0 | 0,003 | 0,0006 |
| Blended biogasoline | 0 | 0,003 | 0,0006 |
| Pure biodiesels | 0 | 0,003 | 0,0006 |
| Blended biodiesels | 0 | 0,003 | 0,0006 |
| Pure bio jet kerosene | 0 | 0,003 | 0,0006 |
| Blended bio jet kerosene | 0 | 0,003 | 0,0006 |
| Other liquid biofuels | 0 | 0,003 | 0,0006 |

**Table 3 — Fuel upstream emission factors [gCO2eq/MJ fuel, NCV], JEC WTW v5.**

| Fuel | Emission factor |
|---|---|
| Hard coal | 15,9 |
| Brown coal | 1,7 |
| Peat | 0 |
| Coal gases | 0 |
| Petroleum Products | 11,6 |
| Natural gas | 12,7 |
| Solid biofuels | 0,7 |
| Liquid biofuels | 46,8 |
| Industrial Waste | 0 |
| Municipal waste | 0 |
| Biogases | 13,7 |
| **Nuclear** | **1,2** |

*(Nuclear at 1,2 gCO2eq/MJ of nuclear heat, combined with the mandated 33 % conversion
efficiency, gives ~3,6 gCO2eq/MJ of nuclear electricity ~ 13 gCO2eq/kWh — which is why
France sits just around the 18 g line.)*

Then **Table A** (2020 country values) as reproduced in section 3.2, with the footnote
"Updated data will be made available by the European Commission on a regular basis."

### 4.5 Part C of DR 2025/2359 — the updated factor set (for comparison)

The low-carbon act rebuilt the same tables with the GWP conversion already folded in (its
Tables 3 and 4 are expressed directly in **gCO2eq/MJ**, not g of substance) and moved
upstream factors into a gas-by-gas Table 1. Worth having if you want a consistent
2025-vintage factor set:

**2025/2359 Table 3 — stationary combustion [gCO2eq/MJ, NCV]** (differences from the 2023
table: CH4/N2O columns now CO2eq; oil shale/oil sands raised from 73,3 to 107,0; other
kerosene 71,9 not 71,5):
Anthracite 98,3 / 0,03 / 0,41; Coking coal 94,6 / 0,03 / 0,41; Other bituminous coal 94,6;
Sub-bituminous 96,1; Lignite 101,0; Patent fuel 97,5; Coke oven coke 107,0; Gas coke 107,0
(N2O 0,03); Coal tar 80,7; Brown coal briquettes 97,5; Gas works gas 44,4 / 0,03 / 0,03;
Coke oven gas 44,4; **Blast furnace gas 260,0**; Other recovered gases 182,0; Peat 106,0;
Oil shale and oil sands **107,0**; Crude oil 73,3 / 0,09 / 0,16; NGL 64,2; Refinery
feedstocks 73,3; Additives/oxygenates 73,3; Other hydrocarbons 73,3; Refinery gas 57,6 /
0,03 / 0,03; Ethane 61,6; LPG 63,1; Motor gasoline 69,3; Aviation gasoline 70,0;
Gasoline-type jet fuel 70,0; Kerosene-type jet fuel 71,5; Other kerosene **71,9**; Naphtha
73,3; Gas oil and diesel oil 74,1; Fuel oil 77,4; White spirit/SBP 73,3; Lubricants 73,3;
Bitumen 80,7; Petroleum coke 97,5; Paraffin waxes 73,3; Other oil products 73,3; **Natural
gas 56,1 / 0,03 / 0,03**; Industrial waste (non-renewable) 143,0 / 0,89 / 1,09;
Non-renewable municipal waste 91,7 / 0,89 / 1,09.

**2025/2359 Table 1 — upstream (lifecycle, excl. use-phase combustion) [g of substance per
MJ of product]:** solid fossil fuels anthracite/coking/bituminous CO2 6,50, CH4 0,390,
N2O 0,00026; sub-bituminous & lignite & brown coal briquettes CO2 1,70, CH4 0, N2O 0;
patent fuel / cokes / coal tar / manufactured gases / blast furnace gas / oil shale /
crude oil / NGL / refinery feedstocks / refinery gas / ethane / LPG / bitumen / petroleum
coke / paraffin waxes / other oil products CO2 5,00, CH4 0,228, N2O 0; peat 0/0/0;
gasolines & kerosenes & naphtha & white spirit CO2 13,40, CH4 1,08 x CH4_crude;
gas oil/diesel and lubricants CO2 15,65, CH4 1,09 x CH4_crude; fuel oil CO2 0,
CH4 1,01 x CH4_crude; **natural gas (excl. LNG liquefaction/shipping/regasification)
CO2 4,90, CH4 0,190, N2O 0,00037**; industrial and municipal waste 0; **nuclear heat
CO2 0,50, CH4 0, N2O 0**. Table 2 (material inputs, gCO2eq/kg) is identical to Part B.2 of
2023/1185 except "SO2" is spelled "Sulphur dioxide (SO2)".

---

## 5. Non-qualifying electricity — the three alternatives (DR 2023/1185, Part A, points 5–6)

Point 5, verbatim:

> "Electricity qualifying as fully renewable according to Article 27(3) of Directive
> (EU)2018/2001, shall be attributed **zero** greenhouse gas emissions."

Point 6, verbatim:

> "**One of the three following alternative methods shall be applied during each calendar
> year** to attribute greenhouse gas emissions values to the electricity taken from the grid
> that does not qualify as fully renewable ... and is used to produce [RFNBOs] and recycled
> carbon fuels:
>
> **(a)** greenhouse gas emissions values shall be attributed **according to part C of this
> Annex**. This is without prejudice to the assessment under State aid rules;
>
> **(b)** greenhouse gas emissions values shall be attributed **depending on the number of
> full load hours** the installation ... is operating. Where the number of full load hours is
> **equal or lower than the number of hours in which the marginal price of electricity was
> set by installations producing renewable electricity or nuclear power plants in the
> preceding calendar year** for which reliable data are available, grid electricity ... shall
> be attributed a greenhouse gas emissions value of **zero gCO2eq/MJ**. Where this number of
> full load hours is exceeded, grid electricity ... shall be attributed a greenhouse gas
> emissions value of **183 gCO2eq/MJ**; or
>
> **(c)** **the greenhouse gas emissions value of the marginal unit generating electricity at
> the time of the production** of the [RFNBO] in the bidding zone **may be used if this
> information is publicly available from the national transmission system operator**.
>
> **If the method set in point (b) is used, it shall also be applied to electricity that is
> used to produce [RFNBOs] and recycled carbon fuels and qualifies as fully renewable**
> according to Article 27(3) of Directive (EU)2018/2001."

Reading the conditions precisely:

- **(a) Country / bidding-zone average.** No preconditions beyond Part C's own rule that
  bidding-zone-level determination is allowed "only if the required data are publicly
  available". In practice: Table A default for EU Member States "until more recent data
  becomes available". This is the workhorse.
- **(b) Full-load-hours test — binary 0 or 183.** The comparator is *not* the number of
  renewable/nuclear hours in the current year but in the **preceding** calendar year "for
  which reliable data are available". The test is on the installation's **annual full load
  hours**, not on which hours it ran — so an electrolyser can run in any hours it likes, and
  the whole year's grid electricity is 0 gCO2eq/MJ so long as annual FLH <= that count.
  **The sting is in the last sentence: choosing (b) forces 0/183 onto *all* electricity,
  including electricity that would otherwise be fully renewable and free.** Above the FLH
  threshold, therefore, (b) is catastrophic: 183 gCO2eq/MJ x ~1.56 = ~285 gCO2eq/MJ of H2,
  ten times the 28.2 limit. (183 gCO2eq/MJ = 658,8 gCO2eq/kWh.)
- **(c) Marginal unit at the time of production.** Two conditions: the value must be that of
  the marginal generating unit at the time of production **in the bidding zone**, and the
  information must be **publicly available from the national TSO**. It is permissive ("may
  be used"), and today almost no TSO publishes an audited marginal-unit emission value, so
  this route is largely theoretical in the EU.
- **Choice is annual and exclusive**: "One of the three ... shall be applied **during each
  calendar year**".

For completeness, there is a *fourth*, different thing in RED III Art. 27(6), 1st subpara,
which is about **how much of the hydrogen is renewable** rather than its emission value:

> "Where electricity is used for the production of renewable fuels of non-biological origin,
> either directly or for the production of intermediate products, **the average share of
> electricity from renewable sources in the country of production, as measured two years
> before the year in question, shall be used to determine the share of renewable energy.**"

Keep these two separate in a model: (i) the *renewable share* of the hydrogen (RED
Art. 27(6) 1st subpara, or 100 % via the 1184 routes); (ii) the *emission factor* attached
to the electricity for the 70 % test (1185 points 5–6).

DR 2025/2359 adds a **fourth** option for low-carbon fuels only: "(b) ... based on the
**hourly average greenhouse gas emissions value of the electricity mix** at the time of
production ... **as forecasted by the transmission system operators for the day-ahead market
... two hours before the market gate closure time**", subject to a harmonised methodology or,
until one exists, approval by the competent authority. Its (c) full-load-hours 0/183 and (d)
marginal-technology options mirror 1185's (b) and (c).

---

## 6. The 70 % threshold — arithmetic confirmed

- **Legal basis:** RED III **Article 29a(1)**: "Energy from renewable fuels of non-biological
  origin shall be counted towards Member States' shares of renewable energy and the targets
  referred to in Articles 3(1), 15a(1), **22a(1)**, 23(1), 24(4) and 25(1) only if the
  greenhouse gas emissions savings from the use of those fuels are **at least 70 %**."
  Art. 29a(2) applies the same 70 % to recycled carbon fuels (Art. 25(1) only). DR 2023/1185
  Art. 2 sets the 70 % for RCFs; Art. 3 makes the Annex methodology binding.
- **Fossil comparator:** 94 gCO2eq/MJ (DR 2023/1185, Annex, Part A, point 2). Confirmed.
- **Maximum intensity:** 94 x 0,30 = **28,2 gCO2eq/MJ** of fuel. Confirmed. (Note the
  delegated act does not state 28,2 itself; it falls out of the savings formula.)
- **Per tonne of H2:** at the RED Annex III LHV of hydrogen, 120 MJ/kg:
  28,2 gCO2eq/MJ x 120 MJ/kg = 3 384 gCO2eq/kg = **3,38 tCO2eq per tonne H2**. Confirmed.
- **Per kWh of H2:** 28,2 x 3,6 = **101,5 gCO2eq/kWh_H2** (secondary sources round to 102).
- **Implied maximum grid intensity** at 1,56 MJ_e/MJ_H2 (~52 kWh/kg, ~64 % LHV):
  28,2 / 1,56 = **18,1 gCO2eq/MJ_e = 65 gCO2eq/kWh_e** — the origin of the 18 g threshold.
  A less efficient electrolyser tightens this proportionally.

---

## 7. Scope — does RFNBO apply to iron-ore reduction, and does Art. 22a bind non-EU producers?

### 7.1 Non-energy industrial use is in scope

Yes. Since CDR (EU) 2024/1408, DR 2023/1184 covers RFNBOs used "in non-energy purposes in
the industrial sector". The Commission's Art. 22a guidance says so explicitly in a footnote:

> "To be noted that in the case of ste[e]l manufacturing through the DRI process, **the
> renewable hydrogen used as reducing agent for the direct reduction of iron would be
> considered an RFNBO used in the industry sector.** The pig iron resulting from the direct
> reduction of iron using renewable hydrogen would **not** be considered an RFNBO."

And in the body (section 2.3): "For the purposes of calculating the numerator, derivatives
are considered products obtained as **direct derivatives of hydrogen**, i.e. resulting from a
chemical link between hydrogen with other molecules. Products containing hydrogen but which
are not direct derivatives of hydrogen (e.g. fertilisers), **or products produced using
hydrogen as a reducing agent (e.g. direct reduced iron) would not be considered RFNBOs**."

So: **the H2 is the RFNBO and counts; the DRI/HBI is not an RFNBO and carries no RED III
status of its own.**

### 7.2 What Article 22a requires

Art. 22a(1), 5th subpara, verbatim:

> "Member States shall ensure that the contribution of renewable fuels of non-biological
> origin used for final energy and non-energy purposes shall be **at least 42 % of the
> hydrogen used for final energy and non-energy purposes in industry by 2030, and 60 % by
> 2035**."

Plus, in Art. 22a(1) 1st subpara, a separate **indicative** target: "an indicative increase
of at least **1,6 percentage points** as an annual average" in the renewable share of energy
used for final energy and non-energy purposes in industry, over 2021–2025 and 2026–2030
(waste heat/cold may count up to 0,4 pp, with a matching uplift of the target).

Definitions and accounting:

- **"Industry"** = RED III Art. 2(18a) = NACE Rev. 2 **sections B (mining and quarrying), C
  (manufacturing), F (construction) and section J division 63 (information service
  activities, i.e. data centres)**. Section D (electricity, gas, steam) is **excluded** —
  so hydrogen burned in a power plant is outside Art. 22a.
- **Denominator** = energy content of hydrogen used for final energy and non-energy purposes
  **in industry**, excluding (i) H2 as intermediate product for conventional transport fuels
  and biofuels, (ii) H2 from decarbonising industrial residual gas used to replace that same
  gas, (iii) **H2 produced as a by-product or derived from by-products** (guidance names
  chloro-alkali, sodium chlorate, cracker by-product H2, styrene/ethylene dehydrogenation,
  and **coke oven gas / blast furnace gas hydrogen in iron and steel making**).
  Denominator counts **hydrogen only, not derivatives**; H2 used to make a derivative is
  counted in the Member State where the derivative is produced.
- **Numerator** = energy content of RFNBOs consumed in industry, excluding RFNBOs used as
  intermediates for conventional transport fuels and biofuels. **It is a consumption
  target**: "RFNBOs are accounted in the numerator of the Member State **where they are
  consumed in their final form** in the industry sector."
- **Energy content** from RED Annex III.
- Anything counted must clear the 70 % GHG saving (Art. 29a(1)).
- RFNBOs used for **non-energy** purposes count towards Art. 22a but **not** towards the
  overall 42,5 % Art. 3 target (guidance section 3). Likewise "renewable electricity used to
  produce RFNBOs will not be counted towards the overall EU renewable energy target".

### 7.3 Article 22b derogation

Art. 22b(1), verbatim:

> "A Member State may **reduce** the contribution of renewable fuels of non-biological origin
> used for final energy and non-energy purposes referred to in Article 22a(1), fifth
> subparagraph, **by 20 % in 2030**, provided that:
> (a) that Member State is **on track towards its national contribution** to the binding
> overall Union target set in Article 3(1) ... at least equivalent to its expected national
> contribution in accordance with the formula referred to in Annex II to Regulation (EU)
> 2018/1999; and
> (b) **the share of hydrogen, or its derivatives, produced from fossil fuels** which is
> consumed in that Member State is **not more than 23 % in 2030 and not more than 20 % in
> 2035**.
> Where any of those conditions are not fulfilled, the reduction ... shall cease to apply."

Guidance section 4: "If these conditions are cumulatively met, the RFNBO target laid down in
Article 22a can be reduced to **33,6 % in 2030 and 48 % in 2035**." (i.e. the 20 % reduction
applies at both dates.) The Art. 22b(1)(b) numerator "includes all hydrogen production
processes that use fossil sources, **including those where CO2 is captured and used or
stored**" — so blue hydrogen counts against the Member State here. For the Art. 22b(1)(b)
denominator, "products produced using hydrogen as a reducing agent (e.g. direct reduced iron)
would not fall under the denominator". Notification to the Commission with the NECP is
required (Art. 22b(2)).

### 7.4 Does Art. 22a bind a non-EU hydrogen producer whose H2 is consumed outside the EU?

**No. Your belief is correct.** Three independent reasons, all textual:

1. **Addressee.** Guidance section 2.2, headed "**Obligation addressed to Member States**":
   "The RFNBO target included in Article 22a ... applies to **Member States**, meaning that it
   is the responsibility of Member States to ensure that the contribution of RFNBOs reaches
   the target ... **The obligation, therefore, does not apply directly to hydrogen
   consumers.**" A directive addressed to Member States cannot, of itself, bind a Brazilian or
   Australian producer.
2. **Consumption trigger.** Guidance section 2.3: "RFNBOs are accounted in the numerator of
   the Member State **where they are consumed in their final form** in the industry sector."
   Hydrogen consumed in Australia or Brazil is consumed in no Member State; it enters neither
   numerator nor denominator anywhere.
3. **Territorial denominator.** Guidance section 2.4: "the denominator includes **only
   hydrogen consumed in the industry sector** as defined in Article 2(18a)" — i.e. in the
   Member State's own industry.

What *is* extraterritorial is DR 2023/1184 itself (Art. 1: "regardless of whether the [RFNBO]
is produced inside or outside the territory of the Union"). But that is a **conditional**
rule: it tells a non-EU producer what it must do **if it wants its hydrogen to be recognised
as an RFNBO in the EU**. It imposes no obligation otherwise. A non-EU producer making
hydrogen for a non-EU DRI plant is under no RED III duty at all; certification is purely
voluntary and only becomes relevant if the *hydrogen or an RFNBO derivative* is itself sold
into the EU.

Caveat worth stating in your model: RED III is not the only EU instrument in this chain.
**CBAM** (Reg. (EU) 2023/956) does reach imported iron and steel and does look at embedded
emissions of the imported good — but it is a separate regime with its own rules, and nothing
in RED III gives an RFNBO certificate a CBAM effect.

---

## 8. Imported hydrogen and imported iron

### 8.1 Can imported hydrogen count toward Art. 22a? — Yes.

Art. 22a is a consumption target, and nothing restricts the origin of the RFNBO. RED III
Art. 22a(3) makes this explicit in the other direction: "Member States shall report the
amount of renewable fuels of non-biological origin that they **expect to import and export**
in their integrated national energy and climate plans ... On the basis of that reporting, the
Commission shall develop a **Union strategy for imported and domestic hydrogen** ..."

Conditions for imported H2 to count: it must meet the RFNBO definition (Art. 2(36)), comply
with DR 2023/1184 (which applies extraterritorially), clear 70 % savings under DR 2023/1185
and Art. 29a(1), be verified under a Member State scheme or a Commission-recognised voluntary
scheme (Art. 30(4)), be traced by **mass balance** (Art. 30), and be entered in the **Union
Database** under Art. 31a. Guidance section 2.3 adds a specific constraint for grid-blended
hydrogen: "it is possible to apply a mass balance system ... **provided the consumer would
physically separate the hydrogen from the mixture of gases**. An allocation of the
sustainability and greenhouse gas emission saving characteristics of hydrogen to natural gas
is **not possible** in the absence of such physical separation."

Note also RED III Art. 7(1): two Member States may, by a notified cooperation agreement,
count RFNBOs consumed in one towards the other's share — a Member-State-to-Member-State
device only; it does not extend to third countries.

### 8.2 Can an EU steelmaker importing HBI made with non-EU hydrogen claim anything under RED III? — No.

- The **HBI is not an RFNBO** (Commission guidance, section 7.1 above: DRI produced using
  hydrogen as a reducing agent is not an RFNBO).
- The **hydrogen was consumed outside the EU**, so it cannot enter any Member State's Art. 22a
  numerator (consumption trigger, section 7.4).
- Art. 22a's denominator likewise never sees it, so importing HBI does not *dilute* the
  obligation either — it simply moves the hydrogen consumption out of the EU statistic
  entirely. That is precisely the leakage the Commission flags in guidance section 2.2:
  national mandatory quotas without support "could lead to carbon leakage and additional
  intra-EU or extra-EU imports of products produced with fossil-based hydrogen".
- The only RED III hooks left to an EU steelmaker buying green HBI are soft: **Art. 22a(2)**
  requires Member States to "promote **voluntary labelling schemes** for industrial products
  that are claimed to be produced with renewable energy and renewable fuels of non-biological
  origin", indicating the percentage of renewable energy/RFNBO used "in the raw material
  acquisition and pre-processing, manufacturing and distribution stage, calculated on the
  basis of ... Commission Recommendation (EU) 2021/2279 or ... ISO 14067:2018". That is a
  labelling claim, not a target contribution.

Corollary for your model: **an EU steelmaker's Art. 22a exposure is driven by where the
hydrogen is consumed, not by where the iron is made.** An EU mill importing HBI has no RFNBO
obligation attached to that HBI; an EU mill running its own H2-DRI shaft sits inside the
denominator and creates the demand the 42 %/60 % target is meant to green.

### 8.3 Certification for a third-country producer

Art. 9 of DR 2023/1184: "**Regardless of whether the [RFNBO] is produced inside or outside
the territory of the Union**, fuel producers may make use of national schemes or
international voluntary schemes recognised by the Commission pursuant to Article 30(4) of
Directive (EU)2018/2001 to demonstrate compliance ... a Member State **shall not require**
the suppliers ... to provide further evidence of compliance."

Recital (14): "**Voluntary schemes and national schemes are expected to play an important
role in the implementation and certification of the rules in third countries** as Member
States are required to accept the evidence obtained from recognised voluntary schemes."

Recognised schemes covering RFNBO/RCF (Commission Implementing Decisions of **19 December
2024**, OJ 20.12.2024) — **note this is from secondary reporting (S&P Global, ISCC, REDcert,
IEA policy database); I was not able to fetch the implementing decisions themselves**:

| Scheme | Implementing Decision (as reported) |
|---|---|
| **CertifHy EU RFNBO Voluntary Scheme** | (EU) 2024/3180 |
| **ISCC EU** | (EU) 2024/3176 |
| **REDcert-EU** | (EU) 2024/3194 |

**TÜV SÜD is not a recognised scheme for RFNBO.** TÜV SÜD operates its own private
"green hydrogen" standards (CMS 70 / GreenHydrogen), which are *not* an Art. 30(4) route into
RED III compliance. TÜV bodies (e.g. TÜV Rheinland) act as **certification bodies/auditors
under** CertifHy/ISCC/REDcert — a different role. **RSB** has published an EU RFNBO/RCF
standard (RSB-STD-11-103 v1.0, Dec 2025) and is presumably seeking recognition; I found no
evidence it has been recognised. Treat the recognised list as: CertifHy, ISCC EU, REDcert-EU.

Verification runs under **Commission Implementing Regulation (EU) 2022/996** (rules to verify
sustainability and GHG criteria), with mass balance under RED Art. 30 and entry into the
Union Database under Art. 31a.

Third-country practicalities established above: bidding-zone equivalence is determined by
the scheme (CertifHy's Brazil decision is the first published example); curtailment
(Route D) is effectively unavailable without a TSO-equivalent publishing redispatch
evidence; GOs or an equivalent tracking certificate must be cancelled in favour of the fuel
producer, and where no such certificate system exists the RSB standard requires "a signed
declaration from the renewable electricity producer confirming that no double selling or
double claiming of the renewable attribute or the associated GHG savings has occurred".

---

## 9. Low-carbon hydrogen

**Status: adopted and in force.** **Commission Delegated Regulation (EU) 2025/2359 of 8 July
2025** (C/2025/4674), published in the OJ (L series, 2025/2359), "specifies the methodology
for calculating the greenhouse gas emissions savings from **low-carbon fuels other than
recycled carbon fuels**" (Art. 1). It supplements the **Hydrogen and Gas Market Directive
(EU) 2024/1788**, not RED III — a different legal base from the RFNBO acts. Press release:
https://ec.europa.eu/commission/presscorner/detail/en/ip_25_1743 ; EP briefing:
https://www.europarl.europa.eu/RegData/etudes/BRIE/2025/777921/EPRS_BRI(2025)777921_EN.pdf

**Threshold:** the same **70 % saving against the same 94 gCO2eq/MJ comparator**, hence the
same **28,2 gCO2eq/MJ** ceiling and **3,38 tCO2eq/t H2**. The 70 % figure comes from
Directive (EU) 2024/1788 (definition of low-carbon hydrogen; Art. 2(13) defines low-carbon
fuels as recycled carbon fuels and synthetic gaseous and liquid fuels whose energy content is
derived from low-carbon hydrogen) — the delegated act supplies only the methodology. Covered
pathways: natural gas with CCS ("blue"), methane pyrolysis, and **electrolysis on grid
electricity that cannot qualify as fully renewable** — the last being the one that matters for
France.

**How it differs from RFNBO:**

| | RFNBO (2023/1184 + 1185) | Low-carbon (2025/2359) |
|---|---|---|
| Legal base | RED III Art. 27(6) / 29a(3) | Dir. (EU) 2024/1788 |
| Electricity must be *renewable* | Yes — additionality, temporal, geographical correlation | **No such criteria at all** |
| Nuclear | Counts only indirectly, via the <18 g zone route | Counts directly, at its computed intensity |
| GHG formula | E = e_i + e_p + e_td + e_u − e_ccs | E = e_i + e_p + e_td + e_u − e_ccs **− e_ccu** (extra term for carbon permanently chemically bound in long-lasting products) |
| Grid-electricity valuation options | 3 (Part C average / FLH 0-or-183 / marginal unit) | **4** — adds TSO **day-ahead forecast hourly average** mix intensity, published two hours before gate closure |
| Grid intensity data | Table A, 2020, generation only | **Table 5, 2019–2023, generation + net imports, any one of the five most recent years may be selected** |
| GWPs | RED Annex V Part C point 4 (CH4 25, N2O 298) | Delegated Regulation (EU) 2020/1044 |
| Counts towards RED III targets (Art. 3, 22a, 25) | Yes | **No** |
| Counts against a Member State under Art. 22b(1)(b) | n/a | Yes if produced from fossil fuels, **including with CCS** |

**The nuclear question is explicitly parked.** Recital (7): "The Commission should, **as soon
as possible, initiate an assessment on the potential introduction of alternative approaches
for recognising low-carbon electricity from nuclear power plants**, based on adequate
criteria. **By 30 June 2026, the Commission should launch a public consultation** on a draft
methodology outlining these criteria. In addition, the Commission should assess the impact
and the implications of **evaluating the greenhouse gas emission intensity of electricity
using average values**." Art. 3 repeats this with a **1 July 2028** assessment deadline,
together with "the introduction of a country- or region-specific approach for standard values
for greenhouse gas emission intensities of inputs as reported in part B", and a commitment
that "When assessing changes to the criteria the Commission shall safeguard existing
projects."

Practical read for France: a French electrolyser on grid power **without** a PPA cannot be an
RFNBO, but at a 2023 grid intensity of 15,4 gCO2eq/MJ it produces hydrogen at roughly
15,4 x 1,56 = ~24 gCO2eq/MJ = ~2,9 tCO2eq/t H2 — comfortably **low-carbon** (70 % saving), and
below the RFNBO ceiling too. What it *cannot* do is count towards Art. 22a. That asymmetry —
a French nuclear-powered electrolyser being "clean enough" but not "renewable enough" — is
the specific thing the recital (7) review is about.

Consistency rule between the two regimes (recital (6) of 2025/2359): "it is appropriate to
set out rules ensuring that **the emission intensity of low-carbon hydrogen and the emission
intensity of renewable hydrogen produced in an electrolyser over the same period are always
the same**, and that the reported energy shares are consistent." Footnote 3 adds: "Where both
renewable fuels of non-biological origin and low-carbon fuels are produced in the same
facility, the period chosen under Regulation (EU)2023/1185 and under this methodology shall
be the same."

---

## 10. Modelling checklist and residual uncertainties

**Country-by-country, as I would set it up:**

| Country | Realistic RFNBO route | Grid EF if not qualifying (gCO2eq/MJ) |
|---|---|---|
| **Germany** | Route E (additionality + correlation) or Route A (direct connection). Art. 11 exemption from Art. 5(a)/(b) if online before 2028. | 99,3 (Table A 2020) / 103,8 (Table 5 2023) |
| **France** | Route C (<18 g) **if** the Table 5 vintage is accepted — PPA + correlation, no additionality; otherwise Route E / Route A. | 19,6 (Table A) / 15,4 (Table 5 2023) |
| **Spain** | Route E or Route A. | 54,1 (Table A) / 47,3 (Table 5 2023) |
| **Australia** | Route A or Route E; bidding zone = NEM region (or SWIS/NT as separate networks). No default table — build from AEMO data via Part C. | build bottom-up |
| **Brazil** | Route B (>90 % RES-E) is worth testing per sub-market; otherwise Route A / Route E. Bidding zones = the four ONS/CCEE sub-markets (CertifHy determination). No default table — build from ONS/EPE data via Part C. | build bottom-up |

**Unit conversions to keep in one place:**
- 18 gCO2eq/MJ_e = 64,8 gCO2eq/kWh_e
- 28,2 gCO2eq/MJ_H2 = 101,5 gCO2eq/kWh_H2 = 3,38 tCO2eq/t H2 (LHV 120 MJ/kg)
- 94 gCO2eq/MJ (comparator) = 11,28 tCO2eq/t H2-equivalent
- 183 gCO2eq/MJ_e = 658,8 gCO2eq/kWh_e

**Open points / conflicts I could not resolve from primary sources:**

1. **Which grid-intensity table governs the Art. 4(2) 18 g test today.** Table A (2020, no
   imports, France 19,6) is the one the RFNBO acts point to; Table 5 (2019–2023, with
   imports, France 15,4–17,8) is newer and built by the same Part C method but was enacted
   for a different regime. The Commission promised annual RFNBO updates (Q&A Q53) and, as
   far as I can find, has never published one. This is the single largest legal uncertainty
   for a French case and should be a sensitivity, not an assumption.
2. **The current version of the Commission Q&A.** I used the 26 July 2023 version via the
   Internet Archive. A 14 March 2024 version exists on CIRCABC and I could not fetch it; some
   answers may have moved.
3. **Whether the 90 %/18 g routes survive the Art. 10 review** due by 1 July 2028; RED III
   Art. 27(6) gives the Commission an explicit power to amend the methodology to "facilitate
   the ramp-up of the hydrogen industry". A 2028–2030 relaxation is a live policy risk in
   both directions.
4. **Australia bidding-zone equivalence** is my inference from the Q&A cascade, not a
   published scheme decision.
5. **The three implementing decisions recognising CertifHy/ISCC/REDcert** are cited from
   secondary reporting; I did not fetch the OJ texts.
