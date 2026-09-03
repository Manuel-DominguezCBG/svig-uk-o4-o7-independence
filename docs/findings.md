# Do O4 and O7 double-count the same evidence?

**Analysis of CancerHotspots v1/v2 against SVIG-UK v1.1 (26/05/2026)**
Prepared for Kevin Baker, Wessex Genomics Laboratory Service (Salisbury).
All figures below are reproducible from `scripts/` and are tabulated in `results/`.

---

## Short answer

1. **No, O4 and O7 are not fully independent** — and the dependency is larger and
   more structural than the shared-cohort argument alone suggests.
2. There are **two separate channels of double counting**, and they need different
   fixes:
   * a **cohort channel** — 47% of the observations underpinning CancerHotspots come
     from MSK-IMPACT, which is also in GENIE; and
   * a **definitional channel** — every cancerhotspots.org route to O7 in SVIG-UK v1.1
     is gated on the count of *the same amino-acid change*, which is the same quantity
     O4 measures.
3. The coupling is now **measured, not assumed**: using COSMIC as an independent
   stand-in for GENIE, **195 of the 198 amino-acid changes that reach O7_strong
   (98.5%) also reach O4_strong**. The +8 combination is not a corner case; at the
   top of the scale it is close to automatic.
4. **A blanket cap of Cancer Hotspots at supporting is not the right fix.** It would
   discard genuine positional evidence at the ~18% of hotspots where the position is
   mutated across many different residues (KRAS G12, NRAS Q61, TP53 R273).
5. **Kevin's conditional cap is the right shape, it can be made operational, and its
   impact is small and targeted.** A "leave-one-variant-out" test on the
   CancerHotspots data decides it objectively. It would cap 55.5% of O7-scoring
   changes, but reduce the score of only 15.1% of changes overall, by a mean of 1.34
   points — and drop just **10 changes** from +8 to +4, **4 of which are already on
   the SVIG-UK canonical (O1) list** and so are unaffected in practice.

---

## 1. The cohort channel: how much of Cancer Hotspots *is* MSK?

The v2 workbook (Chang et al., 2018) reports, for every hotspot position, the split
of its mutation count between MSK-IMPACT (`n_MSK`) and retrospective/public cohorts
(`n_Retro`). Across the 1,024 protein-coding SNV hotspot positions:

| Metric | Value |
|---|---|
| Total mutations underpinning the hotspot calls | 25,814 |
| Contributed by MSK-IMPACT | **12,087 (46.8%)** |
| Contributed by retrospective/public cohorts | 13,727 (53.2%) |
| Median per-position MSK fraction | **44.4%** |
| Positions where MSK contributes ≥50% of the evidence | 432 / 1,024 (42.2%) |
| Positions where MSK contributes ≥75% | 88 / 1,024 (8.6%) |
| Positions where MSK contributes 100% | 22 / 1,024 (2.1%) |

MSK-IMPACT is a GENIE contributing centre. So for a typical hotspot, roughly **half
of the observations that made it statistically significant are the same patient
observations an analyst will then count in GENIE for O4.**

The overlap is not uniform, and this matters clinically. Two of the most frequently
assessed lung variants sit at the extreme:

* **EGFR T790M** — 60 mutations at p.790, of which **58 (96.7%) are MSK**.
* **EGFR L858R** — 144 mutations at p.858, of which **114 (79.2%) are MSK**.

Conversely, `n_Retro` (largely TCGA and other public retrospective series) is *not*
in GENIE, so the non-MSK half of the hotspot evidence is genuinely independent of an
O4 GENIE count. This is why a flat downgrade of O7 is too crude: the degree of
overlap is measurable per position, not constant.

*Tables: `results/06_msk_overlap_summary.tsv`, `results/06b_msk_fraction_distribution.tsv`,
MSK fraction per position in `results/03_gene_position_table.tsv`.*

---

## 2. The definitional channel: O7's thresholds are variant-level, not positional

This is the finding I would put in front of the group first, because it does not
depend on any cohort-overlap argument at all.

Read the O7 tiers in Supplementary Table 1 as an analyst applies them:

| O7 strength | Requirement |
|---|---|
| Strong [+4] | ≥50 at the **position** *and* **≥10 for the same amino-acid change** |
| Moderate [+2] | <50 at the position *and* **≥10 for the same amino-acid change** |
| Supporting [+1] | **2–9 for the same amino-acid change** |

Every tier is gated on the recurrence count of the *exact same amino-acid change*.
The supporting tier is defined on nothing else. The positional count only ever
arbitrates between strong and moderate — it never, on its own, admits O7.

O4's missense thresholds are >10 (strong) and 5–10 (moderate) entries for the variant.
So the ≥10 same-change requirement for O7_strong/moderate sits almost exactly on
O4's strong threshold. **A variant that scores O7_strong from cancerhotspots.org has,
by construction, satisfied a recurrence test nearly identical to the one that gave it
O4_strong — using a patient set that is ~47% shared.** The +8 is substantially one
observation counted twice.

Section 4 below tests this empirically against a second database.

Conceptually, O4 and O7 *can* be independent — O7 is meant to capture positional
selection, O4 variant-level recurrence. But as currently specified, the cancerhotspots
route to O7 does not measure the positional signal in isolation.

---

## 3. How mutation-specific are hotspots in practice?

If O7 is to carry genuinely positional information, positions need to be mutated to
several different residues. Across the 1,024 positions:

| Distinct amino-acid changes at the position | Positions | % |
|---|---|---|
| 1 | 213 | **20.8%** |
| 2 | 293 | 28.6% |
| 3 | 231 | 22.6% |
| 4+ | 287 | 28.0% |

Mean 2.85 distinct changes per position; median 3.

The more informative measure is the **fraction of a position's mutations accounted for
by its most common substitution**:

| Metric | Value |
|---|---|
| Mean top-substitution fraction | 0.703 |
| Median top-substitution fraction | 0.692 |
| Positions ≥80% from one change (**mutation-specific**) | 417 / 1,024 (40.7%) |
| Positions 50–80% from one change (**mixed**) | 418 / 1,024 (40.8%) |
| Positions <50% from one change (**position-specific**) | 189 / 1,024 (18.5%) |

So **one hotspot position in five is a single amino-acid change**, and two in five are
≥80% dominated by one change. For those, "the position is a hotspot" and "this variant
is recurrent" are the same statement, and O7 adds nothing to O4.

But the remaining 18.5% are the opposite case. KRAS p.G12 carries **10 distinct
changes over 2,175 mutations**, the most common (G12D) accounting for only 34.8%.
NRAS p.Q61 carries 6 changes over 422 mutations (top 48.3%). TP53 p.R273 carries 6
changes over 609 mutations (top 41.2%). At these positions the hotspot designation
survives comfortably without the variant under assessment, and O7 is telling you
something O4 is not.

*Tables: `results/02_distribution_unique_changes.tsv`,
`results/03_gene_position_table.tsv` (the complete gene × position table, with the
actual substitutions and their individual counts).*

---

## 4. Testing the coupling against a second, largely independent database

The argument so far is structural. It can be tested. SVIG-UK permits COSMIC as an
alternative to GENIE for O4, and the COSMIC Cancer Mutation Census (CMC v104) gives a
per-substitution sample count that is a direct O4-style measure — drawn from a cohort
that is largely *not* the CancerHotspots cohort. GENIE itself was not available, so
COSMIC is used here as a stand-in, to measure how often O4 and O7 fire on the same
variant rather than to make classification calls.

2,570 of the 2,918 hotspot amino-acid changes (88.1%) match a COSMIC CMC substitution.
Across those:

| Metric | Value |
|---|---|
| Spearman ρ, CancerHotspots count vs COSMIC count | **0.826** |
| Pearson r on log10 counts | 0.854 |
| Spearman ρ restricted to changes reaching O7 moderate or strong | 0.757 |

**The recurrence signal that drives O7 and the recurrence signal that would drive O4
are strongly correlated even across two largely different cohorts.** That is the point
worth making to the group: removing the MSK overlap entirely would not make these two
codes independent, because they are measuring the same underlying phenomenon —
positive selection on a specific residue change.

The cross-tabulation makes the practical consequence concrete. Applying O4's
Supplementary Table 1 missense thresholds to the COSMIC counts:

| | O7 Strong | O7 Moderate | O7 Supporting | O7 Not met |
|---|---|---|---|---|
| **O4 Strong** | **195** | 192 | 663 | 76 |
| O4 Moderate | 1 | 1 | 334 | 131 |
| O4 Not met | 0 | 0 | 200 | 538 |
| No COSMIC match | 2 | 12 | 136 | 190 |

**195 of the 198 changes reaching O7_strong (98.5%) also reach O4_strong**, and 192 of
205 O7_moderate changes (93.7%) do as well. At the top of the scale the two codes fire
together almost deterministically.

Because the guideline warns that COSMIC thresholds "should be much higher" than
GENIE's without specifying them, the analysis repeats everything with a conservative
COSMIC-adjusted scheme (≥51 strong, 20–50 moderate). The coupling holds: 185 of 198
O7_strong changes (93.4%) still reach O4_strong.

*Tables: `results/11_cosmic_o4_o7_crosstab.tsv`, `results/12_o4_o7_correlation.tsv`,
`results/14_per_change_o4_o7_points.tsv`.*

---

## 5. What the cap would actually cost

A cap is only worth adopting if its effect is proportionate. Modelling the proposed
rule over all 2,918 changes, with COSMIC standing in for GENIE:

| | O4 thresholds as written | COSMIC-adjusted thresholds |
|---|---|---|
| Changes scoring both O4 and O7 | 1,386 | 765 |
| Changes currently reaching the full +8 | 195 | 185 |
| Changes reduced at all by the cap | 440 (15.1%) | 104 (3.6%) |
| Mean points lost where capped | 1.34 | 1.91 |
| **Changes dropping from +8 to +4** | **10** | **10** |

The cap is therefore **surgical, not disruptive**. It leaves the score of ~85% of
changes untouched, and only ten amino-acid changes lose the full four points.

Those ten, and their status on the SVIG-UK Canonical Variants List (O1 — standalone
evidence, so the cap cannot change their classification):

| Variant | Change / position total | Distinct changes | COSMIC samples | MSK fraction | On O1 list? |
|---|---|---|---|---|---|
| AKT1 p.E17K | 148 / 149 | 2 | 895 | 63.8% | No |
| EGFR p.L858R | 144 / 144 | 1 | 2,872 | 79.2% | **Yes** |
| FGFR3 p.S249C | 114 / 114 | 1 | 1,582 | 64.0% | No |
| U2AF1 p.S34F | 81 / 85 | 2 | 340 | 50.6% | **Yes** |
| PIK3CA p.R88Q | 75 / 75 | 1 | 357 | 38.7% | No |
| EGFR p.T790M | 59 / 60 | 2 | 472 | 96.7% | **Yes** |
| GNA11 p.Q209L | 57 / 59 | 3 | 446 | 32.2% | No |
| SF3B1 p.K700E | 57 / 59 | 3 | 811 | 33.9% | **Yes** |
| TP53 p.M237I | 55 / 64 | 3 | 348 | 43.8% | No |
| PIK3CA p.E726K | 53 / 55 | 2 | 167 | 49.1% | No |

**Four of the ten are already covered by O1**, where the cap is moot. That leaves six
variants materially affected, each of which the group can review individually — and
each of which is exactly the case the proposal is aimed at: a hotspot that is one
amino-acid change, scored twice.

A 4-point reduction is not cosmetic — it can move a variant from Likely oncogenic
(6–9) to VUS (0–5) if the rest of the workup is thin. That is the trade-off the group
needs to take a view on, and it is the reason for recommending the permissive
calibration of the independence test rather than the strict one.

*Tables: `results/13_points_impact_of_cap.tsv`,
`results/16_variants_materially_affected.tsv`, `results/15_canonical_list_overlap.tsv`.*

---

## 6. Proposed operational rule

This is Kevin's suggestion, with the "demonstrably supported by multiple independent
amino acid changes" clause turned into a test that can be applied at the bench.

> **Proposed amendment.** O4 and O7 may be combined where the evidence is independent.
> Where O7 is supported by cancerhotspots.org and O4 is supported by recurrence of the
> **same amino-acid change**, the combined O4 + O7 contribution is **capped at +4
> points**, unless the hotspot position retains independent recurrence after the
> variant under assessment's own amino-acid change is removed.
>
> **Leave-one-variant-out test.** Subtract the count for the VUA's own amino-acid
> change from the position total. The position provides independent evidence if the
> residual is ≥10 mutations **and** at least one other amino-acid change at that
> position is itself recurrent.

Two calibrations were tested (`results/07_cap_rule_sensitivity.tsv`):

| Test | Definition of "other change is recurrent" | O7-scoring changes left uncapped |
|---|---|---|
| **Permissive** (recommended) | ≥2 mutations | 773 / 1,736 (44.5%) |
| **Strict** | ≥10 mutations — i.e. another change would independently meet SVIG-UK's own bar for O7 above supporting | 599 / 1,736 (34.5%) |

The permissive test is recommended: the residual ≥10 threshold already mirrors the
count SVIG-UK uses to lift O7 above supporting, and the strict test additionally
penalises positions such as KIT p.D816 where the positional signal is real but the
absolute counts are small.

The residual is a conservative construction — the CancerHotspots q-value cannot be
recomputed from the published tables, so "would this still be a significant hotspot
without the VUA?" is approximated by residual recurrence rather than re-run
statistically. Any position passing the test would very likely retain significance.

### Worked examples

| Variant | Count / position total | Distinct changes at position | Residual | MSK fraction | O7 | Recommendation |
|---|---|---|---|---|---|---|
| EGFR p.L858R | 144 / 144 | 1 | 0 | 79.2% | Strong | **Cap at +4** |
| EGFR p.T790M | 59 / 60 | 2 | 1 | 96.7% | Strong | **Cap at +4** |
| FGFR3 p.S249C | 114 / 114 | 1 | 0 | 64.0% | Strong | **Cap at +4** |
| MYD88 p.L265P | 37 / 37 | 1 | 0 | 35.1% | Moderate | **Cap at +4** |
| BRAF p.V600E | 833 / 897 | 6 | 64 (V600M 29, V600K 24) | 32.0% | Strong | Combine (max +8) |
| PIK3CA p.H1047R | 537 / 647 | 4 | 110 | 52.7% | Strong | Combine (max +8) |
| IDH1 p.R132H | 570 / 766 | 6 | 196 | 32.3% | Strong | Combine (max +8) |
| KRAS p.G12D | 757 / 2,175 | 10 | 1,418 | 57.9% | Strong | Combine (max +8) |
| KRAS p.G12C | 384 / 2,175 | 10 | 1,791 | 57.9% | Strong | Combine (max +8) |
| NRAS p.Q61R | 204 / 422 | 6 | 218 | 39.1% | Strong | Combine (max +8) |
| TP53 p.R273H | 251 / 609 | 6 | 358 | 47.0% | Strong | Combine (max +8) |

BRAF V600E is the instructive borderline: the position is ≥80% dominated by V600E, so
it is "mutation-specific" as a position — yet 64 mutations across five other changes
remain after removing V600E, which is itself hotspot-level recurrence. The rule
correctly lets it combine. EGFR L858R is the clean opposite: one change, 144 of 144,
four-fifths of it MSK. There is no positional evidence there that is not L858R's own
recurrence.

Overall, the rule caps **963 of 1,736 (55.5%)** of the amino-acid changes that
currently score any O7 from cancerhotspots.org.

*Table: `results/05_o4_o7_double_counting.tsv` gives the recommendation for every
change; `results/09_worked_examples.tsv` for the subset above.*

---

## 7. Two things already in the guideline worth raising

**(a) The precedent for a cap already exists.** Supplementary Figure 1, note 6 states:

> "Restriction: Sum of points awarding from combining O6 and O7 is limited to 4 pts
> (strong) to avoid double counting key attributes of critical domains / residues
> which are also often measured by in silico tools."

The proposed O4 + O7 cap is the same construction applied to a more strongly
correlated pair. The O4 × O7 cell of that matrix currently carries no footnote at
all — it is permitted without restriction, up to +8.

**(b) The intended mitigation is already in the text but is not usable.** The main
text O7 section ends:

> "MSK and the TCGA should be excluded from this count to enable the use of Cancer
> Hotspots (O7) without double counting evidence."

The intent is exactly the one under discussion. But the sentence sits at the end of a
paragraph about in silico protein modelling and DECIPHER plots, so *"this count"* has
no clear antecedent; and read as an instruction to exclude MSK and TCGA cases from the
O4 count it is not practicable — TCGA is not in GENIE, and the GENIE portal offers no
per-centre view aligned to the CancerHotspots cohort. If a cap is adopted, this
sentence should be revised or replaced rather than left standing alongside it.

---

## 8. Recommendation

**Adopt the conditional cap (Kevin's proposal), with the leave-one-variant-out test as
the operational definition of independence.** Specifically:

1. Add a footnote to the O4 × O7 cell of Supplementary Figure 1, mirroring note 6:
   *"Restriction: where O7 is supported by cancerhotspots.org recurrence of the same
   amino-acid change used to support O4, the sum of points from O4 and O7 is limited to
   4 pts (strong). O4 and O7 may be combined in full where the hotspot position is
   supported by ≥10 further mutations arising from other amino-acid changes."*
2. Do **not** cap Cancer Hotspots at supporting across the board — it would
   under-credit the ~18% of hotspots that carry real positional evidence, and would
   weaken classification of genuinely position-driven drivers such as KRAS G12.
   The conditional cap achieves the same protection against double counting while
   changing the score of only 15% of hotspot changes, and dropping only ten from +8
   to +4 — four of which are already covered by O1.
3. Revise the "MSK and the TCGA should be excluded from this count" sentence, which
   is ambiguous and not practicable as written.
4. Consider stating explicitly that O7 via cancerhotspots.org is a **variant-level**
   test as currently written, and that the positional (O7-proper) argument requires
   evidence of multiple distinct changes at the residue.

---

## 9. Caveats and open points

* **Version currency.** CancerHotspots v2 (2018, 24,592 tumours) predates current
  GENIE releases. The MSK fraction measured here is the MSK share *within the hotspot
  resource*, which is the quantity that matters for whether the hotspot call and the
  GENIE count rest on the same patients — but the absolute overlap with a present-day
  GENIE query will differ. Quantifying that directly requires the GENIE dataset,
  which was not available for this analysis.
* **v1 vs v2 differ substantially.** 459 positions in v1 (2016) vs 1,024 SNV positions
  in v2 (2018); only 252 are shared (`results/08_version_comparison.tsv`). The
  cancerhotspots.org website serves v2. Whichever is used should be version-stamped in
  the report, as the guideline already requires for GENIE.
* **Haematological coverage is a real gap.** JAK2 p.V617F, MPL p.W515, CALR exon 9,
  NPM1 exon 12, ASXL1 and SETBP1 are **absent from both the SNV and indel tables** of
  CancerHotspots v2. The guideline flags limited haematological coverage in general
  terms; `results/10_haematology_coverage_check.tsv` quantifies it for named targets.
  O7 is simply unavailable via this resource for several of the most important
  haematological drivers, which is worth stating explicitly in the guidance.
* **Splice hotspots excluded.** 86 of the 1,110 v2 SNV positions are splice-site
  hotspots; O7 is not applicable to them (Supplementary Figure 1, note 22), so they are
  excluded from the analysis set.
* **Indel hotspots not modelled.** The 55 v2 indel regions are ranges (e.g. EGFR
  745–759) rather than single positions, so the "distinct changes at a position"
  construction does not transfer. They are retained in the merged table for
  completeness but excluded from the O7 tier analysis.
* **">50" vs "≥50".** The main text and Supplementary Table 1 differ at exactly 50
  entries for O7_strong. This analysis uses ≥50 (the supplementary table). Worth
  correcting in the next revision.
* **COSMIC is a stand-in for GENIE, not a substitute.** GENIE was not available, so
  the O4 side of the cross-tabulation is modelled on COSMIC CMC v104. COSMIC and GENIE
  differ in composition, ascertainment and duplicate handling, and the guideline warns
  that COSMIC thresholds should be higher. The COSMIC analysis is therefore evidence
  about *how tightly O4-style and O7-style recurrence track each other*, not a
  prediction of the exact GENIE counts any individual variant would return. Repeating
  section 4 against GENIE is the obvious next step if the dataset can be obtained.
* **Protein-level joins are imperfect.** 11.9% of hotspot changes did not match a
  COSMIC substitution, largely because the two resources use different reference
  transcripts. MYD88 p.L265P appears as p.L273P in COSMIC, and GNAS p.R201 (the
  commonest such case here) does not join at all. This affects the join rate, not the
  correlation among matched pairs.
* **Nonsense and stop-loss changes were excluded from O7 scoring.** 247 of the 2,918
  changes introduce or remove a stop codon. O7 is reserved for missense and small
  in-frame variants, and O2 + O7 is not a permitted combination (Supplementary
  Figure 1, notes 4 and 13), so these are marked ineligible rather than scored. Their
  counts still contribute to the position totals, as they do on the website.
* **Third file.** Two workbooks (four sheets) were supplied against three referenced in
  the covering email. If a further CancerHotspots export exists, the merge should be
  re-run to include it.
* **The workbook was validated against the live resource.** A cancerhotspots.org API
  export (1,165 records: 1,110 single residues + 55 in-frame indel regions) was
  cross-checked against the parsed v2 workbook: **all 3,004 amino-acid changes match,
  with 100% agreement on both per-change counts and position totals**
  (`results/17_api_export_validation.tsv`). The numbers modelled here are the numbers
  an analyst sees on the website.
