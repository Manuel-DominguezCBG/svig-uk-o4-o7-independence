# Are O4 and O7 independent lines of evidence?

`svig-uk-o4-o7-independence`

Analysis supporting a proposed amendment to the **ACGS/SVIG-UK guidelines for the
classification of oncogenicity of somatic variants** (v1.1, 26 May 2026), concerning the
combined use of **O4** (enrichment in a somatic variant database, e.g. GENIE) and **O7**
(mutational hotspot, e.g. cancerhotspots.org).

> **The question.** GENIE case counts support O4 and Cancer Hotspots supports O7, but both
> draw a substantial proportion of their evidence from the same underlying MSK cohort. Are
> the two lines of evidence independent, and should the combined O4 + O7 contribution be
> capped?

Requested by **Kevin Baker** (Principal Clinical Scientist, Oncology; Wessex Genomics
Laboratory Service, Salisbury). Analysis by **Manuel Dominguez Becerra**.

📄 **Read the report:** [**web version**](https://manuel-dominguezcbg.github.io/svig-uk-o4-o7-independence/)
· [docs/report.html](docs/report.html) (source) · [docs/findings.md](docs/findings.md)
(the long-form written answer, with every table referenced)

📦 **Deliverables for Kevin:** [deliverables/](deliverables/) — one dated folder per delivery,
with a handover note, figures and the shareable tables. Latest:
[2026-09-21 — CancerHotspots v3 added](deliverables/2026-09-21_cancerhotspots-v3/).

---

## The short answer

O4 and O7 are **not** fully independent — but the dependency is not primarily the shared
MSK cohort, and a blanket cap would be the wrong remedy.

| | |
|---|---|
| **46.8%** | of the mutation observations underpinning CancerHotspots come from MSK-IMPACT, a GENIE contributing centre. Real overlap — but 53.2% is retrospective/public data, largely TCGA, which is *not* in GENIE. |
| **98.5%** | of amino-acid changes reaching O7_strong also reach O4_strong **on an entirely different database** (COSMIC). The two codes fire together almost deterministically at the top of the scale. |
| **20.8%** | of hotspot residues are a *single* amino-acid change. There, "this position is a hotspot" and "this variant is recurrent" are the same statement, and O7 adds nothing to O4. |
| **18.5%** | of residues are genuinely position-driven (KRAS p.G12 carries 10 distinct changes over 2,175 mutations). There, O7 *is* telling you something O4 is not. |

The real dependency is structural, not cohort-based: **O7 as written is a variant-level
test**, scored from the count of the specific amino-acid change, which is the same
quantity O4 measures. Removing the MSK overlap entirely would not make the two codes
independent, because both are measuring positive selection on the same residue change —
they correlate at Spearman ρ = 0.83 *across two largely different cohorts*.

### What is recommended

A **conditional cap**, not a blanket one:

> Where O7 is supported by cancerhotspots.org recurrence of the *same amino-acid change*
> used to support O4, the combined O4 + O7 contribution is **capped at +4 points** —
> unless the hotspot position retains independent recurrence once the variant under
> assessment's own change is removed.
>
> **Leave-one-variant-out test.** Subtract the count for the VUA's own amino-acid change
> from the residue total. The position provides independent evidence if the residual is
> ≥ 10 mutations *and* at least one other amino-acid change at that residue is itself
> recurrent (≥ 2).

This is surgical rather than disruptive. Modelled over all 2,918 hotspot amino-acid
changes it reduces the score of **15.1%**, and only **10** changes lose the full four
points — **four of which are already covered by O1** (Canonical Variants List), where the
cap is moot. EGFR p.L858R is capped (144 of 144 mutations at the residue are L858R itself);
BRAF p.V600E is not (64 mutations across five other changes remain after removing V600E).

The precedent already exists in the guideline: Supplementary Figure 1 note 6 caps O6 + O7
at 4 points for exactly this reason. The O4 × O7 cell currently carries no footnote at all.

### Confirmed against GENIE v20

Kevin supplied AACR Project GENIE v20 (289,869 samples, 242,866 patients), so the COSMIC
stand-in has been replaced with the database SVIG-UK actually names for O4 — with every
count **de-duplicated to one entry per patient**, as the guideline requires.

| | |
|---|---|
| **197 / 198** | O7_strong changes also reach O4_strong on GENIE (99.5%) — tighter than COSMIC's 98.5%. The one exception is in a gene GENIE panels do not sequence. |
| **197 / 198** | still reach O4_strong **with every MSK-IMPACT patient removed**. Excluding MSK — the mitigation the guideline already proposes — does not prevent a single +8. |
| **ρ = 0.68** | CancerHotspots' retrospective (largely TCGA) counts against GENIE's non-MSK counts, per residue: two cohorts with no patients in common, and the recurrence signal still tracks. |
| **10** | changes drop from +8 to +4 under the cap — the same ten COSMIC predicted, the same four protected by O1. |

The leave-MSK-out result is the decisive one. Remove the shared patients entirely and the
same 197 changes still score +8, so the +8 is not produced by the shared cohort. The double
counting is definitional — O4 and O7 both reward recurrence of the same amino-acid change —
so no cohort exclusion can remove it, and only a cap addresses it.

**Haematological representation.** GENIE v20 is 88.7% solid tumour by patient: **4.9%
myeloid** (12,018 patients) and **5.4% lymphoid** (13,131). Haematological patients
contribute more samples each — 1.98 per myeloid patient against 1.12 per solid-tumour
patient, and 3.1 per myeloid patient at MSK — so de-duplication matters most exactly where
haem-onc applies O4: JAK2 p.V617F falls from 3,393 samples to 1,947 patients. MSK-IMPACT
supplies 74% of GENIE's lymphoid patients and 43% of its myeloid ones, so excluding MSK
costs haematology far more than it costs solid tumours. Counted on haematological patients
only — one of the "on-target" options the guideline gives for AML — 83 of the 197 O7_strong
changes still reach O4_strong.

---

## What was done

1. **Merged and validated the source data.** Both CancerHotspots releases (Chang 2016 and
   Chang 2018, four sheets) parsed into one tidy per-allele table, with assertions that
   per-change counts sum to the published residue total and that `n_MSK + n_Retro`
   reconstitutes it. The analysis aborts if either check fails.
2. **Measured the MSK overlap directly.** The v2 workbook splits every residue's count into
   MSK-IMPACT and retrospective/public contributions — the split that makes the O4/O7
   overlap measurable at all.
3. **Characterised how positional hotspots actually are.** Distinct amino-acid changes per
   residue, dominance of the commonest change, and a mutation-specific / mixed /
   position-specific classification for all 1,024 residues.
4. **Assigned O7 tiers** to every amino-acid change under SVIG-UK v1.1 thresholds, excluding
   splice residues (O7 not applicable, note 22), nonsense/stop-loss changes (O2 + O7 not
   permitted, notes 4 and 13) and synonymous changes (no amino acid changes, so a hotspot
   at the residue says nothing about them).
5. **Tested the coupling against a second database.** COSMIC Cancer Mutation Census v104 as
   an independent stand-in for the O4 side, giving the correlation and the O4 × O7
   cross-tabulation under two threshold schemes.
6. **Modelled the proposed cap** and its sensitivity to how "independent evidence" is
   defined, then bounded its clinical impact against the SVIG-UK Canonical Variants List.
7. **Validated the whole thing against the live resource.** The parsed workbook was
   cross-checked against a cancerhotspots.org API export: **all 3,004 amino-acid changes
   match, 100% agreement on both per-change counts and residue totals**. The numbers
   modelled here are the numbers an analyst sees on the website.
8. **Repeated the O4 side against GENIE v20**, counting unique patients, with a
   leave-MSK-out recount, per-lineage ("on-target") counts, and the haematological
   composition of the cohort. GENIE annotates eight genes on a different isoform from
   CancerHotspots (EZH2 p.Y646 is p.Y641 in GENIE); those are aligned by a per-gene residue
   offset, each one listed for review in `results/18c` (generated locally).

The primary analysis set is the v2 SNV table restricted to protein-coding (non-splice)
residues: **1,024 residues, 2,918 distinct amino-acid changes**, of which 1,649 score O7 at
some strength.

---

## What could be done next, and what data it would need

Full detail in [docs/findings.md § 10](docs/findings.md). None of it is required to act on
the recommendation — these would strengthen it or close gaps it leaves open.

| Analysis | What it would settle | Data needed | Available? |
|---|---|---|---|
| ~~GENIE leave-MSK-out recount~~ | **Done** — 197 of 198 O7_strong changes stay O4_strong without MSK. See the GENIE section above | GENIE v20, supplied by Kevin Baker | ✓ |
| Would the residue still be a hotspot without MSK? | Whether a position was called *only* because of MSK cases | Approximation possible **now** from the published `n_Retro` column; full re-derivation needs the per-sample input to Chang 2018 + its background rate model | Partly in hand; the full version is a large piece of work |
| Haematological hotspot supplement | Restores O7 for JAK2 p.V617F, MPL p.W515, CALR, NPM1, ASXL1, SETBP1 — **none of which appear in CancerHotspots v2 at all** | COSMIC restricted to haematopoietic/lymphoid tissue, the GENIE haem subset, or BeatAML / MDS-CHIP series | GENIE v20 now in hand, and `results/25` already counts the named drivers by lineage |
| Per-database O4 threshold calibration | What COSMIC count is equivalent to GENIE's ≥10 and ≥50 — the guideline says only "much higher" | GENIE + COSMIC joined per substitution | Both datasets now in hand — the cheapest item left |
| Indel hotspot regions | Whether the same double counting applies to the 55 in-frame indel regions, excluded here because they are ranges not residues | Per-sample indel calls, plus an SVIG-UK decision on what independence means across a range | Needs a definitional call as much as data |
| Transcript harmonisation | Recovers the 11.9% of changes that failed to join to COSMIC on transcript choice (MYD88 p.L265P is p.L273P in COSMIC) | MANE Select + VEP or Mutalyzer | Self-contained; the easiest item on the list |
| **Impact on real reported cases** | How many previously classified variants change SVIG-UK class under the cap — a different order of evidence from modelling over the resource | De-identified retrospective extract of variants classified under v1.1 at Wessex GLH, with **per-code point assignments**, not just final classes | Held in-service; needs local IG sign-off, not an external application |
| Currency refresh | Whether tier assignments hold against a modern cohort (v2 is 2018, 24,592 tumours) | A future cancerhotspots.org release, or GENIE | Re-running the pipeline against it is one command |

---

## Repository layout

```
data/raw/        Source workbooks exactly as supplied (unmodified)
data/interim/    Merged, tidied per-allele table produced by script 01
data/reference/  Pinned OncoTree release (2025-10-03) used to assign myeloid / lymphoid / solid
docs/            The report (report.html, built to index.html for Pages), the written
                 findings, and verbatim guideline extracts
scripts/         Analysis pipeline (numbered, run in order)
results/         Deliverable tables (TSV) + a single combined Excel workbook
figures/         Static figures (PNG and SVG), each with a TSV of the values plotted
deliverables/    Dated snapshots of what was sent to Kevin, and what went by email instead
analyses/v3/     The same pipeline with the 528 changes new in CancerHotspots v3 added
```

### Parallel analyses

`scripts/paths.py` decides where every script reads and writes. By default that is
`data/interim/`, `results/` and `figures/`, exactly as before. With `ANALYSIS=v3`, the whole
pipeline runs on CancerHotspots v2 + v3 and writes to `analyses/v3/` instead, leaving the v2
analysis untouched. See [analyses/v3/README.md](analyses/v3/README.md) for the results.

## Source data

| File | Publication | Sheets | Content |
|---|---|---|---|
| `data/raw/cancerhotspots_v1_chang2016.xls` | Chang et al., *Nat Biotechnol* 2016 | `Per Residue`, `Per Allele` | 459 hotspot residues (1,170 residue/change rows) |
| `data/raw/cancerhotspots_v2_chang2018.xls` | Chang et al., *Cancer Discov* 2018 | `SNV-hotspots`, `INDEL-hotspots` | 1,110 SNV residues (incl. 86 splice) + 55 indel regions; 24,592 tumours |
| `data/raw/hotspots_v3.xlsx` | CancerHotspots v3 export | `SNV_Variants`, `INDEL_Variants`, `Hotspot_Residues` | The v2 changes plus 528 changes new in v3, and per-cohort counts (MSK / non-MSK GENIE / TCGA) for the 164 new residues. Used only by the v3 analysis |
| `data/raw/svig_uk_canonical_variants.tsv` | SVIG-UK Supplementary Table 3 | — | 158 canonical (O1) variants, reproduced for reproducibility of script 06 |

CancerHotspots data are made available under the ODC Open Database License (ODbL); see
<https://www.cancerhotspots.org/>.

### Not held in this repository

Two inputs are used by the analysis but are **not redistributable**, and are gitignored:

| | Why | How to obtain |
|---|---|---|
| **COSMIC Cancer Mutation Census v104** and **everything derived from it** — the per-substitution counts (`data/interim/cosmic_cmc_hotspot_genes.tsv`), `results/11`–`14` and `results/cosmic_o4_analysis.xlsx` | COSMIC's academic licence permits use but not redistribution | Download CMC under your own licence and run scripts 04–06, which regenerate every COSMIC table locally. The author holds the files and can supply them on request |
| **SVIG-UK v1.1 guideline PDFs** (`docs/svig-uk/`) | ACGS documents; not republished here | Download from the ACGS best-practice guidelines page. Every passage the analysis depends on is quoted verbatim, with page references, in [docs/svig-uk_o4_o7_extracts.md](docs/svig-uk_o4_o7_extracts.md) |

| **AACR Project GENIE v20** (`data_mutations_extended.txt`, `data_clinical_sample.txt`) and **everything derived from them** — the per-change counts, `results/18`–`28` and `results/genie_o4_analysis.xlsx` | Released under a data-use agreement that forbids redistribution | Register at Synapse (`syn7222066`), accept the terms, point `GENIE_DIR` at the release and run scripts 09–10, which regenerate every GENIE table locally |

Everything else needed to reproduce scripts 01–03 is in the repository; scripts 05–06 also
need COSMIC, and scripts 09–10 need GENIE.

## Running the analysis

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

.venv/bin/python scripts/01_merge_hotspots.py               # -> data/interim/
.venv/bin/python scripts/02_position_analysis.py            # -> results/01..10
.venv/bin/python scripts/04_extract_cosmic_cmc.py           # -> data/interim/  (needs COSMIC)
.venv/bin/python scripts/05_o4_proxy_crosstab.py            # -> results/11..14 (gitignored)
.venv/bin/python scripts/06_canonical_list_overlap.py       # -> results/15..16
.venv/bin/python scripts/07_validate_against_api_export.py  # -> results/17     (needs the API export)
.venv/bin/python scripts/09_extract_genie.py                # -> results/18, 25 (needs GENIE; gitignored)
.venv/bin/python scripts/10_genie_o4_o7.py                  # -> results/19..28 (gitignored)
.venv/bin/python scripts/12_genie_points_per_change.py     # -> results/29..32 (gitignored)
.venv/bin/python scripts/11_figures.py                      # -> figures/ (GENIE figures need 09-10)
.venv/bin/python scripts/03_export_workbook.py              # -> results/cancerhotspots_o7_analysis.xlsx
```

`scripts/08_build_pages.py` regenerates `docs/index.html` from `docs/report.html`; run it
after editing the report. The report is authored as a fragment, and the build step wraps it
in the document skeleton that GitHub Pages needs but the Artifact publisher supplies itself.

Scripts 04, 07 and 09 read large reference files held outside the repository; their
locations can be overridden with environment variables:

| Script | Reference file | Override |
|---|---|---|
| `04_extract_cosmic_cmc.py` | COSMIC CMC v104 (`CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar`) | `COSMIC_CMC_TAR` |
| `07_validate_against_api_export.py` | cancerhotspots.org API export (`cancerhotspots_counts.json`) | `CANCERHOTSPOTS_JSON` |
| `09_extract_genie.py` | AACR Project GENIE release directory | `GENIE_DIR` (default `data/external/genie_v20/`, gitignored) |

## Deliverables (`results/`)

| File | Contents |
|---|---|
| `01_headline_summary.tsv` | Headline metrics for the whole analysis |
| `02_distribution_unique_changes.tsv` | Distribution of hotspots by 1 / 2 / 3 / 4+ distinct amino-acid changes |
| `02b_distribution_unique_changes_exact.tsv` | Same, unbanded |
| `03_gene_position_table.tsv` | **Complete gene × residue table**: substitutions, unique changes, count per change, dominance, MSK split |
| `04_per_allele_o7_tiers.tsv` | **Every amino-acid change with the O7 strength and points it achieves under SVIG-UK v1.1**, carrying its residue's positional columns (distinct changes, the full substitution breakdown, dominance of the commonest change, mutation-specific / mixed / position-specific character) and the MSK split. Self-contained: the Section 5 analysis and the O7 weighting in one table |
| `05_o4_o7_double_counting.tsv` | Per-change independence assessment and cap recommendation |
| `06_msk_overlap_summary.tsv` | MSK-IMPACT contribution to the hotspot evidence base |
| `06b_msk_fraction_distribution.tsv` | Residues banded by MSK fraction |
| `07_cap_rule_sensitivity.tsv` | Sensitivity of the proposed cap rule to threshold choice |
| `08_version_comparison.tsv` | CancerHotspots v1 vs v2 residue overlap |
| `09_worked_examples.tsv` | Familiar variants worked through the proposed rule |
| `10_haematology_coverage_check.tsv` | Coverage of well-known haematological hotspots |
| `15_canonical_list_overlap.tsv` | Every change, flagged against the SVIG-UK Canonical Variants List (O1) (COSMIC columns kept local) |
| `16_variants_materially_affected.tsv` | The changes that drop from +8 to +4, and whether O1 already protects them |
| `17_api_export_validation.tsv` | Workbook vs live cancerhotspots.org API export (100% agreement, 3,004/3,004) |

Also provided as a single combined Excel workbook,
`results/cancerhotspots_o7_analysis.xlsx`.

### COSMIC tables — generated locally, not in the repository

COSMIC's licence permits use but not redistribution, and these tables are derived from it,
so they are gitignored. Anyone with a COSMIC licence regenerates them with scripts 04–05.
Tables 15 and 16 above stay public because they are written without COSMIC columns.

| File | Contents |
|---|---|
| `11_cosmic_o4_o7_crosstab.tsv` | O4 (COSMIC proxy) × O7 cross-tabulation, two threshold schemes |
| `12_o4_o7_correlation.tsv` | Correlation between the CancerHotspots and COSMIC recurrence counts |
| `13_points_impact_of_cap.tsv` | What the proposed cap costs, in changes affected and points lost |
| `14_per_change_o4_o7_points.tsv` | Per-change O4 and O7 points, capped and uncapped |
| `cosmic_o4_analysis.xlsx` | All of the above as one workbook |

### GENIE tables — generated locally, not in the repository

AACR Project GENIE's data-use agreement forbids redistribution, and these tables are
derived from it, so they are gitignored. Anyone with GENIE access regenerates all of them
with `scripts/09_extract_genie.py` and `scripts/10_genie_o4_o7.py`. The findings and the
report quote summary figures only.

| File | Contents |
|---|---|
| `18_genie_cohort_composition.tsv` | **GENIE v20 by lineage (myeloid / lymphoid / solid) and centre**, samples and patients |
| `18b_genie_haem_cancer_types.tsv` | The haematological cancer types in GENIE, in detail |
| `18c_genie_isoform_offsets.tsv` | Genes GENIE annotates on a different isoform, and the residue offset applied |
| `19_genie_headline_summary.tsv` | Headline metrics for the GENIE analysis |
| `20_genie_o4_o7_crosstab.tsv` | O4 (GENIE, unique patients) × O7, with and without MSK-IMPACT |
| `21_genie_points_impact_of_cap.tsv` | What the cap costs with GENIE as the O4 source |
| `22_genie_variants_dropping_8_to_4.tsv` | The changes that drop from +8 to +4, with their GENIE counts |
| `23_genie_leave_msk_out_transitions.tsv` | **The leave-MSK-out recount**: O4 tier on all patients vs non-MSK patients, per O7 tier |
| `24_genie_correlations.tsv` | CancerHotspots vs GENIE recurrence, including the cohort-disjoint comparison |
| `25_genie_haem_named_variants.tsv` | The haem drivers CancerHotspots lacks, counted in GENIE by lineage |
| `26_genie_per_change_o4_o7_points.tsv` | **Per change: GENIE samples, patients, non-MSK and per-lineage patients, O4 tier, capped and uncapped points** |
| `27_genie_vs_cosmic_o4_tiers.tsv` | Agreement between the GENIE and COSMIC O4 tiers |
| `28_genie_on_target_lineage_o4.tsv` | O4 counted on solid, haematological, myeloid or lymphoid patients only |
| `genie_o4_analysis.xlsx` | Tables 18–28 as one workbook |
| `29_genie_points_per_change.tsv` | **Points per change on GENIE**: O7, O4 (unique patients), and the O4 / O7 handling with resulting points under both the permissive and strict criteria |
| `30_genie_cap_impact_by_criterion.tsv` | What the cap costs under each criterion |
| `31_genie_o1_canonical_impact.tsv` | Changes on the O1 canonical list, with both outcomes |
| `32_genie_materially_affected.tsv` | Changes dropping from +8 to +4 under either criterion |
| `genie_points_per_change.xlsx` | Tables 02, 02b, 03 and 29–32 as one workbook |

## Correction (11 September 2026)

The first version of this analysis scored **synonymous** changes for O7. CancerHotspots lists
them at hotspot residues (KRAS p.G60G, 19 mutations), and 87 of the 207 had been given O7 at
supporting or moderate. O7 is reserved for missense and small in-frame variants, so they are
now excluded, as nonsense and stop-loss changes already were. That changed the O7 tier counts
(moderate 205 → 203, supporting 1,333 → 1,248) and the denominators of the independence test
(1,736 → 1,649 O7-scoring changes; 897 rather than 963 capped). It changed **none** of the 198
O7_strong changes, none of the cap-impact figures, and none of the ten changes that drop from
+8 to +4. `results/04_per_allele_o7_tiers.tsv` carries the corrected tiers.

## Limitations

Stated in full in [docs/findings.md § 9](docs/findings.md). GENIE is no longer one of them: sections 4–5 of the findings model O4 on COSMIC, and
section 11 repeats them on GENIE v20, where every conclusion holds. The GENIE counts are
pan-cancer; O4 restricted to a single tumour type will be lower, and per-lineage counts are
provided for laboratories that apply the guideline's "on-target" option.

## Acknowledgement (AACR Project GENIE)

This analysis uses AACR Project GENIE v20.0-public. The terms of GENIE access require that
any resulting publication or presentation cite AACR
Project GENIE Consortium, *Cancer Discov* 2017;7(8):818–31 with the dataset version used,
carry the acknowledgement *"The authors would like to acknowledge the American Association
for Cancer Research and its financial and material support in the development of the AACR
Project GENIE registry, as well as members of the consortium for their commitment to data
sharing. Interpretations are the responsibility of the study authors."*, and display the
AACR Project GENIE logo on all posters and presentations.
