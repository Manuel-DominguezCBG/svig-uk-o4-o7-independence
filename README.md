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
   splice residues (O7 not applicable, note 22) and nonsense/stop-loss changes (O2 + O7 not
   permitted, notes 4 and 13).
5. **Tested the coupling against a second database.** COSMIC Cancer Mutation Census v104 as
   an independent stand-in for the O4 side, giving the correlation and the O4 × O7
   cross-tabulation under two threshold schemes.
6. **Modelled the proposed cap** and its sensitivity to how "independent evidence" is
   defined, then bounded its clinical impact against the SVIG-UK Canonical Variants List.
7. **Validated the whole thing against the live resource.** The parsed workbook was
   cross-checked against a cancerhotspots.org API export: **all 3,004 amino-acid changes
   match, 100% agreement on both per-change counts and residue totals**. The numbers
   modelled here are the numbers an analyst sees on the website.

The primary analysis set is the v2 SNV table restricted to protein-coding (non-splice)
residues: **1,024 residues, 2,918 distinct amino-acid changes**, of which 1,736 score O7 at
some strength.

---

## What could be done next, and what data it would need

Full detail in [docs/findings.md § 10](docs/findings.md). None of it is required to act on
the recommendation — these would strengthen it or close gaps it leaves open.

| Analysis | What it would settle | Data needed | Available? |
|---|---|---|---|
| **GENIE leave-MSK-out recount** | How many variants actually change O4 tier when the shared patients are removed — converts a correlation argument into a count | AACR Project GENIE release (`data_mutations_extended.txt`, `data_clinical_sample.txt`, which carry `CENTER` per sample) | Free via Synapse `syn7222066`; needs registration + data-use agreement. ~1 day's work once obtained |
| Would the residue still be a hotspot without MSK? | Whether a position was called *only* because of MSK cases | Approximation possible **now** from the published `n_Retro` column; full re-derivation needs the per-sample input to Chang 2018 + its background rate model | Partly in hand; the full version is a large piece of work |
| Haematological hotspot supplement | Restores O7 for JAK2 p.V617F, MPL p.W515, CALR, NPM1, ASXL1, SETBP1 — **none of which appear in CancerHotspots v2 at all** | COSMIC restricted to haematopoietic/lymphoid tissue, the GENIE haem subset, or BeatAML / MDS-CHIP series | COSMIC already licensed locally |
| Per-database O4 threshold calibration | What COSMIC count is equivalent to GENIE's ≥10 and ≥50 — the guideline says only "much higher" | GENIE + COSMIC joined per substitution | Falls out cheaply once GENIE is loaded |
| Indel hotspot regions | Whether the same double counting applies to the 55 in-frame indel regions, excluded here because they are ranges not residues | Per-sample indel calls, plus an SVIG-UK decision on what independence means across a range | Needs a definitional call as much as data |
| Transcript harmonisation | Recovers the 11.9% of changes that failed to join to COSMIC on transcript choice (MYD88 p.L265P is p.L273P in COSMIC) | MANE Select + VEP or Mutalyzer | Self-contained; the easiest item on the list |
| **Impact on real reported cases** | How many previously classified variants change SVIG-UK class under the cap — a different order of evidence from modelling over the resource | De-identified retrospective extract of variants classified under v1.1 at Wessex GLH, with **per-code point assignments**, not just final classes | Held in-service; needs local IG sign-off, not an external application |
| Currency refresh | Whether tier assignments hold against a modern cohort (v2 is 2018, 24,592 tumours) | A future cancerhotspots.org release, or GENIE | Re-running the pipeline against it is one command |

---

## Repository layout

```
data/raw/        Source workbooks exactly as supplied (unmodified)
data/interim/    Merged, tidied per-allele table produced by script 01
docs/            The report (report.html, built to index.html for Pages), the written
                 findings, and verbatim guideline extracts
scripts/         Analysis pipeline (numbered, run in order)
results/         Deliverable tables (TSV) + a single combined Excel workbook
```

## Source data

| File | Publication | Sheets | Content |
|---|---|---|---|
| `data/raw/cancerhotspots_v1_chang2016.xls` | Chang et al., *Nat Biotechnol* 2016 | `Per Residue`, `Per Allele` | 459 hotspot residues (1,170 residue/change rows) |
| `data/raw/cancerhotspots_v2_chang2018.xls` | Chang et al., *Cancer Discov* 2018 | `SNV-hotspots`, `INDEL-hotspots` | 1,110 SNV residues (incl. 86 splice) + 55 indel regions; 24,592 tumours |
| `data/raw/svig_uk_canonical_variants.tsv` | SVIG-UK Supplementary Table 3 | — | 158 canonical (O1) variants, reproduced for reproducibility of script 06 |

CancerHotspots data are made available under the ODC Open Database License (ODbL); see
<https://www.cancerhotspots.org/>.

### Not held in this repository

Two inputs are used by the analysis but are **not redistributable**, and are gitignored:

| | Why | How to obtain |
|---|---|---|
| **COSMIC Cancer Mutation Census v104** and the derived per-substitution counts (`data/interim/cosmic_cmc_hotspot_genes.tsv`) | COSMIC's academic licence permits use but not redistribution | Download CMC under your own licence and run `scripts/04_extract_cosmic_cmc.py`, which regenerates the derived table. The author holds the file and can supply it directly on request |
| **SVIG-UK v1.1 guideline PDFs** (`docs/svig-uk/`) | ACGS documents; not republished here | Download from the ACGS best-practice guidelines page. Every passage the analysis depends on is quoted verbatim, with page references, in [docs/svig-uk_o4_o7_extracts.md](docs/svig-uk_o4_o7_extracts.md) |

Everything else needed to reproduce scripts 01–03 and 06 is in the repository.

## Running the analysis

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

.venv/bin/python scripts/01_merge_hotspots.py               # -> data/interim/
.venv/bin/python scripts/02_position_analysis.py            # -> results/01..10
.venv/bin/python scripts/04_extract_cosmic_cmc.py           # -> data/interim/  (needs COSMIC)
.venv/bin/python scripts/05_o4_proxy_crosstab.py            # -> results/11..14
.venv/bin/python scripts/06_canonical_list_overlap.py       # -> results/15..16
.venv/bin/python scripts/07_validate_against_api_export.py  # -> results/17     (needs the API export)
.venv/bin/python scripts/03_export_workbook.py              # -> results/cancerhotspots_o7_analysis.xlsx
```

`scripts/08_build_pages.py` regenerates `docs/index.html` from `docs/report.html`; run it
after editing the report. The report is authored as a fragment, and the build step wraps it
in the document skeleton that GitHub Pages needs but the Artifact publisher supplies itself.

Scripts 04 and 07 read two large reference files held outside the repository; their
locations can be overridden with environment variables:

| Script | Reference file | Override |
|---|---|---|
| `04_extract_cosmic_cmc.py` | COSMIC CMC v104 (`CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar`) | `COSMIC_CMC_TAR` |
| `07_validate_against_api_export.py` | cancerhotspots.org API export (`cancerhotspots_counts.json`) | `CANCERHOTSPOTS_JSON` |

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
| `11_cosmic_o4_o7_crosstab.tsv` | O4 (COSMIC proxy) × O7 cross-tabulation, two threshold schemes |
| `12_o4_o7_correlation.tsv` | Correlation between the CancerHotspots and COSMIC recurrence counts |
| `13_points_impact_of_cap.tsv` | What the proposed cap costs, in changes affected and points lost |
| `14_per_change_o4_o7_points.tsv` | Per-change O4 and O7 points, capped and uncapped |
| `15_canonical_list_overlap.tsv` | Every change, flagged against the SVIG-UK Canonical Variants List (O1) |
| `16_variants_materially_affected.tsv` | The changes that drop from +8 to +4, and whether O1 already protects them |
| `17_api_export_validation.tsv` | Workbook vs live cancerhotspots.org API export (100% agreement, 3,004/3,004) |

Also provided as a single combined Excel workbook,
`results/cancerhotspots_o7_analysis.xlsx`.

## Limitations

Stated in full in [docs/findings.md § 9](docs/findings.md). The most important: **GENIE was
not available**, so the O4 side of the cross-tabulation is modelled on COSMIC. That is
evidence about how tightly O4-style and O7-style recurrence track each other, *not* a
prediction of the exact GENIE count any individual variant would return.

## Acknowledgement required if extended to GENIE

The terms of GENIE access require that any resulting publication or presentation cite AACR
Project GENIE Consortium, *Cancer Discov* 2017;7(8):818–31 with the dataset version used,
carry the acknowledgement *"The authors would like to acknowledge the American Association
for Cancer Research and its financial and material support in the development of the AACR
Project GENIE registry, as well as members of the consortium for their commitment to data
sharing. Interpretations are the responsibility of the study authors."*, and display the
AACR Project GENIE logo on all posters and presentations.
