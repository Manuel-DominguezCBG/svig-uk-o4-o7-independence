# Are O4 and O7 independent lines of evidence?

`svig-uk-o4-o7-independence`

Analysis supporting a proposed amendment to the ACGS/SVIG-UK guidelines for the
classification of oncogenicity of somatic variants (v1.1, 26/05/2026), concerning
the combined use of **O4** (enrichment in a somatic variant database, e.g. GENIE)
and **O7** (mutational hotspot, e.g. cancerhotspots.org).

**Question addressed.** GENIE case counts support O4 and Cancer Hotspots supports
O7, but both draw a substantial proportion of their evidence from the same
underlying MSK cohort. Are the two lines of evidence independent, and should the
combined O4 + O7 contribution be capped?

Requested by Kevin Baker (Principal Clinical Scientist, Oncology; Wessex Genomics
Laboratory Service, Salisbury). Analysis by Manuel Dominguez Becerra.

The written report for SVIG-UK is [docs/report.html](docs/report.html), published at
<https://claude.ai/code/artifact/a7011889-489d-4d9d-8bf0-258d0c4030be>.

---

## Repository layout

```
data/raw/        Source workbooks exactly as supplied (unmodified)
data/interim/    Merged, tidied per-allele table produced by script 01
docs/svig-uk/    SVIG-UK v1.1 main guideline + supplementary material (PDF and extracted text)
docs/            Written findings and verbatim guideline extracts
scripts/         Analysis pipeline (numbered, run in order)
results/         Deliverable tables (TSV) + a single combined Excel workbook
```

## Source data

| File | Publication | Sheets | Content |
|---|---|---|---|
| `data/raw/cancerhotspots_v1_chang2016.xls` | Chang et al., *Nat Biotechnol* 2016 | `Per Residue`, `Per Allele` | 459 hotspot residues (1,170 residue/change rows) |
| `data/raw/cancerhotspots_v2_chang2018.xls` | Chang et al., *Cancer Discov* 2018 | `SNV-hotspots`, `INDEL-hotspots` | 1,110 SNV positions (incl. 86 splice) + 55 indel regions; 24,592 tumours |
| `data/raw/svig_uk_canonical_variants.tsv` | SVIG-UK Supplementary Table 3 | — | 158 canonical (O1) variants, from the local reference set |

CancerHotspots data are made available under the ODC Open Database License (ODbL);
see <https://www.cancerhotspots.org/>.

The v2 workbook carries `n_MSK` and `n_Retro` columns, which split every hotspot
position's mutation count into the contribution from **MSK-IMPACT** (a GENIE
contributing centre) and from **retrospective/public cohorts** (largely TCGA, which
is *not* in GENIE). This split is what makes the O4/O7 overlap directly measurable,
and it is the basis of the central result.

> **Note.** Two workbooks (four sheets) were supplied; the covering email referred
> to three files. If a third CancerHotspots export exists (e.g. the 3D hotspots
> table), it has not been included here.

## Running the analysis

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

.venv/bin/python scripts/01_merge_hotspots.py             # -> data/interim/
.venv/bin/python scripts/02_position_analysis.py          # -> results/01..10
.venv/bin/python scripts/04_extract_cosmic_cmc.py         # -> data/interim/ (needs COSMIC, see below)
.venv/bin/python scripts/05_o4_proxy_crosstab.py          # -> results/11..14
.venv/bin/python scripts/06_canonical_list_overlap.py     # -> results/15..16
.venv/bin/python scripts/07_validate_against_api_export.py  # -> results/17 (needs the API export)
.venv/bin/python scripts/03_export_workbook.py            # -> results/cancerhotspots_o7_analysis.xlsx
```

Scripts 01–03 and 06 need only what is in this repository. Scripts 04 and 07 read two
large reference files held outside it; their locations can be overridden with
environment variables:

| Script | Reference file | Override |
|---|---|---|
| `04_extract_cosmic_cmc.py` | COSMIC Cancer Mutation Census v104 (`CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar`) | `COSMIC_CMC_TAR` |
| `07_validate_against_api_export.py` | cancerhotspots.org API export (`cancerhotspots_counts.json`) | `CANCERHOTSPOTS_JSON` |

Both default to the paths in the local reference set at
`/Users/monkiky/Documents/external/refs/`. COSMIC requires a licence and is not
redistributed here; only the derived per-substitution counts for the 230 hotspot
genes are written to `data/interim/`.

Script 01 asserts that, for every v2 SNV position, the per-allele counts sum to the
stated position total and that `n_MSK + n_Retro` reconstitutes that total. The
analysis aborts if either check fails.

## Deliverables (`results/`)

| File | Contents |
|---|---|
| `01_headline_summary.tsv` | Headline metrics for the whole analysis |
| `02_distribution_unique_changes.tsv` | Distribution of hotspots by 1 / 2 / 3 / 4+ distinct amino-acid changes |
| `02b_distribution_unique_changes_exact.tsv` | Same, unbanded |
| `03_gene_position_table.tsv` | **Complete gene x position table**: substitutions, unique changes, count per change, dominance, MSK split |
| `04_per_allele_o7_tiers.tsv` | Every amino-acid change with the O7 strength it achieves under SVIG-UK v1.1 |
| `05_o4_o7_double_counting.tsv` | Per-change independence assessment and cap recommendation |
| `06_msk_overlap_summary.tsv` | MSK-IMPACT contribution to the hotspot evidence base |
| `06b_msk_fraction_distribution.tsv` | Positions banded by MSK fraction |
| `07_cap_rule_sensitivity.tsv` | Sensitivity of the proposed cap rule to threshold choice |
| `08_version_comparison.tsv` | CancerHotspots v1 vs v2 position overlap |
| `09_worked_examples.tsv` | Familiar variants worked through the proposed rule |
| `10_haematology_coverage_check.tsv` | Coverage of well-known haematological hotspots |
| `11_cosmic_o4_o7_crosstab.tsv` | O4 (COSMIC proxy) x O7 cross-tabulation, two threshold schemes |
| `12_o4_o7_correlation.tsv` | Correlation between the CancerHotspots and COSMIC recurrence counts |
| `13_points_impact_of_cap.tsv` | What the proposed cap costs, in changes affected and points lost |
| `14_per_change_o4_o7_points.tsv` | Per-change O4 and O7 points, capped and uncapped |
| `15_canonical_list_overlap.tsv` | Every change, flagged against the SVIG-UK Canonical Variants List (O1) |
| `16_variants_materially_affected.tsv` | The changes that drop from +8 to +4, and whether O1 already protects them |
| `17_api_export_validation.tsv` | Workbook vs live cancerhotspots.org API export (100% agreement, 3,004/3,004) |

The primary analysis set is the v2 SNV table restricted to protein-coding
(non-splice) positions: 1,024 positions and 2,918 distinct amino-acid changes.
Splice hotspots are excluded because O7 is not applicable to them (SVIG-UK
Supplementary Figure 1, note 22), and within the analysis set a further 247 nonsense
or stop-loss changes are marked ineligible for O7 scoring (notes 4 and 13), leaving
1,736 changes that score O7 at some strength from cancerhotspots.org.

GENIE was not available. COSMIC Cancer Mutation Census v104 is used as an independent
stand-in for the O4 side, which the guideline explicitly permits; see
`docs/findings.md` section 4 for what that does and does not support.

## Findings

See **[docs/findings.md](docs/findings.md)** for the written answer to Kevin's
questions, and **[docs/svig-uk_o4_o7_extracts.md](docs/svig-uk_o4_o7_extracts.md)**
for the verbatim guideline text the analysis is anchored to.
