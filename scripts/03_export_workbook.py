#!/usr/bin/env python3
"""
03_export_workbook.py
---------------------
Bundle the TSV deliverables into a single Excel workbook for circulation.
The TSVs remain the canonical output; this is a convenience view.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "cancerhotspots_o7_analysis.xlsx"

# Excel sheet names are capped at 31 characters.
SHEET_NAMES = {
    "01_headline_summary": "01 Headline summary",
    "02_distribution_unique_changes": "02 Changes per hotspot",
    "02b_distribution_unique_changes_exact": "02b Changes per hotspot exact",
    "03_gene_position_table": "03 Gene x position",
    "04_per_allele_o7_tiers": "04 O7 tiers per change",
    "05_o4_o7_double_counting": "05 O4-O7 independence",
    "06_msk_overlap_summary": "06 MSK overlap",
    "06b_msk_fraction_distribution": "06b MSK fraction bands",
    "07_cap_rule_sensitivity": "07 Cap rule sensitivity",
    "08_version_comparison": "08 v1 vs v2",
    "09_worked_examples": "09 Worked examples",
    "10_haematology_coverage_check": "10 Haem coverage",
    "11_cosmic_o4_o7_crosstab": "11 O4 x O7 crosstab",
    "12_o4_o7_correlation": "12 O4-O7 correlation",
    "13_points_impact_of_cap": "13 Points impact of cap",
    "14_per_change_o4_o7_points": "14 Points per change",
    "15_canonical_list_overlap": "15 O1 canonical overlap",
    "16_variants_materially_affected": "16 Materially affected",
    "17_api_export_validation": "17 API export validation",
}

# GENIE-derived tables go to a separate workbook that is never committed: the
# GENIE data-use agreement forbids redistribution (see .gitignore).
GENIE_OUT = RESULTS / "genie_o4_analysis.xlsx"
GENIE_SHEET_NAMES = {
    "18_genie_cohort_composition": "18 GENIE cohort composition",
    "18b_genie_haem_cancer_types": "18b GENIE haem cancer types",
    "18c_genie_isoform_offsets": "18c GENIE isoform offsets",
    "19_genie_headline_summary": "19 GENIE headline",
    "20_genie_o4_o7_crosstab": "20 GENIE O4 x O7",
    "21_genie_points_impact_of_cap": "21 GENIE cap impact",
    "22_genie_variants_dropping_8_to_4": "22 GENIE 8 to 4",
    "23_genie_leave_msk_out_transitions": "23 GENIE leave-MSK-out",
    "24_genie_correlations": "24 GENIE correlations",
    "25_genie_haem_named_variants": "25 GENIE haem variants",
    "26_genie_per_change_o4_o7_points": "26 GENIE points per change",
    "27_genie_vs_cosmic_o4_tiers": "27 GENIE vs COSMIC O4",
    "28_genie_on_target_lineage_o4": "28 GENIE on-target lineage O4",
}


def write_workbook(out, sheet_names):
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for stem, sheet in sheet_names.items():
            path = RESULTS / f"{stem}.tsv"
            if not path.exists():
                print(f"  skipped (missing): {path.name}")
                continue
            df = pd.read_csv(path, sep="\t")
            df.to_excel(writer, sheet_name=sheet, index=False)
            # Rough auto-width so the workbook is readable without fiddling.
            ws = writer.sheets[sheet]
            for i, col in enumerate(df.columns, start=1):
                width = max(len(str(col)), df[col].astype(str).str.len().max() if len(df) else 0)
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(width + 2, 60)
            ws.freeze_panes = "A2"
    print(f"wrote {out}")


def main():
    write_workbook(OUT, SHEET_NAMES)
    if any((RESULTS / f"{stem}.tsv").exists() for stem in GENIE_SHEET_NAMES):
        write_workbook(GENIE_OUT, GENIE_SHEET_NAMES)


if __name__ == "__main__":
    main()
