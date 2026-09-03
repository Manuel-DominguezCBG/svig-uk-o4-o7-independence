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
}


def main():
    with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
        for stem, sheet in SHEET_NAMES.items():
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
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
