#!/usr/bin/env python3
"""
12_genie_points_per_change.py
-----------------------------
The per-change points table ("14 Points per change") rebuilt on GENIE, with the
proposed O4 + O7 handling applied under BOTH calibrations of the leave-one-variant-out
test (report Section 9):

  permissive  the residue retains >= 10 mutations after removing the change under
              assessment, and at least one other change there has >= 2 mutations;
  strict      as above, but the other change must have >= 10 mutations.

Where the residue passes the test, O4 and O7 combine in full (max +8); where it fails,
the combined O4 + O7 contribution is capped at +4.

O4 is Supplementary Table 1 applied to unique GENIE patients (script 10); O7 is
cancerhotspots.org as scored in script 02.

Everything here is derived from GENIE at the level of individual variants, so every
output is gitignored and goes to Kevin Baker by email, not through the repository.

Inputs:
  results/04_per_allele_o7_tiers.tsv          positional data per change
  results/05_o4_o7_double_counting.tsv        independence test, both calibrations
  results/26_genie_per_change_o4_o7_points.tsv  GENIE counts and O4 tier (script 10)
  results/02_distribution_unique_changes.tsv, 02b_..._exact.tsv, 03_gene_position_table.tsv

Outputs (all local only):
  results/29_genie_points_per_change.tsv      one row per change, both calibrations
  results/30_genie_cap_impact_by_criterion.tsv
  results/31_genie_o1_canonical_impact.tsv    changes on the O1 list
  results/32_genie_materially_affected.tsv    changes dropping from +8 to +4, either test
  results/genie_points_per_change.xlsx        all of the above, with tables 02, 02b, 03
"""

from pathlib import Path

import pandas as pd

import paths

ROOT = Path(__file__).resolve().parents[1]
RESULTS = paths.RESULTS
KEYS = ["hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa"]
VERSION = ["hotspot_version"] if paths.ANALYSIS == "v3" else []
CRITERIA = {"permissive": "", "strict": "_strict"}
COMBINE = "O4 + O7 may be combined (max +8)"
CAP = "Cap combined O4 + O7 at +4"


def truthy(s):
    return s.astype(str).eq("True")


def load():
    dtype = {k: str for k in KEYS}
    pos = pd.read_csv(RESULTS / "04_per_allele_o7_tiers.tsv", sep="\t", dtype=dtype,
                      usecols=KEYS + VERSION + ["same_change_fraction", "substitutions", "top_change",
                                      "top_change_count", "top_change_fraction"])
    test = pd.read_csv(RESULTS / "05_o4_o7_double_counting.tsv", sep="\t", dtype=dtype,
                       usecols=KEYS + ["residual_count", "residual_n_changes",
                                       "residual_max_change_count",
                                       "independent_positional_evidence",
                                       "independent_positional_evidence_strict"])
    genie = pd.read_csv(RESULTS / "26_genie_per_change_o4_o7_points.tsv", sep="\t", dtype=dtype)
    df = (genie.drop(columns=["independent_positional_evidence", "recommended_o4_o7_handling"])
          .merge(test, on=KEYS, how="left", validate="one_to_one")
          .merge(pos, on=KEYS, how="left", validate="one_to_one"))
    assert df["residual_count"].notna().all() and df["substitutions"].notna().all()
    for c in ["gene_in_genie", "on_svig_uk_canonical_list",
              "independent_positional_evidence", "independent_positional_evidence_strict"]:
        df[c] = truthy(df[c])
    return df


def score(df):
    df["o4_o7_uncapped"] = df["o4_o7_uncapped__genie_patients"]
    for name, suffix in CRITERIA.items():
        independent = df[f"independent_positional_evidence{suffix}"]
        df[f"residue_independent__{name}"] = independent
        df[f"o4_o7_handling__{name}"] = independent.map({True: COMBINE, False: CAP})
        df[f"o4_o7_points__{name}"] = df["o4_o7_uncapped"].where(
            independent, df["o4_o7_uncapped"].clip(upper=4))
        df[f"points_lost__{name}"] = df["o4_o7_uncapped"] - df[f"o4_o7_points__{name}"]
    # The handling only changes a score where both codes are applied.
    both = (df["o4_points__genie_patients"] > 0) & (df["o7_points"] > 0)
    for name in CRITERIA:
        df.loc[~both, f"o4_o7_handling__{name}"] = (
            "Not applicable (O4 and O7 not both applied)")
    return df


PER_CHANGE = KEYS + VERSION + [
    # CancerHotspots positional data
    "substitutions", "n_unique_changes", "change_count", "position_total_count",
    "same_change_fraction", "top_change", "top_change_count", "top_change_fraction",
    "hotspot_character", "msk_fraction",
    # O7
    "o7_strength", "o7_points",
    # O4 on GENIE
    "gene_in_genie", "genie_samples", "genie_patients", "genie_patients_non_msk",
    "genie_patients_myeloid", "genie_patients_lymphoid", "genie_patients_solid",
    "o4_strength__genie_patients", "o4_points__genie_patients",
    # combined, and the proposed handling under both tests
    "o4_o7_uncapped",
    "residual_count", "residual_n_changes", "residual_max_change_count",
    "residue_independent__permissive", "o4_o7_handling__permissive",
    "o4_o7_points__permissive", "points_lost__permissive",
    "residue_independent__strict", "o4_o7_handling__strict",
    "o4_o7_points__strict", "points_lost__strict",
    # O1
    "svig_uk_assessment", "on_svig_uk_canonical_list",
]


def impact(df):
    rows = []
    both = (df["o4_points__genie_patients"] > 0) & (df["o7_points"] > 0)
    for name in CRITERIA:
        lost = df[f"points_lost__{name}"]
        drop = (df["o4_o7_uncapped"] == 8) & (df[f"o4_o7_points__{name}"] == 4)
        o1 = df["on_svig_uk_canonical_list"]
        rows.append({
            "criterion": name,
            "changes_assessed": len(df),
            "changes_scoring_both_o4_and_o7": int(both.sum()),
            "changes_capped_where_both_apply": int((both & ~df[f"residue_independent__{name}"]).sum()),
            "changes_whose_score_falls": int((lost > 0).sum()),
            "pct_of_all_changes": round(100 * (lost > 0).mean(), 1),
            "mean_points_lost_where_it_falls": round(lost[lost > 0].mean(), 2),
            "changes_dropping_from_8_to_4": int(drop.sum()),
            "of_which_on_o1_canonical_list": int((drop & o1).sum()),
            "materially_affected_not_on_o1": int((drop & ~o1).sum()),
            "o1_listed_changes_whose_score_falls": int(((lost > 0) & o1).sum()),
        })
    return pd.DataFrame(rows)


def main():
    df = score(load())
    per_change = df[PER_CHANGE].sort_values(["o4_o7_uncapped", "change_count"],
                                            ascending=False)

    o1 = per_change[per_change["on_svig_uk_canonical_list"]].copy()
    o1["cap_moot_because_on_o1"] = (o1["points_lost__permissive"] > 0) | (o1["points_lost__strict"] > 0)

    drop_any = per_change[(per_change["o4_o7_uncapped"] == 8)
                          & ((per_change["o4_o7_points__permissive"] == 4)
                             | (per_change["o4_o7_points__strict"] == 4))].copy()
    drop_any["drops_8_to_4__permissive"] = drop_any["o4_o7_points__permissive"] == 4
    drop_any["drops_8_to_4__strict"] = drop_any["o4_o7_points__strict"] == 4
    drop_any["protected_by_o1"] = drop_any["on_svig_uk_canonical_list"]
    drop_cols = KEYS + VERSION + ["substitutions", "n_unique_changes", "change_count",
                        "position_total_count", "residual_count", "residual_max_change_count",
                        "genie_patients", "genie_patients_non_msk", "msk_fraction",
                        "drops_8_to_4__permissive", "drops_8_to_4__strict",
                        "svig_uk_assessment", "protected_by_o1"]
    drop_any = drop_any[drop_cols]

    summary = impact(df)

    per_change.to_csv(RESULTS / "29_genie_points_per_change.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "30_genie_cap_impact_by_criterion.tsv", sep="\t", index=False)
    o1.to_csv(RESULTS / "31_genie_o1_canonical_impact.tsv", sep="\t", index=False)
    drop_any.to_csv(RESULTS / "32_genie_materially_affected.tsv", sep="\t", index=False)

    about = pd.DataFrame({"": [
        "O4 / O7 points per change on GENIE v20 -- prepared for Kevin Baker, 16 September 2026.",
        "",
        "Sheets",
        "  02 Changes per hotspot / 02b exact -- distribution of hotspot residues by number of distinct changes.",
        "  03 Gene x position -- every residue: substitutions with their counts, unique changes, dominance, MSK split.",
        "  14 Points per change (GENIE) -- every amino-acid change: CancerHotspots counts, O7 (SVIG-UK v1.1),",
        "      O4 on unique GENIE patients (SVIG-UK v1.1), combined points uncapped, and the proposed O4 / O7",
        "      handling with the resulting points under the permissive AND the strict criterion.",
        "  Cap impact by criterion -- what each criterion costs, overall and on the O1 list.",
        "  O1 canonical impact -- every change on the SVIG-UK Canonical Variants List, with both outcomes.",
        "  Materially affected -- changes dropping from +8 to +4 under either criterion, and whether O1 protects them.",
        "",
        "Criteria (report Section 9, leave-one-variant-out test)",
        "  Remove the change under assessment from the residue total. The residue is independent evidence if the",
        "  residual is >= 10 mutations and at least one other change there is itself recurrent:",
        "    permissive: that other change has >= 2 mutations;  strict: >= 10 mutations.",
        "  Independent: O4 + O7 combine in full (max +8). Not independent: O4 + O7 capped at +4.",
        "",
        "O7: cancerhotspots.org v2, Supplementary Table 1. Nonsense, stop-loss and synonymous changes are not applicable.",
        "O4: Supplementary Table 1 thresholds on unique GENIE patients (multiple samples from one patient count once);",
        "    missense > 10 strong, 5-10 moderate; nonsense > 50 strong, 20-50 moderate, 10-19 supporting.",
        "",
        "Contains AACR Project GENIE v20.0-public data at the level of individual variants: not for redistribution.",
    ]})

    sheets = [
        ("About", about, False),
        ("02 Changes per hotspot", pd.read_csv(RESULTS / "02_distribution_unique_changes.tsv", sep="\t"), True),
        ("02b Changes per hotspot exact", pd.read_csv(RESULTS / "02b_distribution_unique_changes_exact.tsv", sep="\t"), True),
        ("03 Gene x position", pd.read_csv(RESULTS / "03_gene_position_table.tsv", sep="\t"), True),
        ("14 Points per change (GENIE)", per_change, True),
        ("Cap impact by criterion", summary, True),
        ("O1 canonical impact", o1, True),
        ("Materially affected", drop_any, True),
    ]
    out = RESULTS / "genie_points_per_change.xlsx"
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for name, frame, header in sheets:
            frame.to_excel(writer, sheet_name=name, index=False, header=header)
            ws = writer.sheets[name]
            if name == "About":
                ws.column_dimensions["A"].width = 120
                continue
            for i, col in enumerate(frame.columns, start=1):
                width = max(len(str(col)), frame[col].astype(str).str.len().max() if len(frame) else 0)
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(width + 2, 50)
            ws.freeze_panes = "E2" if "hugo_symbol" in frame.columns else "A2"
            ws.auto_filter.ref = ws.dimensions

    print(summary.T.to_string(header=False))
    print("\nMATERIALLY AFFECTED (either criterion)")
    print(drop_any.drop(columns=["substitutions"]).to_string(index=False))
    print(f"\nO1-listed changes in the analysis set: {len(o1)}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
