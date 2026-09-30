#!/usr/bin/env python3
"""
14_o7_weighting_options.py
--------------------------
Kevin Baker's two requests of 30 September 2026. Run with ANALYSIS=v3.

REQUEST 1 -- one expanded per-change table.
Every column of the two tables sent before (26_genie_per_change_o4_o7_points and
29_genie_points_per_change), for the CancerHotspots v2 + v3 analysis set, with the
SVIG-UK canonical (O1) overlap, a v2 / v3 comparison, and both the permissive and the
strict handling. Nothing is dropped; the earlier tables are left untouched.

REQUEST 2 -- how to weight O7 on changes the cap does not apply to.
Where the leave-one-variant-out test finds independent positional evidence, O4 and O7
combine, and the question is what O7 is worth. Two options:

  Option 1, total-count:    apply the SVIG-UK O7 thresholds as written -- the variant
                            under assessment's own change count against the full residue
                            total.
  Option 2, residual-count: apply the same thresholds to what is left once the VUA is
                            removed -- the largest OTHER change at the residue against
                            the residual residue total.

Both are scored for every uncapped change, each with the rule it passed or failed, so the
difference can be attributed to the same-change threshold, the residue-total threshold, or
neither.

O7 thresholds (Supplementary Table 1): same change >= 10 and residue total >= 50 -> Strong
(+4); same change >= 10 with residue total < 50 -> Moderate (+2); same change 2-9 ->
Supporting (+1); otherwise not met.

Inputs (from the analysis's own results/):
  29_genie_points_per_change.tsv, 26_genie_per_change_o4_o7_points.tsv,
  15_canonical_list_overlap.tsv, 05_o4_o7_double_counting.tsv

Outputs (all carry GENIE data per variant, so all are gitignored):
  36_genie_points_per_change_expanded.tsv
  37_genie_o7_weighting_options.tsv
  38_genie_o7_weighting_options_summary.tsv
  genie_o4_o7_expanded.xlsx
"""

import pandas as pd

import paths

RESULTS = paths.RESULTS
KEYS = ["hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa"]
O7_CHANGE_HIGH, O7_CHANGE_SUPPORTING, O7_POSITION_STRONG = 10, 2, 50
POINTS = {"Strong": 4, "Moderate": 2, "Supporting": 1, "Not met": 0}


def o7_rule(change_count, position_total, change_label, total_label):
    """The SVIG-UK O7 tier, with the threshold it turned on spelled out."""
    c, t = int(change_count), int(position_total)
    if c >= O7_CHANGE_HIGH:
        if t >= O7_POSITION_STRONG:
            return ("Strong", 4,
                    f"{change_label} {c} >= 10 AND {total_label} {t} >= 50 -> Strong (+4)")
        return ("Moderate", 2,
                f"{change_label} {c} >= 10 BUT {total_label} {t} < 50 -> Moderate (+2)")
    if c >= O7_CHANGE_SUPPORTING:
        return ("Supporting", 1,
                f"{change_label} {c} is 2-9 (< 10) -> Supporting (+1)")
    return ("Not met", 0, f"{change_label} {c} < 2 -> O7 not met (0)")


def difference_category(c_vua, t_full, c_other, t_residual):
    """Which thresholds the two options part company on, as a groupable label."""
    parts = []
    if (c_vua >= O7_CHANGE_HIGH) != (c_other >= O7_CHANGE_HIGH):
        parts.append("same-change threshold (>= 10)")
    if (t_full >= O7_POSITION_STRONG) != (t_residual >= O7_POSITION_STRONG):
        parts.append("residue-total threshold (>= 50)")
    if (c_vua >= O7_CHANGE_SUPPORTING) != (c_other >= O7_CHANGE_SUPPORTING):
        parts.append("supporting threshold (>= 2)")
    return " + ".join(parts) if parts else "no threshold differs"


def difference_reason(c_vua, t_full, c_other, t_residual):
    """The same, with the counts spelled out for auditing one change."""
    reasons = []
    if (c_vua >= O7_CHANGE_HIGH) != (c_other >= O7_CHANGE_HIGH):
        reasons.append(f"same-change threshold (>= 10): VUA {int(c_vua)} vs "
                       f"largest other change {int(c_other)}")
    if (t_full >= O7_POSITION_STRONG) != (t_residual >= O7_POSITION_STRONG):
        reasons.append(f"residue-total threshold (>= 50): total {int(t_full)} vs "
                       f"residual {int(t_residual)}")
    if (c_vua >= O7_CHANGE_SUPPORTING) != (c_other >= O7_CHANGE_SUPPORTING):
        reasons.append(f"supporting threshold (>= 2): VUA {int(c_vua)} vs "
                       f"largest other change {int(c_other)}")
    return "; ".join(reasons) if reasons else "no threshold differs - same O7 tier"


def truthy(s):
    return s.astype(str).eq("True")


def load():
    dtype = {k: str for k in KEYS}
    new = pd.read_csv(RESULTS / "29_genie_points_per_change.tsv", sep="\t", dtype=dtype)
    old = pd.read_csv(RESULTS / "26_genie_per_change_o4_o7_points.tsv", sep="\t", dtype=dtype)
    o1 = pd.read_csv(RESULTS / "15_canonical_list_overlap.tsv", sep="\t", dtype=dtype,
                     usecols=KEYS + ["transcript", "HGVSp_Short"])
    # Columns the earlier table had that the newer one does not.
    extra = [c for c in old.columns if c not in new.columns]
    df = (new.merge(old[KEYS + extra], on=KEYS, how="left", validate="one_to_one")
          .merge(o1, on=KEYS, how="left", validate="one_to_one"))
    assert len(df) == len(new), "merge changed the row count"
    for c in ["gene_in_genie", "on_svig_uk_canonical_list",
              "residue_independent__permissive", "residue_independent__strict"]:
        df[c] = truthy(df[c])
    return df


def add_version_comparison(df):
    """Per change and per residue, which CancerHotspots release it comes from."""
    if "hotspot_version" not in df:          # the v2-only analysis
        df["hotspot_version"] = "v2"
    df["in_cancerhotspots_v2"] = df["hotspot_version"] == "v2"
    df["in_cancerhotspots_v3"] = df["hotspot_version"] == "v3"
    residue = (df.groupby(["hugo_symbol", "amino_acid_position"])["hotspot_version"]
               .agg(lambda s: "+".join(sorted(set(s)))).rename("residue_hotspot_versions"))
    return df.merge(residue, on=["hugo_symbol", "amino_acid_position"], how="left")


EXPANDED_ORDER = [
    # identity and release
    *KEYS, "transcript", "HGVSp_Short", "hotspot_version", "in_cancerhotspots_v2",
    "in_cancerhotspots_v3", "residue_hotspot_versions",
    # CancerHotspots evidence at the residue
    "substitutions", "n_unique_changes", "change_count", "position_total_count",
    "same_change_fraction", "top_change", "top_change_count", "top_change_fraction",
    "hotspot_character", "n_msk", "n_retro", "msk_fraction",
    # O7
    "o7_strength", "o7_points",
    # O4 on GENIE v20
    "gene_in_genie", "genie_rows", "genie_samples", "genie_patients", "genie_patients_msk",
    "genie_patients_non_msk", "genie_patients_myeloid", "genie_patients_lymphoid",
    "genie_patients_solid", "genie_patients_unknown",
    "o4_strength__genie_samples", "o4_strength__genie_patients", "o4_points__genie_patients",
    "o4_o7_uncapped__genie_patients", "o4_o7_capped__genie_patients",
    "o4_strength__genie_patients_excluding_msk", "o4_points__genie_patients_excluding_msk",
    "o4_o7_uncapped__genie_patients_excluding_msk", "o4_o7_capped__genie_patients_excluding_msk",
    "o4_strength__cosmic",
    # the independence test and both calibrations
    "o4_o7_uncapped", "residual_count", "residual_n_changes", "residual_max_change_count",
    "independent_positional_evidence", "recommended_o4_o7_handling",
    "residue_independent__permissive", "o4_o7_handling__permissive",
    "o4_o7_points__permissive", "points_lost__permissive",
    "residue_independent__strict", "o4_o7_handling__strict",
    "o4_o7_points__strict", "points_lost__strict",
    # O1
    "svig_uk_assessment", "on_svig_uk_canonical_list",
]


def weighting_options(df):
    """Request 2: Option 1 and Option 2 for every change the cap does not apply to."""
    uncapped = df[(df["o7_points"] > 0)
                  & (df["o4_points__genie_patients"] > 0)
                  & (df["residue_independent__permissive"]
                     | df["residue_independent__strict"])].copy()

    opt1 = [o7_rule(c, t, "VUA change count", "residue total")
            for c, t in zip(uncapped["change_count"], uncapped["position_total_count"])]
    opt2 = [o7_rule(m, r, "largest other change", "residual residue total")
            for m, r in zip(uncapped["residual_max_change_count"], uncapped["residual_count"])]
    uncapped["option1_o7_strength"] = [t for t, _, _ in opt1]
    uncapped["option1_o7_points"] = [p for _, p, _ in opt1]
    uncapped["option1_o7_rule"] = [r for _, _, r in opt1]
    uncapped["option2_o7_strength"] = [t for t, _, _ in opt2]
    uncapped["option2_o7_points"] = [p for _, p, _ in opt2]
    uncapped["option2_o7_rule"] = [r for _, _, r in opt2]
    # Two alternatives, for comparison (see the delivery note).
    opt3 = [min(a, b) for a, b in zip(uncapped["option1_o7_points"], [p for _, p, _ in opt2])]
    uncapped["option3_lesser_of_both_o7_points"] = opt3
    uncapped["option3_rule"] = [
        f"lesser of Option 1 ({a}) and Option 2 ({b}) -> +{min(a, b)}"
        for a, b in zip(uncapped["option1_o7_points"], uncapped["option2_o7_points"])]
    opt4 = [o7_rule(c, r, "VUA change count", "residual residue total")
            for c, r in zip(uncapped["change_count"], uncapped["residual_count"])]
    uncapped["option4_vua_count_residual_total_o7_strength"] = [t for t, _, _ in opt4]
    uncapped["option4_vua_count_residual_total_o7_points"] = [p for _, p, _ in opt4]
    uncapped["option4_rule"] = [r for _, _, r in opt4]

    uncapped["difference_category"] = [
        difference_category(c, t, m, r) for c, t, m, r in
        zip(uncapped["change_count"], uncapped["position_total_count"],
            uncapped["residual_max_change_count"], uncapped["residual_count"])]
    uncapped["difference_reason"] = [
        difference_reason(c, t, m, r) for c, t, m, r in
        zip(uncapped["change_count"], uncapped["position_total_count"],
            uncapped["residual_max_change_count"], uncapped["residual_count"])]

    o4 = uncapped["o4_points__genie_patients"]
    uncapped["option3_total_points"] = o4 + uncapped["option3_lesser_of_both_o7_points"]
    uncapped["option4_total_points"] = o4 + uncapped["option4_vua_count_residual_total_o7_points"]
    uncapped["option1_total_points"] = o4 + uncapped["option1_o7_points"]
    uncapped["option2_total_points"] = o4 + uncapped["option2_o7_points"]
    uncapped["points_difference_option1_minus_option2"] = (
        uncapped["option1_total_points"] - uncapped["option2_total_points"])
    uncapped["o7_tier_changes_between_options"] = (
        uncapped["option1_o7_strength"] != uncapped["option2_o7_strength"])

    cols = KEYS + ["hotspot_version", "substitutions", "n_unique_changes", "change_count",
                   "position_total_count", "residual_count", "residual_max_change_count",
                   "residue_independent__permissive", "residue_independent__strict",
                   "o7_strength", "o7_points",
                   "genie_patients", "o4_strength__genie_patients", "o4_points__genie_patients",
                   "option1_o7_strength", "option1_o7_points", "option1_o7_rule",
                   "option1_total_points",
                   "option2_o7_strength", "option2_o7_points", "option2_o7_rule",
                   "option2_total_points", "points_difference_option1_minus_option2",
                   "option3_lesser_of_both_o7_points", "option3_rule", "option3_total_points",
                   "option4_vua_count_residual_total_o7_strength",
                   "option4_vua_count_residual_total_o7_points", "option4_rule",
                   "option4_total_points",
                   "o7_tier_changes_between_options", "difference_category", "difference_reason",
                   "svig_uk_assessment", "on_svig_uk_canonical_list"]
    return uncapped[cols].sort_values(
        ["points_difference_option1_minus_option2", "change_count"], ascending=False)


def options_summary(opts, df):
    def block(sub, label):
        rows = [("Changes in this set", len(sub))]
        for opt in ("option1", "option2", "option4_vua_count_residual_total"):
            counts = sub[f"{opt}_o7_strength"].value_counts()
            for tier in ("Strong", "Moderate", "Supporting", "Not met"):
                rows.append((f"{opt}: O7 {tier}", int(counts.get(tier, 0))))
            rows.append((f"{opt}: mean O7 points", round(sub[f"{opt}_o7_points"].mean(), 2)))
            short = opt.split("_")[0]
            rows.append((f"{short}: changes reaching the full +8",
                         int((sub[f"{short}_total_points"] >= 8).sum())))
        rows.append(("option3 (lesser of 1 and 2): mean O7 points",
                     round(sub["option3_lesser_of_both_o7_points"].mean(), 2)))
        rows.append(("option3 (lesser of 1 and 2): changes reaching the full +8",
                     int((sub["option3_total_points"] >= 8).sum())))
        rows += [
            ("O7 tier differs between the options", int(sub["o7_tier_changes_between_options"].sum())),
            ("Mean points lost by using Option 2 instead of Option 1",
             round(sub["points_difference_option1_minus_option2"].mean(), 2)),
            ("Changes losing 3 points under Option 2 (Strong -> Supporting)",
             int((sub["points_difference_option1_minus_option2"] == 3).sum())),
            ("Changes losing 2 points under Option 2", int((sub["points_difference_option1_minus_option2"] == 2).sum())),
            ("Changes losing 1 point under Option 2", int((sub["points_difference_option1_minus_option2"] == 1).sum())),
            ("Changes unchanged between the options", int((sub["points_difference_option1_minus_option2"] == 0).sum())),
            ("Changes scoring MORE under Option 2", int((sub["points_difference_option1_minus_option2"] < 0).sum())),
        ]
        return pd.DataFrame(rows, columns=["metric", "value"]).assign(set=label)

    parts = [block(opts[opts["residue_independent__permissive"]], "uncapped under the permissive criterion"),
             block(opts[opts["residue_independent__strict"]], "uncapped under the strict criterion")]
    reasons = (opts[opts["o7_tier_changes_between_options"]]
               .groupby("difference_category").size().sort_values(ascending=False)
               .reset_index(name="value").rename(columns={"difference_category": "metric"}))
    reasons["metric"] = "Tier difference driven by: " + reasons["metric"]
    reasons["set"] = "all uncapped changes"
    context = pd.DataFrame([
        ("Changes in the analysis set", len(df)),
        ("Changes where both O4 and O7 apply",
         int(((df["o4_points__genie_patients"] > 0) & (df["o7_points"] > 0)).sum())),
        ("... uncapped under the permissive criterion",
         int((((df["o4_points__genie_patients"] > 0) & (df["o7_points"] > 0))
              & df["residue_independent__permissive"]).sum())),
        ("... uncapped under the strict criterion",
         int((((df["o4_points__genie_patients"] > 0) & (df["o7_points"] > 0))
              & df["residue_independent__strict"]).sum())),
    ], columns=["metric", "value"]).assign(set="context")
    return pd.concat([context, *parts, reasons])[["set", "metric", "value"]]


def main():
    df = add_version_comparison(load())
    missing = [c for c in EXPANDED_ORDER if c not in df.columns]
    assert not missing, f"expanded table is missing columns: {missing}"
    expanded = df[EXPANDED_ORDER].sort_values(
        ["o4_o7_uncapped", "change_count"], ascending=False)
    # Whole-number columns stay whole numbers even where the v3 rows are empty,
    # so the v2 rows read exactly as in the table sent before.
    for c in ("n_msk", "n_retro", "genie_rows", "genie_samples", "genie_patients"):
        if c in expanded:
            expanded[c] = expanded[c].astype("Int64")
    expanded.to_csv(RESULTS / "36_genie_points_per_change_expanded.tsv", sep="\t", index=False)

    opts = weighting_options(df)
    opts.to_csv(RESULTS / "37_genie_o7_weighting_options.tsv", sep="\t", index=False)
    summary = options_summary(opts, df)
    summary.to_csv(RESULTS / "38_genie_o7_weighting_options_summary.tsv", sep="\t", index=False)

    about = pd.DataFrame({"": [
        f"O4 / O7 per change, CancerHotspots {'v2 + v3' if paths.ANALYSIS == 'v3' else 'v2'} "
        "with GENIE v20 -- prepared for Kevin Baker, 30 September 2026.",
        "",
        "Sheets",
        "  Points per change (expanded) -- every column of the two tables sent before, plus the",
        "      canonical (O1) overlap, the v2 / v3 comparison, and both handlings. Nothing removed.",
        "  O7 weighting options -- for each change the cap does not apply to, O7 scored under",
        "      Option 1 (total-count) and Option 2 (residual-count), each with the rule it passed",
        "      or failed and the reason they differ.",
        "  O7 options summary -- the same, counted.",
        "",
        "O7 thresholds (SVIG-UK v1.1 Supplementary Table 1)",
        "  same change >= 10 AND residue total >= 50 -> Strong (+4)",
        "  same change >= 10 AND residue total < 50  -> Moderate (+2)",
        "  same change 2-9                            -> Supporting (+1)",
        "Option 1 applies these to the VUA's own change count and the full residue total.",
        "Option 2 applies them to the largest OTHER change and the residue total minus the VUA.",
        "",
        "O4: SVIG-UK thresholds on unique GENIE v20 patients (missense > 10 strong, 5-10 moderate).",
        "Uncapped = the leave-one-variant-out test found independent positional evidence, so O4 and",
        "O7 may be combined; permissive and strict differ in how big the largest other change must be.",
        "",
        "Contains AACR Project GENIE v20.0-public data per variant: not for redistribution.",
    ]})

    out = RESULTS / "genie_o4_o7_expanded.xlsx"
    sheets = [("About", about, False),
              ("Points per change (expanded)", expanded, True),
              ("O7 weighting options", opts, True),
              ("O7 options summary", summary, True)]
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        for name, frame, header in sheets:
            frame.to_excel(writer, sheet_name=name, index=False, header=header)
            ws = writer.sheets[name]
            if name == "About":
                ws.column_dimensions["A"].width = 110
                continue
            for i, col in enumerate(frame.columns, start=1):
                width = max(len(str(col)),
                            frame[col].astype(str).str.len().max() if len(frame) else 0)
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(width + 2, 55)
            ws.freeze_panes = "E2"
            ws.auto_filter.ref = ws.dimensions

    print(f"expanded table: {len(expanded):,} changes x {len(expanded.columns)} columns")
    print(f"uncapped changes scored under both options: {len(opts):,}")
    print(summary.to_string(index=False))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
