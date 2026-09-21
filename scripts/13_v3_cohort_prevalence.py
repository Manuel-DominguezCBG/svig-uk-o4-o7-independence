#!/usr/bin/env python3
"""
13_v3_cohort_prevalence.py
--------------------------
The cohort breakdown CancerHotspots v3 publishes for its new residues, and what it
says about O4 / O7 overlap. Run with ANALYSIS=v3.

For each of the 164 residues new in v3, the export's residue sheet gives the number
of mutations in three cohorts -- MSK-IMPACT, GENIE excluding MSK, and TCGA -- with
each cohort's sample count and tests of whether prevalence differs between them.

These are NOT a breakdown of the mutations behind each hotspot call, as n_MSK /
n_Retro were in v2: at the residue level they total 6,986 mutations, against 2,221
for the per-change counts that O7 is scored on. The file does not say which cohort
the v3 calls were made on. What the columns do show is where these residues are
observed, and MSK-IMPACT and non-MSK GENIE are both GENIE -- the database O4 is
counted in.

Inputs:
  data/raw/hotspots_v3.xlsx                    sheet Hotspot_Residues
  results/05_o4_o7_double_counting.tsv         (analyses/v3) the v3 changes
  results/29_genie_points_per_change.tsv       (analyses/v3, local) GENIE O4 per change

Outputs (analyses/v3/results/):
  33_v3_cohort_prevalence_by_residue.tsv       per residue: counts, prevalence, tests
  34_v3_cohort_prevalence_summary.tsv          totals and headline shares
  35_genie_v3_residue_o4.tsv                   per residue, GENIE O4 (local only)
"""

import pandas as pd

import paths

if paths.ANALYSIS != "v3":
    raise SystemExit("run with ANALYSIS=v3")

RESULTS = paths.RESULTS
V3 = paths.RAW / "hotspots_v3.xlsx"
KEYS = ["hugo_symbol", "amino_acid_position"]
COHORTS = {
    "msk": ("# mut in MSK", "# total samples in MSK"),
    "non_msk_genie": ("# mut in non-MSK-GENIE", "# total samples in non-MSK-GENIE"),
    "tcga": ("# mut in TCGA", "# total samples in TCGA"),
}


def main():
    res = pd.read_excel(V3, sheet_name="Hotspot_Residues")
    out = pd.DataFrame({
        "hugo_symbol": res["Hugo_Symbol"],
        "amino_acid_position": res["Codon_Position"].astype(int).astype(str),
        "codon": res["Codon"],
        "qvalue": res["Q value"],
    })
    for name, (mut, total) in COHORTS.items():
        out[f"mutations_{name}"] = res[mut]
        out[f"samples_{name}"] = res[total]
        out[f"prevalence_{name}_per_10k"] = (1e4 * res[mut] / res[total]).round(3)
    genie = out["mutations_msk"] + out["mutations_non_msk_genie"]
    everything = genie + out["mutations_tcga"]
    out["share_of_cohort_mutations_in_genie"] = (genie / everything).round(4)
    out["share_of_cohort_mutations_in_msk"] = (out["mutations_msk"] / everything).round(4)
    for c in ["pval_MSK_vs_GENIE", "adj_pval_MSK_vs_GENIE", "pval_MSK_vs_TCGA", "adj_pval_MSK_vs_TCGA"]:
        out[c] = res[c]

    # The per-change counts O7 is scored on, summed to the residue.
    dc = pd.read_csv(RESULTS / "05_o4_o7_double_counting.tsv", sep="\t", dtype={"amino_acid_position": str})
    v3 = dc[dc["hotspot_version"] == "v3"]
    per_res = (v3.groupby(KEYS)
               .agg(hotspot_mutations=("change_count", "sum"),
                    n_changes=("variant_aa", "size"),
                    best_o7=("o7_points", "max"))
               .reset_index())
    out = out.merge(per_res, on=KEYS, how="left", validate="one_to_one")
    assert out["hotspot_mutations"].notna().all(), "residue sheet and SNV sheet disagree"
    out = out.sort_values("hotspot_mutations", ascending=False)
    out.to_csv(RESULTS / "33_v3_cohort_prevalence_by_residue.tsv", sep="\t", index=False)

    tot = {n: int(out[f"mutations_{n}"].sum()) for n in COHORTS}
    all_cohorts = sum(tot.values())
    summary = pd.DataFrame([
        ("Residues new in CancerHotspots v3", len(out)),
        ("Genes", out["hugo_symbol"].nunique()),
        ("Mutations behind the v3 hotspot calls (per-change counts, as scored for O7)",
         int(out["hotspot_mutations"].sum())),
        ("Mutations at these residues in MSK-IMPACT", tot["msk"]),
        ("Mutations at these residues in GENIE excluding MSK", tot["non_msk_genie"]),
        ("Mutations at these residues in TCGA", tot["tcga"]),
        ("Share of the three cohorts' mutations that are in GENIE (MSK + non-MSK)",
         round((tot["msk"] + tot["non_msk_genie"]) / all_cohorts, 4)),
        ("Share that are in MSK-IMPACT", round(tot["msk"] / all_cohorts, 4)),
        ("Share that are in TCGA (not in GENIE)", round(tot["tcga"] / all_cohorts, 4)),
        ("Median per-residue share in GENIE", round(out["share_of_cohort_mutations_in_genie"].median(), 4)),
        ("Residues where MSK and non-MSK GENIE prevalence differ (adjusted p < 0.05)",
         int((out["adj_pval_MSK_vs_GENIE"] < 0.05).sum())),
        ("Residues where MSK and TCGA prevalence differ (adjusted p < 0.05)",
         int((out["adj_pval_MSK_vs_TCGA"] < 0.05).sum())),
        ("Residues with a change reaching O7 moderate or strong", int((out["best_o7"] >= 2).sum())),
    ], columns=["metric", "value"])
    summary.to_csv(RESULTS / "34_v3_cohort_prevalence_summary.tsv", sep="\t", index=False)

    # Local only: the same residues with the GENIE O4 of their changes.
    points = RESULTS / "29_genie_points_per_change.tsv"
    if points.exists():
        pc = pd.read_csv(points, sep="\t", dtype={"amino_acid_position": str})
        pc = pc[pc["hotspot_version"] == "v3"]
        g = (pc.groupby(KEYS)
             .agg(genie_patients=("genie_patients", "sum"),
                  genie_patients_non_msk=("genie_patients_non_msk", "sum"),
                  changes_o4_strong=("o4_strength__genie_patients", lambda s: int((s == "Strong").sum())),
                  changes_scoring_both=("o4_o7_uncapped", lambda s: int((s > 4).sum())))
             .reset_index())
        local = out[KEYS + ["codon", "hotspot_mutations", "mutations_msk", "mutations_non_msk_genie",
                            "mutations_tcga"]].merge(g, on=KEYS, how="left")
        local.to_csv(RESULTS / "35_genie_v3_residue_o4.tsv", sep="\t", index=False)
        rho = local["mutations_non_msk_genie"].corr(local["genie_patients_non_msk"], method="spearman")
        print(f"Spearman rho, residue sheet non-MSK GENIE mutations vs GENIE v20 non-MSK patients: {rho:.3f}")

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
