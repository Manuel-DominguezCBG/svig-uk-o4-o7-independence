#!/usr/bin/env python3
"""
10_genie_o4_o7.py
-----------------
Repeat the O4 x O7 analysis with GENIE -- the database SVIG-UK names for O4 --
in place of the COSMIC stand-in, with O4 counted as unique patients.

Three questions, each of which COSMIC could only approximate:

  1. How tightly do O4 and O7 fire together when O4 is counted the way the
     guideline says it should be (GENIE, one entry per patient)?
  2. The leave-MSK-out recount. GENIE records the contributing centre for every
     sample, so the O4 count can be taken with and without MSK-IMPACT -- the cohort
     GENIE shares with CancerHotspots. Once MSK is removed, O4 and O7 rest on
     disjoint patients (CancerHotspots' non-MSK half is TCGA and other public
     cohorts, which are not in GENIE). Whatever coupling survives that is not a
     shared-cohort artefact.
  3. What the proposed conditional cap costs when O4 comes from GENIE.

O4 thresholds are Supplementary Table 1 as written: missense >10 strong, 5-10
moderate; nonsense >50 strong, 20-50 moderate, 10-19 supporting. Stop-loss changes
use the missense thresholds. They and nonsense changes are ineligible for O7, so
they never reach the cap.

Inputs:
  results/05_o4_o7_double_counting.tsv     O7 tier and independence test per change
  results/15_canonical_list_overlap.tsv    O1 status and the COSMIC O4 tier
  data/interim/genie_hotspot_counts.tsv    from script 09

Outputs:
  results/19_genie_headline_summary.tsv
  results/20_genie_o4_o7_crosstab.tsv
  results/21_genie_points_impact_of_cap.tsv
  results/22_genie_variants_dropping_8_to_4.tsv
  results/23_genie_leave_msk_out_transitions.tsv
  results/24_genie_correlations.tsv
  results/26_genie_per_change_o4_o7_points.tsv
  results/27_genie_vs_cosmic_o4_tiers.tsv
  results/28_genie_on_target_lineage_o4.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"

KEYS = ["hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa"]
O7_ORDER = ["Strong", "Moderate", "Supporting", "Not met",
            "Not applicable (nonsense/stop-loss/synonymous)"]
O4_ORDER = ["Strong", "Moderate", "Supporting", "Not met", "Gene not in GENIE"]
POINTS = {"Strong": 4, "Moderate": 2, "Supporting": 1}

# O4 count source -> column. "patients" is the guideline's O4 count; "non_msk"
# is the leave-MSK-out recount; "samples" exists only to show what de-duplication
# changes.
SCHEMES = {
    "genie_patients": "genie_patients",
    "genie_patients_excluding_msk": "genie_patients_non_msk",
}


def o4_tier(count, variant_aa, gene_in_genie):
    if not gene_in_genie:
        return "Gene not in GENIE"
    if variant_aa == "*":                      # nonsense
        if count > 50:
            return "Strong"
        if count >= 20:
            return "Moderate"
        if count >= 10:
            return "Supporting"
        return "Not met"
    if count > 10:                             # missense (and stop-loss)
        return "Strong"
    if count >= 5:
        return "Moderate"
    return "Not met"


def load():
    dtype = {k: str for k in KEYS}
    o7 = pd.read_csv(RESULTS / "05_o4_o7_double_counting.tsv", sep="\t", dtype=dtype)
    o1 = pd.read_csv(RESULTS / "15_canonical_list_overlap.tsv", sep="\t", dtype=dtype,
                     usecols=KEYS + ["cosmic_samples_mutated", "o4_strength__genie_literal",
                                     "svig_uk_assessment", "on_svig_uk_canonical_list"])
    o1 = o1.rename(columns={"o4_strength__genie_literal": "o4_strength__cosmic"})
    genie = pd.read_csv(INTERIM / "genie_hotspot_counts.tsv", sep="\t", dtype=dtype)

    df = o7.merge(genie, on=KEYS, how="left", validate="one_to_one")
    df = df.merge(o1, on=KEYS, how="left", validate="one_to_one")
    assert df["genie_patients"].notna().all(), "analysis-set change missing GENIE counts"
    df["gene_in_genie"] = df["gene_in_genie"].astype(str).eq("True")
    df["independent_positional_evidence"] = (
        df["independent_positional_evidence"].astype(str).eq("True"))
    df["on_svig_uk_canonical_list"] = df["on_svig_uk_canonical_list"].astype(str).eq("True")
    return df


def score(df):
    df["o4_strength__genie_samples"] = [
        o4_tier(c, v, g) for c, v, g in
        zip(df["genie_samples"], df["variant_aa"], df["gene_in_genie"])]
    for name, col in SCHEMES.items():
        tier = [o4_tier(c, v, g) for c, v, g in
                zip(df[col], df["variant_aa"], df["gene_in_genie"])]
        df[f"o4_strength__{name}"] = tier
        df[f"o4_points__{name}"] = [POINTS.get(t, 0) for t in tier]
        uncapped = df[f"o4_points__{name}"] + df["o7_points"]
        capped = uncapped.where(df["independent_positional_evidence"], uncapped.clip(upper=4))
        df[f"o4_o7_uncapped__{name}"] = uncapped
        df[f"o4_o7_capped__{name}"] = capped
    return df


def crosstabs(df):
    out = []
    for name in SCHEMES:
        ct = (pd.crosstab(df[f"o4_strength__{name}"], df["o7_strength"])
              .reindex(index=O4_ORDER, columns=O7_ORDER).fillna(0).astype(int)
              .reset_index().rename(columns={f"o4_strength__{name}": "o4_strength"}))
        ct.insert(0, "o4_count_source", name)
        ct.columns = [c if c in ("o4_count_source", "o4_strength") else f"o7_{c}"
                      for c in ct.columns]
        out.append(ct)
    return pd.concat(out)


def cap_impact(df):
    rows = []
    for name in SCHEMES:
        unc, cap = df[f"o4_o7_uncapped__{name}"], df[f"o4_o7_capped__{name}"]
        affected = unc != cap
        o7s = df["o7_strength"] == "Strong"
        both_strong = o7s & (df[f"o4_strength__{name}"] == "Strong")
        drop = (unc == 8) & (cap == 4)
        rows.append({
            "o4_count_source": name,
            "changes_assessed": len(df),
            "changes_scoring_both_o4_and_o7": int(((df[f"o4_points__{name}"] > 0)
                                                   & (df["o7_points"] > 0)).sum()),
            "o7_strong_changes": int(o7s.sum()),
            "o7_strong_also_o4_strong": int(both_strong.sum()),
            "pct_o7_strong_also_o4_strong": round(100 * both_strong.sum() / o7s.sum(), 1),
            "changes_scoring_the_full_8_points": int((unc >= 8).sum()),
            "changes_reduced_by_the_cap": int(affected.sum()),
            "pct_reduced_by_the_cap": round(100 * affected.mean(), 1),
            "mean_points_lost_where_capped": (round((unc - cap)[affected].mean(), 2)
                                              if affected.any() else 0.0),
            "changes_dropping_from_8_to_4": int(drop.sum()),
            "of_which_on_o1_canonical_list": int((drop & df["on_svig_uk_canonical_list"]).sum()),
        })
    return pd.DataFrame(rows)


def leave_msk_out(df):
    """O4 tier on all patients vs on non-MSK patients, within each O7 tier."""
    out = []
    eligible = df[df["o7_strength"].isin(["Strong", "Moderate", "Supporting"])]
    for o7, sub in eligible.groupby("o7_strength"):
        ct = (pd.crosstab(sub["o4_strength__genie_patients"],
                          sub["o4_strength__genie_patients_excluding_msk"])
              .reindex(index=O4_ORDER, columns=O4_ORDER).fillna(0).astype(int))
        ct = ct.loc[ct.sum(axis=1) > 0]
        ct.index.name = "o4_all_patients"
        ct = ct.reset_index()
        ct.insert(0, "o7_strength", o7)
        ct.columns = [c if c in ("o7_strength", "o4_all_patients") else f"o4_excl_msk_{c}"
                      for c in ct.columns]
        out.append(ct)
    return pd.concat(out).sort_values("o7_strength", key=lambda s: s.map(
        {"Strong": 0, "Moderate": 1, "Supporting": 2}))


def spearman(a, b):
    return round(pd.Series(a).corr(pd.Series(b), method="spearman"), 3)


def correlations(df):
    found = df[df["genie_patients"] > 0]
    o7_mod_plus = found[found["o7_points"] >= 2]

    # Residue level. CancerHotspots publishes its MSK / retrospective split per
    # residue, not per change, so the cohort-disjoint comparison is made there.
    # GENIE residue totals are sums over changes (a patient carrying two changes at
    # one residue is counted twice -- rare, and immaterial to a rank correlation).
    res = (df.groupby(["hugo_symbol", "amino_acid_position"])
           .agg(n_msk=("n_msk", "first"), n_retro=("n_retro", "first"),
                genie_msk=("genie_patients_msk", "sum"),
                genie_non_msk=("genie_patients_non_msk", "sum"),
                gene_in_genie=("gene_in_genie", "first"))
           .reset_index())
    res = res[res["gene_in_genie"]]

    rows = [
        ("Changes in the analysis set", len(df)),
        ("Changes observed in GENIE (>= 1 patient)", len(found)),
        ("Spearman rho: CancerHotspots count vs GENIE patients (all centres)",
         spearman(found["change_count"], found["genie_patients"])),
        ("Pearson r on log10: CancerHotspots count vs GENIE patients (all centres)",
         round(np.log10(found["change_count"]).corr(np.log10(found["genie_patients"])), 3)),
        ("Spearman rho: CancerHotspots count vs GENIE patients excluding MSK",
         spearman(found["change_count"], found["genie_patients_non_msk"])),
        ("Spearman rho, changes reaching O7 moderate or strong, vs GENIE patients",
         spearman(o7_mod_plus["change_count"], o7_mod_plus["genie_patients"])),
        ("Residues with the gene sequenced in GENIE", len(res)),
        ("Spearman rho, residue: CancerHotspots MSK count vs GENIE MSK patients "
         "(same cohort; sanity check)",
         spearman(res["n_msk"], res["genie_msk"])),
        ("Spearman rho, residue: CancerHotspots retrospective count vs GENIE "
         "non-MSK patients (disjoint cohorts)",
         spearman(res["n_retro"], res["genie_non_msk"])),
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def genie_vs_cosmic(df):
    eligible = df[df["o7_points"] > 0]
    ct = (pd.crosstab(eligible["o4_strength__genie_patients"],
                      eligible["o4_strength__cosmic"].fillna("Absent/not counted"))
          .reindex(index=O4_ORDER).dropna(how="all").fillna(0).astype(int))
    ct.index.name = "o4_genie_patients"
    return ct.reset_index()


def headline(df, impact, comp):
    elig = df["o7_points"] > 0
    all_c = comp[comp["cohort"] == "All centres"].set_index("lineage")
    dedup_drop = elig & (df["o4_strength__genie_samples"] != df["o4_strength__genie_patients"])
    o7s = df["o7_strength"] == "Strong"
    strong_all = o7s & (df["o4_strength__genie_patients"] == "Strong")
    strong_both = strong_all & (df["o4_strength__genie_patients_excluding_msk"] == "Strong")
    imp = impact.set_index("o4_count_source")

    def pct(n, d):
        return round(100 * n / d, 1) if d else 0.0

    rows = [
        ("GENIE release", "20.0-public"),
        ("GENIE samples", int(all_c.loc["Total", "samples"])),
        ("GENIE patients", int(all_c.loc["Total", "patients"])),
        ("Myeloid patients (% of GENIE)", f'{all_c.loc["Myeloid", "patients"]:,} '
                                          f'({all_c.loc["Myeloid", "pct_patients"]}%)'),
        ("Lymphoid patients (% of GENIE)", f'{all_c.loc["Lymphoid", "patients"]:,} '
                                           f'({all_c.loc["Lymphoid", "pct_patients"]}%)'),
        ("Solid tumour patients (% of GENIE)", f'{all_c.loc["Solid", "patients"]:,} '
                                               f'({all_c.loc["Solid", "pct_patients"]}%)'),
        ("Hotspot changes observed in GENIE", f'{(df.genie_patients > 0).sum():,} / {len(df):,}'),
        ("Hotspot-change observations: MAF rows", int(df["genie_rows"].sum())),
        ("Hotspot-change observations: distinct samples", int(df["genie_samples"].sum())),
        ("Hotspot-change observations: distinct patients", int(df["genie_patients"].sum())),
        ("Changes whose count falls on de-duplication to patients",
         int((df["genie_samples"] > df["genie_patients"]).sum())),
        ("O7-scoring changes whose O4 tier falls on de-duplication", int(dedup_drop.sum())),
        ("O7_strong changes also O4_strong (GENIE, all patients)",
         f"{strong_all.sum()} / {o7s.sum()} ({pct(strong_all.sum(), o7s.sum())}%)"),
        ("O7_strong changes still O4_strong with MSK excluded",
         f"{strong_both.sum()} / {o7s.sum()} ({pct(strong_both.sum(), o7s.sum())}%)"),
        ("Changes scoring the full +8 (all patients)",
         int(imp.loc["genie_patients", "changes_scoring_the_full_8_points"])),
        ("Changes scoring the full +8 (MSK excluded)",
         int(imp.loc["genie_patients_excluding_msk", "changes_scoring_the_full_8_points"])),
        ("Changes reduced by the cap (all patients)",
         f'{imp.loc["genie_patients", "changes_reduced_by_the_cap"]} '
         f'({imp.loc["genie_patients", "pct_reduced_by_the_cap"]}%)'),
        ("Changes dropping from +8 to +4 under the cap (all patients)",
         int(imp.loc["genie_patients", "changes_dropping_from_8_to_4"])),
        ("... of which on the O1 canonical list",
         int(imp.loc["genie_patients", "of_which_on_o1_canonical_list"])),
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def on_target_lineage(df):
    """O4 counted only on 'on-target' lineages.

    Supplementary Table 1 lets case-counting be limited to 'on-target' tumour
    types, and for AML gives the options AML only, myeloid disorders, or all
    haematopoietic disorders. CancerHotspots is pan-cancer and solid-dominated, so
    restricting O4 to a haematological denominator is where O4 and O7 stop
    measuring the same cohort -- this shows how much of O4 survives it.
    """
    el = df[(df["o7_points"] > 0) & df["gene_in_genie"]]
    sources = [
        ("All patients", el["genie_patients"]),
        ("Solid tumours only", el["genie_patients_solid"]),
        ("Haematological (myeloid + lymphoid)",
         el["genie_patients_myeloid"] + el["genie_patients_lymphoid"]),
        ("Myeloid only", el["genie_patients_myeloid"]),
        ("Lymphoid only", el["genie_patients_lymphoid"]),
    ]
    o7s = el["o7_strength"] == "Strong"
    rows = []
    for label, counts in sources:
        tier = pd.Series([o4_tier(c, v, True) for c, v in zip(counts, el["variant_aa"])],
                         index=el.index)
        pts = tier.map(POINTS).fillna(0)
        rows.append({
            "o4_counted_on": label,
            "o7_scoring_changes": len(el),
            "o4_strong": int((tier == "Strong").sum()),
            "o4_moderate": int((tier == "Moderate").sum()),
            "o4_not_met": int((tier == "Not met").sum()),
            "o7_strong_also_o4_strong": f"{int((o7s & (tier == 'Strong')).sum())} / {int(o7s.sum())}",
            "changes_scoring_the_full_8_points": int(((pts + el["o7_points"]) >= 8).sum()),
        })
    return pd.DataFrame(rows)


def main():
    df = score(load())
    comp = pd.read_csv(RESULTS / "18_genie_cohort_composition.tsv", sep="\t")

    ct = crosstabs(df)
    impact = cap_impact(df)
    lmo = leave_msk_out(df)
    corr = correlations(df)
    vs_cosmic = genie_vs_cosmic(df)
    on_target = on_target_lineage(df)
    head = headline(df, impact, comp)

    drop = df[(df["o4_o7_uncapped__genie_patients"] == 8)
              & (df["o4_o7_capped__genie_patients"] == 4)]
    drop_cols = KEYS + ["change_count", "position_total_count", "n_unique_changes",
                        "residual_count", "genie_patients", "genie_patients_non_msk",
                        "o4_strength__genie_patients_excluding_msk", "msk_fraction",
                        "svig_uk_assessment", "on_svig_uk_canonical_list"]

    per_change_cols = KEYS + [
        "change_count", "position_total_count", "n_unique_changes", "hotspot_character",
        "o7_strength", "o7_points", "independent_positional_evidence",
        "recommended_o4_o7_handling",
        "gene_in_genie", "genie_rows", "genie_samples", "genie_patients",
        "genie_patients_msk", "genie_patients_non_msk",
        "genie_patients_myeloid", "genie_patients_lymphoid", "genie_patients_solid",
        "genie_patients_unknown",
        "o4_strength__genie_samples",
        "o4_strength__genie_patients", "o4_points__genie_patients",
        "o4_o7_uncapped__genie_patients", "o4_o7_capped__genie_patients",
        "o4_strength__genie_patients_excluding_msk",
        "o4_points__genie_patients_excluding_msk",
        "o4_o7_uncapped__genie_patients_excluding_msk",
        "o4_o7_capped__genie_patients_excluding_msk",
        "o4_strength__cosmic", "msk_fraction", "n_msk", "n_retro",
        "svig_uk_assessment", "on_svig_uk_canonical_list",
    ]

    head.to_csv(RESULTS / "19_genie_headline_summary.tsv", sep="\t", index=False)
    ct.to_csv(RESULTS / "20_genie_o4_o7_crosstab.tsv", sep="\t", index=False)
    impact.to_csv(RESULTS / "21_genie_points_impact_of_cap.tsv", sep="\t", index=False)
    (drop[drop_cols].sort_values("change_count", ascending=False)
     .to_csv(RESULTS / "22_genie_variants_dropping_8_to_4.tsv", sep="\t", index=False))
    lmo.to_csv(RESULTS / "23_genie_leave_msk_out_transitions.tsv", sep="\t", index=False)
    corr.to_csv(RESULTS / "24_genie_correlations.tsv", sep="\t", index=False)
    (df[per_change_cols].sort_values("genie_patients", ascending=False)
     .to_csv(RESULTS / "26_genie_per_change_o4_o7_points.tsv", sep="\t", index=False))
    vs_cosmic.to_csv(RESULTS / "27_genie_vs_cosmic_o4_tiers.tsv", sep="\t", index=False)
    on_target.to_csv(RESULTS / "28_genie_on_target_lineage_o4.tsv", sep="\t", index=False)

    for title, table in [("HEADLINE", head), ("O4 x O7", ct), ("CAP IMPACT", impact.T),
                         ("DROPPING 8 -> 4", drop[drop_cols]), ("LEAVE-MSK-OUT", lmo),
                         ("CORRELATIONS", corr), ("GENIE vs COSMIC O4 tier", vs_cosmic),
                         ("ON-TARGET LINEAGE O4", on_target)]:
        print(f"\n=== {title} ===")
        print(table.to_string(index=title in ("CAP IMPACT",)))


if __name__ == "__main__":
    main()
