#!/usr/bin/env python3
"""
02_position_analysis.py
-----------------------
Turn the merged CancerHotspots per-allele table into the deliverables:

  results/01_headline_summary.tsv          one-line-per-metric overview
  results/02_distribution_unique_changes.tsv   1 / 2 / 3 / 4+ changes per hotspot
  results/03_gene_position_table.tsv       complete gene x position table
  results/04_per_allele_o7_tiers.tsv       every amino-acid change with its SVIG-UK O7 tier
  results/05_o4_o7_double_counting.tsv     per-change independence assessment + cap recommendation
  results/06_msk_overlap_summary.tsv       MSK-IMPACT contribution to the hotspot evidence
  results/07_cap_rule_sensitivity.tsv      sensitivity of the cap rule to threshold choice
  results/08_version_comparison.tsv        v1 (2016) vs v2 (2018) consistency check

The primary analysis set is the v2 (Chang et al., 2018) SNV hotspot table, restricted
to protein-coding (non-splice) positions, because that is the resource an analyst
queries at cancerhotspots.org and the only one carrying the MSK / retrospective split.

SVIG-UK v1.1 (26/05/2026), Supplementary Table 1, O7 - cancerhotspots.org tiers:
  Strong     [+4]  >= 50 mutations at the same amino-acid position
                   AND >= 10 mutations for the same amino-acid change
  Moderate   [+2]  <  50 mutations at the same amino-acid position
                   AND >= 10 mutations for the same amino-acid change
  Supporting [+1]  2-9 mutations for the same amino-acid change
  (a change seen only once yields no O7 from this resource)
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"

# --- SVIG-UK O7 thresholds (Supplementary Table 1) -------------------------
O7_POSITION_STRONG = 50     # mutations at the same amino-acid position
O7_CHANGE_HIGH = 10         # mutations for the same amino-acid change
O7_CHANGE_SUPPORTING = 2    # lower bound for supporting

# --- Thresholds for the proposed independence test -------------------------
# "Residual" = the evidence at the position that does NOT come from the variant
# under assessment's own amino-acid change. RESIDUAL_MIN mirrors the >=10
# requirement SVIG-UK already uses to lift O7 above supporting.
RESIDUAL_MIN = 10           # residual mutations required for independent O7
RESIDUAL_RECURRENT_MIN = 2  # at least one other change must itself be recurrent
DOMINANCE_HIGH = 0.80       # >= this fraction from one change => mutation-specific
DOMINANCE_LOW = 0.50        # <  this fraction from one change => position-specific


def o7_applicable(reference_aa, variant_aa):
    """O7 is reserved for missense and small in-frame insertion/deletion variants
    (SVIG-UK Supplementary Figure 1, notes 4 and 13). A nonsense change is O2
    territory and O2 + O7 is not a permitted combination; a stop-loss change is
    O9 territory. Both are excluded from O7 scoring here, though their counts
    still contribute to the published position total, as they do on the website."""
    return variant_aa != "*" and reference_aa != "*"


def o7_tier(change_count, position_total):
    """SVIG-UK O7 strength/points achievable from cancerhotspots.org alone."""
    if change_count >= O7_CHANGE_HIGH:
        if position_total >= O7_POSITION_STRONG:
            return "Strong", 4
        return "Moderate", 2
    if change_count >= O7_CHANGE_SUPPORTING:
        return "Supporting", 1
    return "Not met", 0


def hotspot_character(fraction, n_changes):
    if n_changes == 1 or fraction >= DOMINANCE_HIGH:
        return "mutation-specific"
    if fraction >= DOMINANCE_LOW:
        return "mixed"
    return "position-specific"


def build_position_table(alleles):
    """One row per (gene, amino-acid position)."""
    alleles = alleles.sort_values(["hugo_symbol", "amino_acid_position", "change_count"],
                                  ascending=[True, True, False])

    def summarise(g):
        total = int(g["position_total_count"].iloc[0])
        detail = "|".join(f"{r.variant_aa}:{int(r.change_count)}" for r in g.itertuples())
        top_count = int(g["change_count"].max())
        n_msk = g["n_msk"].iloc[0]
        return pd.Series(
            {
                "reference_aa": g["reference_aa"].iloc[0],
                "n_unique_changes": int(g["variant_aa"].nunique()),
                "position_total_count": total,
                "substitutions": detail,
                "top_change": g["variant_aa"].iloc[0],
                "top_change_count": top_count,
                "top_change_fraction": round(top_count / total, 4),
                "n_changes_ge10": int((g["change_count"] >= O7_CHANGE_HIGH).sum()),
                "n_changes_ge2": int((g["change_count"] >= O7_CHANGE_SUPPORTING).sum()),
                "n_singleton_changes": int((g["change_count"] == 1).sum()),
                "n_msk": int(n_msk),
                "n_retro": int(g["n_retro"].iloc[0]),
                "msk_fraction": round(n_msk / total, 4),
                "qvalue": g["qvalue"].iloc[0],
            }
        )

    pos = (
        alleles.groupby(["hugo_symbol", "amino_acid_position"], sort=False)
        .apply(summarise, include_groups=False)
        .reset_index()
    )
    pos["hotspot_character"] = [
        hotspot_character(f, n)
        for f, n in zip(pos["top_change_fraction"], pos["n_unique_changes"])
    ]
    pos["position_reaches_o7_strong_threshold"] = pos["position_total_count"] >= O7_POSITION_STRONG
    return pos


def build_allele_table(alleles, pos):
    """One row per (gene, position, amino-acid change) = one candidate variant
    under assessment, with its O7 tier and its independence assessment."""
    keys = ["hugo_symbol", "amino_acid_position"]
    df = alleles.merge(
        pos[keys + ["n_unique_changes", "top_change", "top_change_count",
                    "top_change_fraction", "msk_fraction", "hotspot_character",
                    "substitutions"]],
        on=keys, how="left",
    )

    df["change_count"] = df["change_count"].astype(int)
    df["position_total_count"] = df["position_total_count"].astype(int)
    df["n_msk"] = df["n_msk"].astype(int)
    df["n_retro"] = df["n_retro"].astype(int)

    df["o7_applicable"] = [
        o7_applicable(r, v) for r, v in zip(df["reference_aa"], df["variant_aa"])
    ]
    tiers = [
        o7_tier(c, t) if ok else ("Not applicable (nonsense/stop-loss)", 0)
        for c, t, ok in zip(df["change_count"], df["position_total_count"],
                            df["o7_applicable"])
    ]
    df["o7_strength"] = [t[0] for t in tiers]
    df["o7_points"] = [t[1] for t in tiers]

    df["same_change_fraction"] = (df["change_count"] / df["position_total_count"]).round(4)
    df["residual_count"] = df["position_total_count"] - df["change_count"]
    df["residual_n_changes"] = df["n_unique_changes"] - 1

    # Largest single alternative amino-acid change at the position, i.e. the
    # strongest recurrence signal that survives removing the variant under
    # assessment's own change ("leave-one-variant-out").
    counts_by_pos = (
        alleles.groupby(keys)["change_count"]
        .apply(lambda s: sorted(int(x) for x in s))
        .to_dict()
    )
    residual_max = []
    for gene, position, own in zip(df["hugo_symbol"], df["amino_acid_position"],
                                   df["change_count"]):
        others = list(counts_by_pos[(gene, position)])
        others.remove(own)          # drop one instance of this change's own count
        residual_max.append(max(others) if others else 0)
    df["residual_max_change_count"] = residual_max

    # Permissive test: after removing the VUA's own change, does the position
    # retain enough recurrence to stand as a hotspot in its own right?
    df["independent_positional_evidence"] = (
        (df["residual_count"] >= RESIDUAL_MIN)
        & (df["residual_max_change_count"] >= RESIDUAL_RECURRENT_MIN)
    )
    # Conservative test: does at least one OTHER amino-acid change on its own
    # meet the >=10 same-change requirement SVIG-UK uses to lift O7 above
    # supporting? This anchors the independence test to the guideline's own bar.
    df["independent_positional_evidence_strict"] = (
        df["residual_max_change_count"] >= O7_CHANGE_HIGH
    )
    df["recommended_o4_o7_handling"] = [
        "O4 + O7 may be combined (max +8)" if ind
        else "Cap combined O4 + O7 at +4"
        for ind in df["independent_positional_evidence"]
    ]
    df["recommended_o4_o7_handling_strict"] = [
        "O4 + O7 may be combined (max +8)" if ind
        else "Cap combined O4 + O7 at +4"
        for ind in df["independent_positional_evidence_strict"]
    ]
    return df


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    merged = pd.read_csv(INTERIM / "hotspots_merged_per_allele.tsv", sep="\t")

    snv = merged[(merged["source"] == "v2:SNV-hotspots")
                 & (merged["variant_class"] == "snv")].copy()

    pos = build_position_table(snv)
    alleles = build_allele_table(snv, pos)

    # ---- 03 complete gene x position table --------------------------------
    pos_out = pos[[
        "hugo_symbol", "amino_acid_position", "reference_aa", "n_unique_changes",
        "substitutions", "position_total_count", "top_change", "top_change_count",
        "top_change_fraction", "hotspot_character", "n_changes_ge10", "n_changes_ge2",
        "n_singleton_changes", "position_reaches_o7_strong_threshold",
        "n_msk", "n_retro", "msk_fraction", "qvalue",
    ]].sort_values(["position_total_count", "hugo_symbol"], ascending=[False, True])
    pos_out.to_csv(RESULTS / "03_gene_position_table.tsv", sep="\t", index=False)

    # ---- 02 distribution of unique changes per hotspot --------------------
    band = pos["n_unique_changes"].clip(upper=4).map({1: "1", 2: "2", 3: "3", 4: "4+"})
    dist = (
        band.value_counts().rename_axis("unique_changes_at_position")
        .reset_index(name="n_positions")
        .sort_values("unique_changes_at_position")
    )
    dist["pct_of_positions"] = (100 * dist["n_positions"] / len(pos)).round(1)
    exact = (
        pos["n_unique_changes"].value_counts().sort_index()
        .rename_axis("exact_n_unique_changes").reset_index(name="n_positions")
    )
    exact["pct_of_positions"] = (100 * exact["n_positions"] / len(pos)).round(1)
    dist.to_csv(RESULTS / "02_distribution_unique_changes.tsv", sep="\t", index=False)
    exact.to_csv(RESULTS / "02b_distribution_unique_changes_exact.tsv", sep="\t", index=False)

    # ---- 04 per-allele O7 tiers -------------------------------------------
    # The positional columns from the gene x position table travel with each
    # change, so this one file answers both "how positional is the residue?"
    # and "what does O7 award this variant?" without a join.
    allele_out = alleles[[
        "hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa",
        "change_count", "position_total_count", "same_change_fraction",
        "n_unique_changes", "substitutions", "top_change", "top_change_count",
        "top_change_fraction", "hotspot_character",
        "o7_applicable", "o7_strength", "o7_points",
        "n_msk", "n_retro", "msk_fraction",
    ]].sort_values(["change_count"], ascending=False)
    allele_out.to_csv(RESULTS / "04_per_allele_o7_tiers.tsv", sep="\t", index=False)

    # ---- 05 double-counting / cap recommendation --------------------------
    dc = alleles[[
        "hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa",
        "change_count", "position_total_count", "same_change_fraction",
        "n_unique_changes", "substitutions", "residual_count", "residual_n_changes",
        "residual_max_change_count", "hotspot_character", "o7_applicable",
        "o7_strength", "o7_points",
        "independent_positional_evidence", "recommended_o4_o7_handling",
        "independent_positional_evidence_strict", "recommended_o4_o7_handling_strict",
        "n_msk", "n_retro", "msk_fraction",
    ]].sort_values(["o7_points", "change_count"], ascending=False)
    dc.to_csv(RESULTS / "05_o4_o7_double_counting.tsv", sep="\t", index=False)

    # ---- 06 MSK overlap ----------------------------------------------------
    msk_rows = [
        ("Hotspot positions analysed (v2 SNV, non-splice)", len(pos)),
        ("Total mutations underpinning those positions", int(pos["position_total_count"].sum())),
        ("Mutations contributed by MSK-IMPACT (n_MSK)", int(pos["n_msk"].sum())),
        ("Mutations contributed by retrospective/public cohorts (n_Retro)", int(pos["n_retro"].sum())),
        ("Overall MSK fraction of hotspot evidence", round(pos["n_msk"].sum() / pos["position_total_count"].sum(), 4)),
        ("Median per-position MSK fraction", round(pos["msk_fraction"].median(), 4)),
        ("Positions where MSK contributes >= 50% of the evidence", int((pos["msk_fraction"] >= 0.5).sum())),
        ("Positions where MSK contributes >= 75% of the evidence", int((pos["msk_fraction"] >= 0.75).sum())),
        ("Positions where MSK contributes 100% of the evidence", int((pos["msk_fraction"] >= 1.0).sum())),
    ]
    # Left-closed bands, so the band counts reconcile exactly with the ">= 50%"
    # and ">= 75%" figures quoted alongside them.
    msk_band = pd.cut(pos["msk_fraction"], [0, .25, .5, .75, 1.0001], right=False,
                      labels=["0-25%", "25-50%", "50-75%", "75-100%"])
    msk_dist = (msk_band.value_counts().sort_index()
                .rename_axis("msk_fraction_band").reset_index(name="n_positions"))
    msk_dist["pct_of_positions"] = (100 * msk_dist["n_positions"] / len(pos)).round(1)
    pd.DataFrame(msk_rows, columns=["metric", "value"]).to_csv(
        RESULTS / "06_msk_overlap_summary.tsv", sep="\t", index=False)
    msk_dist.to_csv(RESULTS / "06b_msk_fraction_distribution.tsv", sep="\t", index=False)

    # ---- 07 sensitivity of the cap rule ------------------------------------
    sens = []
    scored = alleles[alleles["o7_points"] > 0]
    for residual_min in (5, 10, 20, 50):
        for recurrent_min in (1, 2, 5):
            ind = ((scored["residual_count"] >= residual_min)
                   & (scored["residual_max_change_count"] >= recurrent_min))
            sens.append({
                "residual_min": residual_min,
                "residual_recurrent_min": recurrent_min,
                "n_changes_scoring_o7": len(scored),
                "n_independent_o7": int(ind.sum()),
                "pct_independent_o7": round(100 * ind.mean(), 1),
                "n_capped": int((~ind).sum()),
                "pct_capped": round(100 * (~ind).mean(), 1),
            })
    pd.DataFrame(sens).to_csv(RESULTS / "07_cap_rule_sensitivity.tsv", sep="\t", index=False)

    # ---- 08 v1 vs v2 comparison --------------------------------------------
    def pos_key(df):
        return set(zip(df["hugo_symbol"], df["amino_acid_position"].astype(str)))
    v1 = merged[merged["source"] == "v1:Per Residue"]
    v1_pos, v2_pos = pos_key(v1), pos_key(snv)
    comparison = pd.DataFrame([
        ("v1 (Chang 2016) hotspot positions", len(v1_pos)),
        ("v2 (Chang 2018) SNV hotspot positions (non-splice)", len(v2_pos)),
        ("Positions present in both versions", len(v1_pos & v2_pos)),
        ("Positions only in v1", len(v1_pos - v2_pos)),
        ("Positions only in v2", len(v2_pos - v1_pos)),
    ], columns=["metric", "value"])
    comparison.to_csv(RESULTS / "08_version_comparison.tsv", sep="\t", index=False)

    # ---- 09 worked examples -------------------------------------------------
    examples = [
        ("BRAF", 600, "E"), ("EGFR", 858, "R"), ("EGFR", 790, "M"),
        ("KRAS", 12, "D"), ("KRAS", 12, "C"), ("KRAS", 61, "H"),
        ("NRAS", 61, "R"), ("IDH1", 132, "H"), ("PIK3CA", 1047, "R"),
        ("PIK3CA", 545, "K"), ("TP53", 273, "H"), ("JAK2", 617, "F"),
        ("FGFR3", 249, "C"), ("KIT", 816, "V"), ("MYD88", 265, "P"),
    ]
    want = pd.DataFrame(examples, columns=["hugo_symbol", "amino_acid_position", "variant_aa"])
    want["amino_acid_position"] = want["amino_acid_position"].astype(str)
    dc = dc.copy()
    dc["amino_acid_position"] = dc["amino_acid_position"].astype(str)
    worked = want.merge(dc, on=["hugo_symbol", "amino_acid_position", "variant_aa"], how="left")
    worked.insert(3, "in_cancerhotspots_v2", worked["change_count"].notna())
    worked.to_csv(RESULTS / "09_worked_examples.tsv", sep="\t", index=False)

    # ---- 10 coverage check for haematological hotspots ----------------------
    # SVIG-UK notes that cancerhotspots.org coverage is limited for
    # haematological malignancies; this quantifies that for well-known targets.
    haem = [
        ("JAK2", "617", "V617F"), ("MPL", "515", "W515L/K"),
        ("CALR", None, "exon 9 indels"), ("NPM1", None, "exon 12 indels"),
        ("ASXL1", None, "G646Wfs*12 and other truncating"),
        ("SETBP1", "870", "D868/G870 hotspots"), ("CBL", "371", "Y371 RING finger"),
        ("SF3B1", "700", "K700E"), ("SRSF2", "95", "P95H/L/R"),
        ("U2AF1", "34", "S34F/Y"), ("DNMT3A", "882", "R882H/C"),
        ("IDH2", "140", "R140Q"), ("IDH2", "172", "R172K"),
        ("FLT3", "835", "D835 TKD"), ("KIT", "816", "D816V"),
        ("MYD88", "265", "L265P"),
    ]
    indel = merged[merged["source"] == "v2:INDEL-hotspots"]
    rows = []
    for gene, position, label in haem:
        sub = pos[pos["hugo_symbol"] == gene]
        gene_positions = len(sub)
        gene_indel_regions = indel.loc[indel["hugo_symbol"] == gene,
                                       "amino_acid_position"].nunique()
        if position is not None:
            sub = sub[sub["amino_acid_position"].astype(str) == position]
        rows.append({
            "hugo_symbol": gene,
            "amino_acid_position": position if position else "(any)",
            "clinically_important_variant": label,
            "present_in_snv_sheet": len(sub) > 0,
            "indel_regions_listed_for_this_gene": gene_indel_regions,
            "positions_listed_for_this_gene": gene_positions,
            "position_total_count": int(sub["position_total_count"].iloc[0]) if len(sub) else 0,
            "n_unique_changes": int(sub["n_unique_changes"].iloc[0]) if len(sub) else 0,
            "msk_fraction": float(sub["msk_fraction"].iloc[0]) if len(sub) else None,
        })
    pd.DataFrame(rows).to_csv(RESULTS / "10_haematology_coverage_check.tsv",
                              sep="\t", index=False)

    # ---- 01 headline summary ------------------------------------------------
    scored_all = alleles
    headline = pd.DataFrame([
        ("Hotspot positions analysed (v2 SNV, non-splice)", len(pos)),
        ("Distinct amino-acid changes across those positions", len(alleles)),
        ("Mean unique amino-acid changes per position", round(pos["n_unique_changes"].mean(), 2)),
        ("Median unique amino-acid changes per position", int(pos["n_unique_changes"].median())),
        ("Positions with a single amino-acid change", int((pos["n_unique_changes"] == 1).sum())),
        ("Positions with a single amino-acid change (%)", round(100 * (pos["n_unique_changes"] == 1).mean(), 1)),
        ("Mean fraction from the most common substitution", round(pos["top_change_fraction"].mean(), 4)),
        ("Median fraction from the most common substitution", round(pos["top_change_fraction"].median(), 4)),
        ("Positions classed mutation-specific (top change >= 80%)", int((pos["hotspot_character"] == "mutation-specific").sum())),
        ("Positions classed mixed (50-80%)", int((pos["hotspot_character"] == "mixed").sum())),
        ("Positions classed position-specific (top change < 50%)", int((pos["hotspot_character"] == "position-specific").sum())),
        ("Positions meeting the >=50 count threshold for O7_strong", int(pos["position_reaches_o7_strong_threshold"].sum())),
        ("Amino-acid changes reaching O7_strong (+4)", int((scored_all["o7_strength"] == "Strong").sum())),
        ("Amino-acid changes reaching O7_moderate (+2)", int((scored_all["o7_strength"] == "Moderate").sum())),
        ("Amino-acid changes reaching O7_supporting (+1)", int((scored_all["o7_strength"] == "Supporting").sum())),
        ("Amino-acid changes not meeting any O7 tier", int((scored_all["o7_strength"] == "Not met").sum())),
        ("Amino-acid changes ineligible for O7 (nonsense/stop-loss)", int((~scored_all["o7_applicable"]).sum())),
        ("O7-scoring changes with independent positional evidence", int(scored_all.loc[scored_all["o7_points"] > 0, "independent_positional_evidence"].sum())),
        ("O7-scoring changes recommended for a +4 cap", int((~scored_all.loc[scored_all["o7_points"] > 0, "independent_positional_evidence"]).sum())),
        ("O7-scoring changes independent under the strict test", int(scored_all.loc[scored_all["o7_points"] > 0, "independent_positional_evidence_strict"].sum())),
        ("O7-scoring changes recommended for a +4 cap (strict test)", int((~scored_all.loc[scored_all["o7_points"] > 0, "independent_positional_evidence_strict"]).sum())),
        ("Overall MSK-IMPACT fraction of hotspot evidence", round(pos["n_msk"].sum() / pos["position_total_count"].sum(), 4)),
    ], columns=["metric", "value"])
    headline.to_csv(RESULTS / "01_headline_summary.tsv", sep="\t", index=False)

    print(headline.to_string(index=False))
    print("\nDistribution of unique changes per hotspot position:")
    print(dist.to_string(index=False))
    print("\nWrote", len(list(RESULTS.glob("*.tsv"))), "tables to results/")


if __name__ == "__main__":
    main()
