#!/usr/bin/env python3
"""
05_o4_proxy_crosstab.py
-----------------------
Use COSMIC CMC recurrence counts as a stand-in for GENIE (which was not
available) to measure, rather than merely argue, how strongly O4 and O7 fire on
the same variants.

Two things are quantified:

  1. Correlation between the per-change count that drives O7 (CancerHotspots)
     and an independent per-change count that would drive O4 (COSMIC CMC).
     These are largely different cohorts, so a high correlation shows the
     dependence is structural, not only a shared-cohort artefact.

  2. The O4 x O7 cross-tabulation: how often the +8 combination actually arises,
     and what the proposed cap would cost in points and in classification band.

SVIG-UK explicitly permits COSMIC for O4 but warns that thresholds "should be
much higher if COSMIC is being used compared to GENIE" without specifying them,
so two threshold schemes are reported side by side.

Outputs:
  results/11_cosmic_o4_o7_crosstab.tsv
  results/12_o4_o7_correlation.tsv
  results/13_points_impact_of_cap.tsv
  results/14_per_change_o4_o7_points.tsv
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"

# SVIG-UK Figure 1C: classification categories and exponent sums.
BANDS = [(10, "Oncogenic"), (6, "Likely oncogenic"), (0, "VUS"),
         (-6, "Likely benign"), (-99, "Benign")]

# O4 missense thresholds. "genie_literal" is Supplementary Table 1 as written;
# "cosmic_adjusted" applies the guideline's instruction to raise thresholds for
# COSMIC, reusing the ratio it already uses for truncating variants (5x).
O4_SCHEMES = {
    "genie_literal": {"strong": 11, "moderate": 5},
    "cosmic_adjusted": {"strong": 51, "moderate": 20},
}


def o4_tier(count, scheme):
    if pd.isna(count):
        return "Absent/not counted", 0
    if count >= scheme["strong"]:
        return "Strong", 4
    if count >= scheme["moderate"]:
        return "Moderate", 2
    return "Not met", 0


def band(points):
    for threshold, name in BANDS:
        if points >= threshold:
            return name
    return "Benign"


def main():
    o7 = pd.read_csv(RESULTS / "05_o4_o7_double_counting.tsv", sep="\t")
    cosmic = pd.read_csv(INTERIM / "cosmic_cmc_hotspot_genes.tsv", sep="\t")
    for df in (o7, cosmic):
        df["amino_acid_position"] = df["amino_acid_position"].astype(str)

    keys = ["hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa"]
    df = o7.merge(
        cosmic[keys + ["cosmic_samples_mutated", "cosmic_significance_tier",
                       "cosmic_onc_tsg"]],
        on=keys, how="left",
    )

    # ---- 12 correlation -----------------------------------------------------
    joined = df[df["cosmic_samples_mutated"].notna()]
    corr = pd.DataFrame([
        ("Amino-acid changes in the CancerHotspots analysis set", len(df)),
        ("Matched to a COSMIC CMC protein-level substitution", len(joined)),
        ("Match rate (%)", round(100 * len(joined) / len(df), 1)),
        ("Spearman rho: CancerHotspots count vs COSMIC count",
         round(joined["change_count"].corr(joined["cosmic_samples_mutated"],
                                           method="spearman"), 3)),
        ("Pearson r on log10 counts",
         round(np.log10(joined["change_count"]).corr(
             np.log10(joined["cosmic_samples_mutated"].clip(lower=1))), 3)),
        ("Spearman rho, restricted to changes reaching O7 moderate or strong",
         round(joined.loc[joined["o7_points"] >= 2, "change_count"].corr(
             joined.loc[joined["o7_points"] >= 2, "cosmic_samples_mutated"],
             method="spearman"), 3)),
    ], columns=["metric", "value"])

    # ---- 11 cross-tabulation and 13/14 points impact ------------------------
    crosstabs, impacts = [], []
    for scheme_name, scheme in O4_SCHEMES.items():
        tiers = [o4_tier(c, scheme) for c in df["cosmic_samples_mutated"]]
        df[f"o4_strength__{scheme_name}"] = [t[0] for t in tiers]
        df[f"o4_points__{scheme_name}"] = [t[1] for t in tiers]

        ct = (
            pd.crosstab(df[f"o4_strength__{scheme_name}"], df["o7_strength"])
            .reindex(index=["Strong", "Moderate", "Not met", "Absent/not counted"],
                     columns=["Strong", "Moderate", "Supporting", "Not met"])
            .fillna(0).astype(int)
        )
        ct = ct.reset_index().rename(columns={f"o4_strength__{scheme_name}": "o4_strength"})
        ct.insert(0, "o4_threshold_scheme", scheme_name)
        ct.columns = [c if c in ("o4_threshold_scheme", "o4_strength") else f"o7_{c}"
                      for c in ct.columns]
        crosstabs.append(ct)

        uncapped = df[f"o4_points__{scheme_name}"] + df["o7_points"]
        capped = uncapped.clip(upper=4).where(~df["independent_positional_evidence"], uncapped)
        df[f"o4_o7_uncapped__{scheme_name}"] = uncapped
        df[f"o4_o7_capped__{scheme_name}"] = capped

        affected = uncapped != capped
        impacts.append({
            "o4_threshold_scheme": scheme_name,
            "changes_assessed": len(df),
            "changes_scoring_both_o4_and_o7": int(((df[f"o4_points__{scheme_name}"] > 0)
                                                   & (df["o7_points"] > 0)).sum()),
            "changes_reaching_o4_strong_and_o7_strong": int(
                ((df[f"o4_strength__{scheme_name}"] == "Strong")
                 & (df["o7_strength"] == "Strong")).sum()),
            "changes_scoring_the_full_8_points": int((uncapped >= 8).sum()),
            "changes_reduced_by_the_cap": int(affected.sum()),
            "pct_reduced_by_the_cap": round(100 * affected.mean(), 1),
            "mean_points_lost_where_capped": round((uncapped - capped)[affected].mean(), 2)
            if affected.any() else 0.0,
            "changes_dropping_from_8_to_4": int(((uncapped == 8) & (capped == 4)).sum()),
        })

    pd.concat(crosstabs).to_csv(RESULTS / "11_cosmic_o4_o7_crosstab.tsv",
                                sep="\t", index=False)
    corr.to_csv(RESULTS / "12_o4_o7_correlation.tsv", sep="\t", index=False)
    pd.DataFrame(impacts).to_csv(RESULTS / "13_points_impact_of_cap.tsv",
                                 sep="\t", index=False)

    out_cols = keys + [
        "change_count", "position_total_count", "n_unique_changes",
        "cosmic_samples_mutated", "cosmic_significance_tier", "cosmic_onc_tsg",
        "o7_strength", "o7_points",
        "o4_strength__genie_literal", "o4_points__genie_literal",
        "o4_o7_uncapped__genie_literal", "o4_o7_capped__genie_literal",
        "o4_strength__cosmic_adjusted", "o4_points__cosmic_adjusted",
        "o4_o7_uncapped__cosmic_adjusted", "o4_o7_capped__cosmic_adjusted",
        "independent_positional_evidence", "recommended_o4_o7_handling",
        "msk_fraction",
    ]
    (df[out_cols].sort_values("change_count", ascending=False)
     .to_csv(RESULTS / "14_per_change_o4_o7_points.tsv", sep="\t", index=False))

    print(corr.to_string(index=False))
    print()
    for ct in crosstabs:
        print(ct.to_string(index=False)); print()
    print(pd.DataFrame(impacts).to_string(index=False))


if __name__ == "__main__":
    main()
