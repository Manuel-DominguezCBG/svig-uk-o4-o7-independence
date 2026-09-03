#!/usr/bin/env python3
"""
06_canonical_list_overlap.py
----------------------------
Cross-reference the analysis against the SVIG-UK Canonical Variants List (O1,
Supplementary Table 3).

This matters for the cap proposal: O1 is standalone evidence, so a variant on
that list is classified Oncogenic irrespective of what O4 and O7 contribute. Any
variant that would lose points under the proposed cap but is already on the
canonical list is unaffected in practice - which bounds the clinical risk of
adopting the cap.

Inputs:
  data/raw/svig_uk_canonical_variants.tsv   (gene, transcript, HGVSp_Short, assessment)
  results/14_per_change_o4_o7_points.tsv

Outputs:
  results/15_canonical_list_overlap.tsv       every hotspot change, with O1 status
  results/16_variants_materially_affected.tsv changes losing the full +4, ranked
"""

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

SUBSTITUTION = re.compile(r"^p\.([A-Z*])(\d+)([A-Z*])$")
KEYS = ["hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa"]


def load_canonical():
    df = pd.read_csv(RAW / "svig_uk_canonical_variants.tsv", sep="\t")
    parsed = df["HGVSp_Short"].astype(str).apply(
        lambda s: SUBSTITUTION.match(s).groups() if SUBSTITUTION.match(s) else (None,) * 3
    )
    df[["reference_aa", "amino_acid_position", "variant_aa"]] = pd.DataFrame(
        parsed.tolist(), index=df.index
    )
    df = df.dropna(subset=["amino_acid_position"]).rename(columns={"gene": "hugo_symbol"})
    return df[KEYS + ["transcript", "HGVSp_Short", "svig_uk_assessment"]]


def main():
    changes = pd.read_csv(RESULTS / "14_per_change_o4_o7_points.tsv", sep="\t")
    changes["amino_acid_position"] = changes["amino_acid_position"].astype(str)
    canonical = load_canonical()

    merged = changes.merge(canonical, on=KEYS, how="left")
    merged["on_svig_uk_canonical_list"] = merged["svig_uk_assessment"].notna()
    merged.to_csv(RESULTS / "15_canonical_list_overlap.tsv", sep="\t", index=False)

    # Changes where the cap removes the full 4 points, i.e. +8 becomes +4.
    affected = merged[
        (merged["o4_o7_uncapped__genie_literal"] == 8)
        & (merged["o4_o7_capped__genie_literal"] == 4)
    ].copy()
    affected["protected_by_O1"] = affected["on_svig_uk_canonical_list"]
    cols = KEYS + [
        "change_count", "position_total_count", "n_unique_changes",
        "cosmic_samples_mutated", "msk_fraction", "o7_strength",
        "o4_strength__genie_literal", "svig_uk_assessment", "protected_by_O1",
    ]
    affected = affected[cols].sort_values("change_count", ascending=False)
    affected.to_csv(RESULTS / "16_variants_materially_affected.tsv", sep="\t", index=False)

    n_canon = int(merged["on_svig_uk_canonical_list"].sum())
    print(f"canonical (O1) substitutions parsed: {len(canonical)}")
    print(f"matched into the hotspot analysis set: {n_canon}")
    print(f"\nchanges dropping from +8 to +4 under the cap: {len(affected)}")
    print(f"  of which already on the O1 canonical list (cap is moot): "
          f"{int(affected['protected_by_O1'].sum())}")
    print(f"  materially affected: {int((~affected['protected_by_O1']).sum())}\n")
    print(affected.to_string(index=False))


if __name__ == "__main__":
    main()
