#!/usr/bin/env python3
"""
07_validate_against_api_export.py
---------------------------------
Cross-check the parsed v2 workbook against an independent cancerhotspots.org API
export held in the local reference set. The workbook is supplementary material to
the publication; the API export is what the website serves. If an analyst applying
O7 at the bench sees different numbers from the ones modelled here, the whole
analysis would need re-basing, so this is checked rather than assumed.

Input (outside the repository):
  $CANCERHOTSPOTS_JSON, default
  /Users/monkiky/Documents/external/refs/CancerHotSpots/cancerhotspots_counts.json

Output:
  results/17_api_export_validation.tsv
"""

import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"

DEFAULT_JSON = ("/Users/monkiky/Documents/external/refs/CancerHotSpots/"
                "cancerhotspots_counts.json")


def main():
    path = os.environ.get("CANCERHOTSPOTS_JSON", DEFAULT_JSON)
    if not Path(path).exists():
        sys.exit(f"API export not found: {path}\nSet CANCERHOTSPOTS_JSON to its location.")

    records = json.load(open(path))
    api_rows = []
    for r in records:
        for variant_aa, count in (r.get("variantAminoAcid") or {}).items():
            api_rows.append({
                "hugo_symbol": r["hugoSymbol"],
                "residue": r["residue"],
                "variant_aa": variant_aa,
                "api_change_count": int(count),
                "api_position_total": int(r["tumorCount"]),
                "api_type": r["type"],
            })
    api = pd.DataFrame(api_rows)

    merged = pd.read_csv(INTERIM / "hotspots_merged_per_allele.tsv", sep="\t")
    snv = merged[(merged["source"] == "v2:SNV-hotspots")].copy()
    # The API keys positions as the residue string, e.g. "Q61" or "X187_splice".
    # The API drops the "_splice" suffix from the residue (X187_splice -> X187)
    # and denotes the splice consequence as "sp" rather than "splice".
    snv["residue"] = [
        str(pos).replace("_splice", "") if str(pos).endswith("_splice") else f"{ref}{pos}"
        for ref, pos in zip(snv["reference_aa"], snv["amino_acid_position"])
    ]
    snv["variant_aa"] = snv["variant_aa"].replace({"splice": "sp"})

    keys = ["hugo_symbol", "residue", "variant_aa"]
    api_single = api[api["api_type"] == "single residue"]
    joined = snv.merge(api_single, on=keys, how="outer", indicator=True)

    both = joined[joined["_merge"] == "both"]
    count_matches = (both["change_count"] == both["api_change_count"]).sum()
    total_matches = (both["position_total_count"] == both["api_position_total"]).sum()

    report = pd.DataFrame([
        ("API export records (positions/regions)", len(records)),
        ("API single-residue amino-acid changes", len(api_single)),
        ("Workbook v2 SNV amino-acid changes", len(snv)),
        ("Changes present in both", len(both)),
        ("Changes only in the workbook", int((joined["_merge"] == "left_only").sum())),
        ("Changes only in the API export", int((joined["_merge"] == "right_only").sum())),
        ("Per-change counts identical", int(count_matches)),
        ("Position totals identical", int(total_matches)),
        ("Per-change agreement (%)", round(100 * count_matches / len(both), 2)),
    ], columns=["metric", "value"])

    RESULTS.mkdir(parents=True, exist_ok=True)
    report.to_csv(RESULTS / "17_api_export_validation.tsv", sep="\t", index=False)
    print(report.to_string(index=False))

    mismatched = both[both["change_count"] != both["api_change_count"]]
    if len(mismatched):
        print("\nMismatched counts (first 20):")
        print(mismatched[keys + ["change_count", "api_change_count"]]
              .head(20).to_string(index=False))


if __name__ == "__main__":
    main()
