#!/usr/bin/env python3
"""
01_merge_hotspots.py
--------------------
Parse and merge the CancerHotspots source workbooks into a single tidy
"per-allele" table (one row per gene / amino-acid position / amino-acid change).

Inputs (data/raw/):
  cancerhotspots_v1_chang2016.xls   sheets: "Per Residue", "Per Allele"
  cancerhotspots_v2_chang2018.xls   sheets: "SNV-hotspots", "INDEL-hotspots"

Output:
  data/interim/hotspots_merged_per_allele.tsv

Notes on the source encoding
  * v2 'Reference_Amino_Acid' is "<refAA>:<total mutations at position>"  e.g. "Q:422"
  * v2 'Variant_Amino_Acid'   is "<varAA>:<mutations for this change>"    e.g. "R:204"
  * v1 'Variant Amino Acid'   is a pipe list "<varAA>:<count>|..."        e.g. "E:520|K:33"
  These mutation-count annotations are stripped into explicit numeric columns.

CancerHotspots data are made available under the ODC Open Database License (ODbL).
"""

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"

V1 = RAW / "cancerhotspots_v1_chang2016.xls"
V2 = RAW / "cancerhotspots_v2_chang2018.xls"

# Columns of the harmonised output, in order.
OUT_COLS = [
    "source",              # which workbook/sheet the row came from
    "variant_class",       # snv | indel | splice
    "hugo_symbol",
    "amino_acid_position",  # kept as text: v2 indels use ranges e.g. "27-42"
    "reference_aa",
    "variant_aa",
    "change_count",         # mutations reported for this exact change
    "position_total_count",  # mutations reported at this position (all changes)
    "qvalue",
    "n_msk",                # MSK-IMPACT samples contributing (v2 only)
    "n_retro",              # retrospective/public samples contributing (v2 only)
    "total_samples",
    "tumour_type_composition",
]


def split_annotated(token):
    """'Q:422' -> ('Q', 422). Returns (value, None) when no count is annotated."""
    if pd.isna(token):
        return None, None
    token = str(token).strip()
    if ":" not in token:
        return token, None
    value, _, count = token.rpartition(":")
    try:
        return value, int(count)
    except ValueError:
        return token, None


def parse_pipe_counts(field):
    """'E:520|K:33|R:4' -> [('E', 520), ('K', 33), ('R', 4)]."""
    if pd.isna(field):
        return []
    out = []
    for token in str(field).split("|"):
        token = token.strip()
        if not token:
            continue
        value, count = split_annotated(token)
        out.append((value, count))
    return out


def classify(position, variant_aa, sheet):
    """Assign a coarse variant class used later to gate SVIG-UK O7 eligibility."""
    if "indel" in sheet.lower():
        return "indel"
    if "splice" in str(position).lower() or str(variant_aa).lower() == "splice":
        return "splice"
    return "snv"


def load_v2(sheet):
    df = pd.read_excel(V2, sheet_name=sheet)
    rows = []
    for r in df.itertuples(index=False):
        ref_aa, pos_total = split_annotated(getattr(r, "Reference_Amino_Acid"))
        var_aa, change_count = split_annotated(getattr(r, "Variant_Amino_Acid"))
        position = str(getattr(r, "Amino_Acid_Position")).strip()
        # Indel rows carry a pipe-list of reference peptides; keep the raw text.
        if "indel" in sheet.lower():
            ref_aa = str(getattr(r, "Reference_Amino_Acid"))
            pos_total = None
        rows.append(
            {
                "source": f"v2:{sheet}",
                "variant_class": classify(position, var_aa, sheet),
                "hugo_symbol": getattr(r, "Hugo_Symbol"),
                "amino_acid_position": position,
                "reference_aa": ref_aa,
                "variant_aa": var_aa,
                "change_count": change_count,
                "position_total_count": pos_total
                if pos_total is not None
                else getattr(r, "Mutation_Count"),
                "qvalue": getattr(r, "qvalue"),
                "n_msk": getattr(r, "n_MSK"),
                "n_retro": getattr(r, "n_Retro"),
                "total_samples": getattr(r, "Total_Samples"),
                "tumour_type_composition": getattr(r, "Detailed_Cancer_Types"),
            }
        )
    return pd.DataFrame(rows)


def load_v1_per_residue():
    """v1 'Per Residue' holds the per-allele counts inside 'Variant Amino Acid'."""
    df = pd.read_excel(V1, sheet_name="Per Residue", header=2)
    df = df.rename(
        columns={
            "Hugo Symbol": "gene",
            "Codon": "codon",
            "Variant Amino Acid": "variants",
            "Q-value": "qvalue",
            "Tumor Count": "tumour_count",
            "Tumor Type Composition": "composition",
        }
    )
    df = df[df["gene"].notna()]
    rows = []
    for _, r in df.iterrows():
        codon = str(r["codon"]).strip()               # e.g. "V600"
        m = re.match(r"^([A-Z*]+)(\d+)$", codon)
        ref_aa, position = (m.group(1), m.group(2)) if m else (None, codon)
        pos_total = r["tumour_count"]
        for var_aa, count in parse_pipe_counts(r["variants"]):
            rows.append(
                {
                    "source": "v1:Per Residue",
                    "variant_class": classify(position, var_aa, "Per Residue"),
                    "hugo_symbol": r["gene"],
                    "amino_acid_position": position,
                    "reference_aa": ref_aa,
                    "variant_aa": var_aa,
                    "change_count": count,
                    "position_total_count": pos_total,
                    "qvalue": r["qvalue"],
                    "n_msk": None,
                    "n_retro": None,
                    "total_samples": None,
                    "tumour_type_composition": r["composition"],
                }
            )
    return pd.DataFrame(rows)


def load_v1_per_allele():
    """v1 'Per Allele' has no explicit count; it is recovered from the tumour-type
    composition string (sum of per-tumour-type counts)."""
    df = pd.read_excel(V1, sheet_name="Per Allele", header=1)
    df = df[df["Hugo_Symbol"].notna()]
    rows = []
    for r in df.itertuples(index=False):
        comp = parse_pipe_counts(getattr(r, "Tumor_Type_Composition"))
        count = sum(c for _, c in comp if c is not None) or None
        position = str(getattr(r, "Amino_Acid_Position")).strip()
        if position.endswith(".0"):
            position = position[:-2]
        rows.append(
            {
                "source": "v1:Per Allele",
                "variant_class": classify(position, getattr(r, "Variant_Amino_Acid"),
                                          "Per Allele"),
                "hugo_symbol": getattr(r, "Hugo_Symbol"),
                "amino_acid_position": position,
                "reference_aa": getattr(r, "Reference_Amino_Acid"),
                "variant_aa": getattr(r, "Variant_Amino_Acid"),
                "change_count": count,
                "position_total_count": None,
                "qvalue": None,
                "n_msk": None,
                "n_retro": None,
                "total_samples": None,
                "tumour_type_composition": getattr(r, "Tumor_Type_Composition"),
            }
        )
    out = pd.DataFrame(rows)
    # Position totals are the sum of the per-allele counts for that position.
    totals = out.groupby(["hugo_symbol", "amino_acid_position"])["change_count"].transform("sum")
    out["position_total_count"] = totals
    return out


def main():
    INTERIM.mkdir(parents=True, exist_ok=True)

    parts = [
        load_v2("SNV-hotspots"),
        load_v2("INDEL-hotspots"),
        load_v1_per_residue(),
        load_v1_per_allele(),
    ]
    # Give the count columns a consistent nullable dtype before concatenating:
    # v1 sheets carry no MSK split, so those columns are all-NA in some parts.
    numeric = ["change_count", "position_total_count", "qvalue",
               "n_msk", "n_retro", "total_samples"]
    parts = [p.astype({c: "float64" for c in numeric if c in p}) for p in parts]
    merged = pd.concat(parts, ignore_index=True)[OUT_COLS]

    # Integrity checks that must hold for the downstream analysis to be meaningful.
    v2snv = merged[merged["source"] == "v2:SNV-hotspots"]
    per_pos = v2snv.groupby(["hugo_symbol", "amino_acid_position"]).agg(
        summed=("change_count", "sum"),
        stated=("position_total_count", "first"),
        msk=("n_msk", "first"),
        retro=("n_retro", "first"),
    )
    assert (per_pos["summed"] == per_pos["stated"]).all(), \
        "v2 SNV: per-allele counts do not sum to the stated position total"
    assert (per_pos["msk"] + per_pos["retro"] == per_pos["stated"]).all(), \
        "v2 SNV: n_MSK + n_Retro does not reconstitute the position total"

    out = INTERIM / "hotspots_merged_per_allele.tsv"
    merged.to_csv(out, sep="\t", index=False)

    print(f"wrote {out}  ({len(merged):,} rows)")
    print(merged.groupby(["source", "variant_class"]).size().to_string())
    print(f"\ndistinct (gene, position) across all sources: "
          f"{merged.groupby(['hugo_symbol', 'amino_acid_position']).ngroups:,}")


if __name__ == "__main__":
    main()
