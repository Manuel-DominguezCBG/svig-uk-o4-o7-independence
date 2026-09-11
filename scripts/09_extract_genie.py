#!/usr/bin/env python3
"""
09_extract_genie.py
-------------------
Extract, from an AACR Project GENIE release, the patient-level recurrence of every
CancerHotspots amino-acid change in the analysis set, plus the composition of the
GENIE cohort by haematological lineage.

This replaces the COSMIC stand-in (script 04) with the database SVIG-UK names for
O4. Three things are counted per change, because each answers a different question:

  * unique patients (the O4 count, de-duplicated as the guideline requires:
    "multiple samples from the same patient are only counted as a single entry");
  * unique patients excluding MSK-IMPACT (the leave-MSK-out count, which removes
    the cohort shared with CancerHotspots);
  * unique patients by lineage -- myeloid, lymphoid, solid -- because the guideline
    lets O4 be limited to 'on-target' tumour types, and AML is its worked example.

Lineage comes from the OncoTree tissue of each sample's ONCOTREE_CODE, pinned to
data/reference/oncotree_2025_10_03.json. Codes retired from that release (old
glioma codes) fall back to CANCER_TYPE.

GENIE is distributed under a data-use agreement that forbids redistribution. The
release files are read from GENIE_DIR and never copied into the repository; the
per-change counts written to data/interim/ are gitignored. Only aggregate tables
are written to results/.

Inputs:
  $GENIE_DIR/data_mutations_extended.txt
  $GENIE_DIR/data_clinical_sample.txt
  results/04_per_allele_o7_tiers.tsv            (the analysis set)
  data/reference/oncotree_2025_10_03.json

Outputs:
  data/interim/genie_hotspot_counts.tsv          per-change counts (gitignored)
  results/18_genie_cohort_composition.tsv        lineage x centre, samples and patients
  results/18b_genie_haem_cancer_types.tsv        haematological cancer types in detail
  results/18c_genie_isoform_offsets.tsv          per-gene residue offsets, for review
  results/25_genie_haem_named_variants.tsv       the haem drivers CancerHotspots lacks
"""

import json
import os
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GENIE_DIR = Path(os.environ.get("GENIE_DIR", ROOT / "data" / "external" / "genie_v20"))
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "results"
ONCOTREE = ROOT / "data" / "reference" / "oncotree_2025_10_03.json"

KEYS = ["hugo_symbol", "amino_acid_position", "reference_aa", "variant_aa"]
LINEAGES = ["Myeloid", "Lymphoid", "Solid", "Unknown"]

# Single-residue protein changes, including stop-loss ("p.*468Yext*8") and
# synonymous ("p.G60=", which CancerHotspots records as G60G).
HGVSP = re.compile(r"^p\.([A-Z*])(\d+)([A-Z*=])(?:ext.*)?$")
SNV_CLASSES = {"Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent"}

# Search window for isoform offsets, in residues, and the number of missense
# hotspot changes an offset must recover before it is trusted.
MAX_OFFSET = 300
MIN_OFFSET_MATCHES = 3

# Cancer types whose OncoTree code sits under an organ tissue but which are
# haematological neoplasms (e.g. primary testicular lymphoma).
HAEM_CANCER_TYPE_OVERRIDE = {
    "Non-Hodgkin Lymphoma": "Lymphoid",
    "Hodgkin Lymphoma": "Lymphoid",
    "Mature B-Cell Neoplasms": "Lymphoid",
    "Mature T and NK Neoplasms": "Lymphoid",
}


def load_clinical():
    clin = pd.read_csv(GENIE_DIR / "data_clinical_sample.txt", sep="\t",
                       comment="#", dtype=str)
    clin["center"] = clin["SAMPLE_ID"].str.split("-").str[1]
    clin["is_msk"] = clin["center"] == "MSK"

    tissue = {t["code"]: t["tissue"]
              for t in json.loads(ONCOTREE.read_text())["tumorTypes"]}
    clin["oncotree_tissue"] = clin["ONCOTREE_CODE"].map(tissue)

    def lineage(row):
        override = HAEM_CANCER_TYPE_OVERRIDE.get(row["CANCER_TYPE"])
        if override:
            return override
        t = row["oncotree_tissue"]
        if t in ("Myeloid", "Lymphoid"):
            return t
        if isinstance(t, str):
            return "Solid"
        # Code not in the pinned OncoTree release: retired glioma codes are solid;
        # UNKNOWN stays unknown.
        if row["CANCER_TYPE"] in (None, "UNKNOWN") or pd.isna(row["CANCER_TYPE"]):
            return "Unknown"
        return "Solid"

    clin["lineage"] = clin.apply(lineage, axis=1)
    return clin


def cohort_composition(clin):
    rows = []
    for label, sub in [("All centres", clin), ("MSK-IMPACT", clin[clin.is_msk]),
                       ("Excluding MSK-IMPACT", clin[~clin.is_msk])]:
        n_samples, n_patients = len(sub), sub["PATIENT_ID"].nunique()
        for lin in LINEAGES:
            s = sub[sub.lineage == lin]
            rows.append({
                "cohort": label, "lineage": lin,
                "samples": len(s),
                "pct_samples": round(100 * len(s) / n_samples, 2),
                "patients": s["PATIENT_ID"].nunique(),
                "pct_patients": round(100 * s["PATIENT_ID"].nunique() / n_patients, 2),
            })
        rows.append({"cohort": label, "lineage": "Total", "samples": n_samples,
                     "pct_samples": 100.0, "patients": n_patients, "pct_patients": 100.0})
    comp = pd.DataFrame(rows)

    multi = clin.groupby("PATIENT_ID")["lineage"].nunique()
    comp.attrs["patients_in_more_than_one_lineage"] = int((multi > 1).sum())

    haem = (clin[clin.lineage.isin(["Myeloid", "Lymphoid"])]
            .groupby(["lineage", "CANCER_TYPE"])
            .agg(samples=("SAMPLE_ID", "size"), patients=("PATIENT_ID", "nunique"))
            .reset_index()
            .sort_values(["lineage", "samples"], ascending=[False, False]))
    return comp, haem


def load_maf():
    maf = pd.read_csv(
        GENIE_DIR / "data_mutations_extended.txt", sep="\t", dtype=str, low_memory=False,
        usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode",
                 "HGVSp_Short"],
    )
    return maf.rename(columns={"Hugo_Symbol": "hugo_symbol",
                               "Tumor_Sample_Barcode": "SAMPLE_ID"})


def parse_snvs(maf):
    snv = maf[maf["Variant_Classification"].isin(SNV_CLASSES)].copy()
    parts = snv["HGVSp_Short"].fillna("").str.extract(HGVSP)
    snv["reference_aa"], snv["amino_acid_position"], snv["variant_aa"] = (
        parts[0], parts[1], parts[2])
    snv = snv.dropna(subset=["amino_acid_position"])
    synonymous = snv["variant_aa"] == "="
    snv.loc[synonymous, "variant_aa"] = snv.loc[synonymous, "reference_aa"]
    snv["amino_acid_position"] = snv["amino_acid_position"].astype(int)
    return snv


def isoform_offsets(snv, analysis_set):
    """Per gene, the residue offset that aligns GENIE's annotation isoform to the
    one CancerHotspots used.

    GENIE annotates a few genes on a different transcript from CancerHotspots, so
    the same residue carries a different number (EZH2 p.Y646 is p.Y641 in GENIE;
    FGFR1 p.N577K is p.N546K). For each gene, every offset in +/-MAX_OFFSET is
    tried and scored by how many of the gene's missense hotspot changes it
    recovers, with reference AND variant amino acid required to match. Synonymous
    and nonsense changes are left out of the scoring: (S,S) or (R,*) recur at so
    many residues that they align by chance.

    An offset is adopted only when it recovers at least MIN_OFFSET_MATCHES missense
    changes and at least twice as many as no offset. Across ~230 genes and 600
    candidate offsets, a chance alignment of two changes is expected a few times;
    of three, essentially never -- which is why two is not enough (it admitted
    CDKN1A +146 on a 164-residue protein). Every candidate is written out for review.
    """
    have = {}
    for gene, pos, ref, var in zip(snv["hugo_symbol"], snv["amino_acid_position"],
                                   snv["reference_aa"], snv["variant_aa"]):
        have.setdefault(gene, set()).add((pos, ref, var))

    rows = []
    for gene, sub in analysis_set.groupby("hugo_symbol"):
        missense = sub[(sub["reference_aa"] != sub["variant_aa"])
                       & (sub["reference_aa"] != "*") & (sub["variant_aa"] != "*")]
        changes = list(zip(missense["amino_acid_position"].astype(int),
                           missense["reference_aa"], missense["variant_aa"]))
        in_genie = have.get(gene, set())

        def matched(delta):
            return sum((p + delta, r, v) in in_genie for p, r, v in changes)

        at_zero = matched(0)
        best_delta, best = max(((d, matched(d)) for d in range(-MAX_OFFSET, MAX_OFFSET + 1)),
                               key=lambda x: (x[1], -abs(x[0])))
        adopt = (best_delta != 0 and best >= MIN_OFFSET_MATCHES
                 and best >= 2 * max(at_zero, 1))
        rows.append({"hugo_symbol": gene, "missense_hotspot_changes": len(changes),
                     "matched_at_offset_0": at_zero, "best_offset": best_delta,
                     "matched_at_best_offset": best,
                     "genie_position_offset": best_delta if adopt else 0})
    return pd.DataFrame(rows)


def hotspot_counts(maf, clin, analysis_set):
    snv = parse_snvs(maf)
    offsets = isoform_offsets(snv, analysis_set)

    # Look each hotspot change up at its GENIE residue number; report it under the
    # CancerHotspots number.
    lookup = analysis_set[KEYS].merge(offsets[["hugo_symbol", "genie_position_offset"]],
                                      on="hugo_symbol", how="left")
    lookup["genie_position"] = (lookup["amino_acid_position"].astype(int)
                                + lookup["genie_position_offset"])
    hits = snv.rename(columns={"amino_acid_position": "genie_position"}).merge(
        lookup, on=["hugo_symbol", "genie_position", "reference_aa", "variant_aa"],
        how="inner")
    hits = hits.merge(clin[["SAMPLE_ID", "PATIENT_ID", "is_msk", "lineage"]],
                      on="SAMPLE_ID", how="left")
    assert hits["PATIENT_ID"].notna().all(), "MAF sample missing from clinical file"

    g = hits.groupby(KEYS)
    counts = pd.DataFrame({
        "genie_rows": g.size(),
        "genie_samples": g["SAMPLE_ID"].nunique(),
        "genie_patients": g["PATIENT_ID"].nunique(),
        "genie_patients_msk": hits[hits.is_msk].groupby(KEYS)["PATIENT_ID"].nunique(),
        "genie_patients_non_msk": hits[~hits.is_msk].groupby(KEYS)["PATIENT_ID"].nunique(),
    })
    for lin in LINEAGES:
        counts[f"genie_patients_{lin.lower()}"] = (
            hits[hits.lineage == lin].groupby(KEYS)["PATIENT_ID"].nunique())
    counts = counts.fillna(0).astype(int).reset_index()

    out = lookup[KEYS + ["genie_position_offset"]].merge(counts, on=KEYS, how="left")
    count_cols = [c for c in out.columns if c.startswith("genie_") and c != "genie_position_offset"]
    out[count_cols] = out[count_cols].fillna(0).astype(int)

    # A zero is only meaningful if GENIE panels sequence the gene at all.
    out["gene_in_genie"] = out["hugo_symbol"].isin(set(maf["hugo_symbol"]))
    return out, offsets


def haem_named_variants(maf, clin):
    """The haematological drivers CancerHotspots v2 omits (results/10), counted in
    GENIE. Indel drivers are matched on gene, class and protein position."""
    m = maf.merge(clin[["SAMPLE_ID", "PATIENT_ID", "lineage", "is_msk"]],
                  on="SAMPLE_ID", how="left")
    hg = m["HGVSp_Short"].fillna("")
    pos = pd.to_numeric(hg.str.extract(r"^p\.[A-Z*](\d+)")[0], errors="coerce")
    fs = m["Variant_Classification"].isin(["Frame_Shift_Del", "Frame_Shift_Ins"])

    def exact(gene, *changes):
        return (m.hugo_symbol == gene) & hg.isin([f"p.{c}" for c in changes])

    specs = [
        ("JAK2", "p.V617F", exact("JAK2", "V617F")),
        ("MPL", "p.W515L/K", exact("MPL", "W515L", "W515K")),
        ("CALR", "exon 9 frameshift (aa 352-417)",
         (m.hugo_symbol == "CALR") & fs & pos.between(352, 417)),
        ("NPM1", "exon 12 frameshift (aa 280-298)",
         (m.hugo_symbol == "NPM1") & fs & pos.between(280, 298)),
        ("ASXL1", "p.G646Wfs*12", (m.hugo_symbol == "ASXL1") & (hg == "p.G646Wfs*12")),
        ("SETBP1", "p.D868N / p.G870S", exact("SETBP1", "D868N", "G870S")),
        ("CBL", "p.Y371 (any missense)",
         (m.hugo_symbol == "CBL") & (pos == 371)
         & (m.Variant_Classification == "Missense_Mutation")),
        ("SF3B1", "p.K700E", exact("SF3B1", "K700E")),
        ("SRSF2", "p.P95H/L/R", exact("SRSF2", "P95H", "P95L", "P95R")),
        ("U2AF1", "p.S34F/Y", exact("U2AF1", "S34F", "S34Y")),
        ("DNMT3A", "p.R882H/C", exact("DNMT3A", "R882H", "R882C")),
        ("IDH2", "p.R140Q", exact("IDH2", "R140Q")),
        ("IDH2", "p.R172K", exact("IDH2", "R172K")),
        ("FLT3", "p.D835 (any missense)",
         (m.hugo_symbol == "FLT3") & (pos == 835)
         & (m.Variant_Classification == "Missense_Mutation")),
        ("KIT", "p.D816V", exact("KIT", "D816V")),
        ("MYD88", "p.L265P", exact("MYD88", "L265P")),
    ]

    rows = []
    for gene, label, mask in specs:
        s = m[mask]
        row = {"hugo_symbol": gene, "variant": label,
               "genie_samples": s["SAMPLE_ID"].nunique(),
               "genie_patients": s["PATIENT_ID"].nunique(),
               "genie_patients_non_msk": s.loc[~s.is_msk, "PATIENT_ID"].nunique()}
        for lin in LINEAGES:
            row[f"genie_patients_{lin.lower()}"] = s.loc[s.lineage == lin, "PATIENT_ID"].nunique()
        haem = row["genie_patients_myeloid"] + row["genie_patients_lymphoid"]
        row["haem_share_of_patients"] = (round(haem / row["genie_patients"], 3)
                                         if row["genie_patients"] else None)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    analysis_set = pd.read_csv(RESULTS / "04_per_allele_o7_tiers.tsv", sep="\t", dtype=str)

    clin = load_clinical()
    comp, haem_types = cohort_composition(clin)
    comp.to_csv(RESULTS / "18_genie_cohort_composition.tsv", sep="\t", index=False)
    haem_types.to_csv(RESULTS / "18b_genie_haem_cancer_types.tsv", sep="\t", index=False)

    maf = load_maf()
    counts, offsets = hotspot_counts(maf, clin, analysis_set)
    counts.to_csv(INTERIM / "genie_hotspot_counts.tsv", sep="\t", index=False)
    (offsets[offsets["best_offset"] != 0]
     .sort_values(["genie_position_offset", "matched_at_best_offset"], ascending=[True, False])
     .to_csv(RESULTS / "18c_genie_isoform_offsets.tsv", sep="\t", index=False))

    haem_named = haem_named_variants(maf, clin)
    haem_named.to_csv(RESULTS / "25_genie_haem_named_variants.tsv", sep="\t", index=False)

    print(f"GENIE: {len(clin):,} samples, {clin.PATIENT_ID.nunique():,} patients; "
          f"MAF {len(maf):,} rows")
    print(f"Patients with samples in more than one lineage: "
          f"{comp.attrs['patients_in_more_than_one_lineage']:,}")
    print(comp.to_string(index=False))
    found = (counts["genie_patients"] > 0).sum()
    print(f"\nHotspot changes found in GENIE: {found:,} / {len(counts):,} "
          f"({100 * found / len(counts):.1f}%); genes absent from GENIE: "
          f"{(~counts.gene_in_genie).sum()} changes")
    adopted = offsets[offsets["genie_position_offset"] != 0]
    print(f"\nIsoform offsets adopted for {len(adopted)} genes:")
    print(adopted.to_string(index=False))
    rejected = offsets[(offsets["best_offset"] != 0) & (offsets["genie_position_offset"] == 0)
                       & (offsets["matched_at_best_offset"] > offsets["matched_at_offset_0"])]
    print(f"\nCandidate offsets rejected (improvement too small to trust): {len(rejected)}")
    print(rejected.to_string(index=False))
    print()
    print(haem_named.to_string(index=False))


if __name__ == "__main__":
    main()
