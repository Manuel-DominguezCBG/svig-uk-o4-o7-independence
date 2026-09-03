#!/usr/bin/env python3
"""
04_extract_cosmic_cmc.py
------------------------
Extract a protein-level recurrence table from the COSMIC Cancer Mutation Census
(CMC v104), restricted to the genes that appear in CancerHotspots.

Why: the analysis needs a *second, independent* recurrence resource to stand in
for GENIE, which was not available. SVIG-UK names COSMIC as an acceptable
alternative for O4:

    "Other databases such as the Catalogue Of Somatic Mutations In Cancer
     (COSMIC) are also available and can be used to assess enrichment"

COSMIC is used here to measure *how often O4 and O7 fire on the same variant*,
not to make classification calls. See docs/findings.md for the caveats.

Input (not in this repository; large reference file held outside it):
  $COSMIC_CMC_TAR, default
  /Users/monkiky/Documents/external/refs/COSMIC/CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar

Output:
  data/interim/cosmic_cmc_hotspot_genes.tsv   one row per (gene, AA position,
                                              ref AA, alt AA) with COSMIC sample counts
"""

import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
INTERIM = ROOT / "data" / "interim"

DEFAULT_TAR = ("/Users/monkiky/Documents/external/refs/COSMIC/"
               "CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar")
MEMBER = "CancerMutationCensus_AllData_v104_GRCh37.tsv.gz"

# 1-based column positions in the CMC export (see the bundled README).
COL = {
    "gene": 1, "onc_tsg": 3, "cgc_tier": 4, "mutation_aa": 8,
    "aa_start": 9, "shared_aa": 11, "aa_wt": 14, "aa_mut": 15,
    "desc_aa": 17, "tested": 22, "mutated": 23, "tier": 58,
}

SUBSTITUTION = re.compile(r"^p\.([A-Z*])(\d+)([A-Z*])$")


def hotspot_genes():
    table = pd.read_csv(RESULTS / "03_gene_position_table.tsv", sep="\t")
    return set(table["hugo_symbol"].astype(str))


def stream_cmc(tar_path):
    """Yield decoded lines of the CMC TSV without materialising it on disk."""
    proc = subprocess.Popen(
        ["tar", "-xOf", tar_path, MEMBER],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    with gzip.open(proc.stdout, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            yield line
    proc.stdout.close()
    proc.wait()


def main():
    tar_path = os.environ.get("COSMIC_CMC_TAR", DEFAULT_TAR)
    if not Path(tar_path).exists():
        sys.exit(f"COSMIC CMC archive not found: {tar_path}\n"
                 f"Set COSMIC_CMC_TAR to its location.")

    genes = hotspot_genes()
    idx = {k: v - 1 for k, v in COL.items()}

    # (gene, position, ref AA, alt AA) -> aggregated COSMIC evidence.
    # Distinct codon changes producing the same substitution are summed, which
    # is the protein-level view SVIG-UK works in.
    agg = defaultdict(lambda: {"mutated": 0, "tested": 0, "n_cds_variants": 0,
                               "tiers": set(), "onc_tsg": "", "cgc_tier": ""})

    total = kept = 0
    for n, line in enumerate(stream_cmc(tar_path)):
        if n == 0:
            continue
        total += 1
        f = line.rstrip("\n").split("\t")
        if len(f) <= idx["tier"]:
            continue
        gene = f[idx["gene"]]
        if gene not in genes:
            continue
        m = SUBSTITUTION.match(f[idx["mutation_aa"]])
        if not m:                       # only single-residue substitutions join to hotspots
            continue
        ref_aa, position, alt_aa = m.group(1), m.group(2), m.group(3)
        try:
            mutated = int(f[idx["mutated"]] or 0)
            tested = int(f[idx["tested"]] or 0)
        except ValueError:
            continue
        rec = agg[(gene, position, ref_aa, alt_aa)]
        rec["mutated"] += mutated
        rec["tested"] = max(rec["tested"], tested)
        rec["n_cds_variants"] += 1
        rec["tiers"].add(f[idx["tier"]])
        rec["onc_tsg"] = f[idx["onc_tsg"]]
        rec["cgc_tier"] = f[idx["cgc_tier"]]
        kept += 1
        if total % 2_000_000 == 0:
            print(f"  scanned {total:,} rows, kept {kept:,}", flush=True)

    rows = [
        {
            "hugo_symbol": gene,
            "amino_acid_position": position,
            "reference_aa": ref_aa,
            "variant_aa": alt_aa,
            "cosmic_samples_mutated": rec["mutated"],
            "cosmic_samples_tested": rec["tested"],
            "cosmic_n_cds_variants": rec["n_cds_variants"],
            "cosmic_significance_tier": min(rec["tiers"], key=str),
            "cosmic_onc_tsg": rec["onc_tsg"],
            "cosmic_cgc_tier": rec["cgc_tier"],
        }
        for (gene, position, ref_aa, alt_aa), rec in agg.items()
    ]
    out = pd.DataFrame(rows).sort_values("cosmic_samples_mutated", ascending=False)

    INTERIM.mkdir(parents=True, exist_ok=True)
    path = INTERIM / "cosmic_cmc_hotspot_genes.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(f"scanned {total:,} CMC rows; wrote {len(out):,} protein-level "
          f"substitutions over {out['hugo_symbol'].nunique()} genes to {path}")


if __name__ == "__main__":
    main()
