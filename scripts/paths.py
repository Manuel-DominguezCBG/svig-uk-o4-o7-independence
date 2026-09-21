"""
Where each analysis reads its intermediate files and writes its results.

The pipeline can be run on more than one CancerHotspots release without one run
overwriting another. Set ANALYSIS before running any script:

  (unset) or v2   the original analysis: CancerHotspots v2 (Chang 2018)
                  -> data/interim/, results/, figures/  -- exactly as before
  v3              v2 plus the 528 changes new in CancerHotspots v3
                  -> analyses/v3/interim/, analyses/v3/results/, analyses/v3/figures/

Source data (data/raw, data/reference, data/external) is shared by all analyses.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

ANALYSIS = os.environ.get("ANALYSIS", "v2")
if ANALYSIS == "v2":
    INTERIM = ROOT / "data" / "interim"
    RESULTS = ROOT / "results"
    FIGURES = ROOT / "figures"
elif ANALYSIS == "v3":
    INTERIM = ROOT / "analyses" / "v3" / "interim"
    RESULTS = ROOT / "analyses" / "v3" / "results"
    FIGURES = ROOT / "analyses" / "v3" / "figures"
else:
    raise SystemExit(f"ANALYSIS must be 'v2' or 'v3', not {ANALYSIS!r}")

# CancerHotspots sources included in the SNV analysis set.
SNV_SOURCES = ["v2:SNV-hotspots"] + (["v3:SNV-hotspots"] if ANALYSIS == "v3" else [])
