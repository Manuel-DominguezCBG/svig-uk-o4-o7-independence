# Analysis with CancerHotspots v3

A parallel run of the whole pipeline with the **528 changes new in CancerHotspots v3** added
to the v2 analysis set. It lives entirely in this folder: the original v2 analysis
(`results/`, `figures/`, `data/interim/`) is untouched, and was verified byte-for-byte after
the scripts were changed to support this run.

Source: `data/raw/hotspots_v3.xlsx`, the CancerHotspots v3 export (sheets `SNV_Variants`,
`INDEL_Variants`, `Hotspot_Residues`), supplied by Kevin Baker, 21 September 2026.

## How to run it

Every script reads and writes through `scripts/paths.py`. Setting `ANALYSIS=v3` sends all
output here instead of to the main folders:

```bash
for s in 01_merge_hotspots 02_position_analysis 04_extract_cosmic_cmc 05_o4_proxy_crosstab \
         06_canonical_list_overlap 07_validate_against_api_export 09_extract_genie \
         10_genie_o4_o7 12_genie_points_per_change 03_export_workbook 11_figures \
         13_v3_cohort_prevalence; do
  ANALYSIS=v3 .venv/bin/python scripts/$s.py
done
```

Without `ANALYSIS` (or with `ANALYSIS=v2`) every script behaves exactly as before.

## What v3 adds

| | v2 | v2 + v3 |
|---|---|---|
| Hotspot residues (SNV, non-splice) | 1,024 | **1,188** (+164) |
| Amino-acid changes | 2,918 | **3,446** (+528) |
| Residues with a single change | 213 (20.8%) | 226 (19.0%) |
| Changes reaching O7_strong | 198 | 200 |
| Changes reaching O7_moderate | 203 | 261 |
| Changes reaching O7_supporting | 1,248 | 1,509 |

- **The v2 part of the export is identical to the v2 workbook.** All 2,918 v2 changes match,
  with the same counts; script 01 asserts it.
- **The 164 new residues are almost all new.** Only one (TP53 224) was in v2, and there only
  as a splice site, which is outside the analysis.
- **The v3 counts are small.** 2,221 mutations over 528 changes; the largest is TP53 p.G199V
  with 32. So most v3 changes reach O7 supporting (261) or moderate (58); only 2 reach strong,
  and 18 are not applicable (nonsense).
- **The v3 export has no MSK / retrospective split per residue**, unlike v2's `n_MSK` /
  `n_Retro`. The MSK statistics in tables 06 and 01 therefore use the v2 residues only, and
  `msk_fraction` is empty for v3 rows. Every per-residue and per-change table has a
  `hotspot_version` column saying which release each row comes from.
- **One residue's counts don't add up.** PTEN p.H123 lists 14 mutations over its 5 changes
  but 15 in the tumour-type total, probably a mutation of another type. As for v2, the residue
  total used is the sum of its changes.

## O4 on GENIE and the cap

**Adding v3 changes nothing for the v2 changes.** All 2,918 keep the same O7 tier, GENIE
patient count, O4 tier and points, under both criteria. The per-gene isoform offsets for GENIE
are also the same.

**The 528 new changes:**

- **507 are observed in GENIE v20** (all 53 genes are sequenced there). 225 reach O4_strong.
- **The two new O7_strong changes both score the full +8 and keep it under both criteria:**
  TP53 p.G199V (32 of 53 at the residue; residual 21, largest other change 15) and p.G199E
  (15 of 53).
- **So no new change drops from +8 to +4.** The totals stay at 10 (permissive) and 22 (strict),
  with the same 4 and 9 on the O1 list.
- **Where the cap does bite, it takes one or two points.** 278 v3 changes score both O4 and O7.
  The cap applies to 201 of them under the permissive criterion (39 go from +6 to +4, 105 from
  +5 to +4) and to 237 under the strict one.

| | v2 analysis | | v2 + v3 analysis | |
|---|---|---|---|---|
| | Permissive | Strict | Permissive | Strict |
| Changes scoring both O4 and O7 | 1,522 | 1,522 | 1,800 | 1,800 |
| Changes whose score falls | 671 (23.0%) | 822 (28.2%) | 815 (23.7%) | 992 (28.8%) |
| **Changes dropping from +8 to +4** | **10** | **22** | **10** | **22** |
| … of which on the O1 list | 4 | 9 | 4 | 9 |
| O7_strong changes also O4_strong (GENIE) | 197 / 198 | | 199 / 200 | |
| … still O4_strong with MSK excluded | 197 / 198 | | 199 / 200 | |

## The cohort breakdown at the v3 residues

For each of its 164 new residues, the v3 export gives mutation counts in three cohorts —
MSK-IMPACT, GENIE excluding MSK, and TCGA — with each cohort's sample count and tests of
whether prevalence differs (script 13, tables 33 and 34).

**These are not a breakdown of the evidence behind each hotspot call**, as `n_MSK` / `n_Retro`
were in v2. They total 6,986 mutations at these residues, against 2,221 in the per-change
counts that O7 is scored on, and the file doesn't say which cohort the v3 calls were made on.
What they show is where these residues are observed:

| Cohort | Mutations at the v3 residues | Share |
|---|---|---|
| MSK-IMPACT (a GENIE centre) | 2,641 | 37.8% |
| GENIE excluding MSK | 3,937 | 56.4% |
| TCGA (not in GENIE) | 408 | 5.8% |

- **94% of the mutations recorded at the v3 residues are in GENIE** (MSK plus the other
  centres), and the median residue is 94% GENIE too. For comparison, in v2 about half the
  hotspot evidence (53%) came from retrospective cohorts, largely TCGA, which is not in GENIE.
- **Prevalence hardly differs between cohorts.** It differs between MSK and non-MSK GENIE at
  only 3 of 164 residues (adjusted p < 0.05), and between MSK and TCGA at 1.
- **The export's non-MSK GENIE counts track our own GENIE v20 extraction closely**
  (Spearman ρ = 0.88 per residue).

If the v3 hotspot calls were made on these cohorts, the O7 evidence for v3 residues would come
almost entirely from GENIE — the same database O4 is counted in — and excluding MSK would
remove even less of the overlap than it does for v2. That can't be confirmed from the file;
confirming it needs the method description for the v3 release.

## Files

Tables 01–32 are the same tables as in the main `results/`, for v2 + v3. Two are new:

| File | Contents |
|---|---|
| `results/33_v3_cohort_prevalence_by_residue.tsv` | Each v3 residue: mutations and prevalence in MSK, non-MSK GENIE and TCGA, the GENIE share, prevalence tests, and the per-change counts behind the call |
| `results/34_v3_cohort_prevalence_summary.tsv` | Totals and shares for the table above |

**Kept local, as in the main analysis:** everything derived from GENIE or COSMIC per variant —
tables 11–14, 18–32 and 35 (`35_genie_v3_residue_o4.tsv`: each v3 residue with the GENIE O4
of its changes), the GENIE and COSMIC workbooks, and the GENIE and COSMIC counts in
`interim/`. They are gitignored. For Kevin, `results/genie_points_per_change.xlsx` has the
per-change table with O7, O4 and the handling under both criteria, now including the v3
changes and a `hotspot_version` column.

Figures: `figures/fig1` (residues by number of changes) and `fig3`–`fig4` (GENIE leave-MSK-out
and lineage), drawn for v2 + v3. There is no figure 2 here: it shows the MSK split, which only
v2 residues have, so it would repeat the v2 figure.
