# Points per change on GENIE, both criteria — 16 September 2026

For Kevin Baker. This answers your request of 16 September for the tables from the COSMIC
round, now rebuilt with O4 taken from GENIE:

1. The complete gene × position table for CancerHotspots, with the actual substitutions,
   the number of unique changes per residue and the mutation count for each change, plus
   the distribution of residues by number of changes (1, 2, 3, 4, …).
2. For every amino-acid change: its O7 application under SVIG-UK v1.1 and its O4 code from
   GENIE (unique patients) — the table previously called "14 Points per change".
3. The "O1 canonical impact" and "materially affected" tables.
4. The recommended O4 / O7 handling outcome for each change under **both** the strict and
   the permissive criterion of the proposed amendment (report Section 9).

## Files

| Where | File | What it is |
|---|---|---|
| This folder | `03_gene_position_table.tsv` | Every hotspot residue (1,024): substitutions with their counts (e.g. `D:757\|V:657\|C:384…`), number of unique changes, commonest change and its share, MSK split |
| This folder | `02_distribution_unique_changes.tsv`, `02b_distribution_unique_changes_exact.tsv` | Residues by number of distinct changes: 1 → 213 (20.8%), 2 → 293, 3 → 231, 4 or more → 287 |
| **Email** | `genie_points_per_change.xlsx` | All of the above, plus the per-change GENIE tables below |

The workbook goes by email because it contains AACR Project GENIE data at the level of
individual variants, which the GENIE terms do not allow us to publish. Its sheets:

| Sheet | Contents |
|---|---|
| About | Definitions: O7 and O4 thresholds, and the two criteria |
| 02 / 02b / 03 | As above |
| **14 Points per change (GENIE)** | All 2,918 changes: substitutions at the residue, change and residue counts, O7 tier and points, GENIE samples and unique patients (also non-MSK, myeloid, lymphoid, solid), O4 tier and points, combined points uncapped, then for **each criterion** the handling outcome and the resulting points |
| Cap impact by criterion | What each criterion costs, overall and on the O1 list |
| O1 canonical impact | The 125 changes on the SVIG-UK Canonical Variants List, with both outcomes |
| Materially affected | Changes dropping from +8 to +4 under either criterion, and whether O1 protects them |

## The two criteria

Leave-one-variant-out test: remove the change under assessment from the residue total. The
residue is independent evidence if the residual is ≥ 10 mutations **and** at least one other
change there is itself recurrent — **permissive:** ≥ 2 mutations; **strict:** ≥ 10.
Independent → O4 and O7 combine in full (max +8). Not independent → combined O4 + O7 capped
at +4.

## What the choice of criterion changes, on GENIE

| | Permissive | Strict |
|---|---|---|
| Changes scoring both O4 and O7 | 1,522 | 1,522 |
| … of which capped | 794 | 954 |
| Changes whose score falls | 671 (23.0%) | 822 (28.2%) |
| **Changes dropping from +8 to +4** | **10** | **22** |
| … of which on the O1 list (cap moot) | 4 | 9 |
| **Materially affected** (+8 → +4, not on O1) | **6** | **13** |

The strict criterion caps residues where the positional signal is real but spread across
changes that are individually small. Its twelve extra +8 → +4 changes are of that kind:
TP53 p.R282W (201 of 219 at the residue; the other 18 spread over five changes, the largest
with 9), PIK3CA p.E542K, TP53 p.Y220C, IDH2 p.R172K. Five of the twelve are on the O1 list.
That is the trade-off behind recommending the permissive criterion in the report, and these
tables let the group see it variant by variant.
