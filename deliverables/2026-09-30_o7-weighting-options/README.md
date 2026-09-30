# O4 / O7 per change, and how to weight O7 when the cap does not apply — 30 September 2026

For Kevin Baker, answering the two requests of 30 September, in order. Both use the
**CancerHotspots v2 + v3** analysis set with **GENIE v20** for O4
([`analyses/v3/`](../../analyses/v3/)). The tables sent before are unchanged.

## Request 1 — one expanded per-change table

**`genie_o4_o7_expanded.xlsx`, sheet "Points per change (expanded)"** (sent by email):
**3,446 changes × 60 columns**.

- **Nothing was dropped.** Every column of the two tables sent before is present, and for the
  2,918 v2 changes every value is identical to the file you attached. I checked column by
  column against both.
- **Added: the v2 / v3 comparison** — `hotspot_version`, `in_cancerhotspots_v2`,
  `in_cancerhotspots_v3`, and `residue_hotspot_versions` (which releases call that residue).
  2,918 rows are v2 and 528 are new in v3.
- **Added: the SVIG-UK canonical overlap** — `transcript`, `HGVSp_Short`,
  `svig_uk_assessment` and `on_svig_uk_canonical_list` (126 changes are on the O1 list).
- **Retained: GENIE v20 O4 in full** — counts (rows, samples, unique patients, MSK, non-MSK,
  myeloid, lymphoid, solid, unknown), the O4 strength and points, and the same again with
  MSK-IMPACT excluded.
- **Retained: both handlings** — `o4_o7_handling__permissive` / `__strict` with their points
  and points lost, alongside the residual columns the test uses.

## Request 2 — weighting O7 on uncapped changes

Where the leave-one-variant-out test finds independent positional evidence, O4 and O7
combine and the question is what O7 is worth. Both options are scored for **every uncapped
change** in sheet **"O7 weighting options"**, each with the rule it passed or failed
(`option1_o7_rule`, `option2_o7_rule`), the reason they differ (`difference_category` for
sorting, `difference_reason` with the counts spelled out), and the resulting totals.

O7 thresholds, as in Supplementary Table 1: same change ≥ 10 **and** residue total ≥ 50 →
Strong (+4); same change ≥ 10 with residue total < 50 → Moderate (+2); same change 2–9 →
Supporting (+1).

- **Option 1 (total-count):** the VUA's own change count against the full residue total.
- **Option 2 (residual-count):** the largest **other** change against the residue total minus
  the VUA.

Your worked example reproduces exactly: p.M167I (total 66; I:55, K:5) gives Option 1 = O7
Strong, +4 + 4 = **+8**; Option 2 = O7 Supporting, +4 + 1 = **+5**.

### The result, on 805 uncapped changes (permissive criterion)

| | Option 1 | Option 2 |
|---|---|---|
| O7 Strong | 189 | 281 |
| O7 Moderate | 88 | 328 |
| O7 Supporting | 528 | 196 |
| Mean O7 points | 1.81 | **2.45** |
| Changes reaching the full +8 | 189 | **278** |

### The finding that matters: Option 2 is not uniformly stricter

**Option 2 scores *higher* than Option 1 for 377 of the 805 changes, and lower for only 103**
(12 lose 3 points, 58 lose 2, 33 lose 1; 325 are unchanged). On average it is 0.64 points
**more** generous, not less.

The reason is the asymmetry of the leave-one-out construction. Option 2 scores the VUA on
somebody else's evidence:

| Variant | Its own count | Largest other change | Option 1 | Option 2 |
|---|---|---|---|---|
| TP53 p.R282G | 9 of 219 | 201 (R282W) | Supporting → **+5** | Strong → **+8** |
| TP53 p.R281H | 9 of 69 | 14 | Supporting → **+5** | Strong → **+8** |
| Your p.M167I | 55 of 66 | 5 | Strong → **+8** | Supporting → **+5** |

So Option 2 is stricter exactly where you intended (a residue dominated by the VUA) but
**more generous for minor changes at a dominated residue** — a variant seen 9 times gets the
+4 earned by a different change seen 201 times. Auditing which threshold drives it: the
same-change threshold (≥ 10) accounts for 398 of the 480 tier differences, the residue-total
threshold (≥ 50) for 58, and both together for 24.

### Two alternatives, also scored in the table

- **Option 3 — the lesser of Options 1 and 2.** Takes the more conservative of the two for
  each change, so it removes the inversion. Mean O7 1.58; 119 changes reach +8.
- **Option 4 — the VUA's own count against the residual total.** Keeps the variant scored on
  its own recurrence (as SVIG-UK does) but requires the residue to reach ≥ 50 **without** it
  for O7_strong. Mean O7 1.64; 119 reach +8. On your example: 55 ≥ 10 but residual 11 < 50 →
  Moderate, **+6**.

**My suggestion.** Option 1 is right for routine use: it is the guideline as written, and the
independence test already did the job of excluding double counting, so weighting does not
need to do it twice. If the group wants the weighting itself to carry positional evidence,
**Option 4 is the coherent version** — it is stricter than Option 1 in exactly the cases that
motivated the proposal, keeps each variant scored on its own recurrence, and cannot give a
rare change the credit earned by a common one. Option 2 as defined should not be adopted
without addressing that inversion; Option 3 fixes it, at the cost of being two rules rather
than one.

A wording point for the group: the main text says "> 50 occurrences at the same position"
while Supplementary Table 1 says "≥ 50 entries". These tables use ≥ 50, the operational
table; at exactly 50 the two readings differ.

## Files

| Where | File | What it is |
|---|---|---|
| This folder | `o7_weighting_options_summary.tsv` | The counts above for both criteria, and what drives each tier difference |
| **Email** | `genie_o4_o7_expanded.xlsx` | Sheet 1 the expanded per-change table (request 1); sheet 2 every uncapped change under Options 1–4 with its rules (request 2); sheet 3 the summary |

The workbook goes by email because it carries AACR Project GENIE data per variant.
