# How the residue-independence test works — 21 September 2026

For Kevin Baker, in answer to his question of 21 September about table 05
(`results/05_o4_o7_double_counting.tsv`): the exact definition of the permissive and strict
criteria, and why TP53 p.M237I is capped under both.

## Short answer

Neither criterion uses whether the commonest change accounts for more than 80% of the
residue. TP53 p.M237I is capped because, once M237I itself is removed, only **9** mutations
remain at the residue — one short of the 10 both criteria need.

## The definition

For the variant under assessment (VUA), at its residue:

- **Residual** (`residual_count`) = mutations at the residue − the VUA's own mutations.
- **Largest other change** (`residual_max_change_count`) = the count of the commonest change
  at that residue other than the VUA.

The residue is **independent** (`independent_positional_evidence`, `_strict`) if:

| Criterion | Condition |
|---|---|
| **Permissive** | Residual **≥ 10** and largest other change **≥ 2** |
| **Strict** | Largest other change **≥ 10** (which means the residual is ≥ 10 too) |

Independent → O4 and O7 may be combined (max +8). Not independent → combined O4 + O7 capped
at +4. The question the test asks is: would the residue still be a hotspot without the
variant being assessed? The strict criterion additionally requires another change that
would reach O7 moderate or strong on its own (≥ 10 mutations of that same change).

The report's wording of the strict criterion ("as above, but the other change must have
≥ 10") is equivalent: if another change has ≥ 10, the residual is automatically ≥ 10, so
the only difference between the two criteria is how large the largest other change must be.

## TP53 p.M237I

Residue M237 has 64 mutations: **I:55 | K:5 | V:4**.

| | Value | Permissive | Strict |
|---|---|---|---|
| Residual (64 − 55) | **9** | needs ≥ 10 → **fails** | — |
| Largest other change (K) | **5** | ≥ 2 → passes | needs ≥ 10 → **fails** |
| Result | | **Capped at +4** | **Capped at +4** |

It is a borderline case: one more mutation at the residue from another change would make it
pass the permissive criterion.

## Why it looks like the 80% rule

M237 is also 85.9% M237I, so table 05 labels it "mutation-specific" (`hotspot_character`).
That label is descriptive only; neither criterion reads it. Here the two happen to coincide.

A counter-example: **BRAF p.V600E** is 92.9% of its residue ("mutation-specific") but
passes both criteria — it leaves a residual of 64, and V600M alone has 29 — so it combines
up to +8.

## An asymmetry to be aware of

At the same residue, **M237K and M237V count as independent**: their residuals are 59 and 60,
because M237I stays in the count. The test is deliberately asymmetric — a minor change at a
residue dominated by another change inherits that change's positional evidence, while the
dominant change does not. The practical effect is limited: M237K and M237V only reach O7
supporting (+1), so their O4 + O7 is at most +5.

## Files

| File | What it is |
|---|---|
| `worked_examples_TP53_M237_BRAF_V600E.tsv` | The table 05 rows for the three changes at TP53 M237 and for BRAF p.V600E |

The same clarification is recorded in [the findings](../../docs/findings.md), section 6.
