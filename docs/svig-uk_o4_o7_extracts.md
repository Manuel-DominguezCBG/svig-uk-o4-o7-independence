# SVIG-UK v1.1 (26/05/2026) — verbatim extracts relevant to O4 / O7

Source documents (in `docs/svig-uk/`):

* `svig-uk_v1.1_main.pdf` — ACGS guidelines for the classification of oncogenicity
  of somatic variants in cancer, recommendations by SVIG-UK, v1.1, 26/05/2026.
* `svig-uk_v1.1_supplementary.pdf` — Supplementary material (Supplementary Tables
  1–5, Supplementary Figures 1–3).

Text was extracted with `pdftotext -layout`; the plain-text renderings are kept
alongside the PDFs. Ligatures (`ti`, `fi`, `ff`) are lost by that extractor, so the
quotations below have been retyped against the PDFs.

---

## O4 — Enriched in a somatic variant database (Supplementary Table 1)

**Missense and splice variants**

| Criterion | Strength |
|---|---|
| >10 'on-target' (and/or 'off-target' with confirmed consistent mechanism of action) entries in an international (e.g. GENIE) or national (multicentre) curated database | Strong [+4] |
| 5–10 'on-target' (and/or 'off-target' …) entries in an international (e.g. GENIE) or national (multicentre) curated database | Moderate [+2] |
| >10 VUS 'on-target' (…) entries in a curated in-house database | Moderate [+2] |
| 5–10 VUS 'on-target' (…) entries in a curated in-house database | Supporting [+1] |

Separate count thresholds are given for frameshift/nonsense variants (>50 / 20–50 /
10–19) and for in-frame deletions/insertions (>50 / 20–50 / 10–19).

Key guidance notes:

> "Application of O3 (absence/rarity in a population database) is a prerequisite for
> applying O4."

> "International, national or local databases may be independently used to demonstrate
> enrichment, but local entries must not be combined with other datasets (i.e.
> GENIE/COSMIC). … Where the variant is present in multiple datasets …, the dataset
> giving the strongest level of evidence should be used. Evidence from different data
> sets should not be combined/stacked."

> "The GENIE database contains multiple samples from a single patient … and therefore
> care should be taken to ensure that multiple samples from the same patient are only
> counted as a single entry."

Main text (Points-Based System for Classification):

> "Each code may only be applied once and complementary evidence within each code
> should not be 'stacked' to enable application of a code at a higher strength. For
> example, entries in different databases cannot be combined to apply O4 at a higher
> strength than permitted by the use of one database alone."

---

## O7 — Located in a mutational hotspot and/or critical and well-established functional domain (Supplementary Table 1)

| Criterion | Strength |
|---|---|
| Variant under assessment (VUA) at a mutational hotspot in cancerhotspots.org with **≥50 entries at the same amino acid position, of which there are at least 10 entries for the same amino acid change** | Strong [+4] |
| VUA at a mutational hotspot in cancerhotspots.org with **<50 entries at the same amino acid position, of which there are at least 10 entries for the same amino acid change** — OR — meets criteria for moderate application according to guidance below | Moderate [+2] |
| **2–9 entries for the same amino acid change as the VUA** in cancerhotspots.org — OR — meets criteria for supporting application according to guidance below | Supporting [+1] |

The non-cancerhotspots routes to O7 (moderate/supporting) rest on absence of local
benign variation, local enrichment of oncogenic variants, well-established functional
domain, protein–protein paralogy, or in silico structural modelling.

> "Since functional studies frequently inform delineation of critical domains/hotspots,
> caution should be applied when using O7 together with O10 to avoid any overlap. …
> it is essential that scientific judgement is applied to determine the strength of
> application and **to prevent double counting evidence across multiple codes**."

Main text (Mutational hotspots and functional domains):

> "The published recommendations from the ClinGen Germline/Somatic Variant
> Subcommittee specified that the PM1 criterion (O7 equivalent code) can be applied to
> somatically detected hotspots with ≥10 occurrences in Cancer Hotspots or downgraded
> to supporting for fewer occurrences. We have extended upon these recommendations to
> permit the use of O7 at up to strong for **>50 occurrences at the same position where
> at least 10 of these are the same amino acid change** (as also proposed in the
> guidelines by Horak et al., 2022)."

> "It should be noted that whilst solid tumours have good coverage within this database
> the coverage is more limited for haematological malignancies, and alternative
> resources may need to be considered when assessing the application of O7."

> "**MSK and the TCGA should be excluded from this count to enable the use of Cancer
> Hotspots (O7) without double counting evidence.**"

---

## Code combination guidance (Supplementary Figure 1)

The O4 x O7 cell of the oncogenicity combination matrix carries **no footnote**, i.e.
the combination is currently *permitted without restriction* — up to +8 points.

By contrast, the matrix does restrict analogous pairs for exactly the reason under
discussion here:

| Note | Pair | Guidance |
|---|---|---|
| 5 | O7 x O5 | "Not permitted: O5 and O7 cannot be applied together; use evidence code which achieves the highest strength." |
| **6** | **O7 x O6** | **"Restriction: Sum of points awarding from combining O6 and O7 is limited to 4 pts (strong) to avoid double counting key attributes of critical domains / residues which are also often measured by in silico tools (for example evolutionary conservation)."** |
| 9 | O8 x O7 | "Not permitted: Avoid double counting evidence for constraint; O8 should not be used in combination with O7." |
| 18 | O10 x O7 | "Restriction: O7 can only be used alongside O10 if the hotspot/critical region can be defined as evidence independent of functional studies." |

Note 6 is a direct precedent for the proposed O4 + O7 cap: the same "limited to 4 pts
to avoid double counting" construction, applied to a different correlated pair.

---

## Two textual points to raise with the group

1. **">50" vs "≥50".** The main text says O7 may reach strong for *">50 occurrences at
   the same position"*; Supplementary Table 1 says *"≥50 entries at the same amino acid
   position"*. These differ at exactly 50. This analysis uses the Supplementary Table
   wording (≥50), which is the operational table.

2. **"MSK and the TCGA should be excluded from this count."** The intent — avoiding
   O4/O7 double counting — is clearly present in the guideline, but the sentence sits
   at the end of a paragraph about in silico protein modelling and DECIPHER plots, so
   *"this count"* has no unambiguous antecedent. Read as an instruction to exclude MSK
   and TCGA cases from the O4 case count, it is also difficult to operationalise: GENIE
   does not contain TCGA, and the GENIE portal does not present a per-centre
   deduplication view aligned to the CancerHotspots cohort. This wording is a candidate
   for revision alongside any cap.
