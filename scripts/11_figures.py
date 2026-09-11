#!/usr/bin/env python3
"""
11_figures.py
-------------
Static figures for the report and for deliveries. Every figure is written as PNG
(for email and slides) and SVG, alongside a TSV of exactly the values plotted, so
each figure has a table twin.

CancerHotspots figures are built from tracked results and regenerate anywhere.
The GENIE figures plot summary figures only -- the same ones the report quotes --
and are built only when the gitignored GENIE tables exist locally (scripts 09-10).

Design: the dataviz reference palette (light), used unchanged -- categorical slots 1
and 2 plus neutral ink. One series is one colour; emphasis picks out the rows the
figure is about and greys the rest. Bars are 18 CSS px thick with a 4 px rounded
data end, square at the baseline. Every bar is labelled, so there is no value axis.

Outputs (figures/):
  fig1_changes_per_residue.{png,svg,tsv}
  fig2_msk_share_per_residue.{png,svg,tsv}
  fig3_genie_leave_msk_out.{png,svg,tsv}     needs GENIE
  fig4_genie_lineage.{png,svg,tsv}           needs GENIE
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch, PathPatch  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

# Reference palette, light mode.
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
BASELINE = "#c3c2b7"
SLOT1, SLOT2 = "#2a78d6", "#eb6834"
DEEMPH = "#c3c2b7"

CSS_PX = 1 / 96          # inches
BAR_PX, RADIUS_PX, GAP_PX = 18, 4, 2
ROW_IN = 0.42            # height of one category row

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "svg.fonttype": "none",
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

CH_SOURCE = "Source: Chang et al., Cancer Discov 2018 (cancerhotspots.org v2)."
GENIE_SOURCE = "Source: AACR Project GENIE v20.0-public (summary figures only); CancerHotspots v2."


# ---------------------------------------------------------------- drawing kit
class Canvas:
    """A figure laid out in inches, so marks can be specified in CSS pixels."""

    def __init__(self, width, height):
        self.w, self.h = width, height
        self.fig = plt.figure(figsize=(width, height))

    def axes(self, left, bottom, width, height):
        return self.fig.add_axes([left / self.w, bottom / self.h,
                                  width / self.w, height / self.h])

    def header(self, title, subtitle, left=0.3):
        self.fig.text(left / self.w, 1 - 0.2 / self.h, title, va="top", ha="left",
                      fontsize=12.5, fontweight="semibold", color=INK)
        self.fig.text(left / self.w, 1 - 0.5 / self.h, subtitle, va="top", ha="left",
                      fontsize=9.5, color=INK2, linespacing=1.45)

    def source(self, text, left=0.3):
        self.fig.text(left / self.w, 0.14 / self.h, text, va="bottom", ha="left",
                      fontsize=8, color=MUTED)

    def save(self, name):
        FIGURES.mkdir(exist_ok=True)
        for ext in ("png", "svg"):
            self.fig.savefig(FIGURES / f"{name}.{ext}", dpi=200)
        plt.close(self.fig)


def px(ax, dx=0.0, dy=0.0):
    """CSS pixels -> data units on ax. Limits must already be set."""
    scale = ax.figure.dpi * CSS_PX
    origin = ax.transData.transform((0, 0))
    moved = ax.transData.inverted().transform(origin + [dx * scale, dy * scale])
    return abs(moved[0]), abs(moved[1])


def hbar(ax, y, value, color, thickness=BAR_PX):
    """Horizontal bar from x=0, 4 px rounded at the data end, square at the baseline."""
    if value <= 0:
        return
    rx, _ = px(ax, dx=RADIUS_PX)
    _, ry = px(ax, dy=RADIUS_PX)
    _, h = px(ax, dy=thickness)
    rx = min(rx, value / 2)
    y0, y1 = y - h / 2, y + h / 2
    verts = [(0, y0), (value - rx, y0), (value, y0), (value, y0 + ry),
             (value, y1 - ry), (value, y1), (value - rx, y1), (0, y1), (0, y0)]
    codes = [MPath.MOVETO, MPath.LINETO, MPath.CURVE3, MPath.CURVE3, MPath.LINETO,
             MPath.CURVE3, MPath.CURVE3, MPath.LINETO, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor="none", zorder=3))


def value_label(ax, y, value, text):
    gap, _ = px(ax, dx=6)
    ax.text(value + gap, y, text, va="center", ha="left", fontsize=9, color=INK2, zorder=5)


def category_axis(ax, n, xmax, labels=None):
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.5, n - 0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.tick_params(axis="y", length=0, labelsize=9.5, colors=INK2, pad=8)
    ys = [n - 1 - i for i in range(n)]
    ax.set_yticks(ys, labels if labels is not None else [""] * n)
    ax.axvline(0, color=BASELINE, linewidth=0.75, zorder=4)   # 1 CSS px
    return ys


def write_twin(name, df):
    FIGURES.mkdir(exist_ok=True)
    df.to_csv(FIGURES / f"{name}.tsv", sep="\t", index=False)


def emphasis_chart(name, cats, values, labels, highlight, title, subtitle, source, xmax):
    """One series, with the rows the figure is about in slot 1 and the rest grey."""
    n = len(cats)
    header_in, bottom_in = 1.0, 0.5
    c = Canvas(7.5, header_in + ROW_IN * n + bottom_in)
    ax = c.axes(1.35, bottom_in, 7.5 - 1.35 - 1.4, ROW_IN * n)
    ys = category_axis(ax, n, xmax, cats)
    for y, v, lab, cat in zip(ys, values, labels, cats):
        hbar(ax, y, v, SLOT1 if cat in highlight else DEEMPH)
        value_label(ax, y, v, lab)
    c.header(title, subtitle)
    c.source(source)
    c.save(name)
    write_twin(name, pd.DataFrame({"category": cats, "value": values, "label": labels,
                                   "highlighted": [cat in highlight for cat in cats]}))


# ---------------------------------------------------------------- figures
def fig1_changes_per_residue():
    pos = pd.read_csv(RESULTS / "03_gene_position_table.tsv", sep="\t")
    cats = ["1 change", "2 changes", "3 changes", "4 or more"]
    counts = (pd.cut(pos["n_unique_changes"], [0, 1, 2, 3, 10 ** 6], labels=cats)
              .value_counts().reindex(cats))
    assert counts.tolist() == [213, 293, 231, 287], counts.tolist()
    total = counts.sum()
    emphasis_chart(
        "fig1_changes_per_residue", cats, counts.tolist(),
        [f"{v:,} ({100 * v / total:.1f}%)" for v in counts],
        {"1 change"},
        "One hotspot residue in five carries a single amino-acid change",
        f"CancerHotspots v2 protein-coding hotspot residues (n = {total:,}), by the number of distinct\n"
        "amino-acid changes seen there. Where there is only one, O7 and O4 measure the same thing.",
        CH_SOURCE + " Table: results/03_gene_position_table.tsv.",
        xmax=330,
    )


def fig2_msk_share_per_residue():
    pos = pd.read_csv(RESULTS / "03_gene_position_table.tsv", sep="\t")
    cats = ["0–25%", "25–50%", "50–75%", "75–100%"]
    counts = (pd.cut(pos["msk_fraction"], [0, 0.25, 0.5, 0.75, 1.000001], right=False,
                     labels=cats).value_counts().reindex(cats))
    assert counts.tolist() == [180, 412, 344, 88], counts.tolist()
    total = counts.sum()
    emphasis_chart(
        "fig2_msk_share_per_residue", cats, counts.tolist(),
        [f"{v:,} ({100 * v / total:.1f}%)" for v in counts],
        {"50–75%", "75–100%"},
        "MSK-IMPACT supplies most of the evidence at 42% of hotspot residues",
        f"Hotspot residues (n = {total:,}) by the share of their mutations contributed by MSK-IMPACT,\n"
        "a GENIE contributing centre. Bands are left-closed; highlighted bands are MSK-majority.",
        CH_SOURCE + " Table: results/03_gene_position_table.tsv.",
        xmax=470,
    )


def fig3_genie_leave_msk_out():
    lmo = pd.read_csv(RESULTS / "23_genie_leave_msk_out_transitions.tsv", sep="\t")
    excl_cols = [c for c in lmo.columns if c.startswith("o4_excl_msk_")]
    tiers = ["Strong", "Moderate", "Supporting"]
    rows = []
    for tier in tiers:
        r = lmo[(lmo["o7_strength"] == tier) & (lmo["o4_all_patients"] == "Strong")].iloc[0]
        rows.append((f"O7 {tier.lower()}", int(r[excl_cols].sum()), int(r["o4_excl_msk_Strong"])))
    assert rows[0][1:] == (197, 197), rows[0]
    cats = [r[0] for r in rows]
    n = len(cats)
    header_in, bottom_in = 1.3, 0.5
    c = Canvas(7.5, header_in + 0.62 * n + bottom_in)
    ax = c.axes(1.35, bottom_in, 7.5 - 1.35 - 1.2, 0.62 * n)
    ys = category_axis(ax, n, 1180, cats)
    bar = 16
    _, h = px(ax, dy=bar)
    _, g = px(ax, dy=GAP_PX)
    off = (h + g) / 2
    for y, (_, all_p, excl) in zip(ys, rows):
        hbar(ax, y + off, all_p, SLOT1, bar)
        value_label(ax, y + off, all_p, f"{all_p:,}")
        hbar(ax, y - off, excl, SLOT2, bar)
        value_label(ax, y - off, excl, f"{excl:,}")
    ax.legend([Patch(facecolor=SLOT1), Patch(facecolor=SLOT2)],
              ["All GENIE patients", "MSK-IMPACT patients excluded"],
              loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, frameon=False,
              fontsize=9, labelcolor=INK2, handlelength=0.9, handleheight=0.9,
              borderaxespad=0, columnspacing=1.6)
    c.header("Excluding MSK-IMPACT does not remove a single O4 + O7 double count",
             "Changes reaching O4_strong in GENIE (unique patients), by their O7 tier, with and without\n"
             "MSK-IMPACT. Changes scoring the full +8: 197 with MSK, 197 without.")
    c.source(GENIE_SOURCE)
    c.save("fig3_genie_leave_msk_out")
    write_twin("fig3_genie_leave_msk_out", pd.DataFrame(
        rows, columns=["o7_tier", "o4_strong_all_patients", "o4_strong_msk_excluded"]))


def fig4_genie_lineage():
    comp = pd.read_csv(RESULTS / "18_genie_cohort_composition.tsv", sep="\t")
    allc = comp[comp["cohort"] == "All centres"].set_index("lineage")
    msk = comp[comp["cohort"] == "MSK-IMPACT"].set_index("lineage")
    total_patients = int(allc.loc["Total", "patients"])
    order = [("Solid", "Solid tumour"), ("Lymphoid", "Lymphoid"),
             ("Myeloid", "Myeloid"), ("Unknown", "Unknown")]
    data = pd.DataFrame([{
        "lineage": label,
        "patients": int(allc.loc[k, "patients"]),
        "pct_patients": 100 * allc.loc[k, "patients"] / total_patients,
        "samples_per_patient": allc.loc[k, "samples"] / allc.loc[k, "patients"],
        "msk_share": msk.loc[k, "patients"] / allc.loc[k, "patients"],
    } for k, label in order])
    assert data.loc[1, "patients"] == 13131 and data.loc[2, "patients"] == 12018
    assert f"{data.loc[2, 'pct_patients']:.1f}" == "4.9", data.loc[2, "pct_patients"]
    highlight = {"Lymphoid", "Myeloid"}
    n = len(data)

    header_in, bottom_in, top_pad = 1.05, 0.5, 0.35
    c = Canvas(10.0, header_in + top_pad + ROW_IN * n + bottom_in)
    panels = [
        ("Patients", "patients", 215540 / 0.52,
         lambda r: f"{r.patients:,} ({r.pct_patients:.1f}%)"),
        ("Samples per patient", "samples_per_patient", 2.75,
         lambda r: f"{r.samples_per_patient:.2f}"),
        ("Share of patients from MSK-IMPACT", "msk_share", 1.3,
         lambda r: f"{100 * r.msk_share:.1f}%"),
    ]
    left, width, gap = 1.3, 2.35, 0.45
    for i, (ptitle, col, xmax, fmt) in enumerate(panels):
        ax = c.axes(left + i * (width + gap), bottom_in, width, ROW_IN * n)
        ys = category_axis(ax, n, xmax, data["lineage"].tolist() if i == 0 else None)
        for y, r in zip(ys, data.itertuples()):
            v = getattr(r, col)
            hbar(ax, y, v, SLOT1 if r.lineage in highlight else DEEMPH)
            value_label(ax, y, v, fmt(r))
        ax.set_title(ptitle, loc="left", fontsize=9.5, fontweight="semibold",
                     color=INK2, pad=10)
    c.header("Haematological neoplasms are a tenth of GENIE, sampled more often, and disproportionately from MSK",
             "AACR Project GENIE v20 by lineage. Samples per patient is why de-duplication matters most in\n"
             "haematology. 3,146 patients have samples in more than one lineage, so rows sum to over 100%.")
    c.source(GENIE_SOURCE + " Lineage from OncoTree tissue.")
    c.save("fig4_genie_lineage")
    out = data.copy()
    out["samples_per_patient"] = out["samples_per_patient"].round(2)
    out["msk_share"] = out["msk_share"].round(3)
    out["pct_patients"] = out["pct_patients"].round(2)
    write_twin("fig4_genie_lineage", out)


def main():
    fig1_changes_per_residue()
    fig2_msk_share_per_residue()
    print("wrote fig1, fig2")
    if (RESULTS / "23_genie_leave_msk_out_transitions.tsv").exists():
        fig3_genie_leave_msk_out()
        fig4_genie_lineage()
        print("wrote fig3, fig4")
    else:
        print("GENIE tables not present (they are gitignored); fig3 and fig4 skipped")


if __name__ == "__main__":
    main()
