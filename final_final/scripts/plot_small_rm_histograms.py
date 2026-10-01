#!/usr/bin/env python3
"""Histograms of class and cargo composition for small R-M pIEs.

Set: n_rm + n_other ≤ 3; integrase and Rep do not count.
JAWRJC dropped (HsdS is the 'integrase' call, not RM cargo).
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.patches import Patch

ART = FF / "artifacts/unknown_le3_rm"
FIG = FF / "figures/unknown_le3_rm"
HITS = ART / "rm_allow_rep.tsv"
ATLAS = FF / "artifacts/attL_complement_blast_window/cohort_6960_atlas_join.tsv"
TAX = REPO / "results_organized/misc/gtdb_taxonomy.tsv"

BLUE = "#0072B2"
SKY = "#56B4E9"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERM = "#D55E00"
PINK = "#CC79A7"
GREY = "#999999"
INK = "#222222"
PURPLE = "#AA4499"

RM_ORDER = ["Type IIG", "Type II", "Type IV", "GmrSD", "MTase only"]
RM_COLOR = {
    "Type IIG": GREEN,
    "Type II": VERM,
    "Type IV": ORANGE,
    "GmrSD": PINK,
    "MTase only": SKY,
}
EXTRA_ORDER = ["XRE", "orphan", "TM", "2nd integrase", "PDC-S15", "HTH", "coil", "disordered"]
EXTRA_COLOR = {
    "XRE": "#44AA99",
    "orphan": "#B8B8B8",
    "TM": "#666666",
    "2nd integrase": BLUE,
    "PDC-S15": ORANGE,
    "HTH": "#44AA99",
    "coil": GREY,
    "disordered": "#DDDDDD",
}
# Refined extras after PADLOC + InterPro (not Bakta product names).
IE_EXTRAS = {
    "CAAGAD010000005.1_18": [],
    "JACKEI010000001.1_834": ["orphan"],
    "JAUHVR010000008.1_157": ["TM"],
    "JAUYQH010000169.1_8": ["XRE"],
    "JBCCXD010000004.1_43": ["XRE"],
    "SVFK01000019.1_20": ["coil"],
    "QKJQ01000022.1_9": ["2nd integrase"],
    "JAPXKQ010000125.1_1": ["2nd integrase"],
    "JADJZA010000006.1_80": ["disordered"],
    "JAEZUR010000054.1_8": ["orphan"],
    "CAMCFJ010000009.1_21": ["orphan"],
    "JBFJJN010000005.1_96": ["PDC-S15", "HTH"],
    "WBSG01000027.1_364": ["orphan"],
}

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def letter(ax, s: str) -> None:
    ax.text(-0.12, 1.08, s, transform=ax.transAxes, fontsize=14, fontweight="bold",
            va="bottom", color=INK)


def save(fig, stem: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)
    print("→", FIG / f"{stem}.png")


def rm_class(r) -> str:
    pad = str(r.padloc or "")
    blob = f"{r.rm_products} {r.other_products}".lower()
    if "RM_type_IIG" in pad:
        return "Type IIG"
    if "RM_type_II" in pad:
        return "Type II"
    if "RM_type_IV" in pad or "mrr" in blob or "type iv" in blob:
        return "Type IV"
    if "gmrsd" in blob:
        return "GmrSD"
    if "restriction endonuclease" in blob and (
        "methyl" in blob or "s-adenosylmethionine" in blob
    ):
        return "Type II"
    if "methyl" in blob and "restrict" not in blob:
        # long fused MmeI-like MTase without PADLOC → IIG
        aa = []
        for part in str(r.rm_products).split(";"):
            if "aa," in part:
                try:
                    aa.append(int(part.split("(")[1].split("aa")[0]))
                except (IndexError, ValueError):
                    pass
        if aa and max(aa) >= 800:
            return "Type IIG"
        return "MTase only"
    if "restrict" in blob:
        return "Type IV"
    return "MTase only"


def recipe(r) -> str:
    cls = r.rm_class
    extra = []
    if int(r.n_rep) > 0:
        extra.append("Rep")
    extras = list(IE_EXTRAS.get(r.integrase_id, []))
    extra.extend(extras)
    return cls if not extra else f"{cls} + {' + '.join(extra)}"


def apply_refined_counts(hits: pd.DataFrame) -> pd.DataFrame:
    """JAPXKQ hyp is a second integrase; WBSG SAM is the Type II MTase."""
    hits = hits.copy()
    m = hits.integrase_id == "JAPXKQ010000125.1_1"
    hits.loc[m, ["n_int", "n_other", "n_budget"]] = [2, 0, 2]
    m = hits.integrase_id == "WBSG01000027.1_364"
    hits.loc[m, ["n_rm", "n_other", "n_budget"]] = [2, 1, 3]
    m = hits.integrase_id == "QKJQ01000022.1_9"
    # HEPN is defense cargo; second CDS is integrase (already n_int=2)
    return hits


def annotate_bars(ax, bars, ys, dy=0.15) -> None:
    ymax = max(ys) if len(ys) else 1
    for b, v in zip(bars, ys):
        if v <= 0:
            continue
        ax.text(b.get_x() + b.get_width() / 2, v + dy, str(int(v)),
                ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_ylim(0, ymax + 1.4)


def classify() -> pd.DataFrame:
    hits = pd.read_csv(HITS, sep="\t")
    hits = hits[
        (hits.integrase_id != "JAWRJC010000145.1_11")
        & (hits.func_class == "unknown")
    ].copy()
    atlas = pd.read_csv(ATLAS, sep="\t", usecols=["integrase_id", "sample_id", "primary_label"])
    tax = pd.read_csv(TAX, sep="\t")
    hits = hits.merge(atlas, on="integrase_id", how="left")
    hits = hits.merge(tax, on="sample_id", how="left")
    hits = apply_refined_counts(hits)
    hits["rm_class"] = hits.apply(rm_class, axis=1)
    hits["recipe"] = hits.apply(recipe, axis=1)
    hits["has_rep"] = hits.n_rep > 0
    hits["phylum"] = hits.phylum.fillna("Unknown")
    dest = ART / "rm_allow_rep_classified.tsv"
    hits.to_csv(dest, sep="\t", index=False)
    print("n=", len(hits), "→", dest)
    return hits


def fig_histograms(df: pd.DataFrame) -> None:
    n = len(df)
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2))

    # A — RM class
    ax = axes[0, 0]
    letter(ax, "A")
    counts = df.rm_class.value_counts().reindex(RM_ORDER, fill_value=0)
    bars = ax.bar(range(len(RM_ORDER)), counts.values,
                  color=[RM_COLOR[k] for k in RM_ORDER], width=0.68)
    annotate_bars(ax, bars, counts.values)
    ax.set_xticks(range(len(RM_ORDER)))
    ax.set_xticklabels(RM_ORDER, rotation=20, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title(f"R-M class  (n = {n})")

    # B — refined extras (PADLOC / InterPro, not Bakta)
    ax = axes[0, 1]
    letter(ax, "B")
    counts = [sum(k in IE_EXTRAS.get(iid, []) for iid in df.integrase_id) for k in EXTRA_ORDER]
    bars = ax.bar(range(len(EXTRA_ORDER)), counts,
                  color=[EXTRA_COLOR[k] for k in EXTRA_ORDER], width=0.68)
    annotate_bars(ax, bars, counts)
    ax.set_xticks(range(len(EXTRA_ORDER)))
    ax.set_xticklabels(EXTRA_ORDER, rotation=25, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title("Bakta-unknown extras, re-annotated")

    # C — extra-gene budget (int and Rep free)
    ax = axes[1, 0]
    letter(ax, "C")
    budgets = [1, 2, 3]
    no = np.array([int(((df.n_budget == b) & ~df.has_rep).sum()) for b in budgets])
    yes = np.array([int(((df.n_budget == b) & df.has_rep).sum()) for b in budgets])
    x = np.arange(len(budgets))
    ax.bar(x, no, color=GREEN, width=0.62, label="no Rep")
    ax.bar(x, yes, bottom=no, color=PURPLE, width=0.62, label="+ Rep")
    for i, (a, b) in enumerate(zip(no, yes)):
        tot = a + b
        if tot:
            ax.text(i, tot + 0.12, str(int(tot)), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(["1 gene", "2 genes", "3 genes"])
    ax.set_ylabel("n IE")
    ax.set_ylim(0, max(no + yes) + 1.5)
    ax.set_title("Extra genes besides int (Rep free)")
    ax.legend(frameon=False, fontsize=8)

    # D — host phylum
    ax = axes[1, 1]
    letter(ax, "D")
    phy = df.phylum.value_counts()
    colors_p = {
        "Bacillota": BLUE,
        "Pseudomonadota": VERM,
        "Myxococcota": GREEN,
        "Actinomycetota": ORANGE,
        "Unknown": GREY,
    }
    bars = ax.barh(range(len(phy))[::-1], phy.values[::-1],
                   color=[colors_p.get(k, GREY) for k in phy.index[::-1]], height=0.62)
    for y, v in zip(range(len(phy))[::-1], phy.values[::-1]):
        ax.text(v + 0.12, y, str(int(v)), va="center", fontsize=9)
    ax.set_yticks(range(len(phy))[::-1])
    ax.set_yticklabels(phy.index[::-1])
    ax.set_xlabel("n IE")
    ax.set_xlim(0, phy.values.max() + 1.8)
    ax.set_title("Host phylum")

    fig.tight_layout(w_pad=2.2, h_pad=2.0)
    save(fig, "Figure_small_rm_histograms")


def fig_budget(df: pd.DataFrame) -> None:
    """Extra gene budget (RM + other; int/Rep free) and length."""
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.7))

    ax = axes[0]
    letter(ax, "A")
    budgets = [1, 2, 3]
    no = [(df.n_budget == b) & ~df.has_rep for b in budgets]
    yes = [(df.n_budget == b) & df.has_rep for b in budgets]
    no = np.array([int(s.sum()) for s in no])
    yes = np.array([int(s.sum()) for s in yes])
    x = np.arange(len(budgets))
    ax.bar(x, no, color=GREEN, width=0.62, label="no Rep")
    ax.bar(x, yes, bottom=no, color=PURPLE, width=0.62, label="+ Rep")
    for i, (a, b) in enumerate(zip(no, yes)):
        tot = a + b
        if tot:
            ax.text(i, tot + 0.12, str(int(tot)), ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(["1 gene", "2 genes", "3 genes"])
    ax.set_ylabel("n IE")
    ax.set_ylim(0, max(no + yes) + 1.5)
    ax.set_title("Extra genes besides int (Rep free)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    letter(ax, "B")
    kb = df.ie_len_nt / 1000
    bins = np.arange(2, 11, 1)
    ax.hist(kb[~df.has_rep], bins=bins, color=GREEN, alpha=0.9, label="no Rep")
    ax.hist(kb[df.has_rep], bins=bins, color=PURPLE, alpha=0.85, label="+ Rep")
    ax.set_xlabel("pIE length (kb)")
    ax.set_ylabel("n IE")
    ax.set_title(f"Length  (median {kb.median():.1f} kb)")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout(w_pad=2.4)
    save(fig, "Figure_small_rm_budget_length")


def main() -> None:
    df = classify()
    fig_histograms(df)
    fig_budget(df)


if __name__ == "__main__":
    main()
