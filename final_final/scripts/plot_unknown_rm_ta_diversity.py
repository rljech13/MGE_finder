#!/usr/bin/env python3
"""Diversity histograms: R-M and T/AT on the full unknown cohort (n=3,738)."""
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

ART_RM = FF / "artifacts/unknown_rm"
ART_TA = FF / "artifacts/unknown_ta"
FIG_RM = FF / "figures/unknown_rm"
FIG_TA = FF / "figures/unknown_ta"
N_UNKNOWN = 3738

BLUE = "#0072B2"
SKY = "#56B4E9"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERM = "#D55E00"
PINK = "#CC79A7"
GREY = "#999999"
INK = "#222222"
TEAL = "#44AA99"
PURPLE = "#AA4499"
WINE = "#882255"

RM_ORDER = ["Type I", "Type II", "Type IIG", "Type III", "Type IV", "GmrSD", "MTase only", "Mixed"]
RM_COLOR = {
    "Type I": BLUE,
    "Type II": VERM,
    "Type IIG": GREEN,
    "Type III": TEAL,
    "Type IV": ORANGE,
    "GmrSD": PINK,
    "MTase only": SKY,
    "Mixed": GREY,
}
COMP_ORDER = ["M+R pair", "MTase only", "REase only", "mixed classes"]
COMP_COLOR = {
    "M+R pair": GREEN,
    "MTase only": SKY,
    "REase only": VERM,
    "mixed classes": GREY,
}

TA_ORDER = [
    "RelBE", "HicAB", "AbiE", "VapBC", "MazEF", "HipAB",
    "Phd-Doc", "HigAB", "MqsRA", "YefM-YoeB", "BrnTA", "CcdAB", "unnamed",
]
TA_COLOR = {
    "RelBE": VERM,
    "HicAB": GREEN,
    "AbiE": PINK,
    "VapBC": SKY,
    "MazEF": ORANGE,
    "HipAB": BLUE,
    "Phd-Doc": TEAL,
    "HigAB": PURPLE,
    "MqsRA": WINE,
    "YefM-YoeB": "#117733",
    "BrnTA": "#88CCEE",
    "CcdAB": "#332288",
    "unnamed": GREY,
}
PAIR_ORDER = ["named pair", "likely pair", "toxin only", "antitoxin only"]
PAIR_COLOR = {
    "named pair": GREEN,
    "likely pair": ORANGE,
    "toxin only": VERM,
    "antitoxin only": SKY,
}

PHY_COLOR = {
    "Bacillota": BLUE,
    "Pseudomonadota": VERM,
    "Acidobacteriota": GREEN,
    "Actinomycetota": ORANGE,
    "Desulfobacterota": TEAL,
    "Cyanobacteriota": SKY,
    "Nitrospirota": PINK,
    "Campylobacterota": PURPLE,
    "Chloroflexota": WINE,
    "Bdellovibrionota": "#44AA99",
    "Planctomycetota": "#117733",
    "Unknown": GREY,
    "Other": GREY,
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


def save(fig, folder: Path, stem: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(folder / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)
    print("→", folder / f"{stem}.png")


def annotate_bars(ax, bars, ys, dy=None) -> None:
    ymax = max(ys) if len(ys) else 1
    if dy is None:
        dy = max(0.4, 0.03 * ymax)
    for b, v in zip(bars, ys):
        if v <= 0:
            continue
        ax.text(b.get_x() + b.get_width() / 2, v + dy, str(int(v)),
                ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_ylim(0, ymax + max(4, 0.18 * ymax))


def phylum_top(df: pd.DataFrame, k: int = 8) -> pd.Series:
    phy = df.phylum.fillna("Unknown").replace("", "Unknown")
    vc = phy.value_counts()
    top = vc.head(k)
    rest = int(vc.iloc[k:].sum()) if len(vc) > k else 0
    if rest:
        top = pd.concat([top, pd.Series({"Other": rest})])
    return top


def barh_phylum(ax, df: pd.DataFrame) -> None:
    phy = phylum_top(df, 8)
    ys = np.arange(len(phy))[::-1]
    ax.barh(ys, phy.values, color=[PHY_COLOR.get(k, GREY) for k in phy.index], height=0.62)
    for y, v in zip(ys, phy.values):
        ax.text(v + max(0.3, 0.02 * phy.values.max()), y, str(int(v)), va="center", fontsize=8)
    ax.set_yticks(ys)
    ax.set_yticklabels(phy.index)
    ax.set_xlabel("n IE")
    ax.set_xlim(0, phy.values.max() * 1.18)
    ax.set_title("Host phylum")


def length_hist(ax, kb: pd.Series, color: str) -> None:
    kb = kb[np.isfinite(kb)]
    hi = min(40.0, max(12.0, float(np.quantile(kb, 0.95)) + 2))
    bins = np.arange(0, hi + 2, 2)
    ax.hist(kb[kb < hi], bins=bins, color=color, edgecolor="white", linewidth=0.4)
    med = float(kb.median())
    ax.axvline(med, color=INK, ls="--", lw=1.0, alpha=0.75)
    n_over = int((kb > hi).sum())
    ax.set_xlabel("pIE length (kb)")
    ax.set_ylabel("n IE")
    tail = f"; {n_over} > {hi:.0f} kb" if n_over else ""
    ax.set_title(f"Length  (median {med:.1f} kb{tail})")


def fig_rm(df: pd.DataFrame) -> None:
    n = len(df)
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.3))

    ax = axes[0, 0]
    letter(ax, "A")
    counts = df.rm_class.value_counts().reindex(RM_ORDER, fill_value=0)
    bars = ax.bar(range(len(RM_ORDER)), counts.values,
                  color=[RM_COLOR[k] for k in RM_ORDER], width=0.68)
    annotate_bars(ax, bars, counts.values)
    ax.set_xticks(range(len(RM_ORDER)))
    ax.set_xticklabels(RM_ORDER, rotation=25, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title(f"R-M class  (n = {n} / {N_UNKNOWN:,})")

    ax = axes[0, 1]
    letter(ax, "B")
    cc = df.rm_composition.value_counts().reindex(COMP_ORDER, fill_value=0)
    bars = ax.bar(range(len(COMP_ORDER)), cc.values,
                  color=[COMP_COLOR[k] for k in COMP_ORDER], width=0.68)
    annotate_bars(ax, bars, cc.values)
    ax.set_xticks(range(len(COMP_ORDER)))
    ax.set_xticklabels(["M+R pair", "MTase only", "REase only", "mixed classes"],
                       rotation=20, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title("Composition")

    ax = axes[1, 0]
    letter(ax, "C")
    length_hist(ax, df.ie_len_nt / 1000, GREEN)

    ax = axes[1, 1]
    letter(ax, "D")
    barh_phylum(ax, df)

    fig.tight_layout(w_pad=2.2, h_pad=2.0)
    save(fig, FIG_RM, "Figure_unknown_rm_diversity")


def fig_ta(df: pd.DataFrame) -> None:
    n = len(df)
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.3))

    ax = axes[0, 0]
    letter(ax, "A")
    present = [k for k in TA_ORDER if int((df.ta_class == k).sum()) > 0]
    counts = df.ta_class.value_counts().reindex(present, fill_value=0)
    bars = ax.bar(range(len(present)), counts.values,
                  color=[TA_COLOR[k] for k in present], width=0.68)
    annotate_bars(ax, bars, counts.values)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=30, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title(f"T/AT family  (n = {n} / {N_UNKNOWN:,})")

    ax = axes[0, 1]
    letter(ax, "B")
    pc = df.ta_pair.value_counts().reindex(PAIR_ORDER, fill_value=0)
    bars = ax.bar(range(len(PAIR_ORDER)), pc.values,
                  color=[PAIR_COLOR[k] for k in PAIR_ORDER], width=0.68)
    annotate_bars(ax, bars, pc.values)
    ax.set_xticks(range(len(PAIR_ORDER)))
    ax.set_xticklabels(PAIR_ORDER, rotation=20, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title("Pair completeness")

    ax = axes[1, 0]
    letter(ax, "C")
    length_hist(ax, df.ie_len_nt / 1000, VERM)

    ax = axes[1, 1]
    letter(ax, "D")
    barh_phylum(ax, df)

    fig.tight_layout(w_pad=2.2, h_pad=2.0)
    save(fig, FIG_TA, "Figure_unknown_ta_diversity")


def main() -> None:
    rm = pd.read_csv(ART_RM / "unknown_rm_hits.tsv", sep="\t")
    ta = pd.read_csv(ART_TA / "unknown_ta_hits.tsv", sep="\t")
    rm = rm[rm.rm_class.notna() & (rm.rm_class != "")].copy()
    ta = ta[ta.ta_class.notna() & (ta.ta_class != "")].copy()
    print("RM", len(rm), "TA", len(ta))
    fig_rm(rm)
    fig_ta(ta)


if __name__ == "__main__":
    main()
