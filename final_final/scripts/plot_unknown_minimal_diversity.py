#!/usr/bin/env python3
"""Diversity histograms of composition-minimal unknown IEs after hypothesis."""
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

ART = FF / "artifacts/unknown_ime_classes"
OUT = ART / "after_hypothesis"
FIG = FF / "figures/unknown_ime_classes"
IE = OUT / "unknown_minimal_ie.tsv"
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

CLASS_ORDER = [
    "replication", "TA", "RM", "PDC", "RT-defense", "Par", "SoFic",
    "CRISPR", "AMR", "Abi", "metals", "CBASS", "Wadjet", "Lamassu",
]
CLASS_COLOR = {
    "replication": PURPLE,
    "TA": VERM,
    "RM": GREEN,
    "PDC": ORANGE,
    "RT-defense": PINK,
    "Par": BLUE,
    "SoFic": "#44AA99",
    "CRISPR": SKY,
    "Abi": "#882255",
    "AMR": "#8c510a",
    "CBASS": "#1a9850",
    "Wadjet": "#4575b4",
    "Lamassu": "#fc8d59",
    "metals": "#666666",
}
PHYLUM_ORDER = [
    "Bacillota", "Pseudomonadota", "Actinomycetota", "Acidobacteriota",
    "Cyanobacteriota", "Desulfobacterota", "Nitrospirota", "Chloroflexota",
]
PHYLUM_COLOR = {
    "Bacillota": "#4575b4",
    "Pseudomonadota": "#fc8d59",
    "Actinomycetota": "#8c510a",
    "Acidobacteriota": "#1a9850",
    "Cyanobacteriota": "#66c2a5",
    "Desulfobacterota": "#e78ac3",
    "Nitrospirota": "#a6d854",
    "Chloroflexota": "#ffd92f",
    "Others": "#80b1d3",
    "Unknown": "#bdbdbd",
}

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "pdf.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load() -> pd.DataFrame:
    ie = pd.read_csv(IE, sep="\t")
    cm = ie[ie.composition_minimal == 1].copy()
    atlas = pd.read_csv(ATLAS, sep="\t", usecols=["integrase_id", "sample_id"])
    tax = pd.read_csv(TAX, sep="\t")
    tax["sample_id"] = tax["sample_id"].astype(str).str.strip()
    atlas["sample_id"] = atlas["sample_id"].astype(str).str.strip()
    cm = cm.merge(atlas, on="integrase_id", how="left")
    cm = cm.merge(
        tax[["sample_id", "domain", "phylum", "class", "order", "family"]].drop_duplicates("sample_id"),
        on="sample_id",
        how="left",
    )
    cm["phylum"] = cm["phylum"].fillna("Unknown").replace({"": "Unknown", "nan": "Unknown"})
    cm["phy_bin"] = cm["phylum"].where(cm["phylum"].isin(PHYLUM_ORDER + ["Unknown"]), "Others")
    cm["len_kb"] = cm["ie_len_nt"] / 1000
    cm["resolved"] = cm.resolved_minimal == 1
    return cm


def save(fig, stem: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{stem}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("->", FIG / f"{stem}.png")


def fig_overview(cm: pd.DataFrame) -> None:
    n = len(cm)
    n_res = int(cm.resolved.sum())
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.4))

    ax = axes[0, 0]
    bins = np.linspace(0, 30, 31)
    ax.hist(
        cm.loc[~cm.resolved, "len_kb"].clip(upper=30),
        bins=bins, color=SKY, edgecolor="white", linewidth=0.4,
        label=f"unresolved dark n={n - n_res}",
    )
    ax.hist(
        cm.loc[cm.resolved, "len_kb"].clip(upper=30),
        bins=bins, color=GREEN, edgecolor="white", linewidth=0.4,
        label=f"resolved n={n_res}",
    )
    ax.axvline(cm.len_kb.median(), color=VERM, ls="--", lw=1.2,
               label=f"median {cm.len_kb.median():.1f} kb")
    ax.set_xlabel("IE length (kb)")
    ax.set_ylabel("n composition-minimal")
    ax.set_title(f"Length ({n} composition-minimal after overlay)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    cds = cm.n_cds.dropna().clip(upper=30)
    ax.hist(cds, bins=np.arange(1, 32) - 0.5, color=GREEN, edgecolor="white", linewidth=0.4)
    ax.axvline(cm.n_cds.median(), color=VERM, ls="--", lw=1.2,
               label=f"median {cm.n_cds.median():.0f} CDS")
    ax.set_xlabel("n CDS")
    ax.set_ylabel("n IE")
    ax.set_title("Gene count")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    present = [c for c in CLASS_ORDER if (cm.char_class == c).any()]
    resolved = np.array([int(((cm.char_class == c) & cm.resolved).sum()) for c in present])
    pend = np.array([int(((cm.char_class == c) & ~cm.resolved).sum()) for c in present])
    y = np.arange(len(present))
    ax.barh(y, resolved, color=GREEN, height=0.72, label="resolved")
    ax.barh(y, pend, left=resolved, color=[CLASS_COLOR[c] for c in present],
            height=0.72, alpha=0.55, label="still dark CDS")
    ax.set_yticks(y)
    ax.set_yticklabels(present)
    ax.invert_yaxis()
    ax.set_xlabel("n composition-minimal")
    ax.set_title("By cargo class")
    for i, (a, b) in enumerate(zip(resolved, pend)):
        tot = a + b
        if tot:
            ax.text(tot + 2, i, int(tot), va="center", fontsize=8, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axes[1, 1]
    phy = cm.phy_bin.value_counts().reindex(PHYLUM_ORDER + ["Others", "Unknown"]).fillna(0)
    pcol = [PHYLUM_COLOR[p] for p in phy.index]
    ax.barh(range(len(phy)), phy.values, color=pcol)
    ax.set_yticks(range(len(phy)))
    ax.set_yticklabels(phy.index)
    ax.invert_yaxis()
    ax.set_xlabel("n composition-minimal")
    ax.set_title("Host phylum")
    for i, v in enumerate(phy.values):
        if v:
            ax.text(v + 1, i, int(v), va="center", fontsize=8, color=INK)

    fig.tight_layout()
    save(fig, "Figure_minimal_diversity")


def fig_hosts_by_class(cm: pd.DataFrame) -> None:
    main = ["replication", "TA", "RM", "PDC", "RT-defense"]
    bins = PHYLUM_ORDER + ["Others", "Unknown"]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    bottom = np.zeros(len(main))
    x = np.arange(len(main))
    for phy in bins:
        vals = []
        for cls in main:
            sub = cm[cm.char_class == cls]
            vals.append(int((sub.phy_bin == phy).sum()))
        ax.bar(x, vals, bottom=bottom, color=PHYLUM_COLOR[phy], width=0.72, label=phy)
        bottom += np.array(vals)
    ax.set_xticks(x)
    ax.set_xticklabels(main)
    ax.set_ylabel("n composition-minimal")
    ax.set_title("Host phylum within the five largest minimal classes")
    ax.legend(frameon=False, ncol=5, fontsize=8, loc="upper right")
    fig.tight_layout()
    save(fig, "Figure_minimal_hosts_by_class")


def fig_length_by_class(cm: pd.DataFrame) -> None:
    present = [c for c in CLASS_ORDER if (cm.char_class == c).any()]
    data = [cm.loc[cm.char_class == c, "len_kb"].clip(upper=30).values for c in present]
    fig, ax = plt.subplots(figsize=(9.4, 4.2))
    bp = ax.boxplot(data, tick_labels=present, patch_artist=True, flierprops={"markersize": 3})
    for patch, cls in zip(bp["boxes"], present):
        patch.set_facecolor(CLASS_COLOR[cls])
        patch.set_alpha(0.75)
    ax.set_ylabel("IE length (kb)")
    ax.set_title("Length spread of composition-minimal IEs (after overlay)")
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    fig.tight_layout()
    save(fig, "Figure_minimal_length_by_class")


def fig_delta(cm: pd.DataFrame) -> None:
    """Was-436 vs now-462 by class, plus resolved."""
    was = {
        "replication": 224, "TA": 103, "RM": 63, "PDC": 18, "RT-defense": 10,
        "Par": 5, "SoFic": 3, "CRISPR": 2, "AMR": 2, "Abi": 2,
        "metals": 1, "CBASS": 1, "Wadjet": 1, "Lamassu": 1,
    }
    present = [c for c in CLASS_ORDER if was.get(c, 0) or (cm.char_class == c).any()]
    now = [int((cm.char_class == c).sum()) for c in present]
    old = [was.get(c, 0) for c in present]
    res = [int(((cm.char_class == c) & cm.resolved).sum()) for c in present]
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    y = np.arange(len(present))
    ax.barh(y - 0.22, old, height=0.4, color=GREY, label="was 436 pool")
    ax.barh(y + 0.22, now, height=0.4, color=BLUE, label="now 462")
    ax.barh(y + 0.22, res, height=0.4, color=GREEN, label="resolved 46")
    ax.set_yticks(y)
    ax.set_yticklabels(present)
    ax.invert_yaxis()
    ax.set_xlabel("n composition-minimal")
    ax.set_title("Overlay vs locked 436 search pool")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    save(fig, "Figure_minimal_overlay_delta")


def main() -> None:
    cm = load()
    dest = OUT / "unknown_minimal_with_hosts.tsv"
    cm.to_csv(dest, sep="\t", index=False)
    print("n", len(cm), "resolved", int(cm.resolved.sum()),
          "phylum known", int((cm.phylum != "Unknown").sum()))
    fig_overview(cm)
    fig_hosts_by_class(cm)
    fig_length_by_class(cm)
    fig_delta(cm)


if __name__ == "__main__":
    main()
