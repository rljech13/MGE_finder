#!/usr/bin/env python3
"""Histograms for the locked minimal set: one class, no leftover, ≤1 unannotated ≥70 aa.
Non-obligatory DUF counts as unannotated. DUF + another unannotated ≥70 aa is not allowed.
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

OUT = FF / "artifacts/unknown_ime_classes/after_hypothesis"
FIG = FF / "figures/unknown_ime_classes"
IE = OUT / "unknown_minimal_le1_long.tsv"

BLUE = "#0072B2"
GREEN = "#009E73"
VERM = "#D55E00"
INK = "#222222"

CLASS_ORDER = [
    "replication", "TA", "RM", "PDC", "RT-defense", "AMR", "Par", "metals",
]

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


def letter(ax, s: str) -> None:
    ax.text(-0.12, 1.08, s, transform=ax.transAxes, fontsize=14, fontweight="bold",
            va="bottom", color=INK)


def main() -> None:
    ie = pd.read_csv(IE, sep="\t")
    m = ie[ie.minimal_le1_long == 1].copy()
    m["len_kb"] = m.ie_len_nt / 1000
    n = len(m)
    col = "n_unannotated_ge70" if "n_unannotated_ge70" in m.columns else "n_unresolved_ge70"
    n0 = int((m[col] == 0).sum())
    n1 = int((m[col] == 1).sum())

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2))

    ax = axes[0, 0]
    letter(ax, "A")
    bins = np.linspace(0, 16, 17)
    ax.hist(m.loc[m[col] == 0, "len_kb"].clip(upper=16),
            bins=bins, color=GREEN, edgecolor="white", linewidth=0.4,
            label=f"0 unannotated ≥70 aa  n={n0}")
    ax.hist(m.loc[m[col] == 1, "len_kb"].clip(upper=16),
            bins=bins, color=BLUE, edgecolor="white", linewidth=0.4, alpha=0.85,
            label=f"1 unannotated ≥70 aa  n={n1}")
    ax.axvline(m.len_kb.median(), color=VERM, ls="--", lw=1.2,
               label=f"median {m.len_kb.median():.1f} kb")
    ax.set_xlabel("IE length (kb)")
    ax.set_ylabel("n IE")
    ax.set_title(f"Length  (n = {n})")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    letter(ax, "B")
    cds = m.n_cds.dropna().clip(upper=16)
    ax.hist(cds, bins=np.arange(1, 18) - 0.5, color=GREEN,
            edgecolor="white", linewidth=0.4)
    ax.axvline(m.n_cds.median(), color=VERM, ls="--", lw=1.2,
               label=f"median {m.n_cds.median():.0f} CDS")
    ax.set_xlabel("n CDS")
    ax.set_ylabel("n IE")
    ax.set_title("Gene count")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    letter(ax, "C")
    present = [c for c in CLASS_ORDER if (m.char_class == c).any()]
    z0 = np.array([int(((m.char_class == c) & (m[col] == 0)).sum()) for c in present])
    z1 = np.array([int(((m.char_class == c) & (m[col] == 1)).sum()) for c in present])
    x = np.arange(len(present))
    ax.bar(x, z0, color=GREEN, width=0.68, label="0 unannotated ≥70 aa")
    ax.bar(x, z1, bottom=z0, color=BLUE, width=0.68, label="1 unannotated ≥70 aa")
    for i, tot in enumerate(z0 + z1):
        ax.text(i, tot + 0.8, str(int(tot)), ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(present, rotation=30, ha="right")
    ax.set_ylabel("n IE")
    ax.set_ylim(0, max(z0 + z1) + 8)
    ax.set_title("By cargo class")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    letter(ax, "D")
    counts = [n0, n1]
    bars = ax.bar([0, 1], counts, color=[GREEN, BLUE], width=0.62)
    for b, v in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.5, str(int(v)),
                ha="center", va="bottom", fontsize=10, color=INK)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["0 unannotated ≥70 aa", "1 unannotated ≥70 aa"])
    ax.set_ylabel("n IE")
    ax.set_ylim(0, max(counts) + 12)
    ax.set_title("Unannotated CDS ≥70 aa")

    fig.tight_layout(w_pad=2.2, h_pad=2.0)
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        dest = FIG / f"Figure_minimal_histograms.{ext}"
        fig.savefig(dest, dpi=200 if ext == "png" else None, bbox_inches="tight")
        print("->", dest)
    plt.close(fig)


if __name__ == "__main__":
    main()
