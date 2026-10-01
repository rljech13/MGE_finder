#!/usr/bin/env python3
"""Per-genome fold when a funnel class is present in that genome's unknown IE.

fold_g = (n_IE / L_IE) / (n_rest / L_rest)  on that genome only.
n_rest == 0 → only in IE; fold_rest empty, fold vs own genome = L_genome / L_IE.
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

OUT = FF / "artifacts/unknown_ime_classes/after_hypothesis/genome_enrichment"
FIG = FF / "figures/unknown_ime_classes"

BLUE = "#0072B2"
GREEN = "#009E73"
VERM = "#D55E00"
INK = "#222222"

LAYER = {
    "replication": "A", "TA": "A", "Par": "A",
    "RM": "B", "CRISPR": "B", "Abi": "B", "RT-defense": "B", "SoFic": "B",
    "CBASS": "B", "Shedu": "B", "Thoeris": "B", "BREX": "B", "Wadjet": "B",
    "Lamassu": "B", "PDC": "B", "Rst_DUF4238": "B",
    "AMR": "C", "metals": "C", "bacteriocins": "C",
}
ORDER = [
    "replication", "TA", "Par",
    "RM", "CRISPR", "Abi", "RT-defense", "SoFic", "CBASS", "Shedu",
    "Thoeris", "BREX", "Wadjet", "Lamassu", "PDC", "Rst_DUF4238",
    "AMR", "metals", "bacteriocins",
]

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})


def main() -> None:
    hits = pd.read_csv(OUT / "hits_assigned.tsv", sep="\t")
    lens = pd.read_csv(OUT / "genome_lengths.tsv", sep="\t")
    hits = hits[hits["class"].isin(ORDER)].copy()
    hits["in_ie"] = hits.comp.eq("unk_ie")
    gcls = (
        hits.groupby(["gca", "class"], as_index=False)
        .agg(n_ie=("in_ie", "sum"), n_all=("in_ie", "size"))
    )
    gcls["n_rest"] = gcls.n_all - gcls.n_ie
    m = gcls.merge(lens, on="gca", how="left")
    m = m[(m.L_unk_ie > 0) & (m.L_rest > 0) & (m.genome_nt > 0) & (m.n_ie > 0)].copy()
    m["layer"] = m["class"].map(LAYER)
    m["rho_ie"] = m.n_ie / (m.L_unk_ie / 1e6)
    m["rho_rest"] = np.where(m.n_rest > 0, m.n_rest / (m.L_rest / 1e6), np.nan)
    m["fold_vs_rest"] = m.rho_ie / m.rho_rest
    m["only_in_ie"] = m.n_rest.eq(0)
    m["fold_vs_own_genome"] = m.rho_ie / (m.n_all / (m.genome_nt / 1e6))
    m["expected_ie"] = m.n_all * (m.L_unk_ie / m.genome_nt)
    long_cols = [
        "gca", "class", "layer", "n_ie", "n_rest", "n_all",
        "L_unk_ie", "L_rest", "genome_nt", "n_unk_ie",
        "rho_ie", "rho_rest", "fold_vs_rest", "fold_vs_own_genome",
        "expected_ie", "only_in_ie",
    ]
    m[long_cols].to_csv(OUT / "per_genome_when_in_ie.tsv", sep="\t", index=False)

    recs = []
    for cls in ORDER:
        sub = m[m["class"] == cls]
        if sub.empty:
            continue
        fin = sub.loc[~sub.only_in_ie, "fold_vs_rest"]
        recs.append({
            "class": cls,
            "layer": LAYER[cls],
            "n_genomes_class_in_IE": int(len(sub)),
            "n_only_in_IE": int(sub.only_in_ie.sum()),
            "n_also_in_rest": int((~sub.only_in_ie).sum()),
            "median_fold_vs_rest": float(fin.median()) if len(fin) else np.nan,
            "q25_fold_vs_rest": float(fin.quantile(0.25)) if len(fin) else np.nan,
            "q75_fold_vs_rest": float(fin.quantile(0.75)) if len(fin) else np.nan,
            "frac_fold_gt_1": float((fin > 1).mean()) if len(fin) else np.nan,
            "median_fold_vs_own_genome": float(sub.fold_vs_own_genome.median()),
            "q25_fold_vs_own_genome": float(sub.fold_vs_own_genome.quantile(0.25)),
            "q75_fold_vs_own_genome": float(sub.fold_vs_own_genome.quantile(0.75)),
            "frac_own_gt_1": float((sub.fold_vs_own_genome > 1).mean()),
        })
    summ = pd.DataFrame(recs)
    summ.to_csv(OUT / "per_genome_when_in_ie_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    plot(m, summ)


def plot(m: pd.DataFrame, summ: pd.DataFrame) -> None:
    show = summ[summ.n_genomes_class_in_IE >= 1].copy()
    show = show.sort_values("median_fold_vs_own_genome")
    colors = [{"A": BLUE, "B": VERM, "C": GREEN}[L] for L in show.layer]
    data, labels = [], []
    for cls in show["class"]:
        vals = m.loc[m["class"] == cls, "fold_vs_own_genome"].to_numpy()
        data.append(vals)
        n = int(show.loc[show["class"] == cls, "n_genomes_class_in_IE"].iloc[0])
        only = int(show.loc[show["class"] == cls, "n_only_in_IE"].iloc[0])
        extra = f"  only IE={only}" if only else ""
        labels.append(f"{cls}  n={n}{extra}")

    fig, ax = plt.subplots(figsize=(8.8, 6.8))
    bp = ax.boxplot(
        data, vert=False, widths=0.62, patch_artist=True,
        medianprops={"color": INK, "lw": 1.2},
        whiskerprops={"color": INK}, capprops={"color": INK},
        flierprops={"marker": "o", "ms": 2.5, "alpha": 0.35, "color": INK},
        showfliers=True,
    )
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
        patch.set_edgecolor(INK)
    ax.axvline(1, color=INK, ls="--", lw=1.0)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("per-genome fold = (hits/Mb in unknown IE) / (hits/Mb in that genome)")
    ax.set_title("Fold on genomes where the class is in the unknown IE")
    ax.plot([], [], color=BLUE, lw=6, label="A maintenance")
    ax.plot([], [], color=VERM, lw=6, label="B defense")
    ax.plot([], [], color=GREEN, lw=6, label="C host-adaptive")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        dest = FIG / f"Figure_unknown_ie_per_genome_enrichment.{ext}"
        fig.savefig(dest, dpi=200 if ext == "png" else None, bbox_inches="tight")
        print("->", dest)
    plt.close(fig)


if __name__ == "__main__":
    main()
