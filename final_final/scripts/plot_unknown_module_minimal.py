#!/usr/bin/env python3
"""Histograms of module-minimal unknown IEs — same plate style as small RM."""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

ART = FF / "artifacts/unknown_ime_classes"
FIG = FF / "figures/unknown_ime_classes"
IE = ART / "after_hypothesis/unknown_minimal_ie.tsv"
PROT = FF / "artifacts/unknown_reannotate/proteins.tsv"
DARK = ART / "after_hypothesis/dark_cds_roles.tsv"
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
    "CRISPR", "Abi", "AMR", "CBASS", "Wadjet", "Lamassu", "metals",
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
ACC_ORDER = ["none", "C/HTH", "ArsR/CadC", "MerR", "TA partner", "Xis"]
ACC_COLOR = {
    "none": GREY,
    "C/HTH": "#44AA99",
    "ArsR/CadC": "#8c510a",
    "MerR": ORANGE,
    "TA partner": VERM,
    "Xis": BLUE,
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


def annotate_bars(ax, bars, ys, dy=0.15) -> None:
    ymax = max(ys) if len(ys) else 1
    for b, v in zip(bars, ys):
        if v <= 0:
            continue
        ax.text(b.get_x() + b.get_width() / 2, v + dy, str(int(v)),
                ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_ylim(0, ymax + max(1.4, 0.12 * ymax))


REG = re.compile(
    r"helix-turn-helix|hth |cro/c1|\bxre\b|c protein|"
    r"transcriptional regulator|dna-binding protein|"
    r"\barsr\b|\bsmtb\b|\bmerr\b|\bcadc\b|metalloregulator",
    re.I,
)
XIS = re.compile(r"excisionase|\bxis\b", re.I)


def acc_kind(prod: str) -> str:
    p = prod or ""
    if XIS.search(p):
        return "Xis"
    if re.search(r"arsr|smtb|cadc|metalloregulator", p, re.I):
        return "ArsR/CadC"
    if re.search(r"merr", p, re.I):
        return "MerR"
    if REG.search(p):
        return "C/HTH"
    return ""


def load() -> pd.DataFrame:
    ie = pd.read_csv(IE, sep="\t")
    prot = pd.read_csv(PROT, sep="\t", usecols=["integrase_id", "bakta_product"])
    kinds: dict[str, set[str]] = {}
    n_acc: dict[str, int] = {}
    for r in prot.itertuples():
        k = acc_kind(r.bakta_product)
        if k:
            kinds.setdefault(r.integrase_id, set()).add(k)
            n_acc[r.integrase_id] = n_acc.get(r.integrase_id, 0) + 1
    if DARK.exists():
        dark = pd.read_csv(DARK, sep="\t")
        acc_src = dark[dark.source.isin(["stack_acc", "leftover_family"])]
        for r in acc_src.itertuples():
            kinds.setdefault(r.integrase_id, set()).add("C/HTH")
    ie["accessory_kinds"] = ie.integrase_id.map(
        lambda i: "|".join(sorted(kinds.get(i, set())))
    )
    ie["n_accessory"] = ie.integrase_id.map(lambda i: n_acc.get(i, 0)).fillna(0)
    ie["n_unresolved"] = ie.n_unresolved_dark.fillna(0)
    ie["n_module"] = (ie.n_char_classes == 1).astype(int)
    ie["n_budget"] = (ie.n_cds.fillna(0) - 1 - ie.has_rep.fillna(0)).clip(lower=0)
    ie["module_minimal"] = ie.resolved_minimal.fillna(0).astype(int)
    ie["module_pending"] = (
        (ie.composition_minimal == 1) & (ie.resolved_minimal != 1)
    ).astype(int)
    ie["compact_minimal"] = (
        (ie.composition_minimal == 1) & (ie.n_budget <= 3)
    ).astype(int)
    atlas = pd.read_csv(ATLAS, sep="\t", usecols=["integrase_id", "sample_id"])
    tax = pd.read_csv(TAX, sep="\t")
    tax["sample_id"] = tax["sample_id"].astype(str).str.strip()
    atlas["sample_id"] = atlas["sample_id"].astype(str).str.strip()
    ie = ie.merge(atlas, on="integrase_id", how="left")
    ie = ie.merge(
        tax[["sample_id", "phylum"]].drop_duplicates("sample_id"),
        on="sample_id",
        how="left",
    )
    ie["phylum"] = ie["phylum"].fillna("Unknown").replace({"": "Unknown"})
    ie["len_kb"] = ie.ie_len_nt / 1000
    ie["proven"] = ie.module_minimal == 1
    ie["pending"] = (ie.compact_minimal == 1) & (ie.module_minimal != 1)
    return ie


def compact(df: pd.DataFrame) -> pd.DataFrame:
    """Same budget as small RM/TA: extras ≤ 3, int and Rep free."""
    return df[df.compact_minimal == 1].copy()


def fig_histograms(df: pd.DataFrame) -> None:
    n = len(df)
    n_p = int(df.proven.sum())
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2))

    ax = axes[0, 0]
    letter(ax, "A")
    present = [c for c in CLASS_ORDER if (df.char_class == c).any()]
    proven = np.array([int(((df.char_class == c) & df.proven).sum()) for c in present])
    pend = np.array([int(((df.char_class == c) & ~df.proven).sum()) for c in present])
    x = np.arange(len(present))
    ax.bar(x, proven, color=[CLASS_COLOR[c] for c in present], width=0.68, label="proven")
    ax.bar(x, pend, bottom=proven, color="#DDDDDD", edgecolor=[CLASS_COLOR[c] for c in present],
           linewidth=0.6, width=0.68, label="pending dark")
    for i, (a, b) in enumerate(zip(proven, pend)):
        tot = a + b
        if tot:
            ax.text(i, tot + 0.15, str(int(tot)), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(present, rotation=25, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title(f"Class  (compact n = {n}, proven {n_p})")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, max(proven + pend) + 2.2)

    ax = axes[0, 1]
    letter(ax, "B")
    counts = []
    for k in ACC_ORDER:
        if k == "none":
            counts.append(int((df.accessory_kinds.fillna("") == "").sum()))
        else:
            counts.append(int(df.accessory_kinds.fillna("").str.contains(k, regex=False).sum()))
    bars = ax.bar(range(len(ACC_ORDER)), counts,
                  color=[ACC_COLOR[k] for k in ACC_ORDER], width=0.68)
    annotate_bars(ax, bars, counts)
    ax.set_xticks(range(len(ACC_ORDER)))
    ax.set_xticklabels(ACC_ORDER, rotation=25, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title("Module accessories (named)")

    ax = axes[1, 0]
    letter(ax, "C")
    budgets = [0, 1, 2, 3]
    no = np.array([int(((df.n_budget.fillna(0) == b) & (df.has_rep == 0)).sum()) for b in budgets])
    yes = np.array([int(((df.n_budget.fillna(0) == b) & (df.has_rep == 1)).sum()) for b in budgets])
    x = np.arange(len(budgets))
    ax.bar(x, no, color=GREEN, width=0.62, label="no Rep")
    ax.bar(x, yes, bottom=no, color=PURPLE, width=0.62, label="+ Rep")
    for i, (a, b) in enumerate(zip(no, yes)):
        tot = a + b
        if tot:
            ax.text(i, tot + 0.12, str(int(tot)), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(["0 extra", "1 gene", "2 genes", "3 genes"])
    ax.set_ylabel("n IE")
    ax.set_ylim(0, max(no + yes) + 1.8)
    ax.set_title("Extra genes besides int (Rep free)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    letter(ax, "D")
    phy = df.phylum.value_counts()
    colors_p = {
        "Bacillota": BLUE, "Pseudomonadota": VERM, "Actinomycetota": ORANGE,
        "Acidobacteriota": GREEN, "Unknown": GREY,
    }
    bars = ax.barh(range(len(phy))[::-1], phy.values[::-1],
                   color=[colors_p.get(k, GREY) for k in phy.index[::-1]], height=0.62)
    for y, v in zip(range(len(phy))[::-1], phy.values[::-1]):
        ax.text(v + 0.12, y, str(int(v)), va="center", fontsize=9)
    ax.set_yticks(range(len(phy))[::-1])
    ax.set_yticklabels(phy.index[::-1])
    ax.set_xlabel("n IE")
    ax.set_xlim(0, phy.values.max() + 2.2)
    ax.set_title("Host phylum")

    fig.tight_layout(w_pad=2.2, h_pad=2.0)
    save(fig, "Figure_module_minimal_histograms")


def fig_budget(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.7))
    extra = (df.n_module.fillna(0) + df.n_accessory.fillna(0) + df.n_unresolved.fillna(0))

    ax = axes[0]
    letter(ax, "A")
    budgets = [1, 2, 3, 4, 5]
    no = np.array([int(((extra == b) & (df.has_rep == 0)).sum()) for b in budgets])
    yes = np.array([int(((extra == b) & (df.has_rep == 1)).sum()) for b in budgets])
    x = np.arange(len(budgets))
    ax.bar(x, no, color=GREEN, width=0.62, label="no Rep")
    ax.bar(x, yes, bottom=no, color=PURPLE, width=0.62, label="+ Rep")
    for i, (a, b) in enumerate(zip(no, yes)):
        tot = a + b
        if tot:
            ax.text(i, tot + 0.12, str(int(tot)), ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{b}" for b in budgets])
    ax.set_xlabel("module + accessory + dark")
    ax.set_ylabel("n IE")
    ax.set_ylim(0, max(no + yes) + 1.8 if len(no) else 1)
    ax.set_title("Genes besides int (Rep free)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    letter(ax, "B")
    kb = df.len_kb.clip(upper=16)
    bins = np.arange(1, 17, 1)
    ax.hist(kb[df.proven], bins=bins, color=GREEN, alpha=0.9, label="proven")
    ax.hist(kb[~df.proven], bins=bins, color=GREY, alpha=0.75, label="pending dark")
    ax.axvline(df.len_kb.median(), color=VERM, ls="--", lw=1.1,
               label=f"median {df.len_kb.median():.1f} kb")
    ax.set_xlabel("IE length (kb)")
    ax.set_ylabel("n IE")
    ax.set_title("Length")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout(w_pad=2.4)
    save(fig, "Figure_module_minimal_budget_length")


def fig_scan_target(df: pd.DataFrame) -> None:
    """All one-class no-leftover IEs: proven vs the scan pile."""
    shaped = df[(df.module_minimal == 1) | (df.module_pending == 1)].copy()
    shaped["proven"] = shaped.module_minimal == 1
    shaped["pending"] = shaped.module_pending == 1
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))

    ax = axes[0]
    letter(ax, "A")
    present = [c for c in CLASS_ORDER if (shaped.char_class == c).any()]
    proven = np.array([int(((shaped.char_class == c) & shaped.proven).sum()) for c in present])
    pend = np.array([int(((shaped.char_class == c) & shaped.pending).sum()) for c in present])
    y = np.arange(len(present))
    ax.barh(y, proven, color=GREEN, height=0.62, label=f"proven n={int(shaped.proven.sum())}")
    ax.barh(y, pend, left=proven, color=GREY, height=0.62, label=f"pending n={int(shaped.pending.sum())}")
    ax.set_yticks(y)
    ax.set_yticklabels(present)
    ax.invert_yaxis()
    ax.set_xlabel("n IE")
    ax.set_title("One class, no named leftover")
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axes[1]
    letter(ax, "B")
    pend = shaped[shaped.pending]
    ax.hist(pend.n_unresolved.clip(upper=20), bins=np.arange(1, 22) - 0.5,
            color=ORANGE, edgecolor="white", linewidth=0.4)
    ax.axvline(pend.n_unresolved.median(), color=VERM, ls="--", lw=1.1,
               label=f"median {pend.n_unresolved.median():.0f} dark")
    ax.set_xlabel("n unresolved CDS")
    ax.set_ylabel("n IE")
    ax.set_title("Scan target: dark leftover on pending")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout(w_pad=2.4)
    save(fig, "Figure_module_minimal_scan_target")


def main() -> None:
    ie = load()
    cap = compact(ie)
    dest = ART / "after_hypothesis/unknown_module_minimal_compact.tsv"
    cap.to_csv(dest, sep="\t", index=False)
    ie.to_csv(ART / "after_hypothesis/unknown_module_minimal.tsv", sep="\t", index=False)
    print("compact", len(cap), "proven", int(cap.proven.sum()), "pending", int((~cap.proven).sum()))
    print("shaped", int((ie.proven | ie.pending).sum()))
    fig_histograms(cap)
    fig_budget(cap)
    fig_scan_target(ie)


if __name__ == "__main__":
    main()
