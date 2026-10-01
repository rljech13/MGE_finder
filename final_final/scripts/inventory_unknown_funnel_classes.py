#!/usr/bin/env python3
"""Presence of literature funnel classes on all 4689 unknown IEs."""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

OUT = FF / "artifacts/unknown_ime_classes/after_hypothesis"
FIG = FF / "figures/unknown_ime_classes"
IE = OUT / "unknown_minimal_ie.tsv"
BY = OUT / "by_class.tsv"

BLUE = "#0072B2"
GREEN = "#009E73"
INK = "#222222"

ORDER = [
    "replication", "TA", "Par",
    "RM", "PDC", "CRISPR", "RT-defense", "SoFic", "Abi", "BREX",
    "Thoeris", "CBASS", "Rst_DUF4238", "Wadjet", "Shedu", "Lamassu",
    "metals", "AMR", "bacteriocins",
]

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})


def parse(s) -> set[str]:
    if pd.isna(s) or not str(s).strip():
        return set()
    return {x for x in str(s).split("|") if x and x != "other_defense"}


def main() -> None:
    ie = pd.read_csv(IE, sep="\t")
    by = pd.read_csv(BY, sep="\t")
    cls = ie.classes.map(parse)
    n = len(ie)
    has = cls.map(bool)
    only_rep = cls.map(lambda s: s == {"replication"})
    nonrep = cls.map(lambda s: bool(s - {"replication"}))
    no = cls.map(lambda s: not s)
    A_nr = cls.map(lambda s: bool(s & {"TA", "Par"}))
    B = {"RM", "CRISPR", "Abi", "RT-defense", "SoFic", "CBASS", "Shedu",
         "Thoeris", "BREX", "Wadjet", "Lamassu", "PDC", "Rst_DUF4238"}
    C = {"AMR", "metals", "bacteriocins"}
    has_B = cls.map(lambda s: bool(s & B))
    has_C = cls.map(lambda s: bool(s & C))

    summary = pd.DataFrame([
        {"metric": "unknown_n", "n": n},
        {"metric": "any_funnel_class", "n": int(has.sum())},
        {"metric": "no_funnel_class", "n": int(no.sum())},
        {"metric": "only_replication", "n": int(only_rep.sum())},
        {"metric": "any_nonrep_class", "n": int(nonrep.sum())},
        {"metric": "maintenance_TA_or_Par", "n": int(A_nr.sum())},
        {"metric": "defense_layer", "n": int(has_B.sum())},
        {"metric": "host_adaptive_layer", "n": int(has_C.sum())},
        {"metric": "exactly_one_char_class", "n": int(ie.class_minimal.sum())},
        {"metric": "two_or_more_char_classes", "n": int((ie.n_char_classes >= 2).sum())},
    ])
    summary.to_csv(OUT / "class_inventory_summary.tsv", sep="\t", index=False)

    only_this = {}
    for c in by["class"]:
        only_this[c] = int(cls.map(lambda s, c=c: s == {c}).sum())
    inv = by[["class", "layer", "n_ie_with_class", "n_class_minimal",
              "n_composition_minimal", "genes_included", "presence_lock"]].copy()
    inv["n_only_this_module"] = inv["class"].map(only_this)
    inv.to_csv(OUT / "class_inventory.tsv", sep="\t", index=False)

    present = [c for c in ORDER if int(by.loc[by["class"] == c, "n_ie_with_class"].iloc[0]) > 0]
    vals = [int(by.loc[by["class"] == c, "n_ie_with_class"].iloc[0]) for c in present]
    solo = [only_this[c] for c in present]

    fig, ax = plt.subplots(figsize=(10.2, 5.4))
    x = range(len(present))
    ax.bar(x, vals, color=BLUE, width=0.72, label="IE with class")
    ax.bar(x, solo, color=GREEN, width=0.42,
           label="only this module (+ unannotated)")
    for i, v in enumerate(vals):
        ax.text(i, v + 18, str(v), ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_xticks(list(x))
    ax.set_xticklabels(present, rotation=40, ha="right")
    ax.set_ylabel("n IE / 4689")
    ax.set_title("Unknown IE — funnel class presence (n = 4689)")
    ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(0, max(vals) + 140)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        dest = FIG / f"Figure_unknown_class_inventory.{ext}"
        fig.savefig(dest, dpi=200 if ext == "png" else None, bbox_inches="tight")
        print("->", dest)
    plt.close(fig)
    print(summary.to_string(index=False))
    print(inv[["class", "n_ie_with_class", "n_only_this_module"]].to_string(index=False))


if __name__ == "__main__":
    main()
