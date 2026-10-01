#!/usr/bin/env python3
"""Histograms of class and cargo for small T/AT pIEs.

Set: n_ta + n_other ≤ 3; integrase and Rep do not count.
Unknown only. Role is taken from product names (not the scan 'role' column:
'toxin-antitoxin' in a toxin name used to look like an antitoxin).
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

ART = FF / "artifacts/unknown_le3_ta"
FIG = FF / "figures/unknown_le3_ta"
HITS = ART / "ta_allow_rep.tsv"
GENES = ART / "ta_allow_rep_genes.tsv"
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

TA_ORDER = ["HicAB", "RelBE", "MazEF", "VapBC", "AbiE", "unnamed"]
TA_COLOR = {
    "HicAB": GREEN,
    "RelBE": VERM,
    "MazEF": ORANGE,
    "VapBC": SKY,
    "AbiE": PINK,
    "unnamed": GREY,
}
PAIR_ORDER = ["named pair", "likely pair", "toxin only", "antitoxin only"]
PAIR_COLOR = {
    "named pair": GREEN,
    "likely pair": ORANGE,
    "toxin only": VERM,
    "antitoxin only": SKY,
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


def product_role(prod: str) -> str:
    p = (prod or "").lower()
    p2 = p.replace("toxin-antitoxin", "TA").replace("toxin antitoxin", "TA")
    if "abieii" in p2 or "abigii" in p2:
        return "toxin"
    if "abiei" in p2:
        return "antitoxin"
    if "antitoxin" in p2 or "addiction module" in p2:
        return "antitoxin"
    if "toxin" in p2 or "vapc" in p2 or "mazf" in p2 or "pemk" in p2:
        return "toxin"
    return "ta"


def ta_family(blob: str) -> str:
    b = blob.lower()
    if "hic" in b:
        return "HicAB"
    if "rel" in b or "pare" in b or "addiction module" in b:
        return "RelBE"
    if "maz" in b or "pemk" in b:
        return "MazEF"
    if "vapc" in b or "vapb" in b:
        return "VapBC"
    if "abie" in b or "abigii" in b:
        return "AbiE"
    return "unnamed"


def pair_status(r, roles: list[str]) -> str:
    has_t = "toxin" in roles
    has_a = "antitoxin" in roles
    oth = (str(r.other_products) or "").lower()
    pad = str(r.padloc or "")
    if has_t and has_a:
        return "named pair"
    if "AbiE" in pad and (has_t or has_a):
        return "named pair"
    # adjacent regulator / small hyp as likely missing partner
    if has_t and (
        "helix-turn-helix" in oth
        or "transcriptional regulator" in oth
        or "hth " in oth
        or "hypothetical protein(" in oth
    ):
        return "likely pair"
    if has_a and not has_t:
        return "antitoxin only"
    if has_t and not has_a:
        return "toxin only"
    return "toxin only"


def architecture(r, genes: pd.DataFrame) -> str:
    """attL-left sketch from native order when int is on the right; else reverse."""
    g = genes[genes.integrase_id == r.integrase_id].copy()
    if g.empty:
        return ""
    g = g.sort_values("start")
    is_int = g["kind"].eq("int") | g["product"].str.contains(
        "integrase|recombinase|HNH", case=False, regex=True
    )
    int_mid = float(g.loc[is_int, "start"].max() if is_int.any() else g["start"].max())
    parts = []
    for x in g.itertuples():
        kind = x.kind
        prod = x.product
        if kind == "int" or ("integrase" in prod.lower() or "recombinase" in prod.lower()):
            lab = "int"
        elif kind == "ta":
            role = product_role(prod)
            fam = ta_family(prod)
            if role == "toxin":
                lab = {"HicAB": "HicA", "RelBE": "RelE", "MazEF": "MazF",
                       "VapBC": "VapC", "AbiE": "AbiEii", "unnamed": "toxin"}[fam]
            elif role == "antitoxin":
                lab = {"HicAB": "HicB", "RelBE": "AT", "MazEF": "MazE",
                       "VapBC": "VapB", "AbiE": "AbiEi", "unnamed": "AT"}[fam]
            else:
                lab = fam if fam != "unnamed" else "TA"
        elif kind == "rep":
            lab = "Rep"
        else:
            pl = prod.lower()
            if "regulator" in pl or "helix-turn-helix" in pl or "hth" in pl:
                lab = "regulator"
            elif "hypothetical" in pl:
                lab = "hyp"
            elif "hnh" in pl:
                lab = "int"
            else:
                lab = prod.split()[0][:10]
        parts.append((int(x.start), lab))
    # if integrase is on the left half, reverse so attL-left ≈ cargo then int
    if int_mid < 0.4 * (g["end"].max() or 1):
        parts = list(reversed(parts))
    # collapse consecutive int if HNH already labeled int
    labels = [p[1] for p in parts]
    out = []
    for lab in labels:
        if out and out[-1] == lab == "int":
            continue
        out.append(lab)
    return "attL — " + " — ".join(out) + " — attR"


def classify() -> pd.DataFrame:
    hits = pd.read_csv(HITS, sep="\t")
    hits = hits[hits.func_class == "unknown"].copy()
    genes = pd.read_csv(GENES, sep="\t")
    atlas = pd.read_csv(ATLAS, sep="\t", usecols=["integrase_id", "sample_id", "primary_label"])
    tax = pd.read_csv(TAX, sep="\t")
    hits = hits.merge(atlas, on="integrase_id", how="left")
    hits = hits.merge(tax, on="sample_id", how="left")

    fams, pairs, roles_col, arch = [], [], [], []
    for r in hits.itertuples():
        sub = genes[genes.integrase_id == r.integrase_id]
        ta = sub[sub.kind == "ta"]
        roles = [product_role(p) for p in ta["product"]]
        blob = f"{r.ta_products} {r.padloc}"
        fams.append(ta_family(blob))
        pairs.append(pair_status(r, roles))
        roles_col.append(",".join(roles))
        arch.append(architecture(r, genes))
    hits["ta_class"] = fams
    hits["pair"] = pairs
    hits["ta_roles"] = roles_col
    hits["architecture_attL_left"] = arch
    hits["has_rep"] = hits.n_rep > 0
    hits["phylum"] = hits.phylum.fillna("Unknown")
    dest = ART / "ta_allow_rep_classified.tsv"
    hits.to_csv(dest, sep="\t", index=False)
    loci = hits[[
        "integrase_id", "ta_class", "pair", "n_budget", "n_rep", "n_ta", "n_other",
        "ie_len_nt", "padloc", "in_unknown_le3", "architecture_attL_left",
        "ta_products", "other_products", "phylum", "family",
    ]]
    loci.to_csv(ART / "small_ta_loci.tsv", sep="\t", index=False)
    print("n=", len(hits), "→", dest)
    print(hits[["integrase_id", "ta_class", "pair", "n_budget", "ie_len_nt",
                "architecture_attL_left"]].to_string(index=False))
    return hits


def fig_histograms(df: pd.DataFrame) -> None:
    n = len(df)
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.2))

    ax = axes[0, 0]
    letter(ax, "A")
    counts = df.ta_class.value_counts().reindex(TA_ORDER, fill_value=0)
    bars = ax.bar(range(len(TA_ORDER)), counts.values,
                  color=[TA_COLOR[k] for k in TA_ORDER], width=0.68)
    ymax = max(counts.values) if len(counts) else 1
    for b, v in zip(bars, counts.values):
        if v:
            ax.text(b.get_x() + b.get_width() / 2, v + 0.12, str(int(v)),
                    ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_ylim(0, ymax + 1.3)
    ax.set_xticks(range(len(TA_ORDER)))
    ax.set_xticklabels(TA_ORDER, rotation=20, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title(f"T/AT family  (n = {n})")

    ax = axes[0, 1]
    letter(ax, "B")
    pc = df.pair.value_counts().reindex(PAIR_ORDER, fill_value=0)
    bars = ax.bar(range(len(PAIR_ORDER)), pc.values,
                  color=[PAIR_COLOR[k] for k in PAIR_ORDER], width=0.68)
    ymax = max(pc.values) if len(pc) else 1
    for b, v in zip(bars, pc.values):
        if v:
            ax.text(b.get_x() + b.get_width() / 2, v + 0.12, str(int(v)),
                    ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_ylim(0, ymax + 1.3)
    ax.set_xticks(range(len(PAIR_ORDER)))
    ax.set_xticklabels(PAIR_ORDER, rotation=20, ha="right")
    ax.set_ylabel("n IE")
    ax.set_title("Pair completeness")

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

    ax = axes[1, 1]
    letter(ax, "D")
    phy = df.phylum.value_counts()
    colors_p = {
        "Bacillota": BLUE,
        "Pseudomonadota": VERM,
        "Myxococcota": GREEN,
        "Actinomycetota": ORANGE,
        "Bdellovibrionota": PINK,
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
    save(fig, "Figure_small_ta_histograms")


def fig_budget(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.7))

    ax = axes[0]
    letter(ax, "A")
    budgets = [1, 2, 3]
    no = np.array([int(((df.n_budget == b) & ~df.has_rep).sum()) for b in budgets])
    yes = np.array([int(((df.n_budget == b) & df.has_rep).sum()) for b in budgets])
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
    bins = np.arange(1.5, 5.0, 0.5)
    ax.hist(kb[~df.has_rep], bins=bins, color=GREEN, alpha=0.9, label="no Rep")
    if df.has_rep.any():
        ax.hist(kb[df.has_rep], bins=bins, color=PURPLE, alpha=0.85, label="+ Rep")
    ax.set_xlabel("pIE length (kb)")
    ax.set_ylabel("n IE")
    ax.set_title(f"Length  (median {kb.median():.1f} kb)")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout(w_pad=2.4)
    save(fig, "Figure_small_ta_budget_length")


def main() -> None:
    df = classify()
    fig_histograms(df)
    fig_budget(df)


if __name__ == "__main__":
    main()
