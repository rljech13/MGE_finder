#!/usr/bin/env python3
"""Integrase diversity panel for Stage 8 unclassified Leiden clusters.

Two complementary views:
  A) Forming-integrase HMM subtype (PF00589.27 vs PF22022.2) — discovery axis
  B) Secondary serine-resolvase Pfam (PF02796) co-occurrence in feature matrix

Output:
  figures/04_cluster_unknown_6960/03b_integrase_diversity_by_cluster.{png,svg}
  artifacts/nb04_cluster_unknown_6960/integrase_diversity_by_cluster_c6960_unk.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from scipy.sparse import load_npz

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import FF, REPO
ART = FF / "artifacts/nb04_cluster_unknown_6960"
FIG = FF / "figures/04_cluster_unknown_6960"
PROJ = REPO / "results_organized/genome_clustering/project/genome_matrix"
PATH_TRUTH = (
    REPO / "results_organized/genome_clustering/project/notebook7_extended"
    / "integrase_truth_table_with_ie_atlas.tsv"
)

HMM_COLORS = {
    "PF00589.27": "#3182bd",   # Phage_integrase / RVE core
    "PF22022.2":  "#e6550d",   # Phage_int_M central domain
}
SER_COLORS = {
    "tyrosine_only": "#6baed6",
    "tyrosine+serine": "#756bb1",
    "serine_only": "#e6550d",
}

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def _load_serine_pfam_mask(n_rows: int) -> np.ndarray:
    sys.path.insert(0, str(FF / "03_clustering" / "lib"))
    from cluster_biology import assign_role

    lab = np.load(PROJ / "genome_feature_labels.npz", allow_pickle=True)
    feature_names = np.array([str(f) for f in lab["features"]])
    Xs = load_npz(ART / "feature_subset_c6960_unk.npz")
    ser_cols = [
        i for i, fn in enumerate(feature_names)
        if fn.startswith("pfam_")
        and assign_role(fn[5:].split(".")[0], "", "") == "integrase_serine"
    ]
    return (np.asarray((Xs > 0).astype(np.float32)[:, ser_cols].sum(axis=1)).flatten() > 0)


def build_per_ie_table() -> pd.DataFrame:
    rows = pd.read_csv(ART / "ie_rep_rows_c6960_unk.tsv", sep="\t")
    leid = pd.read_csv(ART / "leiden_clusters_c6960_unk.tsv", sep="\t")
    truth = pd.read_csv(PATH_TRUTH, sep="\t", usecols=["integrase_id_lookup", "integrase_model"])
    truth = truth.drop_duplicates("integrase_id_lookup").rename(
        columns={"integrase_id_lookup": "integrase_id"}
    )

    df = rows.merge(leid, on=["integrase_id", "rep_id"], how="left")
    df = df.merge(truth, on="integrase_id", how="left")
    df["has_serine_pfam"] = _load_serine_pfam_mask(len(df))

    def serine_class(row):
        if row["has_serine_pfam"]:
            return "tyrosine+serine" if str(row["integrase_model"]).startswith(("PF00589", "PF22022")) else "serine_only"
        return "tyrosine_only"

    df["serine_class"] = df.apply(serine_class, axis=1)
    return df


def cluster_summary(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for cl, g in df.groupby("leiden_cluster"):
        if pd.isna(cl) or int(cl) < 0:
            continue
        n = len(g)
        if n < 5:
            continue
        records.append({
            "cluster": int(cl),
            "n_ie": n,
            "n_PF00589": int((g.integrase_model == "PF00589.27").sum()),
            "n_PF22022": int((g.integrase_model == "PF22022.2").sum()),
            "frac_PF00589": round((g.integrase_model == "PF00589.27").mean(), 3),
            "frac_PF22022": round((g.integrase_model == "PF22022.2").mean(), 3),
            "n_serine_pfam": int(g.has_serine_pfam.sum()),
            "frac_serine_pfam": round(g.has_serine_pfam.mean(), 3),
            "n_tyr_only": int((g.serine_class == "tyrosine_only").sum()),
            "n_tyr_plus_ser": int((g.serine_class == "tyrosine+serine").sum()),
        })
    return pd.DataFrame(records).sort_values("cluster")


def _stacked_bar(ax, mat: pd.DataFrame, colors: dict[str, str], *, title: str, ylabel: str):
    x = np.arange(len(mat))
    bottom = np.zeros(len(mat))
    for col in mat.columns:
        vals = mat[col].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=colors[col], edgecolor="white", linewidth=0.6, label=col)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([f"C{int(c)}" for c in mat.index], fontsize=9)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10, pad=8)
    ax.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
    ax.set_axisbelow(True)


def render_figure(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    clusters = summary.sort_values("cluster")["cluster"].tolist()
    sz = summary.set_index("cluster").loc[clusters, "n_ie"]

    hmm_mat = summary.set_index("cluster").loc[clusters, ["frac_PF00589", "frac_PF22022"]]
    hmm_mat.columns = ["PF00589.27", "PF22022.2"]

    ser_mat = summary.set_index("cluster").loc[clusters, ["n_tyr_only", "n_tyr_plus_ser"]].astype(float)
    ser_mat = ser_mat.div(sz.values, axis=0)
    ser_mat.columns = ["tyrosine_only", "tyrosine+serine"]

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.2), gridspec_kw={"wspace": 0.28})

    _stacked_bar(
        axes[0], hmm_mat, HMM_COLORS,
        title="A. Forming-integrase HMM subtype (Stage 1)",
        ylabel="IE fraction",
    )
    _stacked_bar(
        axes[1], ser_mat, SER_COLORS,
        title="B. Serine-resolvase Pfam (PF02796) co-occurrence",
        ylabel="IE fraction",
    )

    fig.suptitle(
        "Integrase diversity · 6960 Unknown (Leiden)",
        fontsize=12, weight="bold", y=1.02,
    )
    fig.text(
        0.5, -0.02,
        "A: HMM model of the forming integrase (every IE is tyrosine-class by Stage 1). "
        "B: an extra serine-resolvase domain in the Pfam feature matrix (~10% of IE; C12 ≈ 100%).",
        ha="center", va="top", fontsize=8, color="#444",
    )
    fig.tight_layout()
    out = FIG / "03b_integrase_diversity_by_cluster.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", pad_inches=0.30)
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.30)
    plt.close(fig)
    print(f"→ {out}")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    df = build_per_ie_table()
    summary = cluster_summary(df)
    out_tsv = ART / "integrase_diversity_by_cluster_c6960_unk.tsv"
    summary.to_csv(out_tsv, sep="\t", index=False)
    print(f"→ {out_tsv}")
    print(summary[["cluster", "n_ie", "frac_PF00589", "frac_PF22022", "frac_serine_pfam"]].to_string(index=False))
    render_figure(df, summary)


if __name__ == "__main__":
    main()
