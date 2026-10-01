#!/usr/bin/env python3
"""PFAM marker heatmap per Leiden cluster — Stage 9 unclassified (same style as Stage 8)."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from scipy.sparse import load_npz, save_npz

from cohort6960_unknown_paths import ART, FIG, PATH_LAB, PATH_X, TAG

PATH_MARKERS = ART / f"cluster_pfam_markers_{TAG}.tsv"
PATH_XS = ART / f"feature_subset_{TAG}.npz"
PATH_LEIDEN = ART / f"leiden_clusters_{TAG}.tsv"
PATH_IE_ROWS = ART / f"ie_rep_rows_{TAG}.tsv"
OUT_PNG = FIG / "03b_heatmap_pfam_markers_fraction.png"

MIN_PIN_LABEL = 0.05

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def _clean_short(name: str, pfam_acc: str) -> str:
    s = str(name or "").strip()
    if s.startswith("(") and "'" in s:
        parts = s.replace("(", "").replace(")", "").replace("'", "").split(",")
        if parts:
            return parts[0].strip()
    if not s or s == "nan":
        return str(pfam_acc)
    return s


def ensure_feature_subset() -> None:
    if PATH_XS.is_file():
        return
    ie = pd.read_csv(PATH_IE_ROWS, sep="\t", usecols=["matrix_row"])
    X = load_npz(PATH_X).astype(np.float32)
    Xs = X[ie["matrix_row"].astype(int).to_numpy()]
    save_npz(PATH_XS, Xs)
    print(f"→ {PATH_XS}  shape={Xs.shape} nnz={Xs.nnz:,}")


def _align_leiden_to_subset() -> np.ndarray:
    leid = pd.read_csv(PATH_LEIDEN, sep="\t")
    ie = pd.read_csv(PATH_IE_ROWS, sep="\t", usecols=["integrase_id"])
    if len(leid) != len(ie):
        leid = leid.set_index("integrase_id").reindex(ie["integrase_id"].astype(str)).reset_index()
    if leid["leiden_cluster"].isna().any():
        raise ValueError("leiden_clusters and ie_rep_rows integrase_id order mismatch")
    return leid["leiden_cluster"].to_numpy()


def _marker_matrix_full() -> tuple[pd.DataFrame, pd.Series, set[tuple[int, str]]]:
    markers = pd.read_csv(PATH_MARKERS, sep="\t")
    lab = np.load(PATH_LAB, allow_pickle=True)
    f2i = {str(f): i for i, f in enumerate(lab["features"])}

    ensure_feature_subset()
    Xs = load_npz(PATH_XS)
    leiden_labels = _align_leiden_to_subset()
    if Xs.shape[0] != len(leiden_labels):
        raise ValueError(f"Xs rows {Xs.shape[0]} != leiden {len(leiden_labels)}")

    clusters = sorted(int(c) for c in np.unique(leiden_labels))
    Xs_bin = (Xs > 0).astype(np.float32)

    marker_feats: list[tuple[str, str]] = []
    seen_feat: set[str] = set()
    for r in markers.itertuples():
        fid = str(r.feature_id)
        if fid not in f2i or fid in seen_feat:
            continue
        seen_feat.add(fid)
        marker_feats.append((_clean_short(r.short_name, r.pfam_acc), fid))

    rank1 = (
        markers.sort_values(["leiden_cluster", "marker_score"], ascending=[True, False])
        .groupby("leiden_cluster")
        .first()
    )
    marker_order: list[str] = []
    seen_m: set[str] = set()
    for cl in clusters:
        if cl not in rank1.index:
            continue
        row = rank1.loc[cl]
        m = _clean_short(row.short_name, row.pfam_acc)
        if m not in seen_m:
            marker_order.append(m)
            seen_m.add(m)
    for m, _ in marker_feats:
        if m not in seen_m:
            marker_order.append(m)
            seen_m.add(m)

    feat_by_marker = {m: fid for m, fid in marker_feats}
    mat = np.zeros((len(marker_order), len(clusters)), dtype=np.float32)
    for j, cl in enumerate(clusters):
        mask = leiden_labels == cl
        n_in = int(mask.sum())
        if n_in == 0:
            continue
        Xc = Xs_bin[mask]
        for i, m in enumerate(marker_order):
            fid = feat_by_marker.get(m)
            if fid is None:
                continue
            col = f2i[fid]
            mat[i, j] = float(Xc[:, col].sum()) / n_in

    pivot = pd.DataFrame(mat, index=marker_order, columns=clusters)
    sizes = pd.Series({cl: int((leiden_labels == cl).sum()) for cl in clusters})

    top5_keys: set[tuple[int, str]] = set()
    for cl, sub in markers.groupby("leiden_cluster"):
        for r in sub.itertuples():
            top5_keys.add((int(cl), _clean_short(r.short_name, r.pfam_acc)))

    return pivot, sizes, top5_keys


def render() -> Path:
    FIG.mkdir(parents=True, exist_ok=True)
    pivot, sizes, top5_keys = _marker_matrix_full()
    cl_hi = int(max(pivot.columns))

    fig_h = max(7.0, 0.32 * pivot.shape[0] + 2.0)
    fig, ax = plt.subplots(figsize=(11, fig_h))
    sns.heatmap(
        pivot.values,
        ax=ax,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "p_in (fraction of IE in the cluster with this Pfam)"},
        xticklabels=[f"C{int(c)}" for c in pivot.columns],
        yticklabels=pivot.index,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
    ax.set_title(
        f"6960 Unknown · Pfam markers × Leiden (C0–C{cl_hi}, n={int(sizes.sum()):,} IE)"
    )

    for i, m in enumerate(pivot.index):
        for j, cl in enumerate(pivot.columns):
            v = pivot.values[i, j]
            if v >= MIN_PIN_LABEL:
                ax.text(
                    j + 0.5, i + 0.5, f"{v:.2f}",
                    ha="center", va="center",
                    color="white" if v > 0.55 else "black", fontsize=7,
                )
            if (int(cl), m) in top5_keys:
                ax.add_patch(plt.Rectangle(
                    (j, i), 1, 1, fill=False, edgecolor="#1a1a1a", linewidth=1.2, zorder=4,
                ))

    size_txt = ", ".join(f"C{int(c)} n={int(sizes[c]):,}" for c in pivot.columns)
    fig.text(
        0.5, 0.01,
        "Rows are the union of the top-5 Pfam hits per cluster (Fisher / marker_score). "
        "Each cell is p_in in that cluster; a box means the marker is in the top 5 for that C. "
        f"{size_txt}",
        ha="center", va="bottom", fontsize=7.5, color="#444",
    )
    fig.subplots_adjust(bottom=0.14)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight", pad_inches=0.30)
    fig.savefig(OUT_PNG.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.30)
    plt.close(fig)
    print("→", OUT_PNG)
    return OUT_PNG


if __name__ == "__main__":
    render()
