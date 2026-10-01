#!/usr/bin/env python3
"""Recompute hallmark composition with HMM imputation for forming integrase.

Problem: `integrase_tyrosine` from Pfam feature matrix alone misses ~1% of IE
even though every cohort member passed Stage-1 HMM (PF00589 / PF22022).

Fix: for each IE, if forming integrase HMM model is PF00589* or PF22022*,
force `integrase_tyrosine` present before cluster aggregation.

Outputs:
  artifacts/nb04_cluster_unknown_6960/cluster_hallmark_composition_c6960_unk.tsv
    (updated; includes frac_ie_with_role_pfam + frac_ie_with_role imputed)
  figures/04_cluster_unknown_6960/03_heatmap_hallmark_fraction.{png,svg}
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
    / "integrase_truth_table_stage9.tsv"
)

ROLES_ORDER = [
    "integrase_tyrosine", "integrase_serine", "transposase", "replication",
    "helicase", "aaa_atpase", "methyltransferase", "restriction_T1",
    "defense_STAND", "toxin_antitoxin", "crispr_cas",
    "regulator_HTH", "tpr_repeat", "zn_ribbon", "phage_virion",
    "smc_chromosome", "peptidase", "iron_sulfur", "transport_abc",
    "coiled_coil", "other",
]

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def _is_tyrosine_hmm(model: str) -> bool:
    m = str(model).strip()
    return m.startswith("PF00589") or m.startswith("PF22022")


def _load_feat_roles(feature_names: np.ndarray):
    sys.path.insert(0, str(FF / "03_clustering" / "lib"))
    from cluster_biology import assign_role

    pfam_desc = pd.read_csv(
        FF / "03_clustering" / "reference" / "pfam_descriptions.tsv",
        sep="\t",
    ).fillna("")
    pfam_cat = pd.read_csv(
        FF / "03_clustering" / "reference" / "pfam_catalog.tsv",
        sep="\t",
    ).fillna("")
    desc_idx = pfam_desc.set_index("pfam_acc")
    cat_idx = pfam_cat.set_index("pfam_acc")

    def get_clan(acc: str) -> str:
        if acc in desc_idx.index and "clan" in desc_idx.columns:
            return str(desc_idx.loc[acc, "clan"])
        return ""

    def get_full(acc: str) -> str:
        if acc in desc_idx.index and "full_name" in desc_idx.columns:
            return str(desc_idx.loc[acc, "full_name"])
        if acc in cat_idx.index and "full_name" in cat_idx.columns:
            return str(cat_idx.loc[acc, "full_name"])
        return ""

    feat_role = np.full(len(feature_names), "other", dtype=object)
    for i, fn in enumerate(feature_names):
        if not fn.startswith("pfam_"):
            continue
        acc = fn[5:].split(".")[0]
        feat_role[i] = assign_role(acc, get_clan(acc), get_full(acc))
    return feat_role


def recompute_hallmark() -> pd.DataFrame:
    lab = np.load(PROJ / "genome_feature_labels.npz", allow_pickle=True)
    feature_names = np.array([str(f) for f in lab["features"]])
    Xs = load_npz(ART / "feature_subset_c6960_unk.npz")
    rows = pd.read_csv(ART / "ie_rep_rows_c6960_unk.tsv", sep="\t")
    leid = pd.read_csv(ART / "leiden_clusters_c6960_unk.tsv", sep="\t")

    leiden_labels = leid["leiden_cluster"].to_numpy()
    sz = pd.Series(leiden_labels).value_counts().sort_index()

    feat_role = _load_feat_roles(feature_names)
    roles_unique = sorted(set(feat_role))

    Xs_bin = (Xs > 0).astype(np.float32).tocsr()
    ie_role_count: dict[str, np.ndarray] = {}
    ie_role_present_pfam: dict[str, np.ndarray] = {}
    for r in roles_unique:
        cols = np.where(feat_role == r)[0]
        if len(cols) == 0:
            ie_role_count[r] = np.zeros(Xs_bin.shape[0], dtype=np.float32)
        else:
            ie_role_count[r] = np.asarray(Xs_bin[:, cols].sum(axis=1)).flatten()
        ie_role_present_pfam[r] = (ie_role_count[r] > 0).astype(np.int8)

    # HMM imputation for forming integrase → integrase_tyrosine
    truth = pd.read_csv(PATH_TRUTH, sep="\t", usecols=["integrase_id_lookup", "integrase_model"])
    truth = truth.drop_duplicates("integrase_id_lookup").rename(
        columns={"integrase_id_lookup": "integrase_id"}
    )
    truth["integrase_id"] = truth["integrase_id"].astype(str)
    truth["hmm_tyrosine"] = truth["integrase_model"].map(_is_tyrosine_hmm)

    rep = rows.copy()
    rep["integrase_id"] = rep["integrase_id"].astype(str)
    rep = rep.merge(truth[["integrase_id", "integrase_model", "hmm_tyrosine"]], on="integrase_id", how="left")
    n_miss_hmm = int(rep["hmm_tyrosine"].isna().sum())
    if n_miss_hmm:
        # Remapped / newly covered IE may lack Stage-1 truth rows — keep Pfam-only signal.
        print(f"WARN: HMM truth missing for {n_miss_hmm:,} IE — skip imputation for those")
        rep["hmm_tyrosine"] = rep["hmm_tyrosine"].fillna(0).astype(np.int8)

    hmm_tyrosine = rep["hmm_tyrosine"].to_numpy(dtype=np.int8)
    n_imputed = int((hmm_tyrosine & ~ie_role_present_pfam["integrase_tyrosine"]).sum())
    print(f"Pfam-only integrase_tyrosine: {ie_role_present_pfam['integrase_tyrosine'].sum():,}/{len(hmm_tyrosine):,}")
    print(f"HMM imputation adds: {n_imputed:,} IE")

    ie_role_present = {r: arr.copy() for r, arr in ie_role_present_pfam.items()}
    ie_role_present["integrase_tyrosine"] = np.maximum(
        ie_role_present_pfam["integrase_tyrosine"], hmm_tyrosine
    ).astype(np.int8)

    roles_order = [r for r in ROLES_ORDER if r in roles_unique]
    records = []
    for cl in sz.index:
        mask = leiden_labels == cl
        n_ie = int(mask.sum())
        if n_ie < 5:
            continue
        for r in roles_order:
            present_pfam = int(ie_role_present_pfam[r][mask].sum())
            present_imp = int(ie_role_present[r][mask].sum())
            hits_in = float(ie_role_count[r][mask].sum())
            n_imp = present_imp - present_pfam if r == "integrase_tyrosine" else 0
            records.append({
                "cluster": int(cl),
                "cluster_size": n_ie,
                "role": r,
                "n_ie_with_role_pfam": present_pfam,
                "frac_ie_with_role_pfam": round(present_pfam / max(n_ie, 1), 3),
                "n_ie_with_role": present_imp,
                "frac_ie_with_role": round(present_imp / max(n_ie, 1), 3),
                "n_imputed_tyrosine": n_imp,
                "mean_pfam_hits_per_ie": round(hits_in / max(n_ie, 1), 2),
                "total_pfam_hits": int(hits_in),
            })

    df = pd.DataFrame(records)
    out_tsv = ART / "cluster_hallmark_composition_c6960_unk.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"→ {out_tsv}  ({len(df):,} rows)")
    return df, sz


def render_heatmap(hallmark_df: pd.DataFrame, sz: pd.Series) -> None:
    wide = hallmark_df.pivot(index="cluster", columns="role", values="frac_ie_with_role").fillna(0)
    show_roles = [
        r for r in ROLES_ORDER
        if r != "other" and r in wide.columns and wide[r].sum() > 0
    ]
    wide = wide[show_roles]
    wide = wide.loc[sz.head(20).index.intersection(wide.index)]
    ylabels = [f"C{int(c)} (n={int(sz.loc[c])})" for c in wide.index]

    fig_h = max(5.5, 0.38 * len(wide) + 2.2)
    fig_w = max(11, 0.55 * len(show_roles) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sns.heatmap(
        wide.values,
        ax=ax,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "frac_ie_with_role (HMM-imputed for integrase_tyrosine)"},
        xticklabels=show_roles,
        yticklabels=ylabels,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right", fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)

    for i in range(wide.shape[0]):
        for j in range(wide.shape[1]):
            v = wide.values[i, j]
            if v >= 0.05:
                ax.text(
                    j + 0.5, i + 0.5, f"{v:.2f}",
                    ha="center", va="center", fontsize=7,
                    color="white" if v > 0.55 else "black",
                )

    ax.set_title(
        "6960 Unknown · hallmark roles × cluster (Leiden)",
        fontsize=11, pad=10,
    )
    fig.text(
        0.5, 0.01,
        "integrase_tyrosine: Pfam in the feature matrix OR the forming-integrase HMM (PF00589/PF22022, Stage 1). "
        "Other roles are Pfam in the matrix only.",
        ha="center", va="bottom", fontsize=8, color="#444", wrap=True,
    )
    fig.subplots_adjust(bottom=0.12)
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    out_png = FIG / "03_heatmap_hallmark_fraction.png"
    fig.savefig(out_png, dpi=160, bbox_inches="tight", pad_inches=0.30)
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.30)
    plt.close(fig)
    print(f"→ {out_png}")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    df, sz = recompute_hallmark()
    tyr = df[df.role == "integrase_tyrosine"]
    print("\nintegrase_tyrosine after imputation:")
    print(tyr[["cluster", "cluster_size", "n_ie_with_role_pfam", "n_ie_with_role", "frac_ie_with_role"]].to_string(index=False))
    render_heatmap(df, sz)


if __name__ == "__main__":
    main()
