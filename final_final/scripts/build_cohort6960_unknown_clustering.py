#!/usr/bin/env python3
"""Diploma-style graph clustering for cohort-6960 atlas Unknown.

Cohort: Figure-4 Unknown (n=4,774) from
  cohort_6960_atlas_classes_conj_split.tsv (partial → unknown; ICE = T4SS only).

Pipeline (mirrors Stage 8/9 unc rebuild):
  cohort → genome_feature_matrix → Jaccard → kNN-15 → Leiden → HDBSCAN
  → PFAM markers → biology → host purity / diversity → figures

Outputs: final_final/artifacts/nb04_cluster_unknown_6960/
Figures:  final_final/figures/04_cluster_unknown_6960/
"""
from __future__ import annotations

import math
import time
import warnings
from collections import Counter
from pathlib import Path

import hdbscan
import igraph as ig
import leidenalg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from scipy.sparse import load_npz, save_npz
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

warnings.filterwarnings("ignore", category=FutureWarning)

from cohort6960_unknown_paths import (  # noqa: E402
    ART,
    FIG,
    GRAPH_K,
    HDBSCAN_MIN_CLUSTER_SIZE,
    HDBSCAN_MIN_SAMPLES,
    IDS,
    LEIDEN_SEED,
    MIN_LEIDEN_CLUSTER_SIZE,
    PATH_CLASSIF,
    PATH_COHORT_JOIN,
    PATH_COHORT_UNKNOWN,
    PATH_EXCLUDE_ATT_ONLY,
    PATH_GTDB,
    PATH_LAB,
    PATH_LEIDEN_OLD,
    PATH_LEIDEN_S8,
    PATH_LEIDEN_S9,
    PATH_PFAM_DESC,
    PATH_SAMPLES_META,
    PATH_TRUTH,
    PATH_X,
    FF,
    REPO,
    TAG,
    ensure_dirs,
)

for _candidate in ("seaborn-v0_8-whitegrid", "ggplot"):
    if _candidate in plt.style.available:
        plt.style.use(_candidate)
        break
rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def load_cohort() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, dict[str, int]]:
    cohort = pd.read_csv(PATH_COHORT_UNKNOWN, sep="\t")
    unk_ids = set(cohort.loc[cohort["func_class"].astype(str) == "unknown", "integrase_id"].astype(str))
    print(f"Figure-4 Unknown (conj_split): {len(unk_ids):,}")

    # Drop attL–attR-only fragments (integrase outside mge_region) — invalid IE
    exclude: set[str] = set()
    if PATH_EXCLUDE_ATT_ONLY.is_file():
        exclude = {
            ln.strip() for ln in PATH_EXCLUDE_ATT_ONLY.read_text().splitlines() if ln.strip()
        }
        before = len(unk_ids)
        unk_ids -= exclude
        print(
            f"  EXCLUDE att-only short IE (integrase outside window): "
            f"{len(exclude):,} → unknown left {len(unk_ids):,} (was {before:,})"
        )

    print("  primary_label breakdown:")
    print(
        cohort.loc[cohort["integrase_id"].astype(str).isin(unk_ids), "primary_label"]
        .value_counts()
        .to_string()
        .replace("\n", "\n    ")
    )

    join = pd.read_csv(PATH_COHORT_JOIN, sep="\t", low_memory=False)
    join = join[join["integrase_id"].astype(str).isin(unk_ids)].copy()

    clf = pd.read_csv(PATH_CLASSIF, sep="\t", low_memory=False)
    clf["integrase_id"] = clf["element_id"].astype(str).str.split(":").str[0]
    clf_cols = [
        c for c in (
            "integrase_id", "element_id", "genome_id", "sample_id",
            "primary_label", "coarse_label", "novel_tier", "IE_core_score", "rep_gbff",
        ) if c in clf.columns
    ]
    clf = clf[clf_cols].drop_duplicates("integrase_id")

    unc = join.merge(clf, on="integrase_id", how="left", suffixes=("", "_clf"))
    # Figure-4 / join labels are authoritative for this cohort
    if "primary_label_clf" in unc.columns:
        unc = unc.drop(columns=["primary_label_clf"])
    if "coarse_label_clf" in unc.columns:
        unc["coarse_label"] = unc.get("coarse_label").fillna(unc["coarse_label_clf"])
        unc = unc.drop(columns=["coarse_label_clf"])
    if "sample_id_clf" in unc.columns:
        unc["sample_id"] = unc["sample_id"].fillna(unc["sample_id_clf"])
        unc = unc.drop(columns=["sample_id_clf"])

    unc["ie_id"] = unc["ie_id"].fillna(unc.get("element_id", "")).astype(str)
    unc["rep_id"] = unc["genome_id"].astype(str)
    unc.loc[unc["rep_id"].isin(["nan", "None", "<NA>"]), "rep_id"] = pd.NA

    n_with_gid = int(unc["genome_id"].notna().sum())
    print(f"  with genome_id from classification: {n_with_gid:,} / {len(unc):,}")

    lab = np.load(PATH_LAB, allow_pickle=True)
    genomes = np.array([str(g) for g in lab["genomes"]])
    feature_names = np.array([str(f) for f in lab["features"]])
    gid_to_i = {g: i for i, g in enumerate(genomes)}

    unc["matrix_row"] = unc["rep_id"].map(lambda r: gid_to_i.get(str(r)) if pd.notna(r) else None)

    # Remap missing via merged_deduplicated/rep_index (numeric rep_id already in matrix)
    n_missing_before = int(unc["matrix_row"].isna().sum())
    if n_missing_before:
        remap_path = ART / "ie_to_matrix_row_missing285.tsv"
        ri_path = REPO / "results_organized/merged_deduplicated/rep_index.tsv"
        if remap_path.is_file():
            remap = pd.read_csv(remap_path, sep="\t")
            remap = remap[remap["matrix_row"].notna() & (pd.to_numeric(remap.get("nnz"), errors="coerce").fillna(0) > 0)]
            lu_row = dict(zip(remap["integrase_id"].astype(str), remap["matrix_row"].astype(int)))
            lu_rep = dict(zip(remap["integrase_id"].astype(str), remap["matrix_rep_id"].astype(str)))
        elif ri_path.is_file():
            ri = pd.read_csv(ri_path, sep="\t", dtype=str)
            lu_row, lu_rep = {}, {}
            miss = unc.loc[unc["matrix_row"].isna(), ["integrase_id", "ie_id"]]
            for _, r in miss.iterrows():
                ie = str(r["ie_id"])
                hits = ri[(ri["cluster_id"] == ie) | ri["cluster_id"].str.startswith(ie + "|", na=False)]
                if hits.empty:
                    continue
                rid = str(hits.iloc[0]["rep_id"])
                if rid in gid_to_i:
                    lu_row[str(r["integrase_id"])] = gid_to_i[rid]
                    lu_rep[str(r["integrase_id"])] = rid
        else:
            lu_row, lu_rep = {}, {}

        if lu_row:
            miss_mask = unc["matrix_row"].isna()
            filled = unc.loc[miss_mask, "integrase_id"].astype(str).map(lu_row)
            unc.loc[miss_mask, "matrix_row"] = filled
            unc.loc[miss_mask, "rep_id"] = unc.loc[miss_mask, "integrase_id"].astype(str).map(lu_rep).fillna(
                unc.loc[miss_mask, "rep_id"]
            )
            print(f"  remapped via rep_index: {int(filled.notna().sum()):,} / {n_missing_before:,}")

    n_missing = int(unc["matrix_row"].isna().sum())
    print(f"  missing matrix_row: {n_missing:,}")

    truth = pd.read_csv(PATH_TRUTH, sep="\t", low_memory=False)
    truth_cols = ["integrase_id_lookup"]
    for col in ("integrase_len_aa", "ie_len_nt", "integrase_model"):
        if col in truth.columns:
            truth_cols.append(col)
    truth = truth[truth_cols].drop_duplicates("integrase_id_lookup")
    truth = truth.rename(columns={"integrase_id_lookup": "integrase_id"})
    unc = unc.merge(truth, on="integrase_id", how="left", suffixes=("", "_truth"))

    if "integrase_len_aa" in unc.columns:
        unc["integrase_len_aa"] = pd.to_numeric(unc["integrase_len_aa"], errors="coerce")
    if "integrase_len_aa" not in unc.columns or unc["integrase_len_aa"].isna().mean() > 0.5:
        s8_join = FF / "artifacts/stage8_atlas_join.tsv"
        if s8_join.is_file():
            s8 = pd.read_csv(
                s8_join, sep="\t", usecols=["integrase_id", "integrase_len_aa"],
                low_memory=False,
            ).drop_duplicates("integrase_id")
            if "integrase_len_aa" in unc.columns:
                unc = unc.drop(columns=["integrase_len_aa"])
            unc = unc.merge(s8, on="integrase_id", how="left")

    # Locked cohort IE lengths (authoritative for 6960)
    ie_len_path = FF / "artifacts/attL_complement_blast_window/cohort_6960_ie_lengths.tsv"
    if ie_len_path.is_file():
        ie_len = pd.read_csv(ie_len_path, sep="\t", usecols=["integrase_id", "ie_len_nt"])
        ie_len["ie_len_nt"] = pd.to_numeric(ie_len["ie_len_nt"], errors="coerce")
        if "ie_len_nt" in unc.columns:
            unc = unc.drop(columns=["ie_len_nt"])
        unc = unc.merge(ie_len, on="integrase_id", how="left")

    if PATH_GTDB.is_file() and "sample_id" in unc.columns:
        gtdb = pd.read_csv(PATH_GTDB, sep="\t", low_memory=False)
        tax_cols = [c for c in ("sample_id", "domain", "phylum", "class", "order", "family") if c in gtdb.columns]
        gtdb = gtdb[tax_cols].drop_duplicates("sample_id")
        unc = unc.merge(gtdb, on="sample_id", how="left")

    if PATH_SAMPLES_META.is_file() and "sample_id" in unc.columns:
        meta = pd.read_csv(PATH_SAMPLES_META, sep="\t", low_memory=False)
        if "Genus" in meta.columns:
            genus_lu = meta[["sample_id", "Genus"]].drop_duplicates("sample_id")
            unc = unc.merge(genus_lu, on="sample_id", how="left")

    covered = unc[unc["matrix_row"].notna()].copy().reset_index(drop=True)
    covered["matrix_row"] = covered["matrix_row"].astype(int)
    covered = covered.sort_values("integrase_id").reset_index(drop=True)
    print(f"Feature coverage: {len(covered):,}")

    no_feat = unc[unc["matrix_row"].isna()].copy()
    return covered, no_feat, genomes, feature_names, gid_to_i


def build_features_and_jaccard(covered: pd.DataFrame) -> tuple:
    t0 = time.time()
    X = load_npz(PATH_X).astype(np.float32)
    print(f"X (full): {X.shape}, nnz={X.nnz:,}  ({time.time() - t0:.1f}s)")

    row_idx = covered["matrix_row"].to_numpy()
    Xs = X[row_idx]
    print(f"Xs ({TAG}): {Xs.shape}, nnz={Xs.nnz:,}")
    save_npz(ART / f"feature_subset_{TAG}.npz", Xs)

    out_cols = [
        "integrase_id", "rep_id", "ie_id", "sample_id", "matrix_row",
        "integrase_len_aa", "ie_len_nt", "phylum", "class", "order", "family", "Genus",
        "novel_tier", "IE_core_score", "primary_label", "coarse_label",
    ]
    out_cols = [c for c in out_cols if c in covered.columns]
    covered[out_cols].to_csv(ART / f"ie_rep_rows_{TAG}.tsv", sep="\t", index=False)
    print(f"→ {ART / f'ie_rep_rows_{TAG}.tsv'}")
    del X

    t0 = time.time()
    sim = (Xs @ Xs.T).toarray().astype(np.float32)
    print(f"  intersections: {time.time() - t0:.1f}s")
    row_sums = np.asarray(Xs.sum(axis=1)).flatten().astype(np.float32)
    denom = row_sums[:, None] + row_sums[None, :] - sim
    mask = denom > 0
    sim_jac = np.where(mask, sim / np.where(mask, denom, 1.0), 0.0).astype(np.float32)
    np.fill_diagonal(sim_jac, 1.0)
    del sim, denom, mask
    off_diag = (sim_jac.sum() - sim_jac.shape[0]) / (sim_jac.shape[0] ** 2 - sim_jac.shape[0])
    print(f"Jaccard: {sim_jac.shape}, mean(off-diag)={off_diag:.4f}")

    np.savez_compressed(
        ART / f"similarity_jaccard_{TAG}.npz",
        sim=sim_jac,
        genomes=covered["rep_id"].astype(str).to_numpy(),
    )
    return Xs, sim_jac


def build_knn_leiden(covered: pd.DataFrame, sim_jac: np.ndarray):
    n = sim_jac.shape[0]
    k_eff = min(GRAPH_K, max(1, n - 1))
    print(f"kNN: k_eff={k_eff}")

    t0 = time.time()
    edge_dict: dict[tuple[int, int], float] = {}
    for i in range(n):
        row = sim_jac[i].copy()
        row[i] = -1.0
        top_k = np.argpartition(row, -k_eff)[-k_eff:]
        for j in top_k:
            w = float(row[j])
            if w <= 0.0:
                continue
            a, b = (i, int(j)) if i < j else (int(j), i)
            prev = edge_dict.get((a, b), 0.0)
            if w > prev:
                edge_dict[(a, b)] = w
    print(f"  edges: {len(edge_dict):,} ({time.time() - t0:.1f}s)")

    edges = list(edge_dict.keys())
    weights = [edge_dict[e] for e in edges]
    G = ig.Graph(n=n, edges=edges, directed=False)
    G.es["weight"] = weights
    G.vs["name"] = [f"n{i}" for i in range(n)]
    comp = G.clusters()
    print(f"  components: {len(comp)}, max={max(len(c) for c in comp)}")

    part = leidenalg.find_partition(
        G, leidenalg.ModularityVertexPartition,
        weights="weight", seed=LEIDEN_SEED,
    )
    leiden_labels = np.asarray(part.membership, dtype=np.int32)
    sz_raw = pd.Series(leiden_labels).value_counts().sort_values(ascending=False)
    keep = set(sz_raw[sz_raw >= MIN_LEIDEN_CLUSTER_SIZE].index.astype(int))
    n_pruned = int((~np.isin(leiden_labels, list(keep))).sum()) if keep else int(len(leiden_labels))
    if n_pruned:
        leiden_labels = leiden_labels.copy()
        leiden_labels[~np.isin(leiden_labels, list(keep))] = -1
        print(
            f"Leiden prune: clusters n<{MIN_LEIDEN_CLUSTER_SIZE} → noise (-1); "
            f"dropped {n_pruned} IE from {int((sz_raw < MIN_LEIDEN_CLUSTER_SIZE).sum())} tiny communities"
        )
    sz = pd.Series(leiden_labels[leiden_labels >= 0]).value_counts().sort_values(ascending=False)
    n_singletons = int((sz == 1).sum()) if len(sz) else 0
    n_noise = int((leiden_labels < 0).sum())
    print(
        f"Leiden: clusters={len(sz)}, singleton={n_singletons}, noise={n_noise}, "
        f"max={int(sz.max()) if len(sz) else 0}, modularity={part.quality():.4f}"
    )

    leiden_df = pd.DataFrame({
        "integrase_id": covered["integrase_id"].astype(str).values,
        "rep_id": covered["rep_id"].astype(str).values,
        "sample_id": covered["sample_id"].astype(str).values,
        "leiden_cluster": leiden_labels,
    })
    leiden_df.to_csv(ART / f"leiden_clusters_{TAG}.tsv", sep="\t", index=False)
    # Keep IE-level graph attributes in sync with pruned labels
    for i, cid in enumerate(leiden_labels):
        G.vs[i]["leiden_cluster"] = int(cid)
    G.write_graphml(str(ART / f"knn_graph_{TAG}.graphml"))
    print(f"→ {ART / f'leiden_clusters_{TAG}.tsv'}")
    print(f"→ {ART / f'knn_graph_{TAG}.graphml'}")
    return G, leiden_labels, sz, part


def run_hdbscan(sim_jac: np.ndarray, covered: pd.DataFrame) -> np.ndarray:
    dist = (1.0 - sim_jac).astype(np.float64)
    np.fill_diagonal(dist, 0.0)
    t0 = time.time()
    h = hdbscan.HDBSCAN(
        metric="precomputed",
        min_cluster_size=HDBSCAN_MIN_CLUSTER_SIZE,
        min_samples=HDBSCAN_MIN_SAMPLES,
    )
    labels = h.fit_predict(dist)
    n_cl = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int((labels == -1).sum())
    print(f"HDBSCAN: clusters={n_cl}, noise={n_noise}, {time.time() - t0:.1f}s")

    pd.DataFrame({
        "integrase_id": covered["integrase_id"].astype(str).values,
        "rep_id": covered["rep_id"].astype(str).values,
        "hdbscan_cluster": labels.astype(int),
    }).to_csv(ART / f"hdbscan_clusters_{TAG}.tsv", sep="\t", index=False)
    return labels


def compute_pfam_markers(Xs, feature_names: np.ndarray, leiden_labels: np.ndarray, sz: pd.Series) -> pd.DataFrame:
    pfam_desc = pd.read_csv(PATH_PFAM_DESC, sep="\t")
    pfam_desc["pfam_acc"] = pfam_desc["pfam_acc"].astype(str)
    pfam_short = dict(zip(pfam_desc["pfam_acc"], pfam_desc["short_name"].astype(str)))
    pfam_long = dict(zip(pfam_desc["pfam_acc"], pfam_desc["full_name"].astype(str)))

    feat2pfam = {}
    for fn in feature_names:
        if fn.startswith("pfam_"):
            feat2pfam[fn] = fn[5:].split(".")[0]

    Xs_bin = (Xs > 0).astype(np.float32)
    n_total = Xs_bin.shape[0]
    feat_freq_global = np.asarray(Xs_bin.sum(axis=0)).flatten()
    eps = 1e-3
    records = []
    for cl in sz.head(20).index:
        mask = leiden_labels == cl
        n_in = int(mask.sum())
        if n_in < 5:
            continue
        in_freq = np.asarray(Xs_bin[mask].sum(axis=0)).flatten()
        out_freq = feat_freq_global - in_freq
        n_out = n_total - n_in
        p_in = in_freq / n_in
        p_out = out_freq / max(n_out, 1)
        enrichment = (p_in + eps) / (p_out + eps)
        score = np.log10(enrichment) * p_in
        for ix in np.argsort(-score)[:5]:
            f = feature_names[ix]
            if not f.startswith("pfam_"):
                continue
            acc = feat2pfam.get(f, "")
            records.append({
                "leiden_cluster": int(cl),
                "cluster_size": n_in,
                "feature_id": f,
                "pfam_acc": acc,
                "short_name": pfam_short.get(acc, ""),
                "full_name": pfam_long.get(acc, ""),
                "count_in": int(in_freq[ix]),
                "p_in": round(float(p_in[ix]), 3),
                "p_out": round(float(p_out[ix]), 3),
                "enrichment": round(float(enrichment[ix]), 2),
                "marker_score": round(float(score[ix]), 3),
            })
    marker_df = pd.DataFrame(records)
    marker_df.to_csv(ART / f"cluster_pfam_markers_{TAG}.tsv", sep="\t", index=False)
    print(f"→ {ART / f'cluster_pfam_markers_{TAG}.tsv'}  ({len(marker_df):,} rows)")
    return marker_df


def compute_biology(covered: pd.DataFrame, leiden_labels: np.ndarray, sz: pd.Series) -> pd.DataFrame:
    merged = covered.copy()
    merged["leiden_cluster"] = leiden_labels
    rows = []
    for cl in sz.head(20).index:
        sub = merged[merged["leiden_cluster"] == int(cl)]
        if len(sub) < 5:
            continue
        top_phy = sub["phylum"].value_counts().head(1) if "phylum" in sub.columns else pd.Series(dtype=int)
        ie_len = pd.to_numeric(sub.get("ie_len_nt"), errors="coerce")
        int_aa = pd.to_numeric(sub.get("integrase_len_aa"), errors="coerce")
        ie_core = pd.to_numeric(sub.get("IE_core_score"), errors="coerce")
        rows.append({
            "leiden_cluster": int(cl),
            "n_IE": len(sub),
            "ie_len_med": int(ie_len.median()) if ie_len.notna().any() else None,
            "ie_len_q25": int(ie_len.quantile(0.25)) if ie_len.notna().any() else None,
            "ie_len_q75": int(ie_len.quantile(0.75)) if ie_len.notna().any() else None,
            "integrase_aa_med": int(int_aa.median()) if int_aa.notna().any() else None,
            "top_phylum": top_phy.index[0] if len(top_phy) else None,
            "top_phylum_frac": round(top_phy.iloc[0] / len(sub), 2) if len(top_phy) else None,
            "n_phyla": sub["phylum"].nunique() if "phylum" in sub.columns else None,
            "frac_high_interest_novel": round(
                (sub.get("novel_tier", pd.Series(dtype=str)) == "high_interest").mean(), 2
            ) if "novel_tier" in sub.columns else None,
            "IE_core_med": round(float(ie_core.median()), 3) if ie_core.notna().any() else None,
        })
    biol_df = pd.DataFrame(rows)
    biol_df.to_csv(ART / f"cluster_biology_{TAG}.tsv", sep="\t", index=False)
    print(f"→ {ART / f'cluster_biology_{TAG}.tsv'}")
    return biol_df, merged


def compute_host_purity_diversity(merged: pd.DataFrame, leiden_labels: np.ndarray, sz: pd.Series) -> None:
    purity_rows = []
    for cl in sz.index:
        mask = merged["leiden_cluster"] == int(cl)
        if mask.sum() < 5:
            continue
        sub = merged.loc[mask, "phylum"].dropna() if "phylum" in merged.columns else pd.Series(dtype=str)
        cnt = Counter(sub.tolist())
        n = sum(cnt.values())
        top_n = max(cnt.values()) if cnt else 0
        purity = top_n / n if n else 0
        h = -sum((c / n) * math.log2(c / n) for c in cnt.values()) if n else 0
        purity_rows.append({
            "cluster": int(cl),
            "n_IE": int(mask.sum()),
            "n_with_phylum": n,
            "purity_phylum": round(purity, 3),
            "shannon_phylum_bits": round(h, 3),
            "effective_phyla": round(2 ** h, 2) if n else 0,
        })
    pd.DataFrame(purity_rows).to_csv(ART / f"cluster_host_purity_{TAG}.tsv", sep="\t", index=False)

    def shannon(values) -> float:
        s = pd.Series(values).value_counts(normalize=True)
        s = s[s > 0]
        return float(-(s * np.log(s)).sum())

    div_rows = []
    for cl, g in merged.groupby("leiden_cluster"):
        if pd.isna(cl) or int(cl) < 0:
            continue
        p_phy = g["phylum"].fillna("unknown") if "phylum" in g.columns else pd.Series(["unknown"] * len(g))
        p_gen = g["Genus"].fillna("unknown") if "Genus" in g.columns else pd.Series(["unknown"] * len(g))
        div_rows.append({
            "leiden_cluster": int(cl),
            "n_IE": len(g),
            "shannon_phylum": shannon(p_phy),
            "shannon_genus": shannon(p_gen),
            "eff_phyla": float(np.exp(shannon(p_phy))),
            "eff_genera": float(np.exp(shannon(p_gen))),
            "top_phylum_frac": float(p_phy.value_counts(normalize=True).iloc[0]),
            "top_genus_frac": float(p_gen.value_counts(normalize=True).iloc[0]) if len(p_gen) else 0.0,
        })
    pd.DataFrame(div_rows).sort_values("n_IE", ascending=False).to_csv(
        ART / f"cluster_diversity_{TAG}.tsv", sep="\t", index=False
    )
    print(f"→ cluster_host_purity / cluster_diversity")


def compare_baselines(leiden_df: pd.DataFrame, covered: pd.DataFrame, sz: pd.Series, part) -> tuple[float, float, pd.DataFrame]:
    old = pd.read_csv(PATH_LEIDEN_OLD, sep="\t")
    old["rep_id"] = old["genome_id"].astype(str)
    merged_old = leiden_df.merge(old[["rep_id", "leiden_hd_cluster"]], on="rep_id", how="left")
    have = merged_old.dropna(subset=["leiden_hd_cluster"]).copy()
    have["leiden_hd_cluster"] = have["leiden_hd_cluster"].astype(int)
    ari_old = adjusted_rand_score(have["leiden_hd_cluster"], have["leiden_cluster"])
    nmi_old = normalized_mutual_info_score(have["leiden_hd_cluster"], have["leiden_cluster"])
    print(f"vs 7,468-IE Leiden (shared {len(have):,}): ARI={ari_old:.4f}, NMI={nmi_old:.4f}")

    ct = pd.crosstab(have["leiden_cluster"], have["leiden_hd_cluster"])

    if PATH_LEIDEN_S8.is_file():
        s8 = pd.read_csv(PATH_LEIDEN_S8, sep="\t")
        m8 = leiden_df.merge(s8[["integrase_id", "leiden_cluster"]].rename(
            columns={"leiden_cluster": "leiden_s8"}
        ), on="integrase_id", how="inner")
        if len(m8):
            ari_s8 = adjusted_rand_score(m8["leiden_s8"], m8["leiden_cluster"])
            nmi_s8 = normalized_mutual_info_score(m8["leiden_s8"], m8["leiden_cluster"])
            print(f"vs Stage8 rebuild (shared {len(m8):,}): ARI={ari_s8:.4f}, NMI={nmi_s8:.4f}")

    if PATH_LEIDEN_S9.is_file():
        s9 = pd.read_csv(PATH_LEIDEN_S9, sep="\t")
        m9 = leiden_df.merge(s9[["integrase_id", "leiden_cluster"]].rename(
            columns={"leiden_cluster": "leiden_s9"}
        ), on="integrase_id", how="inner")
        if len(m9):
            ari_s9 = adjusted_rand_score(m9["leiden_s9"], m9["leiden_cluster"])
            nmi_s9 = normalized_mutual_info_score(m9["leiden_s9"], m9["leiden_cluster"])
            print(f"vs Stage9 unc rebuild (shared {len(m9):,}): ARI={ari_s9:.4f}, NMI={nmi_s9:.4f}")

    return ari_old, nmi_old, ct


def render_figures(
    sz: pd.Series, ct: pd.DataFrame, ari: float, nmi: float,
    n_singletons: int, part, hdbscan_labels: np.ndarray,
    sim_jac: np.ndarray, leiden_labels: np.ndarray,
    integrase_ids: np.ndarray,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    sz_top = sz.head(15)
    ax.barh([f"C{int(c)}" for c in sz_top.index][::-1], sz_top.values[::-1], color="#4f7fc8")
    for y, v in enumerate(sz_top.values[::-1]):
        ax.text(v + max(sz_top.values) * 0.01, y, f" {v:,}", va="center", fontsize=9)
    ax.set_xlabel("IE")
    ax.set_title(f"Leiden ({TAG}) — top 15; {len(sz)} clusters, singleton={n_singletons}")
    ax.set_xlim(0, max(sz_top.values) * 1.18)

    ax = axes[1]
    if ct.shape[0] and ct.shape[1]:
        ct_top = ct.loc[sz.head(15).index.intersection(ct.index)]
        old_order = ct_top.sum(axis=0).sort_values(ascending=False).head(15).index
        ct_top = ct_top[old_order]
        ax.imshow(
            ct_top.values / ct_top.values.sum(axis=1, keepdims=True).clip(1),
            aspect="auto", cmap="Blues",
        )
        ax.set_xticks(range(len(old_order)))
        ax.set_xticklabels([f"C{int(c)}" for c in old_order], fontsize=8)
        ax.set_yticks(range(len(ct_top)))
        ax.set_yticklabels([f"new C{int(c)}" for c in ct_top.index], fontsize=8)
        ax.set_title(f"vs old 7,468-IE Leiden  ARI={ari:.3f}, NMI={nmi:.3f}")
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(
            FIG / f"01_clusters_and_contingency.{ext}",
            dpi=140, bbox_inches="tight", pad_inches=0.30,
        )
    plt.close(fig)
    print(f"→ {FIG / '01_clusters_and_contingency.png'}")

    try:
        import umap
        t0 = time.time()
        reducer = umap.UMAP(
            metric="precomputed", n_neighbors=15, min_dist=0.05,
            random_state=42, n_components=2, verbose=False,
        )
        emb = reducer.fit_transform(1.0 - sim_jac)
        print(f"UMAP: {time.time() - t0:.1f}s")
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 6))
        cmap = plt.get_cmap("tab20")
        top_clusters = sz.head(15).index.tolist()
        ax = axes[0]
        for k, cl in enumerate(top_clusters):
            m = leiden_labels == cl
            ax.scatter(emb[m, 0], emb[m, 1], s=5, alpha=0.6, color=cmap(k % 20), label=f"C{int(cl)} ({m.sum()})")
        other = ~np.isin(leiden_labels, top_clusters)
        if other.any():
            ax.scatter(emb[other, 0], emb[other, 1], s=3, alpha=0.3, color="lightgray")
        ax.set_title("UMAP · Leiden (6960 Unknown)")
        ax.legend(fontsize=6, ncol=2, loc="upper right", framealpha=0.7)
        ax = axes[1]
        h_top = pd.Series(hdbscan_labels[hdbscan_labels >= 0]).value_counts().head(15).index.tolist()
        for k, cl in enumerate(h_top):
            m = hdbscan_labels == cl
            ax.scatter(emb[m, 0], emb[m, 1], s=5, alpha=0.6, color=cmap(k % 20), label=f"C{int(cl)} ({m.sum()})")
        other_h = (hdbscan_labels == -1) | (~np.isin(hdbscan_labels, h_top))
        ax.scatter(emb[other_h, 0], emb[other_h, 1], s=3, alpha=0.25, color="lightgray")
        ax.set_title("UMAP · HDBSCAN")
        ax.legend(fontsize=6, ncol=2, loc="upper right", framealpha=0.7)
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
        plt.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(
                FIG / f"02_umap_clusters.{ext}",
                dpi=140, bbox_inches="tight", pad_inches=0.30,
            )
        plt.close(fig)
        np.savez_compressed(
            ART / f"umap_embedding_{TAG}.npz",
            embedding=emb,
            integrase_id=integrase_ids,
        )
    except Exception as e:
        print(f"UMAP skipped: {e}")


def save_per_ie_and_ids(merged: pd.DataFrame, biol_df: pd.DataFrame, no_feat: pd.DataFrame, sz: pd.Series) -> None:
    out = merged.merge(
        biol_df[["leiden_cluster", "top_phylum"]].rename(columns={"top_phylum": "cluster_top_phylum"}),
        on="leiden_cluster", how="left",
    )
    keep = [
        "integrase_id", "rep_id", "sample_id", "ie_id",
        "integrase_len_aa", "ie_len_nt", "phylum", "class", "order", "family", "Genus",
        "novel_tier", "IE_core_score", "primary_label", "coarse_label",
        "leiden_cluster", "cluster_top_phylum",
    ]
    keep = [c for c in keep if c in out.columns]
    out = out[keep]
    out_path = FF / "artifacts" / f"{TAG}_clusters_rebuilt.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"→ {out_path}  ({len(out):,} rows)")

    for cl in sz.head(20).index:
        if int(sz.loc[cl]) < 5:
            continue
        ids = sorted(out.loc[out["leiden_cluster"] == int(cl), "integrase_id"].astype(str).unique())
        p = IDS / f"{TAG}_rebuilt_C{int(cl):02d}.txt"
        p.write_text("\n".join(ids) + "\n")

    no_ids = sorted(no_feat["integrase_id"].astype(str).unique()) if len(no_feat) else []
    (IDS / f"{TAG}_rebuilt_no_features.txt").write_text("\n".join(no_ids) + "\n")
    print(f"  id_lists: {len(sz)} cluster files, {len(no_ids)} no-features")


def write_report(
    n_atlas: int, n_covered: int, n_no_feat: int, sz: pd.Series,
    n_singletons: int, part, n_h_clusters: int, n_h_noise: int,
    ari: float, nmi: float, biol_df: pd.DataFrame,
) -> None:
    lines = [
        f"# Cohort-6960 Unknown · Leiden clustering (diploma-style)\n",
        f"**Cohort**: Figure-4 Unknown (`cohort_6960_atlas_classes_conj_split.tsv`) = {n_atlas:,} IE",
        f" (with features: {n_covered:,}; no matrix row: {n_no_feat:,})\n",
        f"**Algorithm**: Jaccard → kNN-{GRAPH_K} → Leiden (seed={LEIDEN_SEED})\n",
        f"* clusters: **{len(sz)}** (singleton {n_singletons})",
        f"* modularity: **{part.quality():.4f}**",
        f"* HDBSCAN: clusters {n_h_clusters}, noise {n_h_noise}",
        f"* ARI vs 7,468-IE Leiden: **{ari:.4f}**",
        f"* NMI vs 7,468-IE Leiden: **{nmi:.4f}**\n",
        "## Top clusters\n",
        "| C | n_IE | top phylum | %phylum | ie_len_med |",
        "| --- | ---: | --- | ---: | ---: |",
    ]
    for r in biol_df.head(15).itertuples():
        frac = (r.top_phylum_frac or 0) * 100
        med = r.ie_len_med if pd.notna(r.ie_len_med) else "NA"
        lines.append(f"| C{int(r.leiden_cluster)} | {r.n_IE:,} | {r.top_phylum} | {frac:.0f}% | {med} |")
    lines.append(f"\nArtifacts: `artifacts/nb04_cluster_unknown_6960/`\n")
    lines.append(f"Figures: `figures/04_cluster_unknown_6960/`\n")
    (FIG / "REPORT.md").write_text("\n".join(lines))
    print(f"→ {FIG / 'REPORT.md'}")


def main() -> None:
    ensure_dirs()
    print(f"=== Cohort-6960 Unknown graph clustering ===")
    print(f"ART: {ART}")
    print(f"FIG: {FIG}")

    covered, no_feat, _genomes, feature_names, _gid = load_cohort()
    Xs, sim_jac = build_features_and_jaccard(covered)
    _G, leiden_labels, sz, part = build_knn_leiden(covered, sim_jac)
    n_singletons = int((sz == 1).sum())
    hdbscan_labels = run_hdbscan(sim_jac, covered)
    n_h_clusters = len(set(hdbscan_labels)) - (1 if -1 in hdbscan_labels else 0)
    n_h_noise = int((hdbscan_labels == -1).sum())

    marker_df = compute_pfam_markers(Xs, feature_names, leiden_labels, sz)
    biol_df, merged = compute_biology(covered, leiden_labels, sz)
    compute_host_purity_diversity(merged, leiden_labels, sz)

    leiden_df = pd.DataFrame({
        "integrase_id": covered["integrase_id"].astype(str).values,
        "rep_id": covered["rep_id"].astype(str).values,
        "sample_id": covered["sample_id"].astype(str).values,
        "leiden_cluster": leiden_labels,
    })
    ari, nmi, ct = compare_baselines(leiden_df, covered, sz, part)
    render_figures(
        sz, ct, ari, nmi, n_singletons, part, hdbscan_labels, sim_jac, leiden_labels,
        covered["integrase_id"].astype(str).to_numpy(),
    )
    save_per_ie_and_ids(merged, biol_df, no_feat, sz)
    write_report(
        n_atlas=len(covered) + len(no_feat),
        n_covered=len(covered),
        n_no_feat=len(no_feat),
        sz=sz,
        n_singletons=n_singletons,
        part=part,
        n_h_clusters=n_h_clusters,
        n_h_noise=n_h_noise,
        ari=ari,
        nmi=nmi,
        biol_df=biol_df,
    )
    print("\n✓ clustering complete")


if __name__ == "__main__":
    main()
