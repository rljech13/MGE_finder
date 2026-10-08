#!/usr/bin/env python3
"""Cluster one class from a classification table the user supplies.

  python cluster_pipeline/cluster_cohort.py \
      --classification classification.tsv \
      --matrix genome_feature_matrix.npz \
      --genomes genome_feature_labels.npz \
      --keep-label unclassified \
      --k 15 \
      --min-cluster-size 5 \
      --seed 42 \
      --out leiden.tsv

``--genomes`` is an npz with an array named ``genomes``. Its order is the row
order of ``--matrix``. ``genome_id`` in the classification table must match
those names. Elements of other labels are not clustered.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, load_npz


def jaccard(X: csr_matrix) -> np.ndarray:
    sim = (X @ X.T).toarray().astype(np.float32)
    row_sums = np.asarray(X.sum(axis=1)).flatten().astype(np.float32)
    denom = row_sums[:, None] + row_sums[None, :] - sim
    mask = denom > 0
    out = np.where(mask, sim / np.where(mask, denom, 1.0), 0.0).astype(np.float32)
    np.fill_diagonal(out, 1.0)
    return out


def leiden_knn(sim: np.ndarray, k: int, min_cluster_size: int, seed: int) -> np.ndarray:
    n = sim.shape[0]
    if n == 0:
        return np.array([], dtype=np.int32)
    if n == 1:
        return np.array([-1], dtype=np.int32)
    k_eff = min(k, n - 1)
    edges: dict[tuple[int, int], float] = {}
    for i in range(n):
        row = sim[i].copy()
        row[i] = -1.0
        for j in np.argpartition(row, -k_eff)[-k_eff:]:
            w = float(row[j])
            if w <= 0.0:
                continue
            a, b = (i, int(j)) if i < j else (int(j), i)
            edges[(a, b)] = max(edges.get((a, b), 0.0), w)
    graph = ig.Graph(n=n, edges=list(edges), directed=False)
    if edges:
        graph.es["weight"] = [edges[e] for e in edges]
        part = leidenalg.find_partition(
            graph, leidenalg.ModularityVertexPartition, weights="weight", seed=seed,
        )
    else:
        part = leidenalg.find_partition(
            graph, leidenalg.ModularityVertexPartition, seed=seed,
        )
    labels = np.asarray(part.membership, dtype=np.int32)
    sizes = pd.Series(labels).value_counts()
    keep = set(sizes[sizes >= min_cluster_size].index.astype(int))
    if keep:
        labels = labels.copy()
        labels[~np.isin(labels, list(keep))] = -1
    else:
        labels = np.full(n, -1, dtype=np.int32)
    return labels


def cluster_label(
    classification: pd.DataFrame,
    matrix: csr_matrix,
    genome_names: np.ndarray,
    *,
    keep_label: str,
    label_column: str,
    id_column: str,
    genome_column: str,
    k: int,
    min_cluster_size: int,
    seed: int,
) -> pd.DataFrame:
    chosen = classification[classification[label_column].astype(str) == keep_label].copy()
    gid_to_i = {str(g): i for i, g in enumerate(genome_names)}
    chosen["matrix_row"] = chosen[genome_column].map(lambda g: gid_to_i.get(str(g)))
    covered = chosen[chosen["matrix_row"].notna()].reset_index(drop=True)
    missing = chosen[chosen["matrix_row"].isna()]
    if covered.empty:
        out = chosen[[id_column]].copy()
        out["leiden_cluster"] = -1
        out["error"] = "genome_not_in_matrix"
        return out
    X = matrix[covered["matrix_row"].astype(int).to_numpy()]
    labels = leiden_knn(jaccard(X), k, min_cluster_size, seed)
    out = covered[[id_column, genome_column, label_column]].copy()
    out["leiden_cluster"] = labels
    out["error"] = ""
    if not missing.empty:
        lost = missing[[id_column, genome_column, label_column]].copy()
        lost["leiden_cluster"] = -1
        lost["error"] = "genome_not_in_matrix"
        out = pd.concat([out, lost], ignore_index=True)
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Jaccard kNN Leiden on one class of a cohort.")
    ap.add_argument("--classification", required=True, type=Path)
    ap.add_argument("--matrix", required=True, type=Path, help="scipy sparse npz, rows are genomes")
    ap.add_argument("--genomes", required=True, type=Path, help="npz with array 'genomes'")
    ap.add_argument("--keep-label", default="unclassified")
    ap.add_argument("--label-column", default="coarse_label")
    ap.add_argument("--id-column", default="integrase_id")
    ap.add_argument("--genome-column", default="genome_id")
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--min-cluster-size", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    classification = pd.read_csv(args.classification, sep="\t", dtype=str)
    for col in (args.id_column, args.genome_column, args.label_column):
        if col not in classification.columns:
            raise SystemExit(f"{args.classification} has no column {col}")
    matrix = load_npz(args.matrix).tocsr()
    genomes = np.load(args.genomes, allow_pickle=True)["genomes"]
    if matrix.shape[0] != len(genomes):
        raise SystemExit(
            f"matrix has {matrix.shape[0]} rows but {args.genomes} lists {len(genomes)} genomes"
        )
    out = cluster_label(
        classification, matrix, genomes,
        keep_label=args.keep_label,
        label_column=args.label_column,
        id_column=args.id_column,
        genome_column=args.genome_column,
        k=args.k,
        min_cluster_size=args.min_cluster_size,
        seed=args.seed,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    n_cl = int(out.loc[out["leiden_cluster"] >= 0, "leiden_cluster"].nunique()) if len(out) else 0
    print(f"wrote {args.out}  rows={len(out)}  clusters={n_cl}")


if __name__ == "__main__":
    main()
