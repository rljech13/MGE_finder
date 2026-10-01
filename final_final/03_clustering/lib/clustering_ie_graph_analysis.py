"""
Cluster-level view of an IE graphml: a mixing Poisson null and label permutations.

Called after ``filtered_ie_knn.graphml`` has been written.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats


def _scalar(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ""
    if isinstance(val, (np.integer, int)):
        return int(val)
    if isinstance(val, (np.floating, float)):
        return float(val)
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    s = str(val)
    return s[:6000] + ("…" if len(s) > 6000 else "")


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    ranked[order] = np.arange(1, n + 1)
    adj = p * n / ranked
    adj_sorted = adj[order]
    for i in range(n - 2, -1, -1):
        adj_sorted[i] = min(adj_sorted[i], adj_sorted[i + 1])
    out = np.empty(n)
    out[order] = np.minimum(adj_sorted, 1.0)
    return out


def _compute_pair_stats(O: int, wsum: float, cap: float, E: float) -> dict:
    out: dict = {}
    dens = (O / cap) if cap > 0 else 0.0
    out["cross_cluster_edge_density"] = round(dens, 12)
    if E > 0 and np.isfinite(E):
        z = float((O - E) / np.sqrt(E))
        log2_oe = float(np.log2((O + 0.5) / (E + 0.5)))
        p_enr = float(stats.poisson.sf(O - 1, E))
        fold = float(O / E) if E > 0 else 0.0
    else:
        z = float("nan")
        log2_oe = float("nan") if not (O > 0 and E == 0) else float("inf")
        p_enr = 1.0 if O == 0 else 0.0
        fold = float("inf") if O > 0 and E == 0 else 0.0
    out["relatedness_z"] = z if np.isfinite(z) else ""
    out["log2_obs_over_exp"] = log2_oe if np.isfinite(log2_oe) else ""
    out["poisson_p_enrichment"] = p_enr
    out["fold_enrichment_vs_mixing_null"] = round(fold, 8) if np.isfinite(fold) else ""
    out["n_genome_edges_expected_mixing_null"] = round(E, 8) if np.isfinite(E) else ""
    return out


def _parse_cluster_attr(d: dict, key: str = "leiden_hd_cluster") -> int | None:
    v = d.get(key, -1)
    if v == "" or v is None:
        return None
    try:
        c = int(float(v))
    except (TypeError, ValueError):
        return None
    if c < 0:
        return None
    return c


def label_map_from_graph(G: nx.Graph, *, label_attr: str = "leiden_hd_cluster") -> dict[str, int]:
    m: dict[str, int] = {}
    for n in G.nodes():
        c = _parse_cluster_attr(dict(G.nodes[n]), key=label_attr)
        if c is not None:
            m[str(n)] = c
    return m


def inter_cluster_pair_bincount(
    edges_u: np.ndarray,
    edges_v: np.ndarray,
    labels: np.ndarray,
    K: int,
) -> np.ndarray:
    cu = labels[edges_u]
    cv = labels[edges_v]
    mk = cu != cv
    if not np.any(mk):
        return np.zeros(K * K, dtype=np.int64)
    cu = cu[mk]
    cv = cv[mk]
    a = np.minimum(cu, cv).astype(np.int64)
    b = np.maximum(cu, cv).astype(np.int64)
    pair_id = a * K + b
    return np.bincount(pair_id, minlength=K * K).astype(np.int64)


def export_ie_cluster_overview(
    *,
    path_graphml: Path,
    out_overview_graphml: Path,
    out_pair_mixing_tsv: Path,
    out_internal_external_tsv: Path,
    label_attr: str = "leiden_hd_cluster",
    min_relatedness_z: float = 2.0,
    max_fdr: float = 0.05,
    include_all_observed_edges: bool = False,
) -> None:
    """Compact cluster-only graphml plus the mixing-null TSV."""
    if not path_graphml.is_file():
        print("[cluster overview skip] missing:", path_graphml)
        return

    print(f"[cluster overview] loading {path_graphml.name}…")
    G = nx.read_graphml(str(path_graphml))
    lid_map = label_map_from_graph(G, label_attr=label_attr)
    lid_counts: dict[int, int] = dict(Counter(lid_map.values()))

    inter: dict[tuple[int, int], dict] = defaultdict(lambda: {"n_edges": 0, "weight_sum": 0.0})
    intra_count: dict[int, int] = defaultdict(int)
    intra_wsum: dict[int, float] = defaultdict(float)

    for u, v, data in G.edges(data=True):
        su, sv = str(u), str(v)
        cu, cv = lid_map.get(su), lid_map.get(sv)
        if cu is None or cv is None:
            continue
        w = data.get("weight")
        try:
            w = float(w) if w is not None and str(w) != "" else 1.0
        except (TypeError, ValueError):
            w = 1.0
        if cu == cv:
            intra_count[cu] += 1
            intra_wsum[cu] += w
            continue
        a, b = (cu, cv) if cu < cv else (cv, cu)
        inter[(a, b)]["n_edges"] += 1
        inter[(a, b)]["weight_sum"] += w

    n_genomes_in_graph: dict[int, int] = defaultdict(int)
    for node in G.nodes():
        cid = lid_map.get(str(node))
        if cid is None:
            continue
        n_genomes_in_graph[cid] += 1

    del G

    def n_genomes_assigned(cid: int) -> int:
        return int(lid_counts.get(cid, 0))

    inter_incident: dict[int, int] = defaultdict(int)
    inter_weight_incident: dict[int, float] = defaultdict(float)
    for (a, b), st in inter.items():
        n_e = st["n_edges"]
        w_s = st["weight_sum"]
        inter_incident[a] += n_e
        inter_incident[b] += n_e
        inter_weight_incident[a] += w_s
        inter_weight_incident[b] += w_s

    all_clusters_lid = sorted(lid_counts.keys())
    T_cross = 0.0
    for i, ci in enumerate(all_clusters_lid):
        ni = lid_counts[ci]
        for j in range(i + 1, len(all_clusters_lid)):
            cj = all_clusters_lid[j]
            T_cross += float(ni) * float(lid_counts[cj])

    M_inter = sum(st["n_edges"] for st in inter.values())

    all_clusters = sorted(
        set(lid_counts.keys()) | set(intra_count.keys()) | set(inter_incident.keys()),
    )

    pair_rows: list[dict] = []
    p_enrich_list: list[float] = []
    stats_by_pair: dict[tuple[int, int], dict] = {}

    for i, ci in enumerate(all_clusters):
        ni = n_genomes_assigned(ci)
        for cj in all_clusters[i + 1 :]:
            nj = n_genomes_assigned(cj)
            a, b = (ci, cj) if ci < cj else (cj, ci)
            st = inter.get((a, b), {"n_edges": 0, "weight_sum": 0.0})
            O = int(st["n_edges"])
            wsum = float(st["weight_sum"])
            cap = float(ni) * float(nj)
            if T_cross > 0:
                E = M_inter * cap / T_cross
            else:
                E = float("nan")

            ps = _compute_pair_stats(O, wsum, cap, E)
            p_enrich_list.append(float(ps["poisson_p_enrichment"]))

            row = {
                "cluster_a": a,
                "cluster_b": b,
                "n_reps_a": ni,
                "n_reps_b": nj,
                "cross_cluster_pair_capacity": int(cap) if cap == int(cap) else cap,
                "n_genome_edges_observed": O,
                "edge_weight_sum_observed": wsum,
                **ps,
            }
            stats_by_pair[(a, b)] = row
            pair_rows.append(row)

    p_arr = np.array(p_enrich_list, dtype=float)
    fdr = _bh_fdr(p_arr)
    for idx, row in enumerate(pair_rows):
        row["poisson_p_enrichment_bh_fdr"] = float(fdr[idx])
        z = row.get("relatedness_z", "")
        f = row["poisson_p_enrichment_bh_fdr"]
        O = int(row["n_genome_edges_observed"])
        Estr = row.get("n_genome_edges_expected_mixing_null", "")
        try:
            E = float(Estr) if Estr != "" else float("nan")
        except (TypeError, ValueError):
            E = float("nan")
        try:
            zf = float(z) if z != "" else float("nan")
        except (TypeError, ValueError):
            zf = float("nan")
        z_ok = np.isfinite(zf) and zf >= min_relatedness_z
        f_ok = f <= max_fdr
        excess = O > E if np.isfinite(E) else O > 0
        row["passes_relatedness_filter"] = "yes" if (O >= 1 and excess and (z_ok or f_ok)) else "no"

    out_pair_mixing_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(pair_rows).to_csv(out_pair_mixing_tsv, sep="\t", index=False)
    print(f" → {out_pair_mixing_tsv} ({len(pair_rows)} pairs)")

    H = nx.Graph()
    H.graph["edge_policy"] = (
        "relatedness_filtered" if not include_all_observed_edges else "all_observed_edges"
    )
    H.graph["min_relatedness_z"] = str(min_relatedness_z)
    H.graph["max_fdr_bh"] = str(max_fdr)
    H.graph["source_graph"] = str(path_graphml.name)
    H.graph["description"] = (
        "IE filtered kNN: nodes are Leiden-HD superclusters; edges are aggregated "
        "between-cluster links kept when they exceed a mixing Poisson null."
    )

    for cid in all_clusters:
        nid = f"c{cid}"
        ng = n_genomes_assigned(cid)
        H.add_node(
            nid,
            cluster_id=int(cid),
            cluster_id_str=str(cid),
            node_kind="leiden_hd_cluster_ie",
            n_reps=int(ng),
            n_reps_in_graph=int(n_genomes_in_graph.get(cid, 0)),
            n_internal_edges=int(intra_count.get(cid, 0)),
            internal_edge_weight_sum=_scalar(intra_wsum.get(cid, 0.0)),
            n_inter_cluster_edges_incident=int(inter_incident.get(cid, 0)),
            inter_cluster_weight_sum_incident=_scalar(inter_weight_incident.get(cid, 0.0)),
            internal_edges_per_rep=round(intra_count.get(cid, 0) / max(ng, 1), 6),
        )

    n_added = 0
    for (a, b), st in inter.items():
        O = int(st["n_edges"])
        if O < 1:
            continue
        row = stats_by_pair.get((a, b))
        if row is None:
            continue
        if row.get("passes_relatedness_filter") != "yes" and not include_all_observed_edges:
            continue

        u, v = f"c{a}", f"c{b}"
        wsum = float(st["weight_sum"])
        mean_w = wsum / O if O else 0.0
        z = row.get("relatedness_z", "")
        log2oe = row.get("log2_obs_over_exp", "")
        E = row.get("n_genome_edges_expected_mixing_null", "")
        fdr_v = row.get("poisson_p_enrichment_bh_fdr", "")
        na = n_genomes_assigned(a)
        nb = n_genomes_assigned(b)
        cap = float(na * nb)
        dens = (O / cap) if cap > 0 else 0.0

        H.add_edge(
            u,
            v,
            interaction="aggregated_ie_knn_filtered",
            n_genome_edges=int(O),
            weight_sum=float(wsum),
            weight_mean=round(mean_w, 6),
            relatedness_z=z if z != "" else "",
            log2_obs_over_exp=log2oe if log2oe != "" else "",
            n_genome_edges_expected_mixing_null=E if E != "" else "",
            n_genome_edges_observed=int(O),
            cross_cluster_edge_density=round(dens, 12),
            poisson_p_enrichment_bh_fdr=fdr_v if fdr_v != "" else "",
            fold_enrichment_vs_mixing_null=row.get("fold_enrichment_vs_mixing_null", ""),
        )
        n_added += 1

    if n_added == 0 and not include_all_observed_edges:
        print(
            "[WARN] overview: no edge passed the relatedness filter. "
            "Retry with include_all_observed_edges=True or relax min_relatedness_z / max_fdr."
        )

    nx.write_graphml(H, str(out_overview_graphml))
    print(f" → {out_overview_graphml} ({H.number_of_nodes()} nodes, {n_added} edges)")

    rows_ie = []
    for cid in sorted(all_clusters):
        ng = n_genomes_assigned(cid)
        ne_in = int(intra_count.get(cid, 0))
        max_pairs = ng * (ng - 1) // 2 if ng > 1 else 0
        density_proxy = (ne_in / max_pairs) if max_pairs else 0.0
        rows_ie.append(
            {
                "cluster_id": cid,
                "n_reps": ng,
                "n_reps_in_graph": int(n_genomes_in_graph.get(cid, 0)),
                "n_internal_edges": ne_in,
                "internal_weight_sum": intra_wsum.get(cid, 0.0),
                "n_inter_cluster_edges_incident": int(inter_incident.get(cid, 0)),
                "inter_cluster_weight_sum_incident": inter_weight_incident.get(cid, 0.0),
                "internal_edges_per_rep": ne_in / ng if ng else 0.0,
                "internal_density_proxy": round(density_proxy, 8),
            }
        )
    pd.DataFrame(rows_ie).sort_values("cluster_id").to_csv(out_internal_external_tsv, sep="\t", index=False)
    print(f" → {out_internal_external_tsv} ({len(rows_ie)} clusters)")


def run_ie_cluster_pair_permutation(
    *,
    path_graphml: Path,
    out_pvalues_tsv: Path,
    out_perm_graphml_full: Path,
    out_perm_graphml_filtered: Path | None,
    label_attr: str = "leiden_hd_cluster",
    n_perm: int = 2000,
    seed: int = 42,
    min_obs: int = 1,
    max_p: float = 0.05,
    suggestive_p: float = 0.10,
    skip_filtered_graphml: bool = False,
) -> None:
    """Permute labels on a fixed topology (same idea as cluster_pair_permutation_test.py)."""
    if not path_graphml.is_file():
        print("[permutation skip] missing:", path_graphml)
        return

    print(f"[permutation] loading {path_graphml.name}…")
    G = nx.read_graphml(str(path_graphml))
    lid_map = label_map_from_graph(G, label_attr=label_attr)

    nodes = sorted([str(n) for n in G.nodes() if str(n) in lid_map], key=lambda x: (len(x), x))
    gid_to_idx = {g: i for i, g in enumerate(nodes)}
    N = len(nodes)
    labels = np.zeros(N, dtype=np.int32)
    for g, i in gid_to_idx.items():
        labels[i] = int(lid_map[g])

    max_c = int(labels.max())
    K = max_c + 1
    if K < 2:
        print("[permutation skip] need at least 2 clusters with non-negative labels.")
        return

    eu, ev, ew = [], [], []
    for u, v, data in G.edges(data=True):
        su, sv = str(u), str(v)
        if su not in gid_to_idx or sv not in gid_to_idx:
            continue
        w = data.get("weight", 1.0)
        try:
            w = float(w) if w is not None and str(w) != "" else 1.0
        except (TypeError, ValueError):
            w = 1.0
        eu.append(gid_to_idx[su])
        ev.append(gid_to_idx[sv])
        ew.append(w)

    edges_u = np.array(eu, dtype=np.int32)
    edges_v = np.array(ev, dtype=np.int32)
    edge_w = np.array(ew, dtype=np.float64)
    E = len(edges_u)
    del G

    obs = inter_cluster_pair_bincount(edges_u, edges_v, labels, K)
    wsum_obs = np.zeros(K * K, dtype=np.float64)
    for i in range(E):
        u, v = int(edges_u[i]), int(edges_v[i])
        cu, cv = int(labels[u]), int(labels[v])
        if cu == cv:
            continue
        a, b = (cu, cv) if cu < cv else (cv, cu)
        wsum_obs[a * K + b] += edge_w[i]

    print(f"… labelled reps: {N}, edges with both ends labelled: {E}, K={K}")

    rng = np.random.default_rng(seed)
    n_ge = np.zeros_like(obs, dtype=np.int64)

    print(f"… {n_perm} permutations…")
    for rep in range(n_perm):
        perm = rng.permutation(labels)
        c = inter_cluster_pair_bincount(edges_u, edges_v, perm, K)
        n_ge += (c >= obs).astype(np.int64)
        if (rep + 1) % max(1, n_perm // 10) == 0:
            print(f"    {rep + 1}/{n_perm}")

    p_emp = (n_ge + 1).astype(float) / (n_perm + 1)

    lid_counts_list = Counter(int(labels[i]) for i in range(N))
    rows = []
    for pid in np.nonzero(obs > 0)[0]:
        O = int(obs[pid])
        if O < min_obs:
            continue
        a = int(pid) // K
        b = int(pid) % K
        if a >= b:
            continue
        rows.append(
            {
                "cluster_a": a,
                "cluster_b": b,
                # same counting as cluster_pair_permutation_test.py (kNN edges between reps)
                "n_genome_edges_observed": O,
                "edge_weight_sum_observed": float(wsum_obs[pid]),
                "p_value_permutation": float(p_emp[pid]),
                "n_perm": n_perm,
            }
        )

    out_pvalues_tsv.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        print("[WARN] permutation: no between-cluster pair meets min_obs.")
        pd.DataFrame().to_csv(out_pvalues_tsv, sep="\t", index=False)
        return

    df = pd.DataFrame(rows).sort_values("p_value_permutation")
    df["p_value_permutation_bh_fdr"] = _bh_fdr(df["p_value_permutation"].values)
    df["neg_log10_p"] = -np.log10(np.maximum(df["p_value_permutation"].values, 1e-300))
    df["significance_tier"] = np.where(
        df["p_value_permutation"] <= max_p,
        "core",
        np.where(df["p_value_permutation"] <= suggestive_p, "suggestive", "ns"),
    )
    df["permutation_rank"] = np.arange(1, len(df) + 1)
    df.to_csv(out_pvalues_tsv, sep="\t", index=False)
    print(f" → {out_pvalues_tsv} ({len(df)} pairs)")

    def n_reps_assigned(cid: int) -> int:
        return int(lid_counts_list.get(cid, 0))

    def _build_graph(include_edge: Callable[[pd.Series], bool]) -> tuple[nx.Graph, int]:
        H = nx.Graph()
        H.graph["null_model"] = "label_permutation_fixed_multiset"
        H.graph["n_permutations"] = str(n_perm)
        H.graph["seed"] = str(seed)
        H.graph["source_graph"] = str(path_graphml.name)
        H.graph["significance_core_p"] = str(max_p)
        H.graph["significance_suggestive_p"] = str(suggestive_p)
        H.graph["description"] = (
            "IE kNN: between-cluster edges; p is an empirical permutation (fixed graph, same multiset of labels)."
        )

        all_c = sorted(set(int(x) for x in labels) | set(lid_counts_list.keys()))
        for cid in all_c:
            nid = f"c{cid}"
            attrs: dict = {
                "cluster_id": int(cid),
                "cluster_id_str": str(cid),
                "node_kind": "leiden_hd_cluster_ie",
                "n_reps": n_reps_assigned(cid),
            }
            H.add_node(nid, **attrs)

        n_e = 0
        for _, r in df.iterrows():
            if not include_edge(r):
                continue
            a, b = int(r["cluster_a"]), int(r["cluster_b"])
            O = int(r["n_genome_edges_observed"])
            wsum = float(r["edge_weight_sum_observed"])
            H.add_edge(
                f"c{a}",
                f"c{b}",
                interaction="inter_cluster_ie_knn_edges",
                n_genome_edges=O,
                weight_sum=round(wsum, 6),
                weight_mean=round(wsum / O, 6) if O else 0.0,
                p_value_permutation=float(r["p_value_permutation"]),
                p_value_permutation_bh_fdr=float(r["p_value_permutation_bh_fdr"]),
                neg_log10_p=round(float(r["neg_log10_p"]), 6),
                significance_tier=str(r["significance_tier"]),
                permutation_rank=int(r["permutation_rank"]),
            )
            n_e += 1
        return H, n_e

    H_full, n_full = _build_graph(lambda _r: True)
    nx.write_graphml(H_full, str(out_perm_graphml_full))
    print(f" → {out_perm_graphml_full} ({H_full.number_of_nodes()} nodes, {n_full} edges)")

    if not skip_filtered_graphml and out_perm_graphml_filtered is not None:
        H_filt, n_filt = _build_graph(lambda r: r["p_value_permutation"] <= max_p)
        nx.write_graphml(H_filt, str(out_perm_graphml_filtered))
        print(f" → {out_perm_graphml_filtered} ({H_filt.number_of_nodes()} nodes, {n_filt} edges, p≤{max_p})")


def merge_ie_cluster_overview_with_permutation(
    *,
    path_overview_graphml: Path,
    path_perm_full_graphml: Path,
    out_graphml: Path,
) -> nx.Graph | None:
    """Merge the Poisson-overview metagraph with the full permutation graph (all pairs and p-values).

    * Edges: the union of both graphs (often 40 mixing edges plus up to K(K-1)/2 observed permutation pairs).
    * On shared edges: overview attributes are kept; the permutation layer adds
      ``p_*``, ``significance_tier``, and so on. ``interaction``, ``n_genome_edges``, and weights
      from the permutation side are written as ``perm_*`` so they do not overwrite the mixing statistics.
    * Nodes: attribute union (overview wins, because it carries more within-cluster fields).

    Returns the built ``Graph``, or ``None`` when an input is missing.
    """
    if not path_overview_graphml.is_file() or not path_perm_full_graphml.is_file():
        print(
            "[merge overview+perm skip] both files are required:",
            path_overview_graphml.name,
            path_perm_full_graphml.name,
        )
        return None

    G_o = nx.read_graphml(str(path_overview_graphml))
    G_p = nx.read_graphml(str(path_perm_full_graphml))

    H = nx.Graph()
    H.graph["merge_of"] = "ie_cluster_overview + ie_cluster_overview_permutation_full"
    H.graph["source_overview"] = str(path_overview_graphml.name)
    H.graph["source_permutation_full"] = str(path_perm_full_graphml.name)
    H.graph["description"] = (
        "Union: Poisson/relatedness edges (overview) plus every observed between-cluster pair from "
        "the permutations; on the overlap, mixing fields plus permutation p-values (perm interaction/weights stored as perm_*)."
    )

    nodes = set(str(n) for n in G_o.nodes()) | set(str(n) for n in G_p.nodes())
    for n in sorted(nodes):
        d: dict[str, Any] = {}
        if n in G_o:
            for k, v in dict(G_o.nodes[n]).items():
                d[k] = _scalar(v)
        if n in G_p:
            for k, v in dict(G_p.nodes[n]).items():
                if k not in d or d[k] == "":
                    d[k] = _scalar(v)
        H.add_node(n, **{k: v for k, v in d.items() if v != ""})

    perm_overlay_keys = {
        "p_value_permutation",
        "p_value_permutation_bh_fdr",
        "neg_log10_p",
        "significance_tier",
        "permutation_rank",
    }
    perm_other_conflict = {"n_genome_edges", "weight_sum", "weight_mean"}

    def merge_edge_attrs(d_ov: dict | None, d_perm: dict | None) -> dict:
        out: dict = {}
        tags: list[str] = []
        if d_ov:
            tags.append("mixing_overview")
            for k, v in d_ov.items():
                out[k] = _scalar(v)
        if d_perm:
            tags.append("permutation_full")
            for k, v in d_perm.items():
                if k == "interaction":
                    out["perm_interaction"] = _scalar(v)
                elif k in perm_other_conflict:
                    out[f"perm_{k}"] = _scalar(v)
                else:
                    sk = _scalar(v)
                    if k not in out or out[k] == "":
                        out[k] = sk
                    elif sk != "" and k in perm_overlay_keys:
                        out[k] = sk
        out["edge_sources"] = "+".join(tags)
        out["has_mixing_stats"] = "yes" if d_ov else "no"
        out["has_permutation_stats"] = "yes" if d_perm else "no"
        return out

    seen: set[tuple[str, str]] = set()

    for u, v in G_o.edges():
        su, sv = str(u), str(v)
        a, b = (su, sv) if su < sv else (sv, su)
        seen.add((a, b))
        ed = merge_edge_attrs(dict(G_o[su][sv]), dict(G_p.get_edge_data(su, sv)) if G_p.has_edge(su, sv) else None)
        H.add_edge(a, b, **{k: v for k, v in ed.items() if v != ""})

    for u, v in G_p.edges():
        su, sv = str(u), str(v)
        a, b = (su, sv) if su < sv else (sv, su)
        if (a, b) in seen:
            continue
        seen.add((a, b))
        ed = merge_edge_attrs(None, dict(G_p[su][sv]))
        H.add_edge(a, b, **{k: v for k, v in ed.items() if v != ""})

    out_graphml.parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(H, str(out_graphml))
    print(
        f" → {out_graphml} ({H.number_of_nodes()} nodes, {H.number_of_edges()} edges; "
        "mixing ∪ permutation)"
    )
    return H
