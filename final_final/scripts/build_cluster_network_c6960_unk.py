#!/usr/bin/env python3
"""Cluster-level network for Stage 9 unclassified Leiden clusters.

Pipeline (mirrors final/clustering_filtered_ie/.../filtered_ie_cluster_network_*.graphml):

1. Enrich the IE-level kNN graphml with Leiden cluster labels and per-cluster
   biology (host phylum, atlas IE_core, hallmark roles, padloc systems).
2. ``export_ie_cluster_overview`` → mixing-null + Poisson FDR overview graph.
3. ``run_ie_cluster_pair_permutation`` → permutation p-values + filtered graph.
4. ``merge_ie_cluster_overview_with_permutation`` → merged graph.
5. Apply Cytoscape-style edge derivations (edge_width, edge_opacity,
   edge_significance, edge_strength_score) and write canonical / strong subsets.
6. Render an inline multi-panel matplotlib summary (network + edge metrics).

Outputs (under ``final_final/artifacts/nb04_cluster_unclassified``):
  - knn_graph_c6960_unk_enriched.graphml
  - cluster_overview_mixing_null.graphml + .tsv (pair stats) + internal/external tsv
  - cluster_overview_permutation.graphml + .tsv (full p-values)
  - cluster_overview_perm_filtered.graphml
  - cluster_overview_mixing_perm_merged.graphml
  - cluster_network_cytoscape.graphml (canonical)
  - cluster_network_strong.graphml (fold ≥ 2)

Figures (under ``final_final/figures/04_cluster_unclassified_stage8``):
  - 13_cluster_network_layout.png          (spring layout, Cytoscape-style)
  - 14_cluster_network_panel.png           (multi-panel: layout + matrices)
"""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.collections import LineCollection

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import REPO
NB_LIB = REPO / "final_final" / "03_clustering" / "lib"
sys.path.insert(0, str(NB_LIB))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fig_cluster_network_layout  # noqa: E402  (local script dir added above)
from clustering_ie_graph_analysis import (  # noqa: E402  (after sys.path)
    export_ie_cluster_overview,
    merge_ie_cluster_overview_with_permutation,
    run_ie_cluster_pair_permutation,
)

ART = REPO / "final_final/artifacts/nb04_cluster_unknown_6960"
FIG = REPO / "final_final/figures/04_cluster_unknown_6960"
FIG.mkdir(parents=True, exist_ok=True)

KNN_IN = ART / "knn_graph_c6960_unk.graphml"
KNN_ENRICHED = ART / "knn_graph_c6960_unk_enriched.graphml"

LEIDEN_TSV = ART / "leiden_clusters_c6960_unk.tsv"
IE_ROWS = ART / "ie_rep_rows_c6960_unk.tsv"
BIO_TSV = ART / "cluster_biology_c6960_unk.tsv"
HALL_TSV = ART / "cluster_hallmark_composition_c6960_unk.tsv"
PAD_TSV = ART / "cluster_padloc_summary_c6960_unk.tsv"
HOST_PUR_TSV = ART / "cluster_host_purity_c6960_unk.tsv"
MARKERS_TSV = ART / "cluster_pfam_markers_c6960_unk.tsv"

OVERVIEW_GRAPHML = ART / "cluster_overview_mixing_null.graphml"
OVERVIEW_TSV = ART / "cluster_pair_mixing_null.tsv"
INT_EXT_TSV = ART / "cluster_internal_external_summary.tsv"

PERM_PVAL_TSV = ART / "cluster_pair_permutation_pvalues.tsv"
PERM_FULL_GRAPHML = ART / "cluster_overview_permutation.graphml"
PERM_FILT_GRAPHML = ART / "cluster_overview_perm_filtered.graphml"
MERGED_GRAPHML = ART / "cluster_overview_mixing_perm_merged.graphml"

CYTOSCAPE_GRAPHML = ART / "cluster_network_cytoscape.graphml"
STRONG_GRAPHML = ART / "cluster_network_strong.graphml"

NB_LABEL_ATTR = "leiden_cluster"


# ---------------------------------------------------------------------------
# 1. Enrich IE-level kNN graphml with cluster labels + per-cluster biology.
# ---------------------------------------------------------------------------
def enrich_knn_graph() -> nx.Graph:
    print(f"[enrich] reading {KNN_IN.name}…")
    G = nx.read_graphml(str(KNN_IN))

    leid = pd.read_csv(LEIDEN_TSV, sep="\t").reset_index(drop=True)
    rows = pd.read_csv(IE_ROWS, sep="\t").reset_index(drop=True)

    leid["node_id"] = "n" + leid.index.astype(str)
    rows["node_id"] = "n" + rows.index.astype(str)
    leid_map = dict(zip(leid["node_id"].astype(str), leid["leiden_cluster"].astype(int)))
    rep_map = dict(zip(rows["node_id"].astype(str), rows.get("rep_id", rows.get("ie_id", "")).astype(str)))
    sample_map = dict(zip(rows["node_id"].astype(str), rows.get("sample_id", "").astype(str)))
    ie_map = dict(zip(rows["node_id"].astype(str), rows.get("ie_id", rows.get("rep_id", "")).astype(str)))

    n_assigned = 0
    n_noise = 0
    for n in G.nodes():
        sn = str(n)
        cid = leid_map.get(sn)
        if cid is None:
            continue
        cid_i = int(cid)
        # Tiny Leiden communities are collapsed to -1; omit from cluster network
        if cid_i < 0:
            G.nodes[n][NB_LABEL_ATTR] = -1
            G.nodes[n]["leiden_hd_cluster"] = -1
            n_noise += 1
        else:
            G.nodes[n][NB_LABEL_ATTR] = cid_i
            G.nodes[n]["leiden_hd_cluster"] = cid_i  # alias for the shared utility
            n_assigned += 1
        if rep_map.get(sn):
            G.nodes[n]["rep_id"] = rep_map[sn]
        if sample_map.get(sn):
            G.nodes[n]["sample_id"] = sample_map[sn]
        if ie_map.get(sn):
            G.nodes[n]["ie_id"] = ie_map[sn]

    print(
        f"[enrich] cluster labels written for {n_assigned}/{G.number_of_nodes()} nodes "
        f"(noise/unassigned={n_noise})"
    )

    nx.write_graphml(G, str(KNN_ENRICHED))
    print(f"[enrich] → {KNN_ENRICHED}")
    return G


# ---------------------------------------------------------------------------
# 2-3-4. Mixing-null overview, permutation p-values, merged graph.
# ---------------------------------------------------------------------------
def build_overview_and_perm() -> None:
    export_ie_cluster_overview(
        path_graphml=KNN_ENRICHED,
        out_overview_graphml=OVERVIEW_GRAPHML,
        out_pair_mixing_tsv=OVERVIEW_TSV,
        out_internal_external_tsv=INT_EXT_TSV,
        label_attr="leiden_hd_cluster",
        min_relatedness_z=2.0,
        max_fdr=0.05,
        include_all_observed_edges=False,
    )

    run_ie_cluster_pair_permutation(
        path_graphml=KNN_ENRICHED,
        out_pvalues_tsv=PERM_PVAL_TSV,
        out_perm_graphml_full=PERM_FULL_GRAPHML,
        out_perm_graphml_filtered=PERM_FILT_GRAPHML,
        label_attr="leiden_hd_cluster",
        n_perm=2000,
        seed=42,
        min_obs=1,
        max_p=0.05,
        suggestive_p=0.10,
        skip_filtered_graphml=False,
    )

    merge_ie_cluster_overview_with_permutation(
        path_overview_graphml=OVERVIEW_GRAPHML,
        path_perm_full_graphml=PERM_FULL_GRAPHML,
        out_graphml=MERGED_GRAPHML,
    )


# ---------------------------------------------------------------------------
# 5. Cytoscape-ready cluster network with derived visual fields.
# ---------------------------------------------------------------------------
DROP_EDGE_ATTRS = {
    "perm_interaction",
    "perm_n_genome_edges",
    "perm_weight_sum",
    "perm_weight_mean",
    "has_mixing_stats",
    "has_permutation_stats",
}


def _to_float(val) -> float:
    try:
        f = float(val)
        if math.isnan(f) or math.isinf(f):
            return float("nan")
        return f
    except (TypeError, ValueError):
        return float("nan")


def _classify(sig_tier: str, p_perm_fdr: float, p_poiss_fdr: float) -> str:
    if isinstance(sig_tier, str) and sig_tier in ("core", "suggestive"):
        return sig_tier
    if not math.isnan(p_perm_fdr) and p_perm_fdr <= 0.05:
        return "perm-significant"
    if not math.isnan(p_poiss_fdr) and p_poiss_fdr <= 0.05:
        return "mixing-significant"
    return "ns"


def _clean_edge_attrs(G: nx.Graph) -> None:
    for u, v, attr in list(G.edges(data=True)):
        for k in DROP_EDGE_ATTRS:
            attr.pop(k, None)
        for k in ("n_genome_edges", "n_genome_edges_observed", "permutation_rank"):
            if k in attr:
                try:
                    attr[k] = int(float(attr[k]))
                except (TypeError, ValueError):
                    pass
        for k in (
            "weight_sum",
            "weight_mean",
            "relatedness_z",
            "log2_obs_over_exp",
            "n_genome_edges_expected_mixing_null",
            "cross_cluster_edge_density",
            "poisson_p_enrichment_bh_fdr",
            "fold_enrichment_vs_mixing_null",
            "p_value_permutation",
            "p_value_permutation_bh_fdr",
            "neg_log10_p",
        ):
            if k in attr:
                try:
                    attr[k] = float(attr[k])
                except (TypeError, ValueError):
                    pass

        z = _to_float(attr.get("relatedness_z"))
        fold = _to_float(attr.get("fold_enrichment_vs_mixing_null"))
        p_poiss = _to_float(attr.get("poisson_p_enrichment_bh_fdr"))
        p_perm = _to_float(attr.get("p_value_permutation_bh_fdr"))

        edge_width = 0.0 if math.isnan(z) else max(0.0, min(z, 25.0))
        if math.isnan(p_poiss) or p_poiss <= 0:
            opacity = 0.0
        else:
            opacity = min(1.0, -math.log10(max(p_poiss, 1e-300)) / 50.0)
        sig = _classify(str(attr.get("significance_tier", "")), p_perm, p_poiss)
        attr["edge_width"] = float(edge_width)
        attr["edge_opacity"] = float(round(opacity, 4))
        attr["edge_significance"] = sig
        attr["edge_strength_score"] = float(
            round(0.6 * edge_width + 4.0 * (fold if not math.isnan(fold) else 0.0), 3)
        )


def _enrich_all_graphmls() -> None:
    """Attach full marker / hallmark / padloc / host metadata to all cluster graphmls."""
    scripts_dir = Path(__file__).resolve().parent
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    from enrich_cluster_overview_c6960_unk import enrich_one

    for name in (
        OVERVIEW_GRAPHML.name,
        PERM_FULL_GRAPHML.name,
        PERM_FILT_GRAPHML.name,
        MERGED_GRAPHML.name,
        CYTOSCAPE_GRAPHML.name,
        STRONG_GRAPHML.name,
    ):
        path = ART / name
        if path.is_file():
            enrich_one(path)


def make_cytoscape_subsets(min_fold: float = 2.0) -> None:
    G = nx.read_graphml(str(MERGED_GRAPHML))
    print(f"[cytoscape] merged graph: nodes={G.number_of_nodes()} edges={G.number_of_edges()}")

    keep_canon: list[tuple[str, str]] = []
    keep_strong: list[tuple[str, str]] = []
    for u, v, attr in G.edges(data=True):
        sources = str(attr.get("edge_sources", ""))
        if "mixing_overview" in sources:
            keep_canon.append((u, v))
            fold = _to_float(attr.get("fold_enrichment_vs_mixing_null"))
            if not math.isnan(fold) and fold >= min_fold:
                keep_strong.append((u, v))

    def _subset(edges: list[tuple[str, str]]) -> nx.Graph:
        H = nx.Graph()
        H.graph.update(G.graph)
        for n, attr in G.nodes(data=True):
            H.add_node(n, **attr)
        for u, v in edges:
            H.add_edge(u, v, **G.edges[u, v])
        return H

    canon = _subset(keep_canon)
    strong = _subset(keep_strong)
    _clean_edge_attrs(canon)
    _clean_edge_attrs(strong)

    canon.graph["description"] = (
        "Cytoscape-ready cluster network for Stage 9 unclassified Leiden clusters. "
        "Edges = cross-cluster IE-IE kNN aggregates that passed mixing-null filter "
        "(Poisson FDR ≤ 0.05 + relatedness_z ≥ 2). Derived visual fields: edge_width "
        "(relatedness_z, ≤ 25); edge_opacity (−log10(Poisson FDR)/50, ≤ 1); "
        "edge_significance ∈ {core, perm-significant, mixing-significant, ns}; "
        "edge_strength_score = 0.6·z + 4·fold. "
        "Node colour: marker_family / archetype_color (top PFAM marker). "
        "Node size: node_size (Cytoscape) / node_size_matplotlib; area ∝ node_weight (= n_reps)."
    )
    strong.graph["description"] = canon.graph["description"] + f" Strong subset: fold ≥ {min_fold}."

    nx.write_graphml(canon, str(CYTOSCAPE_GRAPHML))
    nx.write_graphml(strong, str(STRONG_GRAPHML))
    print(f"[canonical] {CYTOSCAPE_GRAPHML} (nodes={canon.number_of_nodes()} edges={canon.number_of_edges()})")
    print(f"[strong]    {STRONG_GRAPHML} (nodes={strong.number_of_nodes()} edges={strong.number_of_edges()})")


# ---------------------------------------------------------------------------
# 6. Inline matplotlib rendering of the cluster network (Cytoscape style).
# ---------------------------------------------------------------------------
CLUSTER_COLORS = {
    0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c", 3: "#9467bd", 4: "#d62728",
    5: "#8c564b", 6: "#e377c2", 7: "#17becf", 8: "#bcbd22", 9: "#aec7e8",
    10: "#98df8a", 11: "#ff9896", 12: "#c5b0d5", 13: "#c49c94",
}


def _cluster_id(n) -> int:
    try:
        return int(str(n)[1:]) if str(n).startswith("c") else int(n)
    except (TypeError, ValueError):
        return -1


def _cluster_color(cid: int) -> str:
    return CLUSTER_COLORS.get(cid, "#7f7f7f")


MARKER_FAMILY_PALETTE = {
    "HTH_regulatory": "#d62728",
    "replication": "#ff7f0e",
    "integrase_SAM": "#1f77b4",
    "integrase_serine": "#e377c2",
    "defense_RM": "#2ca02c",
    "toxin_antitoxin": "#9467bd",
    "transposase": "#8c564b",
    "AAA_motor": "#bcbd22",
    "TPR_scaffold": "#17becf",
    "coiled_coil": "#aec7e8",
    "zn_ribbon": "#ff9896",
    "minimal_integrase": "#5fa9d6",
    "other": "#7f7f7f",
}


def _marker_color(attr: dict) -> str:
    for key in ("marker_color", "archetype_color"):
        c = attr.get(key)
        if isinstance(c, str) and c.startswith("#"):
            return c
    fam = str(attr.get("marker_family", attr.get("archetype_family", "other")))
    return MARKER_FAMILY_PALETTE.get(fam, MARKER_FAMILY_PALETTE["other"])


def _node_n_ie(attr: dict) -> float:
    for key in ("node_size_matplotlib", "node_weight", "bio_n_ie", "bio_n_IE", "n_reps"):
        val = attr.get(key)
        if val is None or val == "":
            continue
        try:
            f = float(val)
            if key == "node_size_matplotlib":
                return f  # already visual size
            if f > 0:
                return f
        except (TypeError, ValueError):
            continue
    return 0.0


def _node_scatter_sizes(G: nx.Graph, *, node_min_size: float = 250,
                         node_max_size: float = 2000) -> np.ndarray:
    """Return matplotlib scatter sizes; prefer precomputed node_size_matplotlib."""
    sizes: list[float] = []
    raw_counts: list[float] = []
    for _, attr in G.nodes(data=True):
        if attr.get("node_size_matplotlib") not in (None, ""):
            try:
                sizes.append(float(attr["node_size_matplotlib"]))
                raw_counts.append(parse_n_ie_fallback(attr))
                continue
            except (TypeError, ValueError):
                pass
        raw_counts.append(_node_n_ie(attr))
        sizes.append(float("nan"))
    raw = np.asarray(raw_counts, dtype=float)
    if np.any(np.isnan(sizes)):
        max_n = float(raw.max()) if raw.size else 1.0
        for i, n_ie in enumerate(raw):
            if math.isnan(sizes[i]):
                frac = n_ie / max(max_n, 1.0)
                sizes[i] = node_min_size + frac * (node_max_size - node_min_size)
    return np.asarray(sizes, dtype=float)


def parse_n_ie_fallback(attr: dict) -> float:
    for key in ("node_weight", "bio_n_ie", "bio_n_IE", "n_reps", "n_reps_in_graph"):
        val = attr.get(key)
        if val is None or val == "":
            continue
        try:
            f = float(val)
            if f > 0:
                return f
        except (TypeError, ValueError):
            continue
    return 0.0


def _size_aware_nudge(G: nx.Graph, pos: dict[str, np.ndarray], n_map: dict[str, float],
                      steps: int = 50, strength: float = 0.07) -> dict[str, np.ndarray]:
    """Push apart nodes so large clusters occupy more layout space."""
    nodes = list(G.nodes())
    out = {n: np.array(p, dtype=float) for n, p in pos.items()}
    for _ in range(steps):
        moved = False
        for i, n in enumerate(nodes):
            for m in nodes[i + 1:]:
                delta = out[m] - out[n]
                dist = float(np.linalg.norm(delta)) + 1e-9
                min_sep = strength * (math.sqrt(n_map.get(n, 1)) + math.sqrt(n_map.get(m, 1)))
                if dist < min_sep:
                    push = (min_sep - dist) / dist * delta * 0.35
                    out[n] -= push
                    out[m] += push
                    moved = True
        if not moved:
            break
    return out


def _layout(G: nx.Graph, seed: int = 7) -> dict[str, np.ndarray]:
    """Spring layout; edge weight = relatedness_z; node separation ∝ sqrt(n_IE)."""
    if G.number_of_edges() == 0:
        return nx.circular_layout(G)
    n_map = {n: parse_n_ie_fallback(a) for n, a in G.nodes(data=True)}
    weights = []
    for _, _, d in G.edges(data=True):
        z = _to_float(d.get("relatedness_z"))
        weights.append(max(0.5, z if not math.isnan(z) else 1.0))
    H = nx.Graph()
    H.add_nodes_from(G.nodes(data=True))
    for (u, v, _), w in zip(G.edges(data=True), weights):
        H.add_edge(u, v, weight=w)
    # Larger k → more space; scale mildly by median cluster size
    med_n = float(np.median(list(n_map.values()))) if n_map else 1.0
    k_base = 2.2 / math.sqrt(max(1, H.number_of_nodes()))
    k = k_base * (1.0 + 0.15 * math.log10(max(med_n, 10)))
    pos = nx.spring_layout(H, seed=seed, k=k, iterations=250, weight="weight")
    pos = _size_aware_nudge(G, pos, n_map)
    iso = [n for n in G.nodes() if G.degree(n) == 0]
    if iso:
        non_iso_xy = np.asarray([pos[n] for n in G.nodes() if n not in iso])
        if non_iso_xy.size:
            cx, cy = non_iso_xy.mean(axis=0)
            r = float(np.linalg.norm(non_iso_xy - np.array([cx, cy]), axis=1).max() * 1.25)
        else:
            cx, cy, r = 0.0, 0.0, 1.0
        n_iso = len(iso)
        for i, n in enumerate(iso):
            theta = 2 * math.pi * i / max(1, n_iso) - math.pi / 2
            pos[n] = np.array([cx + r * math.cos(theta), cy + r * math.sin(theta)])
    return pos


def _draw_network(ax, G: nx.Graph, pos: dict, *, title: str, edge_min_width: float = 0.4,
                  edge_max_width: float = 4.0, node_min_size: float = 250,
                  node_max_size: float = 2000) -> None:
    colors_node = []
    labels = {}
    marker_labels = {}
    for n, attr in G.nodes(data=True):
        cid = _cluster_id(n)
        colors_node.append(_cluster_color(cid))
        labels[n] = str(cid)
        marker_labels[n] = str(
            attr.get("marker_label_short", attr.get("archetype_short", "")) or cid
        )
    node_sizes = _node_scatter_sizes(G, node_min_size=node_min_size, node_max_size=node_max_size)

    edge_segments = []
    edge_colors = []
    edge_widths = []
    for u, v, d in G.edges(data=True):
        if u not in pos or v not in pos:
            continue
        edge_segments.append([pos[u], pos[v]])
        z = _to_float(d.get("relatedness_z"))
        ew = float(d.get("edge_width", z if not math.isnan(z) else 1.0))
        edge_widths.append(edge_min_width + min(ew, 25.0) / 25.0 * (edge_max_width - edge_min_width))
        sig = str(d.get("edge_significance", "ns"))
        opacity = float(d.get("edge_opacity", 0.6))
        if sig == "core":
            base = (0.10, 0.10, 0.10)
        elif sig in ("perm-significant", "mixing-significant"):
            base = (0.18, 0.45, 0.78)
        elif sig == "suggestive":
            base = (0.55, 0.55, 0.55)
        else:
            base = (0.75, 0.75, 0.75)
        a = max(0.18, min(1.0, 0.35 + 0.7 * opacity))
        edge_colors.append((*base, a))

    if edge_segments:
        lc = LineCollection(edge_segments, colors=edge_colors, linewidths=edge_widths, zorder=1)
        ax.add_collection(lc)

    nodes_xy = np.asarray([pos[n] for n in G.nodes()])
    ax.scatter(nodes_xy[:, 0], nodes_xy[:, 1], s=node_sizes, c=colors_node,
               edgecolors="white", linewidths=1.2, zorder=2)
    for n, (x, y) in zip(G.nodes(), nodes_xy):
        ax.text(x, y, labels[n], ha="center", va="center", fontsize=8.5, weight="bold",
                color="white", zorder=3)
        fam = str(G.nodes[n].get("marker_family", G.nodes[n].get("archetype_family", "")))
        short = marker_labels[n]
        if fam or short:
            ax.text(
                x, y - 0.06, short, ha="center", va="top", fontsize=6.5,
                color="#333", zorder=3, clip_on=False,
            )

    ax.set_axis_off()
    ax.set_title(title, fontsize=11, weight="bold")


def render_figures() -> None:
    rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})

    canon = nx.read_graphml(str(CYTOSCAPE_GRAPHML))
    strong = nx.read_graphml(str(STRONG_GRAPHML))
    if canon.number_of_nodes() == 0:
        print("[render] empty canon graph — skip figures")
        return

    pos = _layout(canon, seed=11)

    fig_cluster_network_layout.render(CYTOSCAPE_GRAPHML)

    # --- Figure 14: multi-panel summary ---
    fig = plt.figure(figsize=(16.5, 11.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.6, 1, 1], height_ratios=[1, 1],
                          hspace=0.28, wspace=0.28)
    ax_net = fig.add_subplot(gs[:, 0])
    _draw_network(ax_net, canon, pos, title="A. Cluster network (mixing-null filter)")
    ax_strong = fig.add_subplot(gs[0, 1])
    _draw_network(ax_strong, strong, _layout(strong, seed=11),
                  title=f"B. Strong subset (fold ≥ 2)  ·  edges = {strong.number_of_edges()}",
                  node_min_size=120, node_max_size=900)

    ax_size = fig.add_subplot(gs[0, 2])
    n_ie_per_cluster = sorted(
        [(int(str(n)[1:]) if str(n).startswith("c") else int(n), parse_n_ie_fallback(a))
         for n, a in canon.nodes(data=True)], key=lambda x: -x[1])
    cids = [c for c, _ in n_ie_per_cluster]
    sizes = [s for _, s in n_ie_per_cluster]
    ax_size.bar(range(len(sizes)), sizes,
                color=[_cluster_color(c) for c, _ in n_ie_per_cluster])
    ax_size.set_xticks(range(len(sizes)))
    ax_size.set_xticklabels([f"c{c}" for c in cids], rotation=0, fontsize=8)
    ax_size.set_ylabel("# IE per cluster")
    ax_size.set_title("C. Cluster sizes (cluster-coloured)", fontsize=10)
    for spine in ("top", "right"):
        ax_size.spines[spine].set_visible(False)

    ax_mat = fig.add_subplot(gs[1, 1])
    cids_sorted = sorted(set(int(str(n)[1:]) if str(n).startswith("c") else int(n) for n in canon.nodes()))
    K = len(cids_sorted)
    M = np.zeros((K, K))
    cid_idx = {c: i for i, c in enumerate(cids_sorted)}
    for u, v, d in canon.edges(data=True):
        cu = int(str(u)[1:]) if str(u).startswith("c") else int(u)
        cv = int(str(v)[1:]) if str(v).startswith("c") else int(v)
        z = _to_float(d.get("relatedness_z"))
        if math.isnan(z):
            continue
        i, j = cid_idx[cu], cid_idx[cv]
        M[i, j] = z
        M[j, i] = z
    im = ax_mat.imshow(M, cmap="magma_r", aspect="auto")
    ax_mat.set_xticks(range(K))
    ax_mat.set_yticks(range(K))
    ax_mat.set_xticklabels([f"c{c}" for c in cids_sorted], fontsize=7, rotation=90)
    ax_mat.set_yticklabels([f"c{c}" for c in cids_sorted], fontsize=7)
    ax_mat.set_title("D. relatedness_z (mixing-null only)", fontsize=10)
    fig.colorbar(im, ax=ax_mat, fraction=0.046, pad=0.04, label="z-score")

    ax_pv = fig.add_subplot(gs[1, 2])
    if OVERVIEW_TSV.exists():
        df_p = pd.read_csv(OVERVIEW_TSV, sep="\t")
        if len(df_p):
            df_p = df_p.dropna(subset=["relatedness_z", "fold_enrichment_vs_mixing_null"]).copy()
            mask_pass = df_p["passes_relatedness_filter"].astype(str) == "yes"
            ax_pv.scatter(df_p.loc[~mask_pass, "relatedness_z"],
                          df_p.loc[~mask_pass, "fold_enrichment_vs_mixing_null"].clip(0, 25),
                          c="#cccccc", s=20, alpha=0.7, label="ns", edgecolors="white", linewidths=0.4)
            ax_pv.scatter(df_p.loc[mask_pass, "relatedness_z"],
                          df_p.loc[mask_pass, "fold_enrichment_vs_mixing_null"].clip(0, 25),
                          c="#1f77b4", s=44, alpha=0.95, label="mixing-significant", edgecolors="white", linewidths=0.6)
            mask_strong = mask_pass & (df_p["fold_enrichment_vs_mixing_null"] >= 2.0)
            ax_pv.scatter(df_p.loc[mask_strong, "relatedness_z"],
                          df_p.loc[mask_strong, "fold_enrichment_vs_mixing_null"].clip(0, 25),
                          c="#d95f02", s=70, alpha=1.0, label="strong (fold ≥ 2)", edgecolors="white", linewidths=0.8)
            ax_pv.axvline(2, ls="--", lw=1, color="#888", alpha=0.6)
            ax_pv.axhline(2, ls="--", lw=1, color="#888", alpha=0.6)
            ax_pv.set_xlabel("relatedness_z (Poisson)")
            ax_pv.set_ylabel("fold enrichment over null  (clip 25)")
            ax_pv.set_title("E. Pair significance landscape", fontsize=10)
            ax_pv.legend(fontsize=8, frameon=False, loc="upper left")
            for spine in ("top", "right"):
                ax_pv.spines[spine].set_visible(False)

    fig.suptitle(
        "Cluster-level network — 6960 Unknown Leiden clusters",
        fontsize=12, weight="bold", y=0.995,
    )
    fig.subplots_adjust(left=0.04, right=0.98, top=0.93, bottom=0.08)
    fig.savefig(FIG / "14_cluster_network_panel.png", dpi=160, bbox_inches="tight", pad_inches=0.35)
    fig.savefig(FIG / "14_cluster_network_panel.svg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)
    print("[render] →", FIG / "14_cluster_network_panel.png")


def main() -> None:
    enrich_knn_graph()
    build_overview_and_perm()
    make_cytoscape_subsets(min_fold=2.0)
    _enrich_all_graphmls()
    render_figures()
    print("[done] all artifacts in", ART)
    print("[done] figures in", FIG)


if __name__ == "__main__":
    main()
