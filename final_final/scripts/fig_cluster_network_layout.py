#!/usr/bin/env python3
"""Cytoscape-style rendering of the cluster-level network (figure 13).

Reads ``cluster_network_cytoscape.graphml`` and draws flat circles laid out by a
force-directed layout: node area ∝ n IE, node colour = marker family, label =
top PFAM short name. No titles, no captions — only the colour legend.

Standalone entry point; also imported by ``build_cluster_network_c6960_unk.py``
so the pipeline and a manual restyle produce the same figure.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

REPO = Path("/home/lam34/MGE_finder")
ART = REPO / "final_final/artifacts/nb04_cluster_unknown_6960"
FIG = REPO / "final_final/figures/04_cluster_unknown_6960"

CYTOSCAPE_GRAPHML = ART / "cluster_network_cytoscape.graphml"
OUT_STEM = "13_cluster_network_layout"

FAMILY_FALLBACK_COLOR = "#7f7f7f"
NODE_MIN_AREA = 320.0
NODE_MAX_AREA = 2600.0
STRONG_FOLD = 2.0


def _f(val, default: float = float("nan")) -> float:
    try:
        out = float(val)
    except (TypeError, ValueError):
        return default
    return default if math.isnan(out) or math.isinf(out) else out


def _n_ie(attr: dict) -> float:
    for key in ("node_weight", "bio_n_ie", "n_reps_in_graph", "n_reps"):
        val = _f(attr.get(key), 0.0)
        if val > 0:
            return val
    return 0.0


def _node_color(attr: dict) -> str:
    for key in ("marker_color", "archetype_color"):
        col = attr.get(key)
        if isinstance(col, str) and col.startswith("#"):
            return col
    return FAMILY_FALLBACK_COLOR


def _node_label(attr: dict, node: str) -> str:
    for key in ("marker_label_short", "archetype_short", "markers_top1_short", "markers_top1_pfam"):
        val = attr.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return str(node)


def _family(attr: dict) -> str:
    for key in ("marker_family_human", "marker_family", "archetype_family"):
        val = attr.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return "other"


def _node_areas(G: nx.Graph) -> dict[str, float]:
    counts = {n: _n_ie(a) for n, a in G.nodes(data=True)}
    lo = min(counts.values()) if counts else 0.0
    hi = max(counts.values()) if counts else 1.0
    span = max(hi - lo, 1e-9)
    return {n: NODE_MIN_AREA + (c - lo) / span * (NODE_MAX_AREA - NODE_MIN_AREA)
            for n, c in counts.items()}


def _separate(pos: dict[str, np.ndarray], areas: dict[str, float],
              steps: int = 300, pad: float = 1.35) -> dict[str, np.ndarray]:
    """Relax overlaps so circles touch at most lightly, as in a Cytoscape layout."""
    nodes = list(pos)
    out = {n: np.asarray(p, dtype=float) for n, p in pos.items()}
    span = max(np.ptp(np.asarray(list(out.values())), axis=0).max(), 1e-9)
    radius = {n: pad * math.sqrt(areas[n] / math.pi) / 260.0 * span for n in nodes}
    for _ in range(steps):
        moved = False
        for i, a in enumerate(nodes):
            for b in nodes[i + 1:]:
                delta = out[b] - out[a]
                dist = float(np.linalg.norm(delta)) + 1e-9
                want = radius[a] + radius[b]
                if dist < want:
                    push = (want - dist) / dist * delta * 0.5
                    out[a] -= push
                    out[b] += push
                    moved = True
        if not moved:
            break
    return out


def layout(G: nx.Graph, seed: int = 11) -> dict[str, np.ndarray]:
    areas = _node_areas(G)
    if G.number_of_edges() == 0:
        return nx.circular_layout(G)
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for u, v, d in G.edges(data=True):
        H.add_edge(u, v, weight=max(0.5, _f(d.get("relatedness_z"), 1.0)))
    pos = nx.spring_layout(H, seed=seed, weight="weight", iterations=600,
                           k=2.6 / math.sqrt(max(1, H.number_of_nodes())))

    isolated = [n for n in G.nodes() if G.degree(n) == 0]
    if isolated:
        linked = np.asarray([pos[n] for n in G.nodes() if n not in isolated])
        lo, hi = linked.min(axis=0), linked.max(axis=0)
        drop = (hi[1] - lo[1]) * 0.22
        for i, n in enumerate(isolated):
            frac = (i + 1) / (len(isolated) + 1)
            pos[n] = np.array([lo[0] + frac * (hi[0] - lo[0]), lo[1] - drop])
    return _separate(pos, areas)


def draw(ax, G: nx.Graph, pos: dict[str, np.ndarray], *, label_size: float = 6.8) -> None:
    areas = _node_areas(G)

    segments, colors, widths, styles = [], [], [], []
    for u, v, d in G.edges(data=True):
        if u not in pos or v not in pos:
            continue
        z = _f(d.get("relatedness_z"), 1.0)
        fold = _f(d.get("fold_enrichment_vs_mixing_null"), 0.0)
        significant = str(d.get("edge_significance", "ns")) != "ns"
        segments.append([pos[u], pos[v]])
        widths.append(0.45 + 3.1 * math.sqrt(min(z, 25.0) / 25.0))
        if not significant:
            colors.append((0.55, 0.55, 0.55, 0.55))
            styles.append((0, (4, 3)))
        elif fold >= STRONG_FOLD:
            colors.append((0.10, 0.10, 0.10, 0.85))
            styles.append("solid")
        else:
            colors.append((0.38, 0.38, 0.38, 0.70))
            styles.append("solid")
    if segments:
        ax.add_collection(LineCollection(segments, colors=colors, linewidths=widths,
                                         linestyles=styles, zorder=1, capstyle="round"))

    xy = np.asarray([pos[n] for n in G.nodes()])
    ax.scatter(xy[:, 0], xy[:, 1],
               s=[areas[n] for n in G.nodes()],
               c=[_node_color(a) for _, a in G.nodes(data=True)],
               linewidths=0, zorder=2)

    for (n, attr), (x, y) in zip(G.nodes(data=True), xy):
        ax.text(x, y, _node_label(attr, n), ha="center", va="center",
                fontsize=label_size, color="#111111", zorder=3, clip_on=False)

    ax.set_axis_off()
    ax.margins(0.10)
    ax.autoscale_view()
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + (hi - lo) * 0.18)  # headroom for the colour legend


def colour_legend(ax, G: nx.Graph, **kwargs) -> None:
    seen: dict[str, str] = {}
    order: dict[str, float] = {}
    for _, attr in G.nodes(data=True):
        fam = _family(attr)
        seen.setdefault(fam, _node_color(attr))
        order[fam] = order.get(fam, 0.0) + _n_ie(attr)
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=7,
                      markerfacecolor=seen[fam], markeredgecolor="none", label=fam)
               for fam in sorted(seen, key=lambda f: -order[f])]
    ax.legend(handles=handles, frameon=False, handletextpad=0.5, borderpad=0.0,
              labelspacing=0.55, fontsize=8, **kwargs)


def render(graphml: Path = CYTOSCAPE_GRAPHML, out_stem: str = OUT_STEM,
           seed: int = 11) -> None:
    rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})
    G = nx.read_graphml(str(graphml))
    if G.number_of_nodes() == 0:
        print(f"[render] {graphml.name} has no nodes — skip")
        return

    fig, ax = plt.subplots(figsize=(10.0, 8.0))
    pos = layout(G, seed=seed)
    draw(ax, G, pos)
    colour_legend(ax, G, loc="upper left", bbox_to_anchor=(-0.02, 1.0))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

    FIG.mkdir(parents=True, exist_ok=True)
    for ext, kw in (("png", {"dpi": 300}), ("svg", {})):
        fig.savefig(FIG / f"{out_stem}.{ext}", bbox_inches="tight", pad_inches=0.15,
                    facecolor="white", **kw)
    plt.close(fig)
    print("[render] →", FIG / f"{out_stem}.png")


if __name__ == "__main__":
    render()
