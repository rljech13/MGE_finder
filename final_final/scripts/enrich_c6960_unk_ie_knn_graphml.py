#!/usr/bin/env python3
"""Attach host taxonomy + IE atlas classification overlays to Stage 9 unc kNN graphml.

Input:  artifacts/nb04_cluster_unknown_6960/knn_graph_c6960_unk_enriched.graphml
        (nodes must carry ie_id, sample_id, leiden_cluster from rebuild)
Output: artifacts/nb04_cluster_unknown_6960/knn_graph_c6960_unk_host_ie_atlas.graphml

Mirrors final/clustering_filtered_ie/enrich_filtered_ie_knn_graphml.py for the Stage 9
rebuild cohort. Adds explicit overlay fields for Cytoscape / Gephi:

  leiden_rebuild_cluster   — new Leiden on Stage 9 unc submatrix
  leiden_hd_cluster_atlas  — prior genome-level Leiden from ie_instance_classification
  integrase_id             — locus key
  overlay_coarse_label / overlay_primary_label — atlas labels (duplicate ie_atlas_* for convenience)

Optional: if artifacts/external_detectors_stage9/external_detectors_stage9_hits.tsv exists,
attach external_detector_* fields per IE (benchmark overlay).
"""
from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import REPO
ENRICH_LIB = REPO / "final_final" / "03_clustering" / "lib"
sys.path.insert(0, str(ENRICH_LIB))

from enrich_filtered_ie_knn_graphml import (  # noqa: E402
    attach_atlas_attrs,
    attach_atlas_missing,
    clear_gtdb_fields,
    load_classification,
    load_gtdb_taxonomy_by_sample,
    load_rep_lut,
    _cap,
)

ART = REPO / "final_final/artifacts/nb04_cluster_unknown_6960"
GRAPH_IN = ART / "knn_graph_c6960_unk_enriched.graphml"
GRAPH_OUT = ART / "knn_graph_c6960_unk_host_ie_atlas.graphml"
REP_INDEX = REPO / "results_organized/merged_deduplicated/rep_index.tsv"
GTDB = REPO / "results_organized/misc/gtdb_taxonomy.tsv"
CLASSIF = (
    REPO / "results_organized/genome_clustering/project/notebook7_extended"
    / "ie_instance_classification.tsv"
)
LEIDEN_REBUILD = ART / "leiden_clusters_c6960_unk.tsv"
IE_ROWS = ART / "ie_rep_rows_c6960_unk.tsv"
EXTERNAL_HITS = REPO / "artifacts/external_detectors_stage9/external_detectors_stage9_hits.tsv"

LEIDEN_COLORS = {
    0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c", 3: "#9467bd", 4: "#d62728",
    5: "#8c564b", 6: "#e377c2", 7: "#17becf", 8: "#bcbd22", 9: "#aec7e8",
    10: "#98df8a", 11: "#c49c94",
}


def load_external_lut(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "locus" not in df.columns and "hit_id" in df.columns:
        df["locus"] = df["hit_id"].astype(str).str.split("__").str[0]
    if "locus" not in df.columns:
        return {}
    out: dict[str, dict[str, str]] = {}
    for locus, sub in df.groupby("locus"):
        tools = sorted(sub["tool"].astype(str).unique()) if "tool" in sub.columns else []
        markers = sorted(sub["marker"].astype(str).unique()) if "marker" in sub.columns else []
        out[str(locus)] = {
            "external_detector_tools": "; ".join(tools),
            "external_detector_markers": "; ".join(markers),
            "external_detector_n_hits": str(len(sub)),
            "external_detector_status": "hit",
        }
    return out


def main() -> None:
    if not GRAPH_IN.is_file():
        raise SystemExit(f"Missing input graphml: {GRAPH_IN}")

    G = nx.read_graphml(str(GRAPH_IN))
    rep_lut = load_rep_lut(REP_INDEX)
    tax_lut = load_gtdb_taxonomy_by_sample(GTDB)
    cls_lut = load_classification(CLASSIF)
    ext_lut = load_external_lut(EXTERNAL_HITS)

    # node_id n{i} → integrase_id / ie_id from ie_rep_rows order
    rows = pd.read_csv(IE_ROWS, sep="\t")
    rows["node_id"] = "n" + rows.index.astype(str)
    node_ie = dict(zip(rows["node_id"], rows["ie_id"].astype(str)))
    node_int = dict(zip(rows["node_id"], rows["integrase_id"].astype(str)))

    leid = pd.read_csv(LEIDEN_REBUILD, sep="\t")
    leid["node_id"] = "n" + leid.index.astype(str)
    node_rebuild = dict(zip(leid["node_id"], leid["leiden_cluster"].astype(int)))

    n_miss_host = n_miss_tax = n_miss_cls = 0
    n_ext = 0

    for node, attr in list(G.nodes(data=True)):
        ie_id = str(attr.get("ie_id") or node_ie.get(str(node), "")).strip()
        if not ie_id:
            ie_id = str(node_ie.get(str(node), "")).strip()
        sample_id = str(attr.get("sample_id", "") or "").strip()
        if not sample_id and ie_id in cls_lut:
            sample_id = str(cls_lut[ie_id].get("sample_id", "") or "").strip()

        G.nodes[node]["ie_id"] = ie_id
        G.nodes[node]["integrase_id"] = str(node_int.get(str(node), ie_id.split(":")[0] if ie_id else ""))

        rebuild_cl = node_rebuild.get(str(node), attr.get("leiden_cluster", ""))
        try:
            rebuild_cl_int = int(rebuild_cl)
        except (TypeError, ValueError):
            rebuild_cl_int = -1
        G.nodes[node]["leiden_rebuild_cluster"] = str(rebuild_cl_int)
        G.nodes[node]["leiden_rebuild_color"] = LEIDEN_COLORS.get(rebuild_cl_int, "#7f7f7f")
        G.nodes[node]["cohort"] = "c6960_unk_rebuild"
        G.nodes[node]["graph_level"] = "ie_knn"

        row = rep_lut.get(ie_id)
        if row is None:
            n_miss_host += 1
            G.nodes[node]["host_genus"] = ""
            G.nodes[node]["host_species"] = ""
            G.nodes[node]["host_organism"] = ""
            G.nodes[node]["host_lineage_source"] = ""
        else:
            G.nodes[node]["host_genus"] = _cap(row.get("Genus"))
            G.nodes[node]["host_species"] = _cap(row.get("Species"))
            G.nodes[node]["host_organism"] = _cap(row.get("Organism"))
            G.nodes[node]["host_lineage_source"] = _cap(row.get("source"))

        tax_row = tax_lut.get(sample_id) if sample_id else None
        if tax_row is None:
            n_miss_tax += 1
            clear_gtdb_fields(G, node)
        else:
            G.nodes[node]["host_gtdb_domain"] = _cap(tax_row.get("domain"))
            G.nodes[node]["host_gtdb_phylum"] = _cap(tax_row.get("phylum"))
            G.nodes[node]["host_gtdb_class"] = _cap(tax_row.get("class"))
            G.nodes[node]["host_gtdb_order"] = _cap(tax_row.get("order"))
            G.nodes[node]["host_gtdb_family"] = _cap(tax_row.get("family"))

        clf = cls_lut.get(ie_id)
        if clf is None:
            n_miss_cls += 1
            attach_atlas_missing(G, node)
            G.nodes[node]["overlay_coarse_label"] = ""
            G.nodes[node]["overlay_primary_label"] = ""
            G.nodes[node]["leiden_hd_cluster_atlas"] = ""
        else:
            attach_atlas_attrs(G, node, clf)
            G.nodes[node]["overlay_coarse_label"] = _cap(clf.get("coarse_label"))
            G.nodes[node]["overlay_primary_label"] = _cap(clf.get("primary_label"))
            G.nodes[node]["leiden_hd_cluster_atlas"] = _cap(clf.get("leiden_hd_cluster"))

        int_id = G.nodes[node]["integrase_id"]
        ext = ext_lut.get(int_id) or ext_lut.get(ie_id)
        if ext:
            n_ext += 1
            for k, v in ext.items():
                G.nodes[node][k] = v
        else:
            G.nodes[node]["external_detector_status"] = "no_hit"
            G.nodes[node]["external_detector_tools"] = ""
            G.nodes[node]["external_detector_markers"] = ""
            G.nodes[node]["external_detector_n_hits"] = "0"

    G.graph["description"] = (
        "Stage 9 unclassified IE kNN graph (Jaccard, k=15, Leiden rebuild). "
        "Overlays: host (rep_index + GTDB), ie_atlas_* + overlay_* from ie_instance_classification.tsv, "
        "leiden_rebuild_cluster from nb04_cluster_unknown_6960 rebuild. "
        "external_detector_* populated when artifacts/external_detectors_stage9/external_detectors_stage9_hits.tsv exists."
    )
    G.graph["enrichment_note"] = G.graph["description"]
    nx.write_graphml(G, str(GRAPH_OUT))
    print(f"Wrote {GRAPH_OUT}")
    print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
    print(
        f"missing_rep_index={n_miss_host} missing_gtdb={n_miss_tax} "
        f"missing_classification={n_miss_cls} external_hits={n_ext}"
    )


if __name__ == "__main__":
    main()
