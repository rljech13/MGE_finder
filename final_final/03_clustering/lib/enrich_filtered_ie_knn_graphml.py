#!/usr/bin/env python3
"""Attach GTDB-derived host taxonomy (rep_index) and IE atlas labels (NB7) to filtered_ie_knn nodes.

Inputs:
  - filtered_ie_knn.graphml — nodes carry sample_id, ie_id (element_id).
  - results_organized/merged_deduplicated/rep_index.tsv — Genus-level fields per IE cluster key.
  - results_organized/misc/gtdb_taxonomy.tsv — GTDB lineage (domain…family, including **phylum**) per sample_id.
  - notebook7_extended/ie_instance_classification.tsv — IE atlas classifier (NB07).

Writes a new graphml (does not overwrite the kNN-only file).

Usage:
  conda run -n MGE_finder python enrich_filtered_ie_knn_graphml.py
  conda run -n MGE_finder python enrich_filtered_ie_knn_graphml.py --input path.graphml --out other.graphml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx
import pandas as pd

GRAPH_DEFAULT = Path(__file__).resolve().parent / "filtered_ie_knn.graphml"
REP_INDEX_DEFAULT = Path(__file__).resolve().parents[2] / "results_organized/merged_deduplicated/rep_index.tsv"
CLASSIF_DEFAULT = (
    Path(__file__).resolve().parents[2]
    / "results_organized/genome_clustering/project/notebook7_extended/ie_instance_classification.tsv"
)
GTDB_TAX_DEFAULT = Path(__file__).resolve().parents[2] / "results_organized/misc/gtdb_taxonomy.tsv"
OUT_DEFAULT = Path(__file__).resolve().parent / "filtered_ie_knn_host_ie_atlas.graphml"

TEXT_CAP = 4000


def normalize_cluster(cid: str) -> str:
    s = str(cid).strip()
    return s.split("|sample=", 1)[0].strip() if "|sample=" in s else s


def _cap(val: object) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if pd.isna(s) or s == "nan":
        return ""
    if len(s) > TEXT_CAP:
        return s[: TEXT_CAP - 1] + "…"
    return s


TEXT_ATLAS_KEYS = [
    "primary_label",
    "coarse_label",
    "subtype",
    "confidence",
    "finder_data_missing",
    "secondary_flags",
    "nonstandard_IE_candidate",
    "novel_tier",
    "novel_rationale",
    "annotation_coverage",
    "runner_up_modules",
    "alternate_label_hypotheses",
    "conflict_notes",
    "n_cds",
    "rep_gbff",
    "error",
    "genome_id",
]

FLOAT_ATLAS_KEYS = [
    "IE_core_score",
    "phage_score",
    "inovirus_score",
    "conjugation_score",
    "AICE_DSS_score",
    "replication_score",
    "helicase_only_score",
    "plasmid_score",
    "host_contam_score",
]


def _series_txt(row: pd.Series, col: str) -> str:
    if col not in row.index:
        return ""
    return _cap(row[col])


def _series_float_str(row: pd.Series, col: str) -> str:
    """Scores as canonical decimal strings — GraphML must not mix floats and \"\" on the same attr."""
    if col not in row.index:
        return ""
    v = row[col]
    if pd.isna(v):
        return ""
    s = str(v).strip()
    if s == "" or s.lower() == "nan":
        return ""
    try:
        return str(float(s))
    except (TypeError, ValueError):
        return _cap(v)


def attach_atlas_attrs(G: "nx.Graph", node: str, row: pd.Series) -> None:
    G.nodes[node]["ie_atlas_status"] = "ok"
    for col in TEXT_ATLAS_KEYS:
        G.nodes[node][f"ie_atlas_{col}"] = _series_txt(row, col)
    G.nodes[node]["ie_atlas_classifier_sample_id"] = _series_txt(row, "sample_id")
    G.nodes[node]["ie_atlas_classifier_leiden_hd_cluster"] = _series_txt(row, "leiden_hd_cluster")
    G.nodes[node]["ie_atlas_trna_type_assigned"] = _series_txt(row, "tRNA_type_assigned")

    for col in FLOAT_ATLAS_KEYS:
        G.nodes[node][f"ie_atlas_{col}"] = _series_float_str(row, col)


def attach_atlas_missing(G: "nx.Graph", node: str) -> None:
    G.nodes[node]["ie_atlas_status"] = "missing_classification_row"
    for col in TEXT_ATLAS_KEYS:
        G.nodes[node][f"ie_atlas_{col}"] = ""
    G.nodes[node]["ie_atlas_classifier_sample_id"] = ""
    G.nodes[node]["ie_atlas_classifier_leiden_hd_cluster"] = ""
    G.nodes[node]["ie_atlas_trna_type_assigned"] = ""
    for col in FLOAT_ATLAS_KEYS:
        G.nodes[node][f"ie_atlas_{col}"] = ""


def load_rep_lut(path: Path) -> dict[str, pd.Series]:
    rep = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    out: dict[str, pd.Series] = {}
    for _, row in rep.iterrows():
        k = normalize_cluster(str(row.get("cluster_id", "") or ""))
        if k:
            out[k] = row
    return out


def load_classification(path: Path) -> dict[str, pd.Series]:
    cls = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    if "element_id" not in cls.columns:
        raise SystemExit(f"Missing element_id column in {path}")
    lut: dict[str, pd.Series] = {}
    for _, row in cls.iterrows():
        eid = str(row.get("element_id", "") or "").strip()
        if eid:
            lut[eid] = row
    return lut


def load_gtdb_taxonomy_by_sample(path: Path) -> dict[str, pd.Series]:
    tax = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    if "sample_id" not in tax.columns:
        raise SystemExit(f"Missing sample_id column in {path}")
    lut: dict[str, pd.Series] = {}
    for _, row in tax.iterrows():
        sid = str(row.get("sample_id", "") or "").strip()
        if sid:
            lut[sid] = row
    return lut


def clear_gtdb_fields(G: "nx.Graph", node: str) -> None:
    for k in ("domain", "phylum", "class", "order", "family"):
        G.nodes[node][f"host_gtdb_{k}"] = ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=GRAPH_DEFAULT)
    ap.add_argument("--rep-index", type=Path, default=REP_INDEX_DEFAULT)
    ap.add_argument("--gtdb-taxonomy", type=Path, default=GTDB_TAX_DEFAULT)
    ap.add_argument("--classification", type=Path, default=CLASSIF_DEFAULT)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()

    if not args.input.is_file():
        raise SystemExit(f"Missing graphml: {args.input}")

    if not args.gtdb_taxonomy.is_file():
        raise SystemExit(f"Missing GTDB taxonomy TSV: {args.gtdb_taxonomy}")

    G = nx.read_graphml(args.input)
    rep_lut = load_rep_lut(args.rep_index)
    tax_lut = load_gtdb_taxonomy_by_sample(args.gtdb_taxonomy)
    cls_lut = load_classification(args.classification)

    n_miss_host = 0
    n_miss_tax = 0
    n_miss_cls = 0

    for node, attr in list(G.nodes(data=True)):
        ie_id = str(attr.get("ie_id", "") or "").strip()
        sample_id = str(attr.get("sample_id", "") or "").strip()

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
        else:
            attach_atlas_attrs(G, node, clf)

    G.graph["enrichment_note"] = (
        "host Genus/Organism from rep_index.tsv; host_gtdb_* from misc/gtdb_taxonomy.tsv (GTDB); "
        "ie_atlas_* from ie_instance_classification.tsv (NB07)."
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(G, str(args.out))
    print(f"Wrote {args.out}")
    print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}")
    print(
        f"missing_rep_index_rows={n_miss_host} "
        f"missing_gtdb_taxonomy_sample={n_miss_tax} "
        f"missing_classification_rows={n_miss_cls}"
    )


if __name__ == "__main__":
    main()
