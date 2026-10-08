from pathlib import Path
import json
from typing import Any, Dict, List, Optional

import pandas as pd

from .loaders import (
    find_att_notes,
    hydrate_genes,
    merge_lookup_for_genome,
    parse_gbff_cds,
    rep_filename_for_lookup,
    truth_cluster_key,
    truth_row_meta,
)
from .marker_match import annotate_genes
from .external_detectors import (
    append_external_evidence,
    external_rows_for_element,
    inject_external_hits,
)
from .txsscan import append_txsscan_evidence, inject_txsscan_hits, txsscan_rows_for_element
from .novel import assess_novel, module_conflict_note
from .rules import assign_primary_label, coarse_label_from_primary, secondary_flags as rules_secondary_flags
from .schema import build_alternate_hypotheses, evidence_blob, public_module_scores
from .scores import annotation_coverage, compute_module_scores, runner_up_modules


_SCORE_KEYS = (
    "IE_core_score",
    "phage_score",
    "inovirus_score",
    "conjugation_score",
    "AICE_DSS_score",
    "replication_score",
    "helicase_only_score",
    "plasmid_score",
    "host_contam_score",
    "pici_score",
    "mobilome_score",
    "island_module_score",
    "competence_score",
    "gta_score",
    "satellite_score",
    "is_score",
)


def classify_from_genes(
    genes: List[Dict[str, Any]],
    meta: Dict[str, Any],
    cfg: Dict[str, Any],
    att: Optional[Dict[str, bool]] = None,
    rep_gbff: str = "",
    txsscan_hits_lut: Optional[Dict[str, List[Dict[str, str]]]] = None,
    external_hits_lut: Optional[Dict[str, List[Dict[str, str]]]] = None,
) -> Dict[str, Any]:
    if att is None:
        att = {}
    genes_c = [dict(g) for g in genes]
    genes_sorted = sorted(genes_c, key=lambda x: (int(x.get("start", 0)), int(x.get("end", 0))))
    for j, g in enumerate(genes_sorted):
        g["idx"] = j

    genes_anno, marker_genes = annotate_genes(genes_sorted, cfg)
    eid = str(meta.get("element_id", ""))
    if txsscan_hits_lut:
        tx_rows = txsscan_rows_for_element(eid, txsscan_hits_lut)
        if tx_rows:
            genes_anno = inject_txsscan_hits(genes_anno, tx_rows, cfg)
            marker_genes = append_txsscan_evidence(genes_anno, marker_genes)
    if external_hits_lut:
        ext_rows = external_rows_for_element(eid, external_hits_lut)
        if ext_rows:
            genes_anno = inject_external_hits(genes_anno, ext_rows, cfg)
            marker_genes = append_external_evidence(genes_anno, marker_genes)
    scores = compute_module_scores(genes_anno, cfg)
    primary, subtype, confidence, conflicts = assign_primary_label(scores, genes_anno, meta, cfg, att)
    coarse = coarse_label_from_primary(primary)
    sec = rules_secondary_flags(scores, cfg)
    ann_cov = annotation_coverage(genes_anno)
    conflicts = list(conflicts)
    mc_note = module_conflict_note(scores, cfg)
    if mc_note:
        conflicts.append(mc_note)

    novel_ok, novel_tier, novel_rationale = assess_novel(primary, scores, cfg, conflicts, ann_cov, genes_anno)

    pub = public_module_scores(scores)
    ev = evidence_blob(marker_genes, scores, cfg, primary, conflicts)
    alt = build_alternate_hypotheses(scores, cfg, primary)
    ru_k = int(cfg.get("signal_balance", {}).get("runner_up_k", 3))
    ru = runner_up_modules(scores, k=ru_k)

    rec: Dict[str, Any] = {
        "element_id": meta.get("element_id", ""),
        "genome_id": meta.get("genome_id", ""),
        "sample_id": meta.get("sample_id", ""),
        "leiden_hd_cluster": meta.get("leiden_hd_cluster", ""),
        "tRNA_type_assigned": meta.get("tRNA_type_assigned", ""),
        "finder_data_missing": meta.get("finder_data_missing", ""),
        "primary_label": primary,
        "coarse_label": coarse,
        "subtype": subtype or "",
        "confidence": confidence,
        "secondary_flags": sec,
        "nonstandard_IE_candidate": novel_ok,
        "novel_tier": novel_tier,
        "novel_rationale": novel_rationale,
        "annotation_coverage": round(ann_cov, 4),
        "runner_up_modules": ru,
        "alternate_label_hypotheses": alt,
        "conflict_notes": conflicts,
        "n_cds": len(genes_anno),
        "rep_gbff": rep_gbff,
        "error": "",
        "evidence_json": json.dumps(ev, ensure_ascii=False),
    }
    for k in _SCORE_KEYS:
        rec[k] = round(float(pub.get(k, 0.0)), 6)
    return rec


def classify_truth_row(
    truth_row: pd.Series,
    cfg: Dict[str, Any],
    merge_df: pd.DataFrame,
    rep_file_by_locus: Dict[str, str],
    gbff_dir: Path,
    txsscan_hits_lut: Optional[Dict[str, List[Dict[str, str]]]] = None,
    external_hits_lut: Optional[Dict[str, List[Dict[str, str]]]] = None,
) -> Dict[str, Any]:
    meta = truth_row_meta(truth_row)
    cid = truth_cluster_key(truth_row)
    fname = rep_filename_for_lookup(rep_file_by_locus, cid)
    if not fname:
        return _error_record(meta, "no_rep_index_match_for_cluster", cfg)
    gbff_path = (gbff_dir / fname).resolve()
    if not gbff_path.exists():
        out = _error_record(meta, f"missing_gbff:{gbff_path.name}", cfg)
        out["rep_gbff"] = str(gbff_path.name)
        return out

    genome_id = meta.get("genome_id") or ""
    lut = merge_lookup_for_genome(merge_df, genome_id) if genome_id else {}
    cds = parse_gbff_cds(gbff_path)
    genes = hydrate_genes(cds, lut)
    att = find_att_notes(gbff_path)
    return classify_from_genes(
        genes,
        meta,
        cfg,
        att=att,
        rep_gbff=str(gbff_path.name),
        txsscan_hits_lut=txsscan_hits_lut,
        external_hits_lut=external_hits_lut,
    )


def _error_record(meta: Dict[str, Any], err: str, cfg: Dict[str, Any]) -> Dict[str, Any]:
    _ = cfg
    rec: Dict[str, Any] = {
        "element_id": meta.get("element_id", ""),
        "genome_id": meta.get("genome_id", ""),
        "sample_id": meta.get("sample_id", ""),
        "leiden_hd_cluster": meta.get("leiden_hd_cluster", ""),
        "tRNA_type_assigned": meta.get("tRNA_type_assigned", ""),
        "finder_data_missing": meta.get("finder_data_missing", ""),
        "primary_label": "other.unclassified_IE",
        "coarse_label": "unclassified",
        "subtype": "",
        "confidence": "low",
        "secondary_flags": [],
        "nonstandard_IE_candidate": False,
        "novel_tier": "none",
        "novel_rationale": "",
        "annotation_coverage": 0.0,
        "runner_up_modules": [],
        "alternate_label_hypotheses": [],
        "conflict_notes": [err],
        "n_cds": 0,
        "rep_gbff": "",
        "error": err,
        "evidence_json": json.dumps({"error": err, "marker_genes": [], "module_scores": {}}, ensure_ascii=False),
    }
    for k in _SCORE_KEYS:
        rec[k] = 0.0
    return rec


__all__ = ["classify_from_genes", "classify_truth_row"]
