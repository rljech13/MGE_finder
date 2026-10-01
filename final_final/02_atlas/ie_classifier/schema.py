import json
from typing import Any, Dict, List, Tuple

from . import scores as scorelib

TSV_COLUMNS: List[str] = [
    "element_id",
    "genome_id",
    "sample_id",
    "leiden_hd_cluster",
    "tRNA_type_assigned",
    "finder_data_missing",
    "primary_label",
    "coarse_label",
    "subtype",
    "confidence",
    "secondary_flags",
    "nonstandard_IE_candidate",
    "novel_tier",
    "novel_rationale",
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
    "annotation_coverage",
    "runner_up_modules",
    "alternate_label_hypotheses",
    "conflict_notes",
    "n_cds",
    "rep_gbff",
    "error",
]


def public_module_scores(scores: Dict[str, float]) -> Dict[str, float]:
    return {k: float(v) for k, v in scores.items() if not k.startswith("_")}


def build_alternate_hypotheses(
    scores: Dict[str, float],
    cfg: Dict[str, Any],
    exclude_primary: str,
) -> List[Dict[str, Any]]:
    k = int(cfg.get("signal_balance", {}).get("runner_up_k", 3))
    ru = scorelib.runner_up_modules(scores, k=k * 2)
    out: List[Dict[str, Any]] = []
    label_hints: List[Tuple[str, str]] = [
        ("phage_score", "phage.cryptic_prophage_footprint"),
        ("gta_score", "phage.GTA_like_or_particle_like"),
        ("inovirus_score", "phage.inovirus_like"),
        ("conjugation_score", "conj.partial_ICE_or_CIME"),
        ("competence_score", "conj.competence_island"),
        ("pici_score", "phage.PICI_like"),
        ("mobilome_score", "mob.integrative_transposon"),
        ("island_module_score", "island.PAI_like"),
        ("AICE_DSS_score", "conj.Thermus_DSS_ICETh1_like"),
        ("replication_score", "rep.AEP_PrimPol_helicase"),
        ("plasmid_score", "island.TA_rich"),
    ]
    for mod, val in ru:
        if val < 0.08:
            continue
        hyp = {"runner_up_module": mod, "score": round(val, 4)}
        for mk, lab in label_hints:
            if mk == mod:
                hyp["suggested_label_if_primary_by_module"] = lab
                break
        out.append(hyp)
        if len(out) >= k:
            break
    return out


def evidence_blob(
    marker_genes: List[Dict[str, Any]],
    scores: Dict[str, float],
    cfg: Dict[str, Any],
    primary: str,
    conflicts: List[str],
) -> Dict[str, Any]:
    keep = bool(cfg.get("signal_balance", {}).get("keep_non_primary_marker_genes", True))
    mg = marker_genes if keep else marker_genes
    return {
        "marker_genes": mg,
        "module_scores": public_module_scores(scores),
        "runner_up_modules": scorelib.runner_up_modules(scores, k=int(cfg.get("signal_balance", {}).get("runner_up_k", 3))),
        "alternate_label_hypotheses": build_alternate_hypotheses(scores, cfg, primary),
        "conflicts": list(conflicts),
    }


def flatten_record(rec: Dict[str, Any]) -> Dict[str, str]:
    """Single TSV row; complex fields JSON-serialized."""
    flat = {c: "" for c in TSV_COLUMNS}
    for k, v in rec.items():
        if k not in TSV_COLUMNS and k != "evidence_json":
            continue
        if k == "evidence_json":
            continue
        if isinstance(v, (list, dict)):
            flat[k] = json.dumps(v, ensure_ascii=False)
        elif v is None:
            flat[k] = ""
        else:
            flat[k] = str(v)
    return flat


def jsonl_line(rec: Dict[str, Any]) -> str:
    out: Dict[str, Any] = {}
    for k, v in rec.items():
        if k == "evidence_json":
            continue
        out[k] = v
    evj = rec.get("evidence_json")
    if evj:
        try:
            out["evidence"] = json.loads(str(evj))
        except json.JSONDecodeError:
            out["evidence"] = {}
    return json.dumps(out, ensure_ascii=False)
