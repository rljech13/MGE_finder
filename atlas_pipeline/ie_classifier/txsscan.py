"""Map TXSScan (macsyfinder) hits to IE atlas conjugation markers."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set

# TXSScan gene names -> atlas module_markers keys (ie_atlas_labels.yaml)
TXSSCAN_GENE_TO_MARKER: Dict[str, str] = {
    "T4SS_virb4": "virb4",
    "T4SS_t4cp1": "t4cp",
    "T4SS_I_traI": "relaxase",
    "T4SS_I_traY": "relaxase",
    "T4SS_MOBB": "relaxase",
    "T4SS_T_virB2": "mpf",
    "T4SS_T_virB6": "mpf",
    "T4SS_T_virB8": "mpf",
    "T4SS_T_virB9": "mpf",
    "T4SS_T_virB10": "mpf",
    "T4SS_T_virB11": "mpf",
    "T4SS_I_trbB": "mpf",
    "ComM_comB": "com_machinery",
    "ComM_comC": "com_machinery",
    "ComM_comGA": "com_machinery",
    "ComM_comGB": "com_machinery",
    "ComM_comGC": "com_machinery",
    "ComM_comGD": "com_ea",
    "ComM_comGE": "com_machinery",
    "ComM_comGF": "com_machinery",
    "ComM_comGG": "com_machinery",
    "ComM_comK": "com_machinery",
    "ComM_comX": "com_machinery",
    "ComM_comEA": "com_ea",
    "ComM_comEB": "com_ea",
    "ComM_comEC": "com_ea",
    "Tad_rcpA": "conj_tad",
    "Tad_flp": "conj_tad",
    "TadA": "conj_tad",
    "T4aP_pilB": "conj_t4ap",
    "T4aP_pilC": "conj_t4ap",
    "T4aP_pilQ": "conj_t4ap",
    "T4aP_pilA": "conj_t4ap",
}

# Systems treated as conjugative T4SS (exclude T4aP/T4bP/flagellum/T6SS etc.)
CONJ_SYSTEM_MARKERS: tuple[str, ...] = (
    "pT4SSi",
    "pT4SSt",
)
COMPETENCE_SYSTEM_MARKERS: tuple[str, ...] = ("ComM",)
NON_T4SS_SYSTEM_MARKERS: tuple[str, ...] = ("Tad", "T4aP")


def safe_element_id(element_id: str) -> str:
    return re.sub(r"[^\w\.-]+", "_", element_id or "")[:180]


def load_txsscan_hits_by_element(path: Path) -> Dict[str, List[Dict[str, str]]]:
    """Load long-format hits TSV (merge_txsscan_stage9_results --out-hits)."""
    out: Dict[str, List[Dict[str, str]]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            eid = str(row.get("element_id", "")).strip()
            if eid:
                out.setdefault(eid, []).append(dict(row))
    return out


def txsscan_rows_for_element(
    element_id: str,
    hits_lut: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, str]]:
    return hits_lut.get(safe_element_id(element_id), hits_lut.get(element_id, []))


def parse_hit_locus(hit_id: str) -> str:
    m = re.search(r"locus=([^|]+)", hit_id or "")
    return m.group(1).strip() if m else ""


def is_conj_system(model_fqn: str) -> bool:
    fqn = model_fqn or ""
    return any(tok in fqn for tok in CONJ_SYSTEM_MARKERS)


def is_competence_system(model_fqn: str) -> bool:
    fqn = model_fqn or ""
    return any(tok in fqn for tok in COMPETENCE_SYSTEM_MARKERS)


def is_non_t4ss_system(model_fqn: str) -> bool:
    fqn = model_fqn or ""
    return any(tok in fqn for tok in NON_T4SS_SYSTEM_MARKERS)


def system_allowed_for_marker(model_fqn: str, marker_key: str) -> bool:
    if marker_key in ("com_ea", "com_machinery"):
        return is_competence_system(model_fqn)
    if marker_key in ("conj_tad", "conj_t4ap"):
        return is_non_t4ss_system(model_fqn)
    return is_conj_system(model_fqn)


def marker_weight(cfg: Dict[str, Any], marker_key: str) -> float:
    spec = (cfg.get("module_markers") or {}).get(marker_key) or {}
    return float(spec.get("weight", 0.5))


def inject_txsscan_hits(
    genes: List[Dict[str, Any]],
    txsscan_rows: Iterable[Dict[str, str]],
    cfg: Dict[str, Any],
    *,
    conj_systems_only: bool = False,
) -> List[Dict[str, Any]]:
    """Add marker_hits from TXSScan best_solution rows (matched by locus_tag)."""
    if not txsscan_rows:
        return genes

    by_locus: Dict[str, Dict[str, Any]] = {
        str(g.get("locus_tag", "")).strip(): g for g in genes if g.get("locus_tag")
    }
    seen: Set[tuple[str, str]] = set()

    for row in txsscan_rows:
        model = str(row.get("model_fqn", ""))
        gene_name = str(row.get("gene_name", "")).strip()
        marker_key = TXSSCAN_GENE_TO_MARKER.get(gene_name)
        if not marker_key:
            continue
        if conj_systems_only and not is_conj_system(model):
            continue
        if not conj_systems_only and not system_allowed_for_marker(model, marker_key):
            continue
        locus = parse_hit_locus(str(row.get("hit_id", "")))
        if not locus or locus not in by_locus:
            continue
        g = by_locus[locus]
        pair = (locus, marker_key)
        if pair in seen:
            continue
        existing = {h["key"] for h in g.get("marker_hits") or []}
        if marker_key in existing:
            seen.add(pair)
            continue
        g.setdefault("marker_hits", []).append(
            {
                "key": marker_key,
                "weight": marker_weight(cfg, marker_key),
                "via": "txsscan",
            }
        )
        seen.add(pair)
    return genes


def append_txsscan_evidence(
    genes: List[Dict[str, Any]],
    marker_genes: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Extend flat marker_genes evidence with TXSScan-only hits."""
    for g in genes:
        for h in g.get("marker_hits") or []:
            if h.get("via") != "txsscan":
                continue
            marker_genes.append(
                {
                    "orf_id": g.get("locus_tag", ""),
                    "gene_name_or_product": g.get("product", ""),
                    "marker_type": h["key"],
                    "domain_hits": [],
                    "protein_family_id": g.get("family_id", ""),
                    "coordinates": [g.get("start", 0), g.get("end", 0)],
                    "strand": g.get("strand", "+"),
                    "match_via": "txsscan",
                }
            )
    return marker_genes


def summarize_txsscan_rows(rows: Iterable[Dict[str, str]]) -> Dict[str, Any]:
    """Per-IE summary from best_solution rows."""
    systems: Set[str] = set()
    genes_hit: Set[str] = set()
    markers: Set[str] = set()
    best_score = 0.0
    best_wholeness = 0.0

    for row in rows:
        model = str(row.get("model_fqn", ""))
        if not (is_conj_system(model) or is_competence_system(model) or is_non_t4ss_system(model)):
            continue
        systems.add(model.rsplit("/", 1)[-1])
        genes_hit.add(str(row.get("gene_name", "")))
        mk = TXSSCAN_GENE_TO_MARKER.get(str(row.get("gene_name", "")).strip())
        if mk:
            markers.add(mk)
        try:
            best_score = max(best_score, float(row.get("sys_score", 0) or 0))
            best_wholeness = max(best_wholeness, float(row.get("sys_wholeness", 0) or 0))
        except (TypeError, ValueError):
            pass

    return {
        "n_conj_systems": len(systems),
        "conj_systems": ",".join(sorted(systems)),
        "n_txsscan_genes": len(genes_hit),
        "txsscan_relaxase": "relaxase" in markers,
        "txsscan_t4cp": "t4cp" in markers,
        "txsscan_virb4": "virb4" in markers,
        "txsscan_mpf": "mpf" in markers,
        "txsscan_competence": bool({"com_ea", "com_machinery"} & markers),
        "txsscan_non_t4ss": bool({"conj_tad", "conj_t4ap"} & markers),
        "txsscan_sys_score": best_score,
        "txsscan_sys_wholeness": best_wholeness,
    }


__all__ = [
    "COMPETENCE_SYSTEM_MARKERS",
    "CONJ_SYSTEM_MARKERS",
    "NON_T4SS_SYSTEM_MARKERS",
    "TXSSCAN_GENE_TO_MARKER",
    "append_txsscan_evidence",
    "inject_txsscan_hits",
    "is_competence_system",
    "is_conj_system",
    "is_non_t4ss_system",
    "load_txsscan_hits_by_element",
    "parse_hit_locus",
    "safe_element_id",
    "summarize_txsscan_rows",
    "system_allowed_for_marker",
    "txsscan_rows_for_element",
]
