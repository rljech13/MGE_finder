"""Inject external detector hits (ISEScan, IntegronFinder, DefenseFinder) into IE atlas markers."""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set

def safe_element_id(element_id: str) -> str:
    return re.sub(r"[^\w\.-]+", "_", element_id or "")[:180]


def parse_locus_from_hit_id(hit_id: str) -> str:
    m = re.search(r"locus=([^|]+)", hit_id or "")
    return m.group(1).strip() if m else ""


def gene_overlaps(gene: dict, start: int, end: int) -> bool:
    gs = int(gene.get("start", 0) or 0)
    ge = int(gene.get("end", 0) or 0)
    if gs > ge:
        gs, ge = ge, gs
    lo = min(start, end)
    hi = max(start, end)
    return not (ge < lo or gs > hi)


def load_external_hits_by_element(path: Path) -> Dict[str, List[Dict[str, str]]]:
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


def external_rows_for_element(
    element_id: str,
    hits_lut: Dict[str, List[Dict[str, str]]],
) -> List[Dict[str, str]]:
    return hits_lut.get(safe_element_id(element_id), hits_lut.get(element_id, []))


def marker_weight(cfg: Dict[str, Any], marker_key: str) -> float:
    spec = (cfg.get("module_markers") or {}).get(marker_key) or {}
    return float(spec.get("weight", 0.5))


_ELEMENT_MARKERS = frozenset({"genomad_provirus", "virsorter_viral"})


def _element_level_targets(
    genes: List[Dict[str, Any]],
    marker_key: str,
) -> List[Dict[str, Any]]:
    """Attach sequence-level phage calls to integrase or phage-marked genes."""
    int_genes = [
        g
        for g in genes
        if any(h["key"] == "integrase" for h in (g.get("marker_hits") or []))
    ]
    if int_genes:
        return int_genes[:1]
    phage_genes = []
    for g in genes:
        for h in g.get("marker_hits") or []:
            if str(h.get("key", "")).startswith(("ter", "portal", "mcp", "hk97", "phage", "tail")):
                phage_genes.append(g)
                break
    if phage_genes:
        return phage_genes[:2]
    return genes[:1] if genes else []


def _resolve_locus(
    row: Dict[str, str],
    genes: List[Dict[str, Any]],
    by_locus: Dict[str, Dict[str, Any]],
) -> str | None:
    locus = str(row.get("locus_tag", "")).strip()
    if locus == "__element__":
        return "__element__"
    if locus and locus in by_locus:
        return locus
    locus = parse_locus_from_hit_id(str(row.get("gene_name", "")))
    if locus and locus in by_locus:
        return locus
    ob = row.get("orf_begin", "")
    oe = row.get("orf_end", "")
    if ob and oe:
        try:
            start = int(float(ob))
            end = int(float(oe))
        except (TypeError, ValueError):
            return None
        for g in genes:
            if gene_overlaps(g, start, end):
                return str(g.get("locus_tag", "")).strip() or None
    return None


def inject_external_hits(
    genes: List[Dict[str, Any]],
    external_rows: Iterable[Dict[str, str]],
    cfg: Dict[str, Any],
) -> List[Dict[str, Any]]:
    if not external_rows:
        return genes
    by_locus: Dict[str, Dict[str, Any]] = {
        str(g.get("locus_tag", "")).strip(): g for g in genes if g.get("locus_tag")
    }
    seen: Set[tuple[str, str, str]] = set()
    for row in external_rows:
        marker_key = str(row.get("marker_key", "")).strip()
        if not marker_key:
            continue
        source = str(row.get("source", "")).strip()
        locus = _resolve_locus(row, genes, by_locus)
        if not locus:
            continue
        targets: List[Dict[str, Any]]
        if locus == "__element__" and marker_key in _ELEMENT_MARKERS:
            targets = _element_level_targets(genes, marker_key)
        elif locus in by_locus:
            targets = [by_locus[locus]]
        else:
            continue
        for g in targets:
            lt = str(g.get("locus_tag", "")).strip()
            trip = (lt, marker_key, source)
            if trip in seen:
                continue
            existing = {h["key"] for h in g.get("marker_hits") or []}
            if marker_key in existing:
                seen.add(trip)
                continue
            g.setdefault("marker_hits", []).append(
                {
                    "key": marker_key,
                    "weight": marker_weight(cfg, marker_key),
                    "via": source or "external",
                }
            )
            seen.add(trip)
    return genes


def append_external_evidence(
    genes: List[Dict[str, Any]],
    marker_genes: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    for g in genes:
        for h in g.get("marker_hits") or []:
            via = h.get("via", "")
            if via not in ("isescan", "integron_finder", "defense_finder", "genomad", "virsorter"):
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
                    "match_via": via,
                }
            )
    return marker_genes


def summarize_external_rows(rows: Iterable[Dict[str, str]]) -> Dict[str, Any]:
    sources: Set[str] = set()
    markers: Set[str] = set()
    n = 0
    for row in rows:
        sources.add(str(row.get("source", "")))
        mk = str(row.get("marker_key", ""))
        if mk:
            markers.add(mk)
        n += 1
    return {
        "n_external_hits": n,
        "external_sources": ",".join(sorted(s for s in sources if s)),
        "external_markers": ",".join(sorted(markers)),
        "external_is": "is_transposase" in markers,
        "external_integron": bool({"inti_integron", "attc_integron"} & markers),
        "external_rm": "rm_island" in markers,
        "external_ta": "ta_rich" in markers,
    }


__all__ = [
    "append_external_evidence",
    "external_rows_for_element",
    "inject_external_hits",
    "load_external_hits_by_element",
    "summarize_external_rows",
]
