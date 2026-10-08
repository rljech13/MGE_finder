import re
from typing import Any, Dict, List, Set, Tuple


def _collect_pfam_tokens(gene: Dict[str, Any]) -> Set[str]:
    out: Set[str] = set()
    raw = (gene.get("pfam_id") or "").strip()
    for part in re.split(r"[;,]", raw):
        part = part.strip()
        m = re.search(r"(PF\d{5}(?:\.\d+)?)", part, re.I)
        if m:
            base = m.group(1).split(".")[0].upper()
            out.add(base)
    emp = (gene.get("emapper_PFAMs") or "").upper()
    for m in re.finditer(r"PF\d{5}", emp):
        out.add(m.group(0))
    return out


def _text_blob(gene: Dict[str, Any]) -> str:
    parts = [
        gene.get("product", ""),
        gene.get("eggNOG_Description", ""),
        gene.get("seed_eggNOG_ortholog", ""),
        gene.get("KEGG_ko", ""),
        gene.get("emapper_PFAMs", ""),
    ]
    return " \n ".join(p for p in parts if p).lower()


def match_gene_against_markers(
    gene: Dict[str, Any],
    module_markers: Dict[str, Any],
) -> List[Tuple[str, float, str]]:
    """Return list of (marker_key, weight, match_type) for hits."""
    hits: List[Tuple[str, float, str]] = []
    pfs = _collect_pfam_tokens(gene)
    text = _text_blob(gene)
    for mk, spec in module_markers.items():
        if not isinstance(spec, dict):
            continue
        weight = float(spec.get("weight", 1.0))
        matched = False
        mt = ""
        for pid in spec.get("pfam") or []:
            p = str(pid).strip().split(".")[0].upper()
            if p in pfs or any(x.startswith(p) for x in pfs):
                matched = True
                mt = "pfam"
                break
        if not matched:
            for rx in spec.get("regex") or []:
                try:
                    if re.search(rx, text, re.I):
                        matched = True
                        mt = "regex"
                        break
                except re.error:
                    continue
        if matched:
            hits.append((mk, weight, mt))
    return hits


def annotate_genes(
    genes: List[Dict[str, Any]],
    cfg: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Attach 'marker_hits' to each gene; return flat evidence list for marker_genes."""
    mm = cfg.get("module_markers") or {}
    flat: List[Dict[str, Any]] = []
    for g in genes:
        hits = match_gene_against_markers(g, mm)
        g["marker_hits"] = [{"key": h[0], "weight": h[1], "via": h[2]} for h in hits]
        for h in hits:
            flat.append(
                {
                    "orf_id": g["locus_tag"],
                    "gene_name_or_product": g.get("product", ""),
                    "marker_type": h[0],
                    "domain_hits": sorted(_collect_pfam_tokens(g)),
                    "protein_family_id": g.get("family_id", ""),
                    "coordinates": [g["start"], g["end"]],
                    "strand": g.get("strand", "+"),
                    "match_via": h[2],
                }
            )
    return genes, flat
