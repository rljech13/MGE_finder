import math
from typing import Any, Dict, List, Set, Tuple


def _marker_categories(cfg: Dict[str, Any]) -> Dict[str, str]:
    mm = cfg.get("module_markers") or {}
    return {k: str(v.get("category", "")).strip() for k, v in mm.items() if isinstance(v, dict)}


def _aggregate_marker_weights(weights: List[float], n_genes: int, cfg: Dict[str, Any]) -> float:
    """Module score from unique marker weights. Default: sum/scale (no penalty for few CDS)."""
    s = sum(weights) if weights else 0.0
    if s <= 0:
        return 0.0
    agg = cfg.get("score_aggregation") or {}
    mode = str(agg.get("mode", "marker_sum"))
    if mode == "gene_damped":
        damp = math.sqrt(1.0 + max(0, n_genes - 1) * float(agg.get("gene_damp_slope", 0.12)))
        div = damp * float(agg.get("gene_damp_divisor", 2.5))
        return min(1.0, s / div)
    scale = float(agg.get("sum_scale", 1.35))
    cap = float(agg.get("max_score_cap", 1.0))
    return min(cap, s / scale)


def _unique_marker_weights(
    genes: List[Dict[str, Any]],
    mcat: Dict[str, str],
    predicate,
) -> List[float]:
    seen: Set[str] = set()
    wts: List[float] = []
    for g in genes:
        for h in g.get("marker_hits") or []:
            key = h["key"]
            if key in seen:
                continue
            cat = mcat.get(key, "")
            if predicate(cat, key):
                seen.add(key)
                wts.append(float(h["weight"]))
    return wts


def compute_module_scores(genes: List[Dict[str, Any]], cfg: Dict[str, Any]) -> Dict[str, float]:
    mcat = _marker_categories(cfg)
    ng = max(1, len(genes))
    agg = cfg.get("score_aggregation") or {}
    ie_boost = float(agg.get("ie_core_boost", 1.12))

    def uw(pred) -> float:
        return _aggregate_marker_weights(_unique_marker_weights(genes, mcat, pred), ng, cfg)

    ie_integrase = uw(lambda c, k: k == "integrase")
    ie_core = min(1.0, ie_integrase * ie_boost)

    # Structural and packaging markers only; lysis and reg do not enter phage_score.
    phage_score = uw(
        lambda c, k: c in ("phage_packaging", "phage_head", "phage_tail", "phage_connector"),
    )
    phage_structural_score = uw(
        lambda c, k: c in ("phage_head", "phage_tail", "phage_connector"),
    )
    inovirus_score = uw(lambda c, k: c == "inovirus")
    conjugation_score = uw(lambda c, k: c == "conjugation")
    aice_raw = uw(lambda c, k: c == "AICE")
    dss_raw = uw(lambda c, k: c == "DSS")
    rep_raw = uw(lambda c, k: c == "replication")
    hel_only_wt = _unique_marker_weights(
        genes,
        mcat,
        lambda c, k: k == "helicase",
    )
    pol_wt = _unique_marker_weights(
        genes,
        mcat,
        lambda c, k: k in ("primpol_aep", "dna_pol_a", "dna_pol_b"),
    )
    plasmid_score = uw(lambda c, k: c == "plasmid")
    host_score = uw(lambda c, k: c == "host")

    mix_aice = float(agg.get("aice_dss_mix_aice", 1.0))
    mix_dss = float(agg.get("aice_dss_mix_dss", 0.32))
    mix_repabc = float(agg.get("aice_dss_mix_repabc", 0.22))
    AICE_DSS_score = min(1.0, mix_aice * aice_raw + mix_dss * dss_raw + mix_repabc * uw(lambda c, k: k == "repABC"))
    replication_score = rep_raw
    helicase_only_score = _aggregate_marker_weights(hel_only_wt, ng, cfg) if not pol_wt else 0.0
    if pol_wt:
        helicase_only_score = 0.0

    pici_score = max(
        uw(lambda c, k: c == "PICI"),
        uw(lambda c, k: k == "phage_reg"),
    )
    mobilome_score = uw(lambda c, k: c == "mobilome")
    is_score = uw(lambda c, k: c == "IS")
    pai_score = uw(lambda c, k: c == "PAI")
    rm_island_score = uw(lambda c, k: c == "RM_island")
    ta_island_score = max(
        uw(lambda c, k: c == "TA_island"),
        uw(lambda c, k: k == "ta_system"),
    )
    competence_score = uw(lambda c, k: c == "competence")
    conj_non_t4ss_score = uw(lambda c, k: c == "conj_non_t4ss")
    satellite_score = uw(lambda c, k: c == "satellite")
    gta_score = max(satellite_score, uw(lambda c, k: k == "gta_marker"))
    island_module_score = min(
        1.0,
        pai_score + rm_island_score * 0.85 + ta_island_score * 0.75,
    )

    return {
        "IE_core_score": ie_core,
        "phage_score": phage_score,
        "_phage_structural_score": phage_structural_score,
        "inovirus_score": inovirus_score,
        "conjugation_score": conjugation_score,
        "mobilization_score": min(1.0, conjugation_score * 0.9),
        "AICE_DSS_score": AICE_DSS_score,
        "replication_score": replication_score,
        "helicase_only_score": helicase_only_score,
        "plasmid_score": plasmid_score,
        "host_contam_score": host_score,
        "pici_score": pici_score,
        "mobilome_score": mobilome_score,
        "is_score": is_score,
        "pai_score": pai_score,
        "rm_island_score": rm_island_score,
        "ta_island_score": ta_island_score,
        "island_module_score": island_module_score,
        "competence_score": competence_score,
        "conj_non_t4ss_score": conj_non_t4ss_score,
        "satellite_score": satellite_score,
        "gta_score": gta_score,
        "_integrase": ie_integrase,
        "_pol_weights_present": bool(pol_wt),
        "_hel_weights_present": bool(hel_only_wt),
    }


def phage_structural_flags(genes: List[Dict[str, Any]], cfg: Dict[str, Any]) -> Dict[str, Any]:
    keys_hit: Set[str] = set()
    for g in genes:
        for h in g.get("marker_hits") or []:
            keys_hit.add(h["key"])
    terl = "terl" in keys_hit
    hk97_mcp = "hk97_mcp" in keys_hit
    portal = "portal" in keys_hit
    mcp = "mcp" in keys_hit
    tail = "tail" in keys_hit
    connector = "connector" in keys_hit
    packaging_only = terl and not hk97_mcp
    return {
        "has_terl": terl,
        "has_hk97_mcp": hk97_mcp,
        "has_portal_or_mcp": portal or mcp or hk97_mcp,
        "has_tail_or_connector": tail or connector,
        "packaging_only_structural": packaging_only,
    }


def conjugation_flags(genes: List[Dict[str, Any]]) -> Dict[str, Any]:
    keys: Set[str] = set()
    for g in genes:
        for h in g.get("marker_hits") or []:
            keys.add(h["key"])
    return {
        "relaxase": "relaxase" in keys,
        "t4cp": "t4cp" in keys,
        "virb4": "virb4" in keys,
        "integrase": "integrase" in keys,
    }


def synteny_primpol_helicase(genes: List[Dict[str, Any]], window: int) -> bool:
    idx_pol: List[int] = []
    idx_hel: List[int] = []
    for g in genes:
        keys = {h["key"] for h in g.get("marker_hits") or []}
        ix = int(g["idx"])
        if {"primpol_aep", "dna_pol_a", "dna_pol_b"} & keys:
            idx_pol.append(ix)
        if "helicase" in keys:
            idx_hel.append(ix)
    for a in idx_pol:
        for b in idx_hel:
            if abs(a - b) <= window:
                return True
    return False


def synteny_scc_pattern2(genes: List[Dict[str, Any]], window: int) -> bool:
    """Ordered ccpol -> mp_partner -> cch2 within span."""
    order = []
    for g in sorted(genes, key=lambda x: x["idx"]):
        keys = {h["key"] for h in g.get("marker_hits") or []}
        tags: List[str] = []
        if "ccpol" in keys:
            tags.append("P")
        if "mp_partner" in keys:
            tags.append("M")
        if "cch2" in keys:
            tags.append("C")
        for t in tags:
            order.append((int(g["idx"]), t))
    if len(order) < 3:
        return False
    for i in range(len(order)):
        if order[i][1] != "P":
            continue
        for j in range(i + 1, len(order)):
            if order[j][1] != "M" or order[j][0] - order[i][0] > window * 3:
                continue
            for k in range(j + 1, len(order)):
                if order[k][1] == "C" and order[k][0] - order[j][0] <= window * 3:
                    return True
    return False


def synteny_scc_pattern1(genes: List[Dict[str, Any]], window: int) -> bool:
    for g in genes:
        keys = {h["key"] for h in g.get("marker_hits") or []}
        if "cch" not in keys and "helicase" not in keys:
            continue
        ix = int(g["idx"])
        for g2 in genes:
            if abs(int(g2["idx"]) - ix) > window:
                continue
            k2 = {h["key"] for h in g2.get("marker_hits") or []}
            if "ssb" in k2:
                return True
    return False


def annotation_coverage(genes: List[Dict[str, Any]]) -> float:
    if not genes:
        return 0.0
    n = 0
    for g in genes:
        if (g.get("eggNOG_Description") or "").strip() or (g.get("pfam_id") or "").strip():
            n += 1
    return n / len(genes)


def marker_key_hits(genes: List[Dict[str, Any]], keys: Set[str]) -> int:
    got: Set[str] = set()
    for g in genes:
        for h in g.get("marker_hits") or []:
            if h["key"] in keys:
                got.add(h["key"])
    return len(got)


def has_marker_keys(genes: List[Dict[str, Any]], keys: Set[str]) -> bool:
    return marker_key_hits(genes, keys) > 0


def count_marker_key_occurrences(genes: List[Dict[str, Any]], key: str) -> int:
    n = 0
    for g in genes:
        for h in g.get("marker_hits") or []:
            if h["key"] == key:
                n += 1
    return n


def runner_up_modules(scores: Dict[str, float], k: int, exclude: Tuple[str, ...] = ()) -> List[Tuple[str, float]]:
    keys = [
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
    ]
    items = [(kk, scores.get(kk, 0.0)) for kk in keys if kk not in exclude]
    items.sort(key=lambda x: -x[1])
    return items[:k]


def distinct_inovirus_marker_count(genes: List[Dict[str, Any]]) -> int:
    want = {"inov_zot", "inov_rcr", "inov_rep", "inov_ssb", "inov_coat", "inov_assembly"}
    got: Set[str] = set()
    for g in genes:
        for h in g.get("marker_hits") or []:
            if h["key"] in want:
                got.add(h["key"])
    return len(got)


def has_inovirus_zot(genes: List[Dict[str, Any]]) -> bool:
    for g in genes:
        for h in g.get("marker_hits") or []:
            if h["key"] == "inov_zot":
                return True
    return False


def ie_length_nt(meta: Dict[str, Any], genes: List[Dict[str, Any]]) -> int:
    raw = meta.get("ie_len_nt")
    if raw not in (None, ""):
        try:
            return int(float(raw))
        except (TypeError, ValueError):
            pass
    if not genes:
        return 0
    return max(int(g.get("end", 0) or 0) for g in genes) - min(int(g.get("start", 0) or 0) for g in genes) + 1


def phage_category_count(genes: List[Dict[str, Any]], cfg: Dict[str, Any]) -> int:
    mcat = _marker_categories(cfg)
    cats: Set[str] = set()
    for g in genes:
        for h in g.get("marker_hits") or []:
            c = mcat.get(h["key"], "")
            if c and c.startswith("phage"):
                cats.add(c)
    return len(cats)
