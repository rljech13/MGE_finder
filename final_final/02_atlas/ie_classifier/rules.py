from typing import Any, Dict, List, Optional, Tuple

from . import scores as scorelib

# Granular primary_label → simplified bucket for reporting / downstream plots.
_COARSE_PHAGE_CAUDO = frozenset({"phage.caudoviral.complete"})
_COARSE_PHAGE_INO = frozenset({"phage.inovirus_like"})
_COARSE_PHAGE_GTA = frozenset({
    "phage.GTA_like_or_particle_like",
    "phage.satellite_or_dependent",
    "phage.cryptic_prophage_footprint",
})
_COARSE_PHAGE_PICI = frozenset({"phage.PICI_like"})
_COARSE_MOB = frozenset({
    "mob.integrative_transposon",
    "mob.IS_heavy_microisland",
    "mob.integron_like",
    "mob.group2_intron_casposon",
})
_COARSE_ISLAND = frozenset({
    "island.PAI_like",
    "island.RM_module",
    "island.TA_rich",
})


def coarse_label_from_primary(primary: str) -> str:
    """
    User-facing coarse taxonomy:
      - phage_caudoviral — strict caudoviruses (HK97+TerL).
      - phage_inovirus — inovirus-like.
      - phage_gta — GTA / satellite / cryptic prophage footprints.
      - phage_pici — PICI-like islands.
      - conj — all conj.* merged.
      - mob_element — transposons, IS islands, integrons, group II / casposons.
      - island_module — PAI, RM, TA-rich islands.
      - rep — replication architectures.
      - unclassified — other.unclassified_IE and unknown.
    """
    p = (primary or "").strip()
    if p in _COARSE_PHAGE_CAUDO:
        return "phage_caudoviral"
    if p in _COARSE_PHAGE_INO:
        return "phage_inovirus"
    if p in _COARSE_PHAGE_GTA:
        return "phage_gta"
    if p in _COARSE_PHAGE_PICI:
        return "phage_pici"
    if p.startswith("conj."):
        return "conj"
    if p in _COARSE_MOB:
        return "mob_element"
    if p in _COARSE_ISLAND:
        return "island_module"
    if p.startswith("rep."):
        return "rep"
    if p == "other.unclassified_IE":
        return "unclassified"
    if p.startswith("phage."):
        return "unclassified"
    return "unclassified"


def _boost_ie_core(meta: Dict[str, Any], genes: List[Dict[str, Any]], base_ie: float) -> float:
    ie = base_ie
    if meta.get("tRNA_type_assigned"):
        ie = min(1.0, ie + 0.15)
    int_lookup = (meta.get("integrase_id_lookup") or "").strip()
    if int_lookup:
        for g in genes:
            if g.get("locus_tag") == int_lookup or str(g.get("locus_tag", "")).startswith(int_lookup.split("_")[0]):
                ie = min(1.0, ie + 0.1)
                break
    return ie


def _iec_ok(rules: Dict[str, Any], iec: float) -> bool:
    if "min_ie_core" not in rules:
        return True
    return iec >= float(rules.get("min_ie_core", 0))


def _ie_len_ok(rules: Dict[str, Any], meta: Dict[str, Any], genes: List[Dict[str, Any]]) -> bool:
    ie_len = scorelib.ie_length_nt(meta, genes)
    if "min_ie_len_nt" in rules and ie_len < int(rules["min_ie_len_nt"]):
        return False
    if "max_ie_len_nt" in rules and ie_len > int(rules["max_ie_len_nt"]):
        return False
    return True


def _caudoviral_core_ok(ps: Dict[str, Any], ant: Dict[str, Any], rules: Dict[str, Any]) -> bool:
    need_mcp = bool(rules.get("require_hk97_mcp", ant.get("any_caudoviral_requires_hk97_mcp", True)))
    need_terl = bool(rules.get("require_terl", ant.get("any_caudoviral_requires_terl", True)))
    if need_mcp and not ps.get("has_hk97_mcp"):
        return False
    if need_terl and not ps.get("has_terl"):
        return False
    return True


def _phage_blocked_by_conjugation(scores: Dict[str, float], cfg: Dict[str, Any], cf: Dict[str, Any]) -> bool:
    ant = cfg.get("anti_phage_overcall") or {}
    if not ant.get("forbid_conjugation_markers", False):
        return False
    thr = float(ant.get("conjugation_block_min_score", 0.11))
    if scores.get("conjugation_score", 0) >= thr:
        return True
    return bool(cf.get("relaxase") or cf.get("t4cp") or cf.get("virb4"))


def _not_full_caudo(ps: Dict[str, Any]) -> bool:
    return not (ps.get("has_terl") and ps.get("has_hk97_mcp"))


def _has_integrase(cf: Dict[str, Any], scores: Dict[str, float]) -> bool:
    return bool(cf.get("integrase")) or float(scores.get("_integrase", 0)) >= 0.2


def _passes(
    label: str,
    scores: Dict[str, float],
    genes: List[Dict[str, Any]],
    meta: Dict[str, Any],
    cfg: Dict[str, Any],
    cf: Dict[str, Any],
    cj: Dict[str, Any],
    ant: Dict[str, Any],
    win: int,
) -> Tuple[bool, List[str]]:
    """Return (ok, notes)."""
    lr_all = cfg.get("label_rules") or {}
    rules = lr_all.get(label) or {}
    notes: List[str] = []
    iec = float(scores.get("IE_core_score", 0.0))

    if label == "phage.caudoviral.complete":
        if not rules:
            return False, notes
        if _phage_blocked_by_conjugation(scores, cfg, cf):
            return False, notes
        ps = scorelib.phage_structural_flags(genes, cfg)
        if not _caudoviral_core_ok(ps, ant, rules):
            return False, notes
        if not _ie_len_ok(rules, meta, genes):
            return False, notes
        if scores.get("phage_score", 0) < float(rules.get("min_phage_score", 0.5)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        min_cat = int(ant.get("min_distinct_phage_marker_categories_for_complete", 1))
        if scorelib.phage_category_count(genes, cfg) < min_cat:
            return False, notes
        return True, notes

    if label == "phage.inovirus_like":
        if _phage_blocked_by_conjugation(scores, cfg, cf):
            return False, notes
        if scores.get("inovirus_score", 0) < float(rules.get("min_inovirus_score", 0.5)):
            return False, notes
        if rules.get("require_zot", True) and not scorelib.has_inovirus_zot(genes):
            return False, notes
        if scorelib.distinct_inovirus_marker_count(genes) < int(rules.get("min_distinct_inov_markers", 1)):
            return False, notes
        if not _ie_len_ok(rules, meta, genes):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "phage.PICI_like":
        if _phage_blocked_by_conjugation(scores, cfg, cf):
            return False, notes
        if not _not_full_caudo(cj):
            return False, notes
        if scores.get("pici_score", 0) < float(rules.get("min_pici_score", 0.22)):
            return False, notes
        if rules.get("require_integrase", True) and not _has_integrase(cf, scores):
            return False, notes
        if rules.get("require_phage_reg_or_cpt", True):
            has_cpt = scorelib.has_marker_keys(genes, {"pici_cpt"})
            n_reg = scorelib.count_marker_key_occurrences(genes, "phage_reg")
            if not (has_cpt or n_reg >= int(rules.get("min_phage_reg_hits", 2))):
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "phage.GTA_like_or_particle_like":
        if _phage_blocked_by_conjugation(scores, cfg, cf):
            return False, notes
        ps = scorelib.phage_structural_flags(genes, cfg)
        if ps["has_terl"] and ps["has_hk97_mcp"]:
            return False, notes
        phage_ok = scores.get("phage_score", 0) >= float(rules.get("min_phage_score", 0.14))
        gta_ok = scores.get("gta_score", 0) >= float(rules.get("min_gta_score", 0.12))
        if not (phage_ok or gta_ok):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "phage.satellite_or_dependent":
        if _phage_blocked_by_conjugation(scores, cfg, cf):
            return False, notes
        if not _not_full_caudo(cj):
            return False, notes
        if scores.get("satellite_score", 0) < float(rules.get("min_satellite_score", 0.22)):
            return False, notes
        if scores.get("phage_score", 0) < float(rules.get("min_phage_score", 0.08)):
            return False, notes
        if not _ie_len_ok(rules, meta, genes):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "phage.cryptic_prophage_footprint":
        if _phage_blocked_by_conjugation(scores, cfg, cf):
            return False, notes
        if not _not_full_caudo(cj):
            return False, notes
        if rules.get("require_integrase", True) and not _has_integrase(cf, scores):
            return False, notes
        if scores.get("phage_score", 0) < float(rules.get("min_phage_score", 0.10)):
            return False, notes
        if "max_phage_categories" in rules:
            if scorelib.phage_category_count(genes, cfg) > int(rules["max_phage_categories"]):
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "conj.T4SS_ICE":
        cfj = scorelib.conjugation_flags(genes)
        if not cfj["relaxase"] or not cfj["t4cp"] or not cfj["virb4"]:
            return False, notes
        if not cfj["integrase"] and scores.get("IE_core_score", 0) < 0.15:
            return False, notes
        if scores.get("conjugation_score", 0) < float(rules.get("min_conjugation_score", 0.45)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "conj.IME_mobilizable":
        cfj = scorelib.conjugation_flags(genes)
        if not cfj["relaxase"] or not cfj["integrase"]:
            return False, notes
        thr = float(rules.get("virb4_strong_threshold", 0.35))
        if cfj["virb4"] and scores.get("conjugation_score", 0) >= thr:
            return False, notes
        if scores.get("conjugation_score", 0) < float(rules.get("min_conjugation_score", 0.25)):
            return False, notes
        if rules.get("defer_if_primpol_helicase_synteny", True):
            pr = scorelib.synteny_primpol_helicase(genes, win)
            has_prim = any(
                any(h["key"] == "primpol_aep" for h in (g.get("marker_hits") or [])) for g in genes
            )
            if pr and has_prim:
                rs = float(scores.get("replication_score", 0))
                if rs >= float(rules.get("min_replication_to_defer_ime", 0.25)):
                    return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "conj.partial_ICE_or_CIME":
        cfj = scorelib.conjugation_flags(genes)
        sc = scores.get("conjugation_score", 0)
        n_markers = sum([cfj["relaxase"], cfj["t4cp"], cfj["virb4"]])
        if sc < float(rules.get("min_conjugation_score", 0.15)):
            return False, notes
        if n_markers >= 2 and not (cfj["relaxase"] and cfj["t4cp"] and cfj["virb4"]):
            return True, notes
        if cfj["virb4"] and not cfj["relaxase"]:
            return True, notes
        return False, notes

    if label == "conj.non_T4SS_machinery":
        cfj = scorelib.conjugation_flags(genes)
        if cfj["relaxase"] and cfj["t4cp"] and cfj["virb4"]:
            return False, notes
        if scores.get("conj_non_t4ss_score", 0) < float(rules.get("min_conj_non_t4ss_score", 0.24)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.55)):
            return False, notes
        if rules.get("require_integrase_or_ie", True):
            if not _has_integrase(cf, scores) and iec < 0.15:
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "conj.competence_island":
        if scores.get("competence_score", 0) < float(rules.get("min_competence_score", 0.28)):
            return False, notes
        if rules.get("require_integrase_or_com", True):
            has_com = scorelib.has_marker_keys(genes, {"com_ea", "com_machinery"})
            if not (_has_integrase(cf, scores) or has_com):
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "conj.Thermus_DSS_ICETh1_like":
        from .marker_match import match_gene_against_markers

        has_tdt_text = any("tdta" in (g.get("product") or "").lower() for g in genes)
        mm = cfg.get("module_markers") or {}
        hit_tdt = any(
            any(h[0] == "tdtA_like" for h in match_gene_against_markers(g, mm)) for g in genes
        )
        hit_nur = any(any(h[0] == "nurA" for h in match_gene_against_markers(g, mm)) for g in genes)
        hit_rm = any(any(h[0] == "rm" for h in match_gene_against_markers(g, mm)) for g in genes)
        dss_sc = float(scores.get("AICE_DSS_score", 0))
        min_dss = float(rules.get("min_dss_score", 0.55))
        weak = float(rules.get("min_dss_score_weak_evidence", 0.45))
        has_dss_anchor = hit_tdt or hit_nur or hit_rm or has_tdt_text
        thr = min_dss if not has_dss_anchor else weak
        if dss_sc < thr:
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "mob.integrative_transposon":
        if scores.get("mobilome_score", 0) < float(rules.get("min_mobilome_score", 0.24)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.45)):
            return False, notes
        if rules.get("require_transposon_signature", True):
            has_tn3 = scorelib.has_marker_keys(genes, {"tn3_resolvase"})
            has_tn = scorelib.has_marker_keys(genes, {"tn_transposase"})
            has_int = _has_integrase(cf, scores)
            n_is = scorelib.count_marker_key_occurrences(genes, "is_transposase")
            if n_is >= int(rules.get("defer_if_is_hits", 2)) and not has_tn3:
                return False, notes
            if not (has_tn3 or (has_tn and has_int)):
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "mob.integron_like":
        if scores.get("mobilome_score", 0) < float(rules.get("min_mobilome_score", 0.20)):
            return False, notes
        if rules.get("require_integron_signature", True):
            has_inti = scorelib.has_marker_keys(genes, {"inti_integron"})
            has_attc = scorelib.has_marker_keys(genes, {"attc_integron"})
            has_int = _has_integrase(cf, scores)
            if not (has_inti or (has_attc and has_int)):
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "mob.IS_heavy_microisland":
        if scorelib.has_marker_keys(genes, {"tn3_resolvase"}):
            return False, notes
        if not _ie_len_ok(rules, meta, genes):
            return False, notes
        if scores.get("is_score", 0) < float(rules.get("min_is_score", 0.22)):
            return False, notes
        if scores.get("replication_score", 0) > float(rules.get("max_replication_score", 0.45)):
            return False, notes
        n_is = scorelib.count_marker_key_occurrences(genes, "is_transposase")
        min_hits = int(rules.get("min_is_hits", 2))
        if n_is < min_hits and not (n_is >= 1 and _has_integrase(cf, scores)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "mob.group2_intron_casposon":
        if scores.get("mobilome_score", 0) < float(rules.get("min_mobilome_score", 0.30)):
            return False, notes
        if rules.get("require_group2_or_casposon", True):
            if not scorelib.has_marker_keys(genes, {"group2_rt", "casposon"}):
                return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "island.PAI_like":
        if not _has_integrase(cf, scores):
            return False, notes
        if scores.get("pai_score", 0) < float(rules.get("min_pai_score", 0.22)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.40)):
            return False, notes
        if scores.get("phage_score", 0) > float(rules.get("max_phage_score", 0.55)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "island.RM_module":
        if not _has_integrase(cf, scores):
            return False, notes
        if scores.get("rm_island_score", 0) < float(rules.get("min_rm_island_score", 0.22)):
            return False, notes
        if scores.get("replication_score", 0) > float(rules.get("max_replication_score", 0.40)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.40)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "island.TA_rich":
        if not _has_integrase(cf, scores):
            return False, notes
        if scores.get("ta_island_score", 0) < float(rules.get("min_ta_island_score", 0.28)):
            return False, notes
        if scores.get("replication_score", 0) > float(rules.get("max_replication_score", 0.35)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "rep.AEP_PrimPol_helicase":
        if scores.get("_phage_structural_score", 0) > float(rules.get("max_phage_structural", 0.5)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.45)):
            return False, notes
        pr = scorelib.synteny_primpol_helicase(genes, win)
        if rules.get("require_helicase_near_primpol_window", True) and not pr:
            return False, notes
        has_prim = any(
            any(h["key"] == "primpol_aep" for h in (g.get("marker_hits") or [])) for g in genes
        )
        if not has_prim:
            return False, notes
        rs = float(scores.get("replication_score", 0))
        if rs < float(rules.get("min_replication_pair_score", 0.4)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "rep.PolA_or_PolB_helicase":
        has_pol_ab = any(
            any(h["key"] in ("dna_pol_a", "dna_pol_b") for h in (g.get("marker_hits") or []))
            for g in genes
        )
        if not has_pol_ab:
            return False, notes
        if scores.get("_phage_structural_score", 0) > float(rules.get("max_phage_structural", 0.5)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.45)):
            return False, notes
        rs = float(scores.get("replication_score", 0))
        thr = float(rules.get("min_pol_helicase_score", 0.4))
        if rs < thr:
            return False, notes
        if not scorelib.synteny_primpol_helicase(genes, win) and rs < 0.55:
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "rep.helicase_only":
        if scores.get("_pol_weights_present"):
            return False, notes
        if scores.get("helicase_only_score", 0) < float(rules.get("min_helicase_only_score", 0.35)):
            return False, notes
        if scores.get("_phage_structural_score", 0) > float(rules.get("max_phage_structural", 0.5)):
            return False, notes
        if scores.get("conjugation_score", 0) > float(rules.get("max_conjugation_score", 0.45)):
            return False, notes
        if not _iec_ok(rules, iec):
            return False, notes
        return True, notes

    if label == "other.unclassified_IE":
        return bool(rules.get("fallback", True)), notes

    return False, notes


def assign_primary_label(
    scores: Dict[str, float],
    genes: List[Dict[str, Any]],
    meta: Dict[str, Any],
    cfg: Dict[str, Any],
    att: Dict[str, bool],
) -> Tuple[str, Optional[str], str, List[str]]:
    """primary, subtype, confidence, conflict_messages."""
    scores = dict(scores)
    scores["IE_core_score"] = _boost_ie_core(meta, genes, float(scores.get("IE_core_score", 0)))
    if att.get("attL_attR_present"):
        scores["IE_core_score"] = min(1.0, float(scores["IE_core_score"]) + 0.08)

    priority = list(cfg.get("priority_order") or [])
    ant = cfg.get("anti_phage_overcall") or {}
    win = int(cfg.get("synteny", {}).get("window_genes", 8))

    cf = scorelib.conjugation_flags(genes)
    cj = scorelib.phage_structural_flags(genes, cfg)
    conflicts: List[str] = []

    if cj.get("packaging_only_structural") and scores.get("phage_score", 0) > 0.4:
        conflicts.append("packaging_only_structural: anti_phage_overcall")

    chosen = "other.unclassified_IE"
    for label in priority:
        ok, _notes = _passes(label, scores, genes, meta, cfg, cf, cj, ant, win)
        if not ok:
            continue
        chosen = label
        break

    subtype = None
    if chosen.startswith("phage.") and scores.get("plasmid_score", 0) > 0.35:
        subtype = "phage_plasmid_like"

    conf = _confidence_for(chosen, scores, cfg, cj, ant)
    chosen, conf = apply_packaging_cap(chosen, conf, scores, genes, cfg, cj)
    return chosen, subtype, conf, conflicts


def _confidence_for(label: str, scores: Dict[str, float], cfg: Dict[str, Any], pj: Dict[str, Any], ant: Dict[str, Any]) -> str:
    cb = cfg.get("confidence") or {}
    high = cb.get("high") or {}
    med = cb.get("medium") or {}
    iec = float(scores.get("IE_core_score", 0))
    top_spec = max(
        float(scores.get("phage_score", 0)),
        float(scores.get("inovirus_score", 0)),
        float(scores.get("conjugation_score", 0)),
        float(scores.get("AICE_DSS_score", 0)),
        float(scores.get("replication_score", 0)),
        float(scores.get("pici_score", 0)),
        float(scores.get("mobilome_score", 0)),
        float(scores.get("island_module_score", 0)),
        float(scores.get("gta_score", 0)),
    )
    if iec >= float(high.get("min_ie_core", 0.35)) and top_spec >= float(high.get("min_specialized_module", 0.45)):
        if label == "phage.caudoviral.complete" and not (pj.get("has_terl") and pj.get("has_hk97_mcp")):
            return "medium"
        return "high"
    if iec >= float(med.get("min_ie_core", 0.22)) and top_spec >= float(med.get("min_specialized_module", 0.25)):
        return "medium"
    return "low"


def secondary_flags(scores: Dict[str, float], cfg: Dict[str, Any]) -> List[str]:
    out = []
    trig = cfg.get("secondary_flag_triggers") or {}
    if scores.get("phage_score", 0) >= float(trig.get("has_phage_markers", {}).get("min_phage_module", 0.2)):
        out.append("has_phage_markers")
    if scores.get("conjugation_score", 0) >= float(trig.get("has_conjugation_markers", {}).get("min_conjugation_module", 0.2)):
        out.append("has_conjugation_markers")
    if scores.get("replication_score", 0) >= float(trig.get("has_replication_module", {}).get("min_replication_module", 0.25)):
        out.append("has_replication_module")
    if scores.get("plasmid_score", 0) >= float(trig.get("has_plasmid_maintenance", {}).get("min_plasmid_module", 0.25)):
        out.append("has_plasmid_maintenance")
    if scores.get("mobilome_score", 0) >= 0.22:
        out.append("has_mobilome_markers")
    if scores.get("island_module_score", 0) >= 0.22:
        out.append("has_island_module_markers")
    return out


def apply_packaging_cap(
    primary: str,
    confidence: str,
    scores: Dict[str, float],
    genes: List[Dict[str, Any]],
    cfg: Dict[str, Any],
    pj: Dict[str, Any],
) -> Tuple[str, str]:
    ant = cfg.get("anti_phage_overcall") or {}
    if primary != "phage.caudoviral.complete":
        return primary, confidence
    if pj.get("packaging_only_structural"):
        cap = str(ant.get("packaging_alone_max_primary", "phage.cryptic_prophage_footprint"))
        return cap, str(ant.get("packaging_alone_cap_confidence", "medium"))
    min_cat = int(ant.get("min_distinct_phage_marker_categories_for_complete", 3))
    if scorelib.phage_category_count(genes, cfg) < min_cat:
        if scores.get("gta_score", 0) >= 0.12:
            return "phage.GTA_like_or_particle_like", "medium"
        return "phage.cryptic_prophage_footprint", "medium"
    return primary, confidence
