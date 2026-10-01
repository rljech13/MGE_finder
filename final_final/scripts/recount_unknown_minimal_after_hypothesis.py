#!/usr/bin/env python3
"""Recount composition-minimal unknown IEs after leftover + hypothesis.

Does NOT overwrite the locked 436 tables
(unknown_minimal_ie.tsv / ime_class_funnel.tsv). Writes
artifacts/unknown_ime_classes/after_hypothesis/.

Locked meaning of "minimal" (workspace):
  integrase + exactly one funnel class
  + accessories of that class
  + optional replicase.

Composition-minimal (436, hyp allowed) remains the search pool word.
This recount overlays stack / leftover-strict / hypothesis / CRISPR-array
on previously dark CDS.

  - HTH / C / TA partner / RM piece / C-box / C-candidate → accessory
  - transposase / phage structural / terminase / second cargo → leftover
  - Pfam E=1 / COG S / foldseek ID-only / bare IPR accession → still dark
  - hypothesis is not a gene name (does not resolve darkness)
  - array / ncRNA / promoter do not make an IE minimal;
    CRISPR array DOES add the CRISPR class (Cas ∨ array)
  - DefenseFinder tokens outside the funnel are ignored (promiscuous),
    same policy as protein-level PADLOC
  - leftover family one-word stems are mapped only when unambiguous

resolved_minimal = composition_minimal AND n_unresolved_dark == 0
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from count_unknown_minimal_ie import (  # noqa: E402
    CARGO,
    ORPHAN_ARRAY,
    RM_EXCLUDE,
    RM_RE,
    RM_REASE,
    gene_tag,
)
from recount_unknown_minimal_annotated import FUNNEL_META, PAD  # noqa: E402

ART = FF / "artifacts/unknown_ime_classes"
OUT = ART / "after_hypothesis"
IE_OLD = ART / "unknown_minimal_ie.tsv"
ANN = FF / "artifacts/unknown_reannotate/unknown_all_annotated.tsv"
STACK = FF / "artifacts/unknown_reannotate/still_dark_all/still_dark_calls.tsv"
LEFTOVER = FF / "artifacts/unknown_reannotate/still_dark_all/leftover/leftover_calls.tsv"
HYP = FF / "artifacts/unknown_reannotate/still_dark_all/hypothesis/hypothesis_calls.tsv"
CIS = FF / "artifacts/unknown_reannotate/still_dark_all/hypothesis/cis_elements.tsv"

CARGO_SET = set(CARGO) - {"other_defense"}

STACK_ACC = re.compile(r"^(C/HTH|HTH/regulator|C-candidate)", re.I)
IPR_BARE = re.compile(r"^ipr:", re.I)
CLEAR_ACC = re.compile(
    r"sigma-?70|copg|ribbon-helix-helix|wyl domain|"
    r"helix-turn-helix|\bhth\b|\bxre\b|\bcro\b|arsr|cadc|merr|"
    r"transcriptional regulator|dna-binding|dna binding domain|"
    r"excisionase|\bxis\b",
    re.I,
)
CLEAR_LEFTOVER = re.compile(
    r"transposase|\btnp[a-z]?\b|insertion sequence|"
    r"terminase|capsid|portal|tape.?measure|baseplate|\bholin\b|endolysin|"
    r"major tail|minor tail|tail fiber|tail sheath|head completion|"
    r"scaffold protein|neck protein|phage portal|"
    r"zona occludens|\bzot\b|zonular occludens|"
    r"relaxase|mobilisation|mobilization|\bmob[abcqtv]\b|"
    r"invertase|resolvase",
    re.I,
)
PHAGE_GENERIC = re.compile(r"(^|:)phage protein(\b|$)", re.I)
AMR_BARE = re.compile(
    r"glyoxylase|beta-lactamase superfamily|metallo-beta-lactamase domain",
    re.I,
)
BACT_FALSE = re.compile(r"barnase|endou-colicin|colicine5", re.I)
PIN_TA = re.compile(r"\bpin domain\b", re.I)
HYP_SKIP = re.compile(
    r"psort location|uncharacterized|hypothetical protein|"
    r"unknown function|domain of unknown",
    re.I,
)
FOLDSEEK_ID = re.compile(r"^[0-9][a-z0-9]{3}-assembly", re.I)
KEGG_ONLY = re.compile(r"^kegg_ko=", re.I)
HYP_REP = re.compile(
    r"replication initiator|replication initiation|"
    r"bacteriophage replication|initiator rep|\brep protein\b|"
    r"primase|primpol|huh.?tag|huh motif",
    re.I,
)
HYP_NOT_REP = re.compile(r"\bdnaa\b|replisome organizer", re.I)

# Truncated leftover family stems. Unlisted one-word stems stay dark.
FAMILY_WORD = {
    "hth": "ok_adj",
    "arsr": "ok_adj",
    "asnc": "ok_adj",
    "sigma-70": "ok_adj",
    "cro": "ok_adj",
    "xre": "ok_adj",
    "kila-n": "ok_adj",
    "homeodomain": "ok_adj",
    "gata-type": "ok_adj",
    "swim-type": "ok_adj",
    "transcriptional": "ok_adj",
    "dna-binding": "ok_adj",
    "helix-turn-helix": "ok_adj",
    "transposase": "leftover",
    "terminase": "leftover",
    "phage": "leftover",
    "bacteriophage": "leftover",
    "prevent-host-death": "TA",
    "antitoxin": "TA",
    "replication": "replication",
    "rep": "replication",
    "restriction": "rm_part",
    "integrase": "int",
    "replication-relaxation": "leftover",
}

DF_MAP = [
    (re.compile(r"^RM_type|^RM__|^RM_", re.I), "RM"),
    (re.compile(r"AbiEii|RosmerTA|ShosTA|ietAS", re.I), "TA"),
    (re.compile(r"\bAbi[LDUZ]\b", re.I), "Abi"),
    (re.compile(r"Abi[AK]|AbiP2|Retron|^DRT", re.I), "RT-defense"),
    (re.compile(r"SoFic", re.I), "SoFic"),
    (re.compile(r"cbass", re.I), "CBASS"),
    (re.compile(r"shedu", re.I), "Shedu"),
    (re.compile(r"Thoeris|Ths[AB]", re.I), "Thoeris"),
    (re.compile(r"BREX|Pgl[XZ]|Brx[A-Z]", re.I), "BREX"),
    (re.compile(r"Wadjet|Jet[A-D]", re.I), "Wadjet"),
    (re.compile(r"Lamassu|Lmu[ABC]", re.I), "Lamassu"),
    (re.compile(r"CRISPR|\bcas\b|Cnu[AB]|Type-[IVX]+-[A-Z]", re.I), "CRISPR"),
    (re.compile(r"PDC-|HEC-|PD-T4|PD-Lambda", re.I), "PDC"),
    (re.compile(r"DUF4238|Rst_DUF", re.I), "Rst_DUF4238"),
    (re.compile(r"DnaG|UG6_DnaG", re.I), "replication"),
    (re.compile(r"^VSPR$|^DMS_other", re.I), None),
]


def parse_classes(s: object) -> set[str]:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return set()
    return {x for x in str(s).split("|") if x and x != "other_defense"}


def defense_role(tok: str) -> str | None:
    t = tok.strip()
    if not t:
        return None
    for rx, cls in DF_MAP:
        if rx.search(t):
            if cls is None:
                return None
            if cls == "RM":
                return "rm_part"
            return cls
    return None


def post_filter(role: str, text: str) -> str:
    if role == "AMR" and AMR_BARE.search(text):
        return "ok_adj"
    if role == "bacteriocins" and (
        BACT_FALSE.search(text) or re.search(r"transport|permease|abc-type", text, re.I)
    ):
        return "ok_adj"
    return role


def is_clear_leftover(text: str) -> bool:
    t = text or ""
    if CLEAR_LEFTOVER.search(t) or PHAGE_GENERIC.search(t):
        return True
    return False


def tag_text(text: str, *, allow_aaa: bool, default: str = "ok_adj") -> str:
    """Classify a named product. Default is accessory/dark, not leftover.

    gene_tag leftover is the catch-all for unnamed enzymes — too broad for
    hypothesis and generic PSC domains. Only explicit leftover marks drop
    composition-minimal.
    """
    t = text or ""
    if STACK_ACC.search(t) or CLEAR_ACC.search(t) or IPR_BARE.search(t):
        return "ok_adj"
    if t.lower().startswith("duf") and re.match(r"duf\d+", t, re.I):
        return "ok_adj"
    if is_clear_leftover(t):
        return "leftover"
    if PIN_TA.search(t):
        return "TA"
    if HYP_REP.search(t) and not HYP_NOT_REP.search(t):
        return "replication"
    role = gene_tag(t, allow_aaa_as_rep=allow_aaa)
    if role in CARGO_SET or role in {"rm_part", "int"}:
        return post_filter(role, t)
    return default


def classify_stack(call: str, *, allow_aaa: bool) -> tuple[str, str, bool]:
    """Return (role, source, resolved)."""
    if not call or call == "still_dark":
        return "ok_adj", "dark", False
    if STACK_ACC.search(call):
        return "ok_adj", "stack_acc", call.startswith("C/HTH") or call.startswith("HTH/")
    if call.startswith("ipr:"):
        return "ok_adj", "ipr_bare", False
    if call.startswith("defense:"):
        # Same policy as protein-level PADLOC: DF HMMs do not add classes.
        return "ok_adj", "df_ignore", False
    prod = call
    for pref in ("psc:", "family:"):
        if call.startswith(pref):
            prod = call[len(pref):]
            break
    return tag_text(prod, allow_aaa=allow_aaa), "stack_psc", True


def leftover_strict_parts(lab: object) -> list[str]:
    if lab is None or (isinstance(lab, float) and pd.isna(lab)):
        return []
    out = []
    for part in re.split(r";\s*", str(lab)):
        if part.startswith("pfamE1:"):
            continue
        if part.strip():
            out.append(part.strip())
    return out


def classify_leftover(parts: list[str], *, allow_aaa: bool) -> tuple[str, str, bool]:
    if not parts:
        return "ok_adj", "none", False
    ranked: list[tuple[int, str, str]] = []
    for part in parts:
        if part.startswith("expert:"):
            ranked.append((0, part[7:], "leftover_expert"))
        elif part.startswith("pscc:"):
            ranked.append((1, part[5:], "leftover_pscc"))
        elif part.startswith("family:"):
            ranked.append((2, part[7:], "leftover_family"))
        elif part.startswith("psc_vs:"):
            continue
        elif part.startswith("KO:") or part.startswith("EC:"):
            continue
        else:
            ranked.append((3, part, "leftover_other"))
    if not ranked:
        return "ok_adj", "leftover_weak", False
    ranked.sort(key=lambda x: x[0])
    _prio, text, src = ranked[0]
    key = text.strip().lower()
    if src == "leftover_family" and " " not in key and "/" not in key:
        if key.startswith("duf"):
            return "ok_adj", src, True
        if key in FAMILY_WORD:
            return FAMILY_WORD[key], src, True
        return "ok_adj", "family_stem_dark", False
    return tag_text(text, allow_aaa=allow_aaa), src, True


def extract_hyp_text(lab: object) -> str:
    if lab is None or (isinstance(lab, float) and pd.isna(lab)):
        return ""
    bits = []
    for part in re.split(r";\s*", str(lab)):
        if part.startswith("eggnog:"):
            for m in re.finditer(r"(?:Preferred_name|Description)=([^;]+)", part):
                val = m.group(1).strip()
                if not val or HYP_SKIP.search(val) or KEGG_ONLY.search(val):
                    continue
                if val.lower() in {"membrane", "cytoplasmic", "cytoplasmicmembrane"}:
                    continue
                if val.lower().startswith("psort"):
                    continue
                bits.append(val)
        elif part.startswith("hh:"):
            if HYP_SKIP.search(part):
                continue
            bits.append(part)
        elif part.startswith("foldseek:"):
            rest = part[9:].strip()
            if FOLDSEEK_ID.match(rest) and " " not in rest:
                continue
            if HYP_SKIP.search(rest):
                continue
            bits.append(rest)
    return " ".join(bits)


def classify_hyp(lab: object, *, allow_aaa: bool) -> tuple[str, str, bool]:
    text = extract_hyp_text(lab)
    if not text:
        return "ok_adj", "hyp_empty", False
    return tag_text(text, allow_aaa=allow_aaa), "hyp", False  # never resolves


def overlay_role(
    stack_call: str,
    leftover_lab: object,
    hyp_lab: object,
    *,
    allow_aaa: bool,
) -> tuple[str, str, bool]:
    role, src, resolved = classify_stack(stack_call, allow_aaa=allow_aaa)
    if src != "dark":
        return role, src, resolved
    parts = leftover_strict_parts(leftover_lab)
    if parts:
        role, src, resolved = classify_leftover(parts, allow_aaa=allow_aaa)
        if src != "leftover_weak" and src != "family_stem_dark" and src != "none":
            return role, src, resolved
        if src == "family_stem_dark":
            # still try hyp
            pass
        else:
            return role, src, resolved
    if hyp_lab and str(hyp_lab).strip():
        return classify_hyp(hyp_lab, allow_aaa=allow_aaa)
    return "ok_adj", "still_dark", False


def char_from_classes(cls: set[str]) -> tuple[int, str, int]:
    has_rep = "replication" in cls
    cargo_no_rep = cls - {"replication"}
    if cargo_no_rep:
        return len(cargo_no_rep), "|".join(sorted(cargo_no_rep)), int(has_rep)
    if has_rep:
        return 1, "replication", 1
    return 0, "", 0


def write_funnel(ie: pd.DataFrame) -> pd.DataFrame:
    pad = pd.read_csv(PAD, sep="\t")
    pad_toks = {}
    for r in pad.itertuples():
        pad_toks[r.integrase_id] = str(getattr(r, "padloc_systems", "") or "").lower()
    rows = []
    for layer, name, genes, not_inc, src, presence in FUNNEL_META:
        if name == "virulence":
            rows.append({
                "class": name, "genes_included": genes, "not_included": not_inc,
                "layer": layer, "n_ie_with_class": 0, "n_composition_minimal": 0,
                "n_comp_min_plus_rep": 0, "n_class_minimal": 0,
                "n_tight_minimal_okadj_le3": 0, "n_comp_min_le6_cds": 0,
                "n_resolved_minimal": 0,
                "literature": src, "presence_lock": presence,
            })
            continue
        if name in {"Gabija", "Septu", "Hachiman", "Zorya", "DISARM"}:
            n_has = sum(1 for iid in ie.integrase_id if name.lower() in pad_toks.get(iid, ""))
            rows.append({
                "class": name, "genes_included": genes, "not_included": not_inc,
                "layer": layer, "n_ie_with_class": n_has, "n_composition_minimal": 0,
                "n_comp_min_plus_rep": 0, "n_class_minimal": 0,
                "n_tight_minimal_okadj_le3": 0, "n_comp_min_le6_cds": 0,
                "n_resolved_minimal": 0,
                "literature": src, "presence_lock": "absent" if n_has == 0 else "yes",
            })
            continue
        has = ie.classes.fillna("").apply(lambda s, n=name: n in set(s.split("|")) if s else False)
        sub = ie[ie.char_class == name]
        rows.append({
            "class": name, "genes_included": genes, "not_included": not_inc,
            "layer": layer, "n_ie_with_class": int(has.sum()),
            "n_composition_minimal": int(sub.composition_minimal.sum()),
            "n_comp_min_plus_rep": int(((sub.composition_minimal == 1) & (sub.has_rep == 1) & (sub.char_class != "replication")).sum()),
            "n_class_minimal": int(sub.class_minimal.sum()),
            "n_tight_minimal_okadj_le3": int(sub.tight_minimal.sum()),
            "n_comp_min_le6_cds": int(sub.comp_min_le6.sum()),
            "n_resolved_minimal": int(sub.resolved_minimal.sum()),
            "literature": src, "presence_lock": presence,
        })
    by = pd.DataFrame(rows)[[
        "class", "genes_included", "not_included", "layer",
        "n_ie_with_class", "n_composition_minimal", "n_comp_min_plus_rep",
        "n_class_minimal", "n_tight_minimal_okadj_le3", "n_comp_min_le6_cds",
        "n_resolved_minimal", "literature", "presence_lock",
    ]]
    by.to_csv(OUT / "by_class.tsv", sep="\t", index=False)
    return by


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(IE_OLD, sep="\t")
    ann = pd.read_csv(ANN, sep="\t", usecols=["protein_id", "integrase_id", "bakta_dark", "pfam_hit"])
    stack = pd.read_csv(STACK, sep="\t", usecols=["protein_id", "integrase_id", "aa_len", "call"])
    lo = pd.read_csv(LEFTOVER, sep="\t", usecols=["protein_id", "leftover_label", "leftover_named"])
    hyp = pd.read_csv(HYP, sep="\t", usecols=["protein_id", "hyp_label", "hyp_named"])
    cis = pd.read_csv(CIS, sep="\t", usecols=["integrase_id", "kind"])

    crispr_ie = set(cis.loc[cis.kind.str.startswith("crispr:"), "integrase_id"])

    dark = stack.merge(lo, on="protein_id", how="left").merge(hyp, on="protein_id", how="left")
    by_ie_dark: dict[str, list[dict]] = {}
    for r in dark.itertuples():
        by_ie_dark.setdefault(r.integrase_id, []).append({
            "protein_id": r.protein_id,
            "aa_len": r.aa_len,
            "call": r.call,
            "leftover_label": r.leftover_label,
            "hyp_label": r.hyp_label,
        })

    named_by_ie = (
        ann[(ann.bakta_dark == 0) | (ann.pfam_hit == 1)]
        .groupby("integrase_id").size().to_dict()
    )
    n_cds_ann = ann.groupby("integrase_id").size().to_dict()

    rows = []
    cds_rows = []
    for r in old.itertuples():
        iid = r.integrase_id
        allow_aaa = bool(r.from_rep)
        cls = parse_classes(r.classes)
        recs = by_ie_dark.get(iid, [])
        n_ok = int(r.n_ok_adj) if pd.notna(r.n_ok_adj) else 0
        n_leftover = int(r.n_leftover) if pd.notna(r.n_leftover) else 0
        leftover_ex = [
            x for x in str(r.leftover_examples or "").split("; ")
            if x and x.lower() != "nan"
        ]
        overlay_leftover = []
        overlay_cargo = []
        n_unresolved = 0
        n_overlay = 0
        has_re_new = False

        for x in recs:
            n_overlay += 1
            role, src, resolved = overlay_role(
                x["call"], x["leftover_label"], x["hyp_label"], allow_aaa=allow_aaa,
            )
            cds_rows.append({
                "integrase_id": iid,
                "protein_id": x["protein_id"],
                "aa_len": x["aa_len"],
                "role": role,
                "source": src,
                "resolved": int(resolved),
            })
            if not resolved:
                n_unresolved += 1
            # Dark CDS were already ok_adj on the 436 count. Re-role them.
            if role == "leftover":
                n_leftover += 1
                n_ok = max(0, n_ok - 1)
                blob = x["call"] if x["call"] and x["call"] != "still_dark" else (
                    x["leftover_label"] or x["hyp_label"] or ""
                )
                overlay_leftover.append(f"{src}:{blob}")
            elif role in CARGO_SET:
                n_ok = max(0, n_ok - 1)
                if role not in cls:
                    overlay_cargo.append(role)
                cls.add(role)
            elif role == "rm_part":
                n_ok = max(0, n_ok - 1)
                blob = " ".join(str(v or "") for v in (x["call"], x["leftover_label"], x["hyp_label"]))
                if not RM_EXCLUDE.search(blob) and RM_REASE.search(blob):
                    has_re_new = True
                    if "RM" not in cls:
                        overlay_cargo.append("RM")
                    cls.add("RM")
            elif role == "int":
                n_ok = max(0, n_ok - 1)

        if iid in ORPHAN_ARRAY or iid in crispr_ie:
            if "CRISPR" not in cls:
                overlay_cargo.append("CRISPR")
            cls.add("CRISPR")
        if has_re_new:
            cls.add("RM")
        cls.discard("other_defense")

        n_char, char, has_rep = char_from_classes(cls)
        class_min = n_char == 1
        comp_min = class_min and n_leftover == 0
        n_cds = r.n_cds if pd.notna(r.n_cds) else n_cds_ann.get(iid)
        resolved_min = int(comp_min and n_unresolved == 0 and iid in n_cds_ann)
        was_cm = int(r.composition_minimal)
        if was_cm and not comp_min:
            if n_char > 1:
                drop = "second_class"
            elif n_leftover > 0:
                drop = "named_leftover"
            else:
                drop = "lost_class"
        elif not was_cm and comp_min:
            drop = "gained"
        else:
            drop = ""

        d = r._asdict()
        d.update({
            "classes": "|".join(sorted(cls)),
            "n_classes": len(cls),
            "char_class": char,
            "n_char_classes": n_char,
            "has_rep": has_rep,
            "n_ok_adj": n_ok,
            "n_leftover": n_leftover,
            "leftover_examples": "; ".join((leftover_ex + overlay_leftover)[:8]),
            "class_minimal": int(class_min),
            "composition_minimal": int(comp_min),
            "compact_le6": int(class_min and n_cds is not None and n_cds <= 6),
            "comp_min_le6": int(comp_min and n_cds is not None and n_cds <= 6),
            "tight_minimal": int(comp_min and n_ok <= 3),
            "n_unresolved_dark": n_unresolved,
            "n_overlay_dark": n_overlay,
            "n_named_bakta_pfam": int(named_by_ie.get(iid, 0)),
            "resolved_minimal": resolved_min,
            "overlay_cargo_added": "|".join(sorted(set(overlay_cargo))),
            "n_overlay_leftover": len(overlay_leftover),
            "crispr_array": int(iid in crispr_ie),
            "was_composition_minimal": was_cm,
            "delta": drop,
        })
        if "Index" in d:
            d.pop("Index")
        rows.append(d)

    ie = pd.DataFrame(rows)
    ie.to_csv(OUT / "unknown_minimal_ie.tsv", sep="\t", index=False)
    pd.DataFrame(cds_rows).to_csv(OUT / "dark_cds_roles.tsv", sep="\t", index=False)

    n = len(ie)
    n_class = int(ie.class_minimal.sum())
    n_comp = int(ie.composition_minimal.sum())
    n_comp_rep = int(((ie.composition_minimal == 1) & (ie.has_rep == 1) & (ie.char_class != "replication")).sum())
    n_comp_norep = int(((ie.composition_minimal == 1) & ((ie.has_rep == 0) | (ie.char_class == "replication"))).sum())
    n_res = int(ie.resolved_minimal.sum())
    old_cm = set(old.loc[old.composition_minimal == 1, "integrase_id"])
    new_cm = set(ie.loc[ie.composition_minimal == 1, "integrase_id"])
    lost = old_cm - new_cm
    gained = new_cm - old_cm
    kept = old_cm & new_cm

    summary = pd.DataFrame([
        {"metric": "unknown_n", "n": n},
        {"metric": "class_minimal", "n": n_class},
        {"metric": "composition_minimal_was", "n": len(old_cm)},
        {"metric": "composition_minimal_now", "n": n_comp},
        {"metric": "comp_min_kept", "n": len(kept)},
        {"metric": "comp_min_lost", "n": len(lost)},
        {"metric": "comp_min_gained", "n": len(gained)},
        {"metric": "comp_min_module_plus_rep", "n": n_comp_rep},
        {"metric": "comp_min_module_or_rep_only", "n": n_comp_norep},
        {"metric": "resolved_minimal", "n": n_res},
        {"metric": "tight_minimal_okadj_le3", "n": int(ie.tight_minimal.sum())},
        {"metric": "comp_min_le6_cds", "n": int(ie.comp_min_le6.sum())},
        {"metric": "class_min_le6_cds", "n": int(ie.compact_le6.sum())},
        {"metric": "no_cargo_class", "n": int((ie.n_char_classes == 0).sum())},
        {"metric": "multi_class", "n": int((ie.n_char_classes > 1).sum())},
        {"metric": "crispr_array_ie", "n": int(ie.crispr_array.sum())},
        {"metric": "comp_min_with_unresolved", "n": int(((ie.composition_minimal == 1) & (ie.n_unresolved_dark > 0)).sum())},
        {"metric": "comp_min_zero_unresolved", "n": int(((ie.composition_minimal == 1) & (ie.n_unresolved_dark == 0)).sum())},
    ])
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    by = write_funnel(ie)

    delta = ie[ie.delta.isin(["second_class", "named_leftover", "lost_class", "gained"])][
        ["integrase_id", "delta", "char_class", "classes", "n_char_classes",
         "n_leftover", "n_unresolved_dark", "overlay_cargo_added",
         "leftover_examples", "was_composition_minimal", "composition_minimal",
         "resolved_minimal", "crispr_array"]
    ]
    delta.to_csv(OUT / "delta_vs_436.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\ncomposition-minimal", n_comp, "was", len(old_cm),
          "kept", len(kept), "lost", len(lost), "gained", len(gained))
    print("resolved_minimal", n_res)
    print("\nby class composition-minimal / resolved:")
    print(by[["class", "n_ie_with_class", "n_class_minimal",
              "n_composition_minimal", "n_resolved_minimal"]].to_string(index=False))
    if len(lost):
        lost_ie = ie[ie.integrase_id.isin(lost)]
        print("\nlost reasons:")
        print(lost_ie.delta.value_counts().to_string())
        print(lost_ie.char_class.value_counts().to_string())
    if len(gained):
        gie = ie[ie.integrase_id.isin(gained)]
        print("\ngained char:")
        print(gie.char_class.value_counts().to_string())
        print("overlay cargo", gie.overlay_cargo_added.value_counts().head(12).to_string())
    print("wrote", OUT)


if __name__ == "__main__":
    main()
