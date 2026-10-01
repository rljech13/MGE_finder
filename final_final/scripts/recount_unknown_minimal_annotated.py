#!/usr/bin/env python3
"""Recount composition-minimal unknown IEs after Pfam on dark CDS.

Locked definition (not proven):
  integrase + exactly one funnel class
  + accessories of that class (HTH/C, ArsR, TA partner, Xis, hyp/DUF)
  + optional replicase
  Named leftover (tnp, Mob, other enzymes) disqualifies.
  Unresolved dark CDS do not.

Pfam may:
  - name a leftover enzyme on a hyp → IE drops out of composition-minimal
  - name a high-confidence cargo hallmark → class is added
    (new solo module, or second class → not minimal)
  - name int / HTH / DUF → still allowed

PADLOC protein hits are not used to add classes (promiscuous HMMs).
IE-level PADLOC already sits in the seed from the Bakta/PADLOC funnel.
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
    gene_tag,
    refine_rm,
)

ART = FF / "artifacts/unknown_ime_classes"
IE_OLD = ART / "unknown_minimal_ie.tsv"
ANN = FF / "artifacts/unknown_reannotate/unknown_all_annotated.tsv"
PAD = FF / "artifacts/nb04_cluster_unclassified/padloc_per_ie_stage8_unc.tsv"

CARGO_SET = set(CARGO) - {"other_defense"}

# High-confidence Pfam → cargo / int / acc / leftover. Anything else stays hyp.
PFAM_INT = re.compile(r"^(Phage_integrase|rve|Phage_int_SAM)$", re.I)
PFAM_ACC = re.compile(
    r"^(HTH_\d+|Cro|Xre|ArsR|MerR|SmtB|Sigma70_r4|Zn_Ribbon_1|TPR_\d+|DUF\d+)$",
    re.I,
)
PFAM_REP = re.compile(
    r"^(Viral_Rep|Rep_trans|RepA|RepSA|Rep_3|DnaB_C|DnaG|REP_ORF2-G2P)$",
    re.I,
)
PFAM_RM_RE = re.compile(r"^(Mrr_cat|GmrSD_N|RE_|Endonuc)", re.I)
PFAM_RM_MT = re.compile(r"^(N6_Mtase|N6_N4_Mtase|MmeI_Mtase|Methylase_S)$", re.I)
PFAM_TA = re.compile(
    r"^(HicA_toxin|HicB|RelE|ParE_toxin|MazE_antitoxin|MazF|PemK_toxin|"
    r"AbiEi|AbiEii|VapC|MqsA_antitoxin)$",
    re.I,
)
PFAM_LEFTOVER = re.compile(
    r"^(Zot|Peptidase_|Trypsin|Transpos|DDE_Tnp|Mob_|Relaxase|"
    r"Pkinase|Flg_|NikA|Resolvase|Invertase|IS[0-9])",
    re.I,
)

FUNNEL_META = [
    ("A", "replication",
     "Rep helicase; Rep / Rep_trans / RepSA; PrimPol (AEP); PolA; PolB; DnaG-like primase",
     "", "Bellanger 2017", "yes"),
    ("A", "TA",
     "HicA+HicB; RelE+RelB / ParE+ParD; VapC+VapB; MazF+MazE / PemK+PemI; AbiEii+AbiEi (type IV)",
     "AbiL/D/U/Z (residual Abi)", "Bellanger 2017", "yes"),
    ("A", "Par",
     "ParA (Soj); ParB",
     "ParE/ParD (that pair is TA)", "Bellanger 2017", "yes"),
    ("B", "RM",
     "I: HsdR+HsdM+HsdS; II: REase+MTase; IIG: fused RM (MmeI-like); IV: Mrr, GmrSD",
     "C/Xre (regulator); MTase-only", "Kobayashi; Botelho 2023", "yes"),
    ("B", "CRISPR",
     "any Cas (Cas1–14, Cmr, Csm, Cse…) or a CRISPR array; either one is enough",
     "a complete operon is not required", "Pinilla-Redondo 2023", "yes"),
    ("B", "Abi",
     "AbiL; AbiD; AbiU; AbiZ",
     "AbiE (that pair is TA); AbiA/K/P2 (RT-defense)", "DefenseFinder", "yes"),
    ("B", "RT-defense",
     "DRT (RT+effector); retron RT; AbiA; AbiK; AbiP2",
     "", "Gao; Millman", "yes"),
    ("B", "SoFic",
     "SoFic (Fic/AMPylase), single-gene",
     "a bare Fic/Fido without the SoFic name", "Millman 2022", "yes"),
    ("B", "CBASS",
     "Cdn/Dnc cyclase + Cap effector (Cap2/4/5, TIR, SAVED…)",
     "", "Botelho 2023", "yes"),
    ("B", "Shedu", "SduA", "", "DefenseFinder", "yes"),
    ("B", "Thoeris", "ThsA (NAD hydrolase); ThsB (TIR)", "", "Botelho 2023", "rare"),
    ("B", "BREX", "BrxA; BrxB; BrxC; PglX; PglZ; BrxL", "", "DefenseFinder", "rare"),
    ("B", "Wadjet", "JetA; JetB; JetC; JetD", "", "DefenseFinder", "rare"),
    ("B", "Lamassu", "LmuA; LmuB; LmuC", "", "DefenseFinder", "rare"),
    ("B", "PDC",
     "PADLOC PDC-S02, S15, S70, M34 and other PDC/HEC HMMs",
     "not a biological family", "PADLOC", "yes"),
    ("B", "Rst_DUF4238", "DUF4238", "", "Rousset 2022", "yes"),
    ("B", "Gabija", "GajA; GajB", "", "DefenseFinder", "absent"),
    ("B", "Septu", "PtuA; PtuB", "", "DefenseFinder", "absent"),
    ("B", "Hachiman", "HamA; HamB", "", "DefenseFinder", "absent"),
    ("B", "Zorya", "ZorA; ZorB; ZorC; ZorD; ZorE", "", "DefenseFinder", "absent"),
    ("B", "DISARM", "DrmA; DrmB; DrmC; DrmD; DrmE", "", "DefenseFinder", "absent"),
    ("C", "AMR",
     "named Bla, Tet, Erm, Aph, Aac, Lnu, Van…",
     "a bare metallo-beta-lactamase domain", "Bellanger 2017", "yes"),
    ("C", "metals",
     "CadA / Cd-ATPase (transporting/translocating/exporting); "
     "MerA/B/C/E/P/T; ArsA/B/C/D/I/N + ACR3; CopA/B/C/D/Z; Cus; Czc; ChrA",
     "ArsR/MerR/CopG; polymerases; HMA / cation transporter / multicopper / cytochrome",
     "Bellanger 2017", "yes"),
    ("C", "bacteriocins",
     "circularin A; uberolysin; bacteriocin structural",
     "", "Bellanger 2017", "yes"),
    ("C", "virulence", "not scanned", "", "Botelho 2023", "not_scanned"),
]


def tag_pfam(hmm: object) -> str | None:
    if hmm is None or (isinstance(hmm, float) and pd.isna(hmm)):
        return None
    h = str(hmm)
    if PFAM_INT.search(h):
        return "int"
    if PFAM_ACC.search(h):
        return "ok_adj"
    if PFAM_REP.search(h):
        return "replication"
    if PFAM_RM_RE.search(h):
        return "rm_part"
    if PFAM_RM_MT.search(h):
        return "rm_part"
    if PFAM_TA.search(h):
        return "TA"
    if PFAM_LEFTOVER.search(h):
        return "leftover"
    return None


def cds_role(bakta: str, hmm: object, allow_aaa: bool) -> tuple[str, str]:
    """Return (role, source). Bakta first; Pfam only overrides hyp/ok_adj."""
    bt = gene_tag(bakta or "", allow_aaa_as_rep=allow_aaa)
    pf = tag_pfam(hmm)
    if bt in CARGO_SET or bt == "rm_part":
        return bt, "bakta"
    if bt == "int":
        return bt, "bakta"
    if bt == "leftover":
        # Pfam can rescue a misnamed leftover into int/cargo/acc
        if pf in CARGO_SET or pf in {"rm_part", "int", "ok_adj"}:
            return pf, "pfam"
        return "leftover", "bakta"
    # bakta ok_adj / hyp: Pfam may name cargo, leftover, int, acc
    if pf is None:
        return "ok_adj", "hyp"
    return pf, "pfam"


def parse_classes(s: object) -> set[str]:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return set()
    return {x for x in str(s).split("|") if x and x != "other_defense"}


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
            "literature": src, "presence_lock": presence,
        })
    by = pd.DataFrame(rows)[[
        "class", "genes_included", "not_included", "layer",
        "n_ie_with_class", "n_composition_minimal", "n_comp_min_plus_rep",
        "n_class_minimal", "n_tight_minimal_okadj_le3", "n_comp_min_le6_cds",
        "literature", "presence_lock",
    ]]
    by.to_csv(ART / "ime_class_funnel.tsv", sep="\t", index=False)
    by.to_csv(ART / "ie_class_funnel.tsv", sep="\t", index=False)
    by.to_csv(ART / "unknown_minimal_by_class.tsv", sep="\t", index=False)
    return by


def main() -> None:
    old = pd.read_csv(IE_OLD, sep="\t")
    ann = pd.read_csv(ANN, sep="\t")
    by_ie: dict[str, list[dict]] = {}
    for r in ann.itertuples():
        by_ie.setdefault(r.integrase_id, []).append({
            "bakta": r.bakta_product,
            "hmm": r.pfam_hmm,
        })

    rows = []
    for r in old.itertuples():
        iid = r.integrase_id
        allow_aaa = bool(r.from_rep)
        seed = parse_classes(r.classes)
        recs = by_ie.get(iid)
        if not recs:
            d = r._asdict()
            d["ann_used"] = 0
            rows.append(d)
            continue
        roles, sources = [], []
        for x in recs:
            role, src = cds_role(x["bakta"], x["hmm"], allow_aaa)
            roles.append(role)
            sources.append(src)
        products = [x["bakta"] for x in recs]
        cls = set(seed)
        for t in roles:
            if t in CARGO_SET:
                cls.add(t)
        if any(t == "replication" for t in roles) or r.from_rep:
            cls.add("replication")
        if iid in ORPHAN_ARRAY:
            cls.add("CRISPR")
        # MTase-only from Pfam is not RM unless the seed already had RM
        if "RM" not in seed:
            has_re = any(
                tag_pfam(x["hmm"]) == "rm_part" and PFAM_RM_RE.search(str(x["hmm"] or ""))
                for x in recs
            )
            has_re = has_re or any(
                gene_tag(p, allow_aaa_as_rep=allow_aaa) == "rm_part"
                and re.search(r"restrict|rease|hsd[rms]|gmrsd|\bmrr\b|endonuclease", p or "", re.I)
                for p in products
            )
            if not has_re:
                cls.discard("RM")
        if "RM" in seed:
            cls.add("RM")
        cls.discard("other_defense")

        n_leftover = sum(t == "leftover" for t in roles)
        n_ok = sum(t == "ok_adj" for t in roles)
        leftover_ex = []
        for x, t, src in zip(recs, roles, sources):
            if t == "leftover":
                leftover_ex.append(f"{x['bakta']}|{x['hmm'] or '-'}:{src}")
        cargo_no_rep = cls - {"replication"}
        has_rep = "replication" in cls
        if cargo_no_rep:
            n_char = len(cargo_no_rep)
            char = "|".join(sorted(cargo_no_rep))
        elif has_rep:
            n_char = 1
            char = "replication"
        else:
            n_char = 0
            char = ""
        class_min = n_char == 1
        comp_min = class_min and n_leftover == 0
        n_cds = r.n_cds if pd.notna(r.n_cds) else len(recs)
        d = r._asdict()
        d.update({
            "ann_used": 1,
            "classes": "|".join(sorted(cls)),
            "n_classes": len(cls),
            "char_class": char,
            "n_char_classes": n_char,
            "has_rep": int(has_rep),
            "n_ok_adj": n_ok,
            "n_leftover": n_leftover,
            "leftover_examples": "; ".join(leftover_ex[:8]),
            "class_minimal": int(class_min),
            "composition_minimal": int(comp_min),
            "compact_le6": int(class_min and n_cds is not None and n_cds <= 6),
            "comp_min_le6": int(comp_min and n_cds is not None and n_cds <= 6),
            "tight_minimal": int(comp_min and n_ok <= 3),
        })
        rows.append(d)

    ie = pd.DataFrame(rows)
    if "Index" in ie.columns:
        ie = ie.drop(columns=["Index"])
    ie.to_csv(ART / "unknown_minimal_ie.tsv", sep="\t", index=False)

    n = len(ie)
    n_class = int(ie.class_minimal.sum())
    n_comp = int(ie.composition_minimal.sum())
    n_comp_rep = int(((ie.composition_minimal == 1) & (ie.has_rep == 1) & (ie.char_class != "replication")).sum())
    n_comp_norep = int(((ie.composition_minimal == 1) & ((ie.has_rep == 0) | (ie.char_class == "replication"))).sum())
    summary = pd.DataFrame([
        {"metric": "unknown_n", "n": n},
        {"metric": "class_minimal", "n": n_class},
        {"metric": "composition_minimal", "n": n_comp},
        {"metric": "comp_min_module_plus_rep", "n": n_comp_rep},
        {"metric": "comp_min_module_or_rep_only", "n": n_comp_norep},
        {"metric": "tight_minimal_okadj_le3", "n": int(ie.tight_minimal.sum())},
        {"metric": "comp_min_le6_cds", "n": int(ie.comp_min_le6.sum())},
        {"metric": "class_min_le6_cds", "n": int(ie.compact_le6.sum())},
        {"metric": "no_cargo_class", "n": int((ie.n_char_classes == 0).sum())},
        {"metric": "multi_class", "n": int((ie.n_char_classes > 1).sum())},
    ])
    summary.to_csv(ART / "unknown_minimal_summary.tsv", sep="\t", index=False)
    by = write_funnel(ie)

    old_cm = set(old.loc[old.composition_minimal == 1, "integrase_id"])
    new_cm = set(ie.loc[ie.composition_minimal == 1, "integrase_id"])
    print(summary.to_string(index=False))
    print("\ncomposition-minimal", n_comp, "was", len(old_cm),
          "kept", len(old_cm & new_cm), "lost", len(old_cm - new_cm),
          "gained", len(new_cm - old_cm))
    print("\nby class composition-minimal:")
    print(by[["class", "n_ie_with_class", "n_class_minimal", "n_composition_minimal", "n_comp_min_plus_rep"]].to_string(index=False))
    lost = ie[ie.integrase_id.isin(old_cm - new_cm)]
    if len(lost):
        print("\nlost", len(lost), "char counts:")
        print(lost.char_class.value_counts().to_string())
        print("sample leftover:")
        for rec in lost.head(15).itertuples():
            print(f"  {rec.integrase_id} {rec.char_class} n_char={rec.n_char_classes} {str(rec.leftover_examples)[:100]}")
    gained = ie[ie.integrase_id.isin(new_cm - old_cm)]
    if len(gained):
        print("\ngained", len(gained))
        print(gained.char_class.value_counts().to_string())


if __name__ == "__main__":
    main()
