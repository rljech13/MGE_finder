#!/usr/bin/env python3
"""Score module-minimal unknown IEs.

Locked definition:
  integrase + exactly one funnel module
  + accessories necessary or possible for that module
  + optional replicase
  Unresolved hyp/DUF and named leftover are not minimal.

The IE already has an integrase (atlas). Bakta often leaves it unnamed
(200–550 aa). PADLOC / the class scan may likewise leave the module
unnamed — one leftover hyp can be that module. Remaining dark CDS stay
unresolved (the scan target).
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
from count_unknown_minimal_ie import gene_tag  # noqa: E402

ART = FF / "artifacts/unknown_ime_classes"
IE = ART / "unknown_minimal_ie.tsv"
PROT = FF / "artifacts/unknown_reannotate/proteins.tsv"
SMALL_RM = FF / "artifacts/unknown_le3_rm/rm_allow_rep.tsv"
SMALL_TA = FF / "artifacts/unknown_le3_ta/ta_allow_rep.tsv"
DROP_RM = {"JAWRJC010000145.1_11"}  # HsdS counted as integrase; not RM cargo

HYP = re.compile(r"hypothetical|uncharacterized|conserved domain protein|duf\d+", re.I)
BARE_TOXIN = re.compile(r"^(putative )?toxin\b", re.I)
REG = re.compile(
    r"helix-turn-helix|hth |cro/c1|\bxre\b|c protein|"
    r"transcriptional regulator|dna-binding protein|"
    r"\barsr\b|\bsmtb\b|\bmerr\b|\bcadc\b|metalloregulator",
    re.I,
)
XIS = re.compile(r"excisionase|\bxis\b", re.I)
CBASS_ACC = re.compile(r"\be1\b|\be2\b|\bjab\b|\bhorma\b|\bquec\b|\bcap[0-9]?\b", re.I)
PDC_CORE = re.compile(r"duf1016|duf5343|pdc-", re.I)
RST_CORE = re.compile(r"duf4238", re.I)
PADLOC_CAN_BE_DARK = {
    "RM", "PDC", "CBASS", "RT-defense", "SoFic", "Wadjet", "Lamassu",
    "Shedu", "Thoeris", "BREX", "CRISPR", "Abi",
}


def raw_role(prod: str, char: str, has_rep_label: bool) -> str:
    p = prod or ""
    if BARE_TOXIN.search(p.strip()):
        return "module" if char == "TA" else "other_class"
    tag = gene_tag(p, allow_aaa_as_rep=has_rep_label)
    if tag == "int":
        return "int"
    if tag == "replication":
        return "rep"
    if tag in {"RM", "rm_part"}:
        return "module" if char == "RM" else "other_class"
    if tag == char:
        return "module"
    if tag in {
        "TA", "Par", "CRISPR", "Abi", "RT-defense", "SoFic", "CBASS",
        "Shedu", "Thoeris", "BREX", "Wadjet", "Lamassu", "PDC",
        "Rst_DUF4238", "AMR", "metals", "bacteriocins",
    }:
        return "other_class"
    if char == "PDC" and PDC_CORE.search(p):
        return "module"
    if char == "Rst_DUF4238" and RST_CORE.search(p):
        return "module"
    if char == "CBASS" and CBASS_ACC.search(p):
        return "accessory"
    if XIS.search(p):
        return "accessory"
    if REG.search(p):
        return "accessory"
    if HYP.search(p):
        return "unresolved"
    return "leftover"


def acc_kind(prod: str) -> str:
    p = prod or ""
    if XIS.search(p):
        return "Xis"
    if re.search(r"arsr|smtb|cadc|metalloregulator", p, re.I):
        return "ArsR/CadC"
    if re.search(r"merr", p, re.I):
        return "MerR"
    if CBASS_ACC.search(p):
        return "CBASS adaptor"
    if REG.search(p):
        return "C/HTH"
    if HYP.search(p):
        return "TA partner"
    return "other acc"


def score_ie(rows: list[dict], char: str, has_rep_label: bool) -> dict:
    roles = [raw_role(r["bakta_product"], char, has_rep_label) for r in rows]
    # Atlas integrase: one Tyr-recombinase-sized hyp if Bakta hid it.
    if roles.count("int") == 0:
        cands = [
            i for i, (r, role) in enumerate(zip(rows, roles))
            if role == "unresolved" and 200 <= int(r["aa_len"]) <= 550
        ]
        if cands:
            # prefer ~350–420 aa (typical Tyr rec here)
            pick = min(cands, key=lambda i: abs(int(rows[i]["aa_len"]) - 380))
            roles[pick] = "int"

    # PADLOC/scan already called the class; one remaining dark ORF can be it.
    if roles.count("module") == 0 and char in PADLOC_CAN_BE_DARK:
        cands = [i for i, role in enumerate(roles) if role == "unresolved"]
        if cands:
            pick = max(cands, key=lambda i: int(rows[i]["aa_len"]))
            roles[pick] = "module"

    # TA: one unnamed neighbour = partner.
    if (
        char == "TA"
        and roles.count("unresolved") == 1
        and roles.count("module") >= 1
        and roles.count("leftover") == 0
        and roles.count("other_class") == 0
    ):
        roles = ["accessory" if r == "unresolved" else r for r in roles]

    n = {k: roles.count(k) for k in (
        "int", "rep", "module", "accessory", "unresolved", "leftover", "other_class"
    )}
    kinds = []
    for r, role in zip(rows, roles):
        if role == "accessory":
            kinds.append(acc_kind(r["bakta_product"]))
    proven = (
        n["leftover"] == 0
        and n["unresolved"] == 0
        and n["other_class"] == 0
        and (n["module"] + n["rep"] > 0)
        and n["int"] >= 1
    )
    pending = (
        n["leftover"] == 0
        and n["other_class"] == 0
        and n["unresolved"] > 0
        and (n["module"] + n["rep"] > 0)
    )
    return {
        "n_int": n["int"],
        "n_rep_cds": n["rep"],
        "n_module": n["module"],
        "n_accessory": n["accessory"],
        "n_unresolved": n["unresolved"],
        "n_named_leftover": n["leftover"],
        "n_other_class": n["other_class"],
        "n_module_plus_acc": n["module"] + n["accessory"],
        "accessory_kinds": "|".join(sorted(set(kinds))),
        "module_minimal": int(proven),
        "module_pending": int(pending),
    }


def main() -> None:
    ie = pd.read_csv(IE, sep="\t")
    prot = pd.read_csv(PROT, sep="\t")
    by: dict[str, list[dict]] = {}
    for r in prot.itertuples():
        by.setdefault(r.integrase_id, []).append({
            "bakta_product": r.bakta_product,
            "aa_len": r.aa_len,
        })
    rows = []
    for r in ie.itertuples():
        recs = by.get(r.integrase_id)
        if not recs:
            rows.append({
                "integrase_id": r.integrase_id,
                "in_protein_table": 0,
                "module_minimal": 0,
                "module_pending": 0,
            })
            continue
        sc = score_ie(recs, str(r.char_class or ""), bool(r.from_rep))
        sc["integrase_id"] = r.integrase_id
        sc["in_protein_table"] = 1
        rows.append(sc)
    sc = pd.DataFrame(rows)
    out = ie.merge(sc, on="integrase_id", how="left")
    out.loc[out.n_char_classes != 1, ["module_minimal", "module_pending"]] = 0

    rm_plate = set(pd.read_csv(SMALL_RM, sep="\t").query("func_class=='unknown'").integrase_id) - DROP_RM
    ta_plate = set(pd.read_csv(SMALL_TA, sep="\t").query("func_class=='unknown'").integrase_id)
    out["on_small_rm_plate"] = out.integrase_id.isin(rm_plate).astype(int)
    out["on_small_ta_plate"] = out.integrase_id.isin(ta_plate).astype(int)
    # Small RM/TA plates win: MTase-only and PDC-as-extra stay in that class.
    out.loc[out.on_small_rm_plate == 1, "char_class"] = "RM"
    out.loc[out.on_small_ta_plate == 1, "char_class"] = "TA"

    int_n = out.n_int.fillna(0)
    rep_n = out.n_rep_cds.fillna(0)
    cds = out.n_cds.fillna(0)
    out["n_budget"] = (cds - int_n - rep_n).clip(lower=0)
    # Same budget as small RM/TA: extras ≤ 3; int and Rep free. budget 0 = int+Rep only.
    on_budget = (out.n_budget <= 3) & (out.in_protein_table == 1) & out.char_class.fillna("").ne("")
    one = out.n_char_classes == 1
    out["compact_minimal"] = (
        ((on_budget & one) | out.on_small_rm_plate.eq(1) | out.on_small_ta_plate.eq(1))
        & ~out.integrase_id.isin(DROP_RM)
    ).astype(int)
    # JAWRJC was force-proven as dark IIG; drop it
    out.loc[out.integrase_id.isin(DROP_RM), ["module_minimal", "compact_minimal"]] = 0

    dest = ART / "unknown_module_minimal.tsv"
    out.to_csv(dest, sep="\t", index=False)
    mm = out[out.module_minimal == 1]
    pend = out[out.module_pending == 1]
    cap = out[out.compact_minimal == 1]
    print("scored", int(out.in_protein_table.sum()), "/", len(out))
    print("module_minimal (named, no dark)", len(mm))
    print(mm.char_class.value_counts().to_string() if len(mm) else "  (none)")
    print("compact_minimal (budget≤3, same as small RM/TA)", len(cap))
    print(cap.char_class.value_counts().to_string() if len(cap) else "  (none)")
    print("  RM plate", int(cap.on_small_rm_plate.sum()), "TA plate", int(cap.on_small_ta_plate.sum()))
    print("wrote", dest)


if __name__ == "__main__":
    main()
