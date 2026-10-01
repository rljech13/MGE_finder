#!/usr/bin/env python3
"""Class funnel + compositional-minimal unknown IEs (4689).

Minimal = integrase + exactly one cargo class from the locked IME funnel.
Replication may sit with that class and still counts as minimal.
Replication alone is a characterizing class.

Leftovers allowed on composition-minimal: hyp/DUF, HTH/Cro/Xre, Xis,
TPR / coiled-coil / Zn-ribbon. Named extra functions (transposase, Mob,
other enzymes) disqualify composition-minimal but not class-minimal.
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from unknown_metals_markers import is_metal_cargo, metal_subsystems

ART = FF / "artifacts/unknown_ime_classes"
V2 = FF / "artifacts/attL_complement_blast_window/cohort_6960_four_class_v2.tsv"
SCAN = FF / "artifacts/unknown_rm/unknown_rm_ta_scan.tsv"
PAD = FF / "artifacts/nb04_cluster_unclassified/padloc_per_ie_stage8_unc.tsv"
IDX = REPO / "results_organized_rep_6960/atlas_rep_index.tsv"
REP_ROOT = REPO / "results_organized_rep_6960"
ORPHAN_ARRAY = {"JAGPMA010000011.1_5"}

CARGO = [
    "replication",
    "TA",
    "Par",
    "RM",
    "CRISPR",
    "Abi",
    "RT-defense",
    "SoFic",
    "CBASS",
    "Shedu",
    "Thoeris",
    "BREX",
    "Wadjet",
    "Lamassu",
    "PDC",
    "Rst_DUF4238",
    "AMR",
    "metals",
    "bacteriocins",
    "other_defense",
]

INT_RE = re.compile(
    r"integrase|recombinase|xerc|xerd|tyrosine-type|tyrosine recombinase",
    re.I,
)
XIS_RE = re.compile(r"excisionase|\bxis\b|dna binding domain-containing protein, excisionase", re.I)
REP_RE = re.compile(
    r"\brep\b|replication initiation|replication initiator|initiator rep|"
    r"rep protein|repa\b|repc\b|repabc|primase|primpol|dnag|"
    r"bacteriophage replication gene a|\bgpa\b|plasmid replication|"
    r"rep_trans|replitron|repsa|dna primase|dna polymerase|"
    r"polymerase nucleotidyl|dna polymerase subunit|"
    r"helicase rep|rep helicase",
    re.I,
)
HEL_REP_RE = re.compile(
    r"superfamily ii dna helicase|dna helicase|helicase.*?aaa|"
    r"aaa family atpase",
    re.I,
)
RM_RE = re.compile(
    r"restrict|rease|hsd[rms]\b|gmrsd|\bmrr\b|mmei|"
    r"type [ivx]+g?\s+restrict|restriction[- ]modification|"
    r"restriction endonuclease|site-specific dna-methyl|"
    r"dna methylase|dna methyltransferase|dna \(cytosine-5|"
    r"n-?6.?adenine|adenine-specific dna|cytosine-specific|"
    r"c-5 cytosine|modification methylase|n6.mtase|dam methyl|\bdcm\b|"
    r"type ii restriction|type i restriction|type iv restriction",
    re.I,
)
RM_EXCLUDE = re.compile(
    r"rrna|trna methyl|ribosom|protein methyl|histone|uroporphyrin|"
    r"cobalamin|chemotaxis|release factor|\bhemk\b|\bprm[abc]\b|"
    r"\bksga\b|\brsm[a-z]\b|\brrm[a-z]\b|\berm[abc]\b|\bmraw\b|"
    r"ubiquinone|menaquinone|biotin|rnam\b|antirestriction|\bard[ac]\b",
    re.I,
)
RM_REASE = re.compile(
    r"restrict|rease|hsd[rms]\b|gmrsd|\bmrr\b|type [ivx]+g?\s+restrict|"
    r"restriction endonuclease|type ii restriction|type i restriction|"
    r"type iv restriction",
    re.I,
)
TA_RE = re.compile(
    r"toxin[- ]antitoxin|antitoxin|addiction module|"
    r"\bhic[ab]\b|\brel[eb]\b|\bpar[ed]\b|\bmaz[ef]\b|\bpem[ik]\b|"
    r"\bvap[bc]\b|\bccd[ab]\b|\bhig[ab]\b|\bmqs[ar]\b|"
    r"abiei|abieii|abigii|\bphd\b|\bdoc\b|\byoeb\b|\byefm\b|"
    r"\bbrn[ta]\b|\bhip[ab]\b|zeta toxin|type ii toxin|"
    r"prevent-host-death|plasmid stabilization protein",
    re.I,
)
TA_EXCLUDE = re.compile(
    r"rtx|botulinum|tetanus|cholera|shiga|leukotoxin|exotoxin|"
    r"endotoxin|enterotoxin|neurotoxin|cytotoxin|pertussis|"
    r"zona occludens|zonular occludens|\bzot\b|holin-like|"
    r"\btc toxin|metallopeptidase toxin|pre-toxin tg",
    re.I,
)
PAR_RE = re.compile(
    r"chromosome partitioning|\bplasmid partition|\bparb\b|\bpara\b|"
    r"\bsoj\b|partitioning protein par|sporulation initiation inhibitor protein soj",
    re.I,
)
PAR_EXCLUDE = re.compile(r"\bpar[ed]\b|pare family|pard family", re.I)
CRISPR_RE = re.compile(
    r"crispr|\bcas[1-9]\b|\bcas1[0-4]\b|\bcas[a-z]\b|"
    r"crispr-associated|\bcmr[0-9]\b|\bcsm[0-9]\b|\bcse[0-9]\b",
    re.I,
)
ABI_RE = re.compile(r"\babil\b|\babid\b|\babiu\b|\babiz\b|abortive infection protein", re.I)
RTDEF_RE = re.compile(
    r"\bdrt\b|retron|reverse transcriptase|group ii intron|"
    r"\babia\b|\babik\b|\babip2\b",
    re.I,
)
SOFIC_RE = re.compile(r"\bsofic\b", re.I)
CBASS_RE = re.compile(r"\bcbass\b|cyclic oligonucleotide|cgas", re.I)
SHEDU_RE = re.compile(r"\bshedu\b", re.I)
THOERIS_RE = re.compile(r"thoeris|\bthsa\b|\bthsb\b", re.I)
BREX_RE = re.compile(r"\bbrex\b", re.I)
WADJET_RE = re.compile(r"\bwadjet\b", re.I)
LAMASSU_RE = re.compile(r"lamassu", re.I)
RST_RE = re.compile(r"duf4238", re.I)
AMR_RE = re.compile(
    r"(?<!metallo-)beta-lactamase(?! domain)|"
    r"chloramphenicol acetyl|tetracycline resistance|"
    r"tetracycline[\s-]efflux|aminoglycoside phospho|"
    r"aminoglycoside resistance|\blnu\(|"
    r"\btet\([a-z]|\berm[abc]\b|"
    r"vancomycin resistance|\bvana\b|\bvanb\b|"
    r"macrolide.*resistance|lincosamide resistance|"
    r"aac\([0-9]|aph\([0-9]",
    re.I,
)
# Metals: unknown_metals_markers.is_metal_cargo (aliases + operon partners;
# regulators and generic Cu/Zn domains stay out).
BACT_RE = re.compile(
    r"bacteriocin|circularin|uberolysin|lantibiotic|microcin|colicin|enterocin",
    re.I,
)
OK_ADJ_RE = re.compile(
    r"hypothetical protein|uncharacterized|conserved domain protein|"
    r"conserved hypothetical|domain of unknown function|duf\d+|"
    r"helix-turn-helix|hth |cro/c1|\bxre\b|transcriptional regulator|"
    r"dna-binding protein|dna binding domain|"
    r"coiled-coil|tetratricopeptide|\btpr\b|zinc.?ribbon|zinc finger|"
    r"excisionase|\bxis\b",
    re.I,
)
TNP_RE = re.compile(
    r"transposase|\btnp[a-z]?\b|insertion sequence|\bis[0-9]|"
    r"invertase|resolvase|dna-invertase",
    re.I,
)
MOB_RE = re.compile(r"relaxase|mobilization protein|mobilisation|\bmob[abcqtv]\b|\bvird2\b", re.I)
BLOB_RE = re.compile(r"^(.+)\((\d+)aa,(\d+)-(\d+)\)$")

PAD_MAP = [
    (re.compile(r"^RM_type|^RM_", re.I), "RM"),
    (re.compile(r"^AbiE$|RosmerTA|ShosTA|^ietAS$", re.I), "TA"),
    (re.compile(r"^Abi[LDUZ]$", re.I), "Abi"),
    (re.compile(r"^Abi[AK]$|^AbiP2$|^DRT|^retron", re.I), "RT-defense"),
    (re.compile(r"^SoFic$", re.I), "SoFic"),
    (re.compile(r"^cbass", re.I), "CBASS"),
    (re.compile(r"^shedu$", re.I), "Shedu"),
    (re.compile(r"^Thoeris", re.I), "Thoeris"),
    (re.compile(r"^BREX", re.I), "BREX"),
    (re.compile(r"^Wadjet", re.I), "Wadjet"),
    (re.compile(r"^Lamassu", re.I), "Lamassu"),
    (re.compile(r"^CRISPR|^cas", re.I), "CRISPR"),
    (re.compile(r"^PDC-|^HEC-|^PD-T4|^PD-Lambda", re.I), "PDC"),
    (re.compile(r"DUF4238|Rst_DUF", re.I), "Rst_DUF4238"),
    (re.compile(r"^VSPR$|^DMS_other$", re.I), None),
]


def padloc_class(tok: str) -> str | None:
    t = tok.strip()
    if not t:
        return None
    for rx, cls in PAD_MAP:
        if rx.search(t):
            return cls
    if re.search(r"septu|gabija|hachiman|zorya|disarm|AVAST|SEFIR|viperin|argonaute|PT_Ssp", t, re.I):
        return "other_defense"
    return "other_defense"


def explode_padloc(s: str) -> list[str]:
    if not s or pd.isna(s):
        return []
    out = []
    for part in re.split(r"[|;,\s]+", str(s)):
        if part:
            out.append(part)
    return out


def parse_blob(blob: str) -> list[str]:
    if not blob or pd.isna(blob):
        return []
    prods = []
    for chunk in str(blob).split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = BLOB_RE.match(chunk)
        prods.append(m.group(1).strip() if m else chunk)
    return prods


def gene_tag(prod: str, *, allow_aaa_as_rep: bool = False) -> str:
    p = prod or ""
    if INT_RE.search(p) and not RM_REASE.search(p):
        return "int"
    if XIS_RE.search(p):
        return "ok_adj"
    if RST_RE.search(p):
        return "Rst_DUF4238"
    if BACT_RE.search(p):
        return "bacteriocins"
    if is_metal_cargo(p):
        return "metals"
    if AMR_RE.search(p):
        return "AMR"
    if CRISPR_RE.search(p):
        return "CRISPR"
    if PAR_RE.search(p) and not PAR_EXCLUDE.search(p):
        return "Par"
    if ABI_RE.search(p):
        return "Abi"
    if not TA_EXCLUDE.search(p) and TA_RE.search(p):
        return "TA"
    if THOERIS_RE.search(p):
        return "Thoeris"
    if SHEDU_RE.search(p):
        return "Shedu"
    if CBASS_RE.search(p):
        return "CBASS"
    if BREX_RE.search(p):
        return "BREX"
    if WADJET_RE.search(p):
        return "Wadjet"
    if LAMASSU_RE.search(p):
        return "Lamassu"
    if SOFIC_RE.search(p):
        return "SoFic"
    if RTDEF_RE.search(p):
        return "RT-defense"
    if not RM_EXCLUDE.search(p) and RM_RE.search(p):
        return "rm_part"
    if REP_RE.search(p) or (allow_aaa_as_rep and HEL_REP_RE.search(p)):
        return "replication"
    if OK_ADJ_RE.search(p):
        return "ok_adj"
    if TNP_RE.search(p):
        return "leftover"
    if MOB_RE.search(p):
        return "leftover"
    return "leftover"


def resolve_gbff(gbff_rel: str, realm: str) -> str | None:
    rel = str(gbff_rel or "")
    realm = str(realm or "")
    if not rel or rel == "nan":
        return None
    p = REP_ROOT / rel
    if p.is_file() and p.stat().st_size > 0:
        return str(p)
    stem = Path(rel).name.replace(".annotation.gbff.merged", "")
    if realm and realm != "nan":
        p2 = REP_ROOT / realm / stem / "annotation.gbff.merged"
        if p2.is_file() and p2.stat().st_size > 0:
            return str(p2)
    return None


def parse_gbff_products(path: str) -> list[str]:
    rec = SeqIO.read(path, "genbank")
    out = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        out.append((f.qualifiers.get("product") or [""])[0])
    return out


def _parse_one(job: dict) -> tuple[str, list[str], bool]:
    path = resolve_gbff(job["gbff_path"], job["source_realm"])
    if not path:
        return job["integrase_id"], [], False
    try:
        return job["integrase_id"], parse_gbff_products(path), True
    except Exception:
        return job["integrase_id"], [], False


def ie_classes(tags: list[str], pad_cls: set[str], has_rep_label: bool, iid: str) -> set[str]:
    cls: set[str] = set(pad_cls)
    rm_parts = [t for t in tags if t == "rm_part"]
    if "RM" in cls:
        pass
    elif rm_parts:
        # MTase-only is not RM unless PADLOC already called it
        blob_has_re = any(True for t in tags if t == "rm_part")
        # decide from products stored alongside? we only have tags here.
        # Conservative: require a restriction-like tag count via later pass.
        cls.add("RM")  # refined below
    for t in tags:
        if t in CARGO:
            cls.add(t)
    if has_rep_label or "replication" in tags:
        cls.add("replication")
    if iid in ORPHAN_ARRAY:
        cls.add("CRISPR")
    return cls


def refine_rm(cls: set[str], products: list[str], pad_cls: set[str]) -> set[str]:
    if "RM" in pad_cls:
        cls.add("RM")
        return cls
    has_re = any(not RM_EXCLUDE.search(p) and RM_REASE.search(p) for p in products)
    has_mt = any(not RM_EXCLUDE.search(p) and RM_RE.search(p) and not RM_REASE.search(p) for p in products)
    if has_re:
        cls.add("RM")
    elif "RM" in cls and not has_re and has_mt:
        cls.discard("RM")
    return cls


def main() -> None:
    ART.mkdir(parents=True, exist_ok=True)
    v2 = pd.read_csv(V2, sep="\t")
    unk = v2[v2.func_class == "unknown"].copy()
    scan = pd.read_csv(SCAN, sep="\t")
    pad = pd.read_csv(PAD, sep="\t")
    idx = pd.read_csv(IDX, sep="\t", usecols=["integrase_id", "gbff_path", "source_realm"])

    pad_by = {}
    for r in pad.itertuples():
        toks = explode_padloc(r.padloc_systems)
        pad_by[r.integrase_id] = {c for t in toks if (c := padloc_class(t))}

    scan_prods: dict[str, list[str]] = {}
    scan_flags: dict[str, dict] = {}
    for r in scan.itertuples():
        prods = []
        for col in ("rm_products", "ta_products", "rep_products", "other_products"):
            prods.extend(parse_blob(getattr(r, col, "")))
        scan_prods[r.integrase_id] = prods
        scan_flags[r.integrase_id] = {
            "has_rm": bool(r.has_rm) and str(getattr(r, "rm_composition", "") or "") != "MTase only",
            "has_ta": bool(r.has_ta),
        }

    missing = unk[~unk.integrase_id.isin(scan_prods)].merge(idx, on="integrase_id", how="left")
    jobs = [
        {
            "integrase_id": r.integrase_id,
            "gbff_path": getattr(r, "gbff_path", "") or "",
            "source_realm": getattr(r, "source_realm", "") or "",
        }
        for r in missing.itertuples()
    ]
    print("parse GBFF for", len(jobs), "rep.* IEs")
    extra: dict[str, list[str]] = {}
    extra_ok = 0
    if jobs:
        with ProcessPoolExecutor(max_workers=8) as ex:
            for i, (iid, prods, ok) in enumerate(ex.map(_parse_one, jobs, chunksize=20), 1):
                extra[iid] = prods
                extra_ok += int(ok)
                if i % 200 == 0:
                    print(f"  {i}/{len(jobs)}")
    print("gbff_ok", extra_ok, "/", len(jobs))

    rows = []
    metal_rows = []
    for r in unk.itertuples():
        iid = r.integrase_id
        prods = scan_prods.get(iid)
        from_scan = prods is not None
        if prods is None:
            prods = extra.get(iid, [])
        has_rep_label = str(r.primary_label).startswith("rep.")
        allow_aaa = has_rep_label
        tags = [gene_tag(p, allow_aaa_as_rep=allow_aaa) for p in prods]
        pad_cls = pad_by.get(iid, set())
        cls = set(pad_cls)
        for t in tags:
            if t in CARGO:
                cls.add(t)
        if has_rep_label or any(t == "replication" for t in tags):
            cls.add("replication")
        if iid in ORPHAN_ARRAY:
            cls.add("CRISPR")
        flags = scan_flags.get(iid, {})
        if flags.get("has_rm"):
            cls.add("RM")
        if flags.get("has_ta"):
            cls.add("TA")
        cls = refine_rm(cls, prods, pad_cls)
        cls.discard("other_defense")

        n_leftover = sum(t == "leftover" for t in tags)
        n_ok = sum(t == "ok_adj" for t in tags)
        leftover_ex = [p for p, t in zip(prods, tags) if t == "leftover"][:8]
        metal_prods = [p for p in prods if is_metal_cargo(p)]
        subs: list[str] = []
        for p in metal_prods:
            for s in metal_subsystems(p):
                if s not in subs:
                    subs.append(s)
        if metal_prods:
            metal_rows.append({
                "integrase_id": iid,
                "n_metal_genes": len(metal_prods),
                "subsystems": "|".join(subs),
                "products": "; ".join(metal_prods[:12]),
            })

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
        n_cds = r.n_cds if pd.notna(r.n_cds) else (len(prods) + 1 if prods else None)

        rows.append({
            "integrase_id": iid,
            "ie_id": r.ie_id,
            "primary_label": r.primary_label,
            "from_rep": int(has_rep_label),
            "from_scan": int(from_scan),
            "ie_len_nt": r.ie_len_nt,
            "n_cds": n_cds,
            "n_products_parsed": len(prods),
            "classes": "|".join(sorted(cls)),
            "n_classes": len(cls),
            "char_class": char,
            "n_char_classes": n_char,
            "has_rep": int(has_rep),
            "n_ok_adj": n_ok,
            "n_leftover": n_leftover,
            "leftover_examples": "; ".join(leftover_ex),
            "class_minimal": int(class_min),
            "composition_minimal": int(comp_min),
            "compact_le6": int(class_min and n_cds is not None and n_cds <= 6),
            "comp_min_le6": int(comp_min and n_cds is not None and n_cds <= 6),
            "tight_minimal": int(comp_min and n_ok <= 3),
        })

    ie = pd.DataFrame(rows)
    ie.to_csv(ART / "unknown_minimal_ie.tsv", sep="\t", index=False)
    pd.DataFrame(metal_rows).to_csv(ART / "metals_calls.tsv", sep="\t", index=False)

    n = len(ie)
    n_class = int(ie.class_minimal.sum())
    n_comp = int(ie.composition_minimal.sum())
    n_comp_rep = int(((ie.composition_minimal == 1) & (ie.has_rep == 1) & (ie.char_class != "replication")).sum())
    n_comp_norep = int(((ie.composition_minimal == 1) & ((ie.has_rep == 0) | (ie.char_class == "replication"))).sum())

    # layer, class, genes_included, not_included, literature, presence
    funnel_meta = [
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
        ("B", "Shedu",
         "SduA",
         "", "DefenseFinder", "yes"),
        ("B", "Thoeris",
         "ThsA (NAD hydrolase); ThsB (TIR)",
         "", "Botelho 2023", "rare"),
        ("B", "BREX",
         "BrxA; BrxB; BrxC; PglX; PglZ; BrxL",
         "", "DefenseFinder", "rare"),
        ("B", "Wadjet",
         "JetA; JetB; JetC; JetD",
         "", "DefenseFinder", "rare"),
        ("B", "Lamassu",
         "LmuA; LmuB; LmuC",
         "", "DefenseFinder", "rare"),
        ("B", "PDC",
         "PADLOC PDC-S02, S15, S70, M34 and other PDC/HEC HMMs",
         "not a biological family", "PADLOC", "yes"),
        ("B", "Rst_DUF4238",
         "DUF4238",
         "", "Rousset 2022", "yes"),
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

    by_rows = []
    for layer, name, genes, not_inc, src, presence in funnel_meta:
        if name == "virulence":
            n_has = 0
            n_cm = n_cmin = n_plusrep = n_le6 = n_tight = 0
        else:
            # presence on IE
            has = ie.classes.fillna("").apply(lambda s, n=name: n in set(s.split("|")) if s else False)
            if name in {"Gabija", "Septu", "Hachiman", "Zorya", "DISARM"}:
                has = ie.classes.fillna("").str.contains(name, case=False, regex=False) | (
                    ie.classes.fillna("").apply(lambda s: "other_defense" in s.split("|") if s else False)
                    & False
                )
                # also search leftover/pad via classes other_defense is too broad; check token in classes only if we mapped them
            n_has = int(has.sum())
            sub = ie[ie.char_class == name]
            n_cm = int(sub.class_minimal.sum())
            n_cmin = int(sub.composition_minimal.sum())
            n_plusrep = int(((sub.composition_minimal == 1) & (sub.has_rep == 1) & (sub.char_class != "replication")).sum())
            n_le6 = int(sub.comp_min_le6.sum())
            n_tight = int(sub.tight_minimal.sum()) if "tight_minimal" in sub.columns else 0
        by_rows.append({
            "layer": layer,
            "class": name,
            "genes_included": genes,
            "not_included": not_inc,
            "literature": src,
            "presence_lock": presence,
            "n_ie_with_class": n_has,
            "n_class_minimal": n_cm,
            "n_composition_minimal": n_cmin,
            "n_comp_min_plus_rep": n_plusrep,
            "n_tight_minimal_okadj_le3": n_tight,
            "n_comp_min_le6_cds": n_le6,
        })

    # fix empty DF families: search padloc-mapped other_defense is not the same.
    # Recompute Gabija family from raw pad if needed — already mapped to other_defense.
    # Leave n_ie_with_class 0 unless we find the token in classes (we map them to other_defense).
    # Add a dedicated recount from pad tokens stored? We didn't store raw tokens.
    # Re-join pad for those five names.
    pad_toks = {}
    for r in pad.itertuples():
        pad_toks[r.integrase_id] = " ".join(explode_padloc(r.padloc_systems)).lower()
    for name in ("Gabija", "Septu", "Hachiman", "Zorya", "DISARM"):
        key = name.lower()
        n_has = sum(1 for iid in ie.integrase_id if key in pad_toks.get(iid, ""))
        for row in by_rows:
            if row["class"] == name:
                row["n_ie_with_class"] = n_has
                row["presence_lock"] = "absent" if n_has == 0 else "yes"

    by = pd.DataFrame(by_rows)
    by = by[[
        "class", "genes_included", "not_included", "layer",
        "n_ie_with_class", "n_composition_minimal", "n_comp_min_plus_rep",
        "n_class_minimal", "n_tight_minimal_okadj_le3", "n_comp_min_le6_cds",
        "literature", "presence_lock",
    ]]
    by.to_csv(ART / "ime_class_funnel.tsv", sep="\t", index=False)
    by.to_csv(ART / "unknown_minimal_by_class.tsv", sep="\t", index=False)

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

    print(summary.to_string(index=False))
    print("\nby class (composition-minimal > 0 or n_has>0):")
    print(by.to_string(index=False))
    print("wrote", ART)


if __name__ == "__main__":
    main()
