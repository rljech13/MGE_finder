#!/usr/bin/env python3
"""Scan the full unknown cohort (n=3,738) for R-M and T/AT cargo.

No size filter. Integrase / Rep are recorded but do not define the set.
R-M: PADLOC RM_* plus DNA-specific Bakta names (not generic methylases).
T/AT: PADLOC AbiE plus named / toxin-antitoxin Bakta products.
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import REPO, FF

import csv
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ART_RM = FF / "artifacts/unknown_rm"
ART_TA = FF / "artifacts/unknown_ta"
FUNC = FF / "artifacts/attL_complement_blast_window/cohort_6960_four_func_classes.tsv"
ATLAS = FF / "artifacts/attL_complement_blast_window/cohort_6960_atlas_join.tsv"
IDX = REPO / "results_organized_rep_6960/atlas_rep_index.tsv"
TAX = REPO / "results_organized/misc/gtdb_taxonomy.tsv"
REP_ROOT = REPO / "results_organized_rep_6960"
N_WORKERS = 8

INT_RE = re.compile(
    r"integrase|recombinase|xerc|xerd|tyrosine-type|tyrosine recombinase",
    re.I,
)
REP_RE = re.compile(
    r"\brep\b|replication initiation|replication initiator|initiator rep|"
    r"rep protein|repa\b|repc\b|repabc|primase|primpol|dnag|"
    r"bacteriophage replication gene a|gpa\b|plasmid replication",
    re.I,
)
HEL_REP_RE = re.compile(r"helicase rep|rep helicase|trw[c]|tr[ai] helicase", re.I)

# DNA R-M only. Generic / rRNA / protein methylases are excluded.
RM_EXCLUDE = re.compile(
    r"rrna|trna methyl|ribosom|protein methyl|histone|uroporphyrin|"
    r"cobalamin|chemotaxis|release factor|\bhemk\b|\bprm[abc]\b|"
    r"\bksga\b|\brsm[a-z]\b|\brrm[a-z]\b|\berm[abc]\b|\bmraw\b|"
    r"ubiquinone|menaquinone|biotin|rnam\b|antirestriction|\bard[ac]\b",
    re.I,
)
RM_RE = re.compile(
    r"restrict|rease|hsd[rms]\b|gmrsd|\bmrr\b|mmei|"
    r"type [ivx]+g?\s+restrict|restriction[- ]modification|"
    r"restriction endonuclease|site-specific dna-methyl|"
    r"dna methylase|dna methyltransferase|dna \(cytosine-5|"
    r"n-?6.?adenine|adenine-specific dna|cytosine-specific|"
    r"c-5 cytosine|modification methylase|n6.mtase|dam methyl|\bdcm\b",
    re.I,
)

TA_EXCLUDE = re.compile(
    r"rtx|botulinum|tetanus|cholera|shiga|leukotoxin|exotoxin|"
    r"endotoxin|enterotoxin|neurotoxin|cytotoxin|pertussis|"
    r"zona occludens|zonular occludens|\bzot\b|holin-like|"
    r"\btc toxin|metallopeptidase toxin|pre-toxin tg",
    re.I,
)
TA_RE = re.compile(
    r"toxin[- ]antitoxin|antitoxin|addiction module|"
    r"\bhic[ab]\b|\brel[eb]\b|\bpar[ed]\b|\bmaz[ef]\b|\bpem[ik]\b|"
    r"\bvap[bc]\b|\bccd[ab]\b|\bhig[ab]\b|\bmqs[ar]\b|"
    r"abiei|abieii|abigii|\bphd\b|\bdoc\b|\byoeb\b|\byefm\b|"
    r"\bbrn[ta]\b|\bhip[ab]\b|zeta toxin|putative toxin|"
    r"type ii toxin|\btoxin\b",
    re.I,
)
PAD_TA_RE = re.compile(r"^AbiE$", re.I)


def resolve_paths(gbff_rel: str, realm: str) -> tuple[str | None, str | None]:
    gbff = pad = None
    rel = str(gbff_rel or "")
    realm = str(realm or "")
    if not rel or rel == "nan":
        return None, None
    p = REP_ROOT / rel
    if p.is_file() and p.stat().st_size > 0:
        gbff = str(p)
    stem = Path(rel).name.replace(".annotation.gbff.merged", "")
    if realm and realm != "nan":
        p2 = REP_ROOT / realm / stem / "annotation.gbff.merged"
        if gbff is None and p2.is_file() and p2.stat().st_size > 0:
            gbff = str(p2)
        pc = REP_ROOT / realm / stem / "padloc.csv"
        if pc.is_file() and pc.stat().st_size > 0:
            pad = str(pc)
    if gbff:
        pc2 = Path(gbff).parent / "padloc.csv"
        if pad is None and pc2.is_file() and pc2.stat().st_size > 0:
            pad = str(pc2)
    return gbff, pad


def kind_of(prod: str) -> str:
    p = prod or ""
    if INT_RE.search(p):
        return "int"
    if not RM_EXCLUDE.search(p) and RM_RE.search(p):
        return "rm"
    if not TA_EXCLUDE.search(p) and TA_RE.search(p):
        return "ta"
    if REP_RE.search(p) or HEL_REP_RE.search(p):
        return "rep"
    return "other"


def ta_role(prod: str) -> str:
    p = (prod or "").lower().replace("toxin-antitoxin", "TA").replace("toxin antitoxin", "TA")
    if "abieii" in p or "abigii" in p:
        return "toxin"
    if "abiei" in p:
        return "antitoxin"
    if "antitoxin" in p or "addiction module" in p:
        return "antitoxin"
    if "toxin" in p or "vapc" in p or "mazf" in p or "pemk" in p:
        return "toxin"
    if TA_RE.search(prod or ""):
        return "ta"
    return ""


def parse_gbff(path: str) -> list[dict]:
    rec = SeqIO.read(path, "genbank")
    rows = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        prod = (f.qualifiers.get("product") or [""])[0]
        rows.append({
            "start": int(f.location.start) + 1,
            "end": int(f.location.end),
            "strand": f.location.strand or 1,
            "product": prod,
            "kind": kind_of(prod),
            "aa": (int(f.location.end) - int(f.location.start)) // 3,
            "role": ta_role(prod),
            "padloc_sys": "",
        })
    return rows


def read_padloc(path: str) -> list[tuple[int, int, str]]:
    out = []
    try:
        with open(path, newline="") as fh:
            rdr = csv.DictReader(fh)
            for r in rdr:
                sys = str(r.get("system") or "")
                if not sys:
                    continue
                is_rm = sys.startswith("RM_") or "RM_type" in sys
                is_ta = bool(PAD_TA_RE.match(sys))
                if not (is_rm or is_ta):
                    continue
                try:
                    a, b = int(r["start"]), int(r["end"])
                except (KeyError, TypeError, ValueError):
                    continue
                if a > b:
                    a, b = b, a
                out.append((a, b, sys))
    except Exception:
        return []
    return out


def overlay_padloc(genes: list[dict], spans: list[tuple[int, int, str]]) -> None:
    for g in genes:
        for a, b, sys in spans:
            ov = min(g["end"], b) - max(g["start"], a)
            if ov <= 0.5 * min(g["end"] - g["start"] + 1, b - a + 1):
                continue
            if sys.startswith("RM_") or "RM_type" in sys:
                g["kind"] = "rm"
                g["padloc_sys"] = sys if not g["padloc_sys"] else f"{g['padloc_sys']}|{sys}"
            elif PAD_TA_RE.match(sys):
                if g["kind"] != "rm":
                    g["kind"] = "ta"
                    if not g["role"]:
                        g["role"] = "ta"
                g["padloc_sys"] = sys if not g["padloc_sys"] else f"{g['padloc_sys']}|{sys}"


def padloc_rm_classes(systems: set[str]) -> set[str]:
    out = set()
    for s in systems:
        if s == "RM_type_I":
            out.add("Type I")
        elif s == "RM_type_II":
            out.add("Type II")
        elif s == "RM_type_IIG":
            out.add("Type IIG")
        elif s == "RM_type_III":
            out.add("Type III")
        elif s in {"RM_type_IV", "RM_type_HNH"}:
            out.add("Type IV")
    return out


def bakta_rm_classes(rm_genes: list[dict]) -> set[str]:
    if not rm_genes:
        return set()
    blob = " ".join(g["product"] for g in rm_genes).lower()
    max_aa = max(g["aa"] for g in rm_genes)
    has_gmr = "gmrsd" in blob
    has_mrr = "mrr" in blob or "type iv" in blob
    has_hsd = bool(re.search(r"hsd[rms]", blob))
    has_re = bool(re.search(r"restrict|rease", blob))
    has_mt = bool(re.search(r"methyl", blob))
    out = set()
    if has_gmr:
        out.add("GmrSD")
    if has_mrr:
        out.add("Type IV")
    if has_hsd or "type i restriction" in blob:
        out.add("Type I")
    if "type iii" in blob:
        out.add("Type III")
    if has_re and has_mt:
        out.add("Type IIG" if max_aa >= 800 else "Type II")
    elif has_mt and not has_re and not has_gmr and not has_mrr:
        out.add("Type IIG" if max_aa >= 800 else "MTase only")
    elif has_re and not has_mt and not has_gmr and not has_mrr:
        out.add("Type II")
    return out


def pick_rm_class(pad_cls: set[str], bak_cls: set[str]) -> str:
    cls = set(pad_cls)
    if not cls:
        cls = set(bak_cls)
    else:
        if "GmrSD" in bak_cls:
            cls.add("GmrSD")
    if not cls:
        return "other"
    if len(cls) == 1:
        return next(iter(cls))
    return "Mixed"


def rm_composition(rm_class: str, rm_genes: list[dict]) -> str:
    if rm_class == "Mixed":
        return "mixed classes"
    if rm_class == "MTase only":
        return "MTase only"
    if rm_class in {"Type IV", "GmrSD"}:
        return "REase only"
    if rm_class in {"Type I", "Type IIG", "Type III"}:
        return "M+R pair"
    blob = " ".join(g["product"] for g in rm_genes).lower()
    has_re = bool(re.search(r"restrict|rease", blob))
    has_mt = bool(re.search(r"methyl", blob))
    if has_re and has_mt:
        return "M+R pair"
    if has_mt:
        return "MTase only"
    if has_re:
        return "REase only"
    return "M+R pair"


def ta_family(blob: str) -> str:
    b = blob.lower()
    if "hic" in b:
        return "HicAB"
    if "hip" in b:
        return "HipAB"
    if "ccd" in b:
        return "CcdAB"
    if "phd" in b or re.search(r"\bdoc\b", b):
        return "Phd-Doc"
    if "mqs" in b:
        return "MqsRA"
    if "hig" in b:
        return "HigAB"
    if "yefm" in b or "yoeb" in b:
        return "YefM-YoeB"
    if "brn" in b:
        return "BrnTA"
    if "rel" in b or "pare" in b or "yafq" in b or "addiction module" in b:
        return "RelBE"
    if "maz" in b or "pemk" in b:
        return "MazEF"
    if "vapc" in b or "vapb" in b or "pin family" in b:
        return "VapBC"
    if "abie" in b or "abigii" in b:
        return "AbiE"
    return "unnamed"


def pair_status(roles: list[str], other: str, pad: str) -> str:
    has_t = "toxin" in roles
    has_a = "antitoxin" in roles
    oth = (other or "").lower()
    if has_t and has_a:
        return "named pair"
    if "AbiE" in (pad or "") and (has_t or has_a):
        return "named pair"
    if has_t and (
        "helix-turn-helix" in oth
        or "transcriptional regulator" in oth
        or "hth " in oth
        or "hypothetical protein(" in oth
    ):
        return "likely pair"
    if has_a and not has_t:
        return "antitoxin only"
    if has_t and not has_a:
        return "toxin only"
    return "toxin only"


def fmt_genes(genes: list[dict], kind: str) -> str:
    return "; ".join(
        f"{g['product']}({g['aa']}aa,{g['start']}-{g['end']})"
        for g in genes if g["kind"] == kind
    )


def scan_one(job: dict) -> dict:
    gbff, pad = resolve_paths(job["gbff_path"], job["source_realm"])
    rec = {
        "integrase_id": job["integrase_id"],
        "ie_id": job["ie_id"],
        "ie_len_nt": job["ie_len_nt"],
        "n_cds_atlas": job["n_cds"],
        "gbff_ok": bool(gbff),
        "padloc_ok": bool(pad),
        "n_cds_gbff": 0,
        "n_int": 0,
        "n_rep": 0,
        "n_rm": 0,
        "n_ta": 0,
        "n_other": 0,
        "padloc": "",
        "rm_products": "",
        "ta_products": "",
        "rep_products": "",
        "other_products": "",
        "rm_class": "",
        "rm_classes": "",
        "rm_composition": "",
        "ta_class": "",
        "ta_pair": "",
        "ta_roles": "",
        "has_rm": False,
        "has_ta": False,
        "genes": [],
    }
    if not gbff:
        return rec
    genes = parse_gbff(gbff)
    pad_systems: set[str] = set()
    if pad:
        spans = read_padloc(pad)
        overlay_padloc(genes, spans)
        pad_systems = {s for _a, _b, s in spans}
    n_int = sum(g["kind"] == "int" for g in genes)
    if n_int == 0 and genes:
        n_int = 1
        others = [g for g in genes if g["kind"] == "other"]
        if others:
            # unnamed finder integrase: leave counts as-is (int missing from names)
            pass
    rec["n_cds_gbff"] = len(genes)
    rec["n_int"] = n_int
    rec["n_rep"] = sum(g["kind"] == "rep" for g in genes)
    rec["n_rm"] = sum(g["kind"] == "rm" for g in genes)
    rec["n_ta"] = sum(g["kind"] == "ta" for g in genes)
    rec["n_other"] = sum(g["kind"] == "other" for g in genes)
    rec["padloc"] = "|".join(sorted(pad_systems))
    rec["rm_products"] = fmt_genes(genes, "rm")
    rec["ta_products"] = fmt_genes(genes, "ta")
    rec["rep_products"] = fmt_genes(genes, "rep")
    rec["other_products"] = fmt_genes(genes, "other")
    rec["has_rm"] = rec["n_rm"] > 0
    rec["has_ta"] = rec["n_ta"] > 0
    rec["genes"] = [
        {"integrase_id": job["integrase_id"], **g} for g in genes
        if g["kind"] in {"rm", "ta"}
    ]
    if rec["has_rm"]:
        rm_genes = [g for g in genes if g["kind"] == "rm"]
        pad_cls = padloc_rm_classes(pad_systems)
        bak_cls = bakta_rm_classes(rm_genes)
        rec["rm_class"] = pick_rm_class(pad_cls, bak_cls)
        rec["rm_classes"] = "|".join(sorted(pad_cls | bak_cls))
        rec["rm_composition"] = rm_composition(rec["rm_class"], rm_genes)
    if rec["has_ta"]:
        ta_genes = [g for g in genes if g["kind"] == "ta"]
        roles = [g["role"] or ta_role(g["product"]) for g in ta_genes]
        blob = f"{rec['ta_products']} {rec['padloc']}"
        rec["ta_class"] = ta_family(blob)
        rec["ta_pair"] = pair_status(roles, rec["other_products"], rec["padloc"])
        rec["ta_roles"] = ",".join(roles)
    return rec


def main() -> None:
    ART_RM.mkdir(parents=True, exist_ok=True)
    ART_TA.mkdir(parents=True, exist_ok=True)
    func = pd.read_csv(FUNC, sep="\t")
    atlas = pd.read_csv(ATLAS, sep="\t", usecols=["integrase_id", "sample_id", "n_cds"])
    idx = pd.read_csv(IDX, sep="\t")
    tax = pd.read_csv(TAX, sep="\t")
    unk = func[func.func_class == "unknown"].copy()
    unk = unk.merge(atlas, on="integrase_id", how="left")
    unk = unk.merge(idx[["integrase_id", "gbff_path", "source_realm"]], on="integrase_id", how="left")
    unk = unk.merge(tax, on="sample_id", how="left")
    jobs = []
    for r in unk.itertuples():
        jobs.append({
            "integrase_id": r.integrase_id,
            "ie_id": r.ie_id,
            "ie_len_nt": int(r.ie_len_nt) if pd.notna(r.ie_len_nt) else 0,
            "n_cds": r.n_cds if pd.notna(r.n_cds) else None,
            "gbff_path": getattr(r, "gbff_path", "") or "",
            "source_realm": getattr(r, "source_realm", "") or "",
        })
    print("unknown", len(jobs), "workers", N_WORKERS)
    rows = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for i, rec in enumerate(ex.map(scan_one, jobs, chunksize=40), 1):
            rows.append(rec)
            if i % 500 == 0:
                print(f"  {i}/{len(jobs)}")
    genes = []
    for rec in rows:
        genes.extend(rec.pop("genes"))
    out = pd.DataFrame(rows)
    meta = unk[[
        "integrase_id", "sample_id", "primary_label", "phylum", "class", "order", "family",
    ]].copy()
    if "phylum" not in meta.columns:
        meta["phylum"] = "Unknown"
    out = out.merge(meta, on="integrase_id", how="left")
    out["phylum"] = out["phylum"].fillna("Unknown")
    dest_all = ART_RM / "unknown_rm_ta_scan.tsv"
    keep = [c for c in out.columns]
    out.to_csv(dest_all, sep="\t", index=False)
    rm = out[out.has_rm].copy()
    ta = out[out.has_ta].copy()
    rm.to_csv(ART_RM / "unknown_rm_hits.tsv", sep="\t", index=False)
    ta.to_csv(ART_TA / "unknown_ta_hits.tsv", sep="\t", index=False)
    pd.DataFrame(genes).to_csv(ART_RM / "unknown_rm_ta_genes.tsv", sep="\t", index=False)
    print("parsed", int(out.gbff_ok.sum()), "missing_gbff", int((~out.gbff_ok).sum()))
    print("padloc.csv", int(out.padloc_ok.sum()))
    print("RM+", len(rm), f"/ {len(out)}")
    if len(rm):
        print(rm.rm_class.value_counts().to_string())
        print(rm.rm_composition.value_counts().to_string())
    print("TA+", len(ta), f"/ {len(out)}")
    if len(ta):
        print(ta.ta_class.value_counts().to_string())
        print(ta.ta_pair.value_counts().to_string())
    print("wrote", dest_all)


if __name__ == "__main__":
    main()
