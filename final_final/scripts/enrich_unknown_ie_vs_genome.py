#!/usr/bin/env python3
"""Size-normalized enrichment of funnel hallmarks: unknown IE vs rest of genome.

Same compact HMM set on Prodigal orfs.faa for both compartments (IE Bakta
names are not used). Cohort = unknown 4689 on genomes that have finder orfs.
Rest includes other atlas MGEs unless the chrom-only contrast is used.

Stages: hmm, intervals, faa, search, summarize, plot.
"""
from __future__ import annotations

import sys
from pathlib import Path as _MgePath
sys.path.insert(0, str(_MgePath(__file__).resolve().parent))
from mge_paths import CAS, DF, FF, HMMFETCH, HMMSEARCH, HMMPRESS, HMMSTAT, PADLOC, PFAM, REPO

import argparse
import math
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from scipy.stats import fisher_exact

ART = FF / "artifacts/unknown_ime_classes/after_hypothesis"
OUT = ART / "genome_enrichment"
FIG = FF / "figures/unknown_ime_classes"
IE_TSV = ART / "unknown_minimal_ie.tsv"
ATLAS = REPO / "results_organized_rep_6960/atlas_rep_index.tsv"
FINDER = REPO / "finder_pipeline/results"

E_CUT = 1e-5
BIT_CUT = 22.0
BATCH_N = 220
N_SEARCH_WORKERS = 4
CPU_PER_SEARCH = 12

BLUE = "#0072B2"
GREEN = "#009E73"
VERM = "#D55E00"
GREY = "#7A7A7A"
INK = "#222222"

LAYER = {
    "replication": "A",
    "TA": "A",
    "Par": "A",
    "RM": "B",
    "CRISPR": "B",
    "Abi": "B",
    "RT-defense": "B",
    "SoFic": "B",
    "CBASS": "B",
    "Shedu": "B",
    "Thoeris": "B",
    "BREX": "B",
    "Wadjet": "B",
    "Lamassu": "B",
    "PDC": "B",
    "Rst_DUF4238": "B",
    "Gabija": "B",
    "Septu": "B",
    "Hachiman": "B",
    "Zorya": "B",
    "DISARM": "B",
    "AMR": "C",
    "metals": "C",
    "bacteriocins": "C",
    "control_PolA": "ctrl",
    "control_DnaG": "ctrl",
}

ORDER = [
    "replication", "TA", "Par",
    "RM", "CRISPR", "Abi", "RT-defense", "SoFic", "CBASS", "Shedu",
    "Thoeris", "BREX", "Wadjet", "Lamassu", "PDC", "Rst_DUF4238",
    "Gabija", "Septu", "Hachiman", "Zorya", "DISARM",
    "AMR", "metals", "bacteriocins",
    "control_PolA", "control_DnaG",
]

PDC_TYPES = [
    "PDC-S01", "PDC-S02", "PDC-S03", "PDC-S04", "PDC-S05", "PDC-S06",
    "PDC-S09", "PDC-S11", "PDC-S13", "PDC-S14", "PDC-S15", "PDC-S16",
    "PDC-S18", "PDC-S19", "PDC-S20", "PDC-S23", "PDC-S24", "PDC-S26",
    "PDC-S27", "PDC-S30", "PDC-S31", "PDC-S32", "PDC-S33", "PDC-S37",
    "PDC-S39", "PDC-S46", "PDC-S47", "PDC-S54", "PDC-S57", "PDC-S60",
    "PDC-S62", "PDC-S63", "PDC-S65", "PDC-S70", "PDC-S71", "PDC-S72",
    "PDC-S73",
    "PDC-M03", "PDC-M04", "PDC-M05", "PDC-M06", "PDC-M11", "PDC-M12",
    "PDC-M13", "PDC-M19", "PDC-M21", "PDC-M30", "PDC-M31", "PDC-M32",
    "PDC-M33", "PDC-M34", "PDC-M37", "PDC-M41", "PDC-M42", "PDC-M45",
    "HEC-02", "HEC-05", "HEC-06", "HEC-09",
    "PD-T4-6", "PD-T7-3", "PD-T7-4", "PD-Lambda-1", "PD-Lambda-5",
]

DF_FILES = [
    # RM — REase / fused / Type IV only (MTase-only is not RM)
    ("RM__Type_I_REases_FAM_0.einsi_trimmed.hmm", "RM", "rease"),
    ("RM__Type_I_REases_FAM_1.einsi_trimmed.hmm", "RM", "rease"),
    ("RM__Type_I_REases_FAM_2.einsi_trimmed.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_1.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_2.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_3.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_4.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_5.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_FAM_0.einsi_trimmed.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_FAM_1.einsi_trimmed.hmm", "RM", "rease"),
    ("RM_Type_IIG__Type_IIG_FAM_2.einsi_trimmed.hmm", "RM", "rease"),
    ("RM_Type_III__Type_III_REases_FAM_0.einsi_trimmed.hmm", "RM", "rease"),
    ("RM_Type_III__Type_III_REases_FAM_1.einsi_trimmed.hmm", "RM", "rease"),
    ("RM_Type_IV__FAM_0.hmm", "RM", "rease"),
    ("RM_Type_IV__FAM_1.hmm", "RM", "rease"),
    ("RM_Type_IV__FAM_2.hmm", "RM", "rease"),
    ("RM_Type_IV__Type_IV_01.hmm", "RM", "rease"),
    ("RM_Type_IV__Type_IV_03.hmm", "RM", "rease"),
    ("RM_Type_IV__Type_IV_05.hmm", "RM", "rease"),
    ("RM_Type_IV__Type_IV_21.hmm", "RM", "rease"),
    ("RM_Type_IV__Type_IV_22.hmm", "RM", "rease"),
    *[(f"RM_Type_II__Type_II_REase{i:02d}.hmm", "RM", "rease") for i in range(1, 39) if i != 36],
    # TA
    ("mEp506_HicAB__HicA.hmm", "TA", "core"),
    ("mEp506_HicAB__HicB.hmm", "TA", "core"),
    ("MazEF__MazE.hmm", "TA", "core"),
    ("MazEF__MazF.hmm", "TA", "core"),
    ("AbiEii__AbiEii.hmm", "TA", "core"),
    ("AbiEii__AbiEi_1.hmm", "TA", "core"),
    ("AbiEii__AbiEi_2.hmm", "TA", "core"),
    ("AbiEii__AbiEi_3.hmm", "TA", "core"),
    ("AbiEii__AbiEi_4.hmm", "TA", "core"),
    # residual Abi
    ("AbiL__AbiLi.hmm", "Abi", "core"),
    ("AbiL__AbiLii.hmm", "Abi", "core"),
    ("AbiD__AbiD.hmm", "Abi", "core"),
    ("AbiU__AbiU.hmm", "Abi", "core"),
    ("AbiZ__AbiZ.hmm", "Abi", "core"),
    # RT-defense cores
    ("AbiA_large__AbiA_large.hmm", "RT-defense", "core"),
    ("AbiK__AbiK.hmm", "RT-defense", "core"),
    ("AbiP2__AbiP2.hmm", "RT-defense", "core"),
    ("DRT2__DRT2.hmm", "RT-defense", "core"),
    ("DRT6__DRT6.hmm", "RT-defense", "core"),
    ("DRT7__DRT7.hmm", "RT-defense", "core"),
    ("DRT8__DRT8.hmm", "RT-defense", "core"),
    ("DRT9__DRT9.hmm", "RT-defense", "core"),
    ("DRT_1__drt1a.hmm", "RT-defense", "core"),
    ("DRT_2__drt2.hmm", "RT-defense", "core"),
    ("DRT_3__drt3a.hmm", "RT-defense", "core"),
    ("DRT_4__drt4.hmm", "RT-defense", "core"),
    ("DRT_5__drt5.hmm", "RT-defense", "core"),
    *[(f"Retron__RT_Tot_{i}.hmm", "RT-defense", "core") for i in range(1, 16)],
    ("SoFic__SoFic.hmm", "SoFic", "core"),
    ("CBASS__A_Cyclase_AGS_C.hmm", "CBASS", "core"),
    ("CBASS__Cyclase_II.hmm", "CBASS", "core"),
    ("CBASS__Cyclase_SMODS.hmm", "CBASS", "core"),
    ("CBASS__Cyclase_new.hmm", "CBASS", "core"),
    ("CBASS__Effector_TIR.hmm", "CBASS", "core"),
    ("CBASS__SAVED_Cap4.hmm", "CBASS", "core"),
    ("CBASS__Jab.hmm", "CBASS", "core"),
    ("Shedu__SduA.hmm", "Shedu", "core"),
    ("Thoeris__ThsA_SIRim.hmm", "Thoeris", "core"),
    ("Thoeris__ThsA_Caps4.hmm", "Thoeris", "core"),
    ("Thoeris__ThsB_1.hmm", "Thoeris", "core"),
    ("BREX__brxA_DUF1819.hmm", "BREX", "core"),
    ("BREX__brxB_DUF1788.hmm", "BREX", "core"),
    ("BREX__brxC.hmm", "BREX", "core"),
    ("BREX__brxL.hmm", "BREX", "core"),
    ("BREX__pglX1.hmm", "BREX", "core"),
    ("BREX__pglZA.hmm", "BREX", "core"),
    ("Wadjet__JetA_I.hmm", "Wadjet", "core"),
    ("Wadjet__JetB_I.hmm", "Wadjet", "core"),
    ("Wadjet__JetC_I.hmm", "Wadjet", "core"),
    ("Wadjet__JetD_I.hmm", "Wadjet", "core"),
    ("Lamassu__LmuA_HNH.hmm", "Lamassu", "core"),
    ("Lamassu__LmuB_Long.hmm", "Lamassu", "core"),
    ("Lamassu__LmuB_Short.hmm", "Lamassu", "core"),
    ("Lamassu__LmuC_Clade_I.hmm", "Lamassu", "core"),
    ("Rst_DUF4238__DUF4238_Pers.hmm", "Rst_DUF4238", "core"),
    ("Gabija__GajA.hmm", "Gabija", "core"),
    ("Gabija__GajB_1.hmm", "Gabija", "core"),
    ("Septu__PtuA.hmm", "Septu", "core"),
    ("Septu__PtuB.hmm", "Septu", "core"),
    ("Hachiman__HamA_1.hmm", "Hachiman", "core"),
    ("Hachiman__HamB.hmm", "Hachiman", "core"),
    ("Zorya__ZorA.hmm", "Zorya", "core"),
    ("Zorya__ZorB.hmm", "Zorya", "core"),
    ("DISARM__drmA.hmm", "DISARM", "core"),
    ("DISARM__drmB.hmm", "DISARM", "core"),
    ("DISARM__drmC.hmm", "DISARM", "core"),
]

CAS_FILES = [
    ("cas1_I_II_III_IV_V_VI_1.hmm", "CRISPR", "core"),
    ("cas2_I_II_III_IV_V_VI_1.hmm", "CRISPR", "core"),
    ("cas3_I_2.hmm", "CRISPR", "core"),
    ("cas9_II-A_II-B_II-C_1.hmm", "CRISPR", "core"),
    ("cas10_III_1.hmm", "CRISPR", "core"),
]

PFAM_FETCH = [
    ("Rep3_N", "replication", "core"),
    ("Rep3_C", "replication", "core"),
    ("Rep-A_N", "replication", "core"),
    ("RepB_primase", "replication", "core"),
    ("RepB_primase_C", "replication", "core"),
    ("Replitron_C", "replication", "core"),
    ("Bac_RepA_C", "replication", "core"),
    ("AEP-TOTE", "replication", "core"),
    ("PriCT_1", "replication", "core"),
    ("PriCT_2", "replication", "core"),
    ("RepB-RCR_reg", "replication", "core"),
    ("RepD-like_N", "replication", "core"),
    ("ParB", "Par", "core"),
    ("HTH_ParB", "Par", "core"),
    ("ParB_C", "Par", "core"),
    ("RelE", "TA", "core"),
    ("PemK_toxin", "TA", "core"),
    ("HicB", "TA", "core"),
    ("GmrSD_N", "RM", "rease"),
    ("GmrSD_C", "RM", "rease"),
    ("Cas3_HD", "CRISPR", "core"),
    ("Cas9_RuvC", "CRISPR", "core"),
    ("Cas10-Cmr2_palm2", "CRISPR", "core"),
    ("Beta-lactamase", "AMR", "core"),
    ("Beta-lactamase2", "AMR", "core"),
    ("ArsA_ATPase", "metals", "core"),
    ("ArsD", "metals", "core"),
    ("Cation_ATPase_N", "metals", "core"),
    ("Bacteriocin_IId", "bacteriocins", "core"),
    ("DNA_pol_A_exo1", "control_PolA", "control"),
    ("DNAG_N", "control_DnaG", "control"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def run(cmd: list[str], **kw):
    return subprocess.run(cmd, check=True, **kw)


def gca_of(rep_gbff: str) -> str | None:
    m = re.search(r"(GCA_\d+\.\d+)", str(rep_gbff or ""))
    return m.group(1) if m else None


def parse_element(eid: str) -> tuple[str, int, int] | None:
    parts = str(eid).split(":")
    if len(parts) != 3:
        return None
    contig, span = parts[1], parts[2]
    if "-" not in span:
        return None
    a, b = span.rsplit("-", 1)
    try:
        lo, hi = int(a), int(b)
    except ValueError:
        return None
    if hi < lo:
        lo, hi = hi, lo
    return contig, lo, hi


def merge_iv(ivs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not ivs:
        return []
    xs = sorted(ivs)
    out = [xs[0]]
    for a, b in xs[1:]:
        if a <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def covered(ivs: list[tuple[int, int]]) -> int:
    return sum(b - a + 1 for a, b in merge_iv(ivs))


def overlaps(lo: int, hi: int, ivs: list[tuple[int, int]]) -> bool:
    mid = (lo + hi) // 2
    for a, b in ivs:
        if a <= mid <= b:
            return True
    return False


def hmm_name_from_text(txt: str) -> str:
    for line in txt.splitlines():
        if line.startswith("NAME"):
            return line.split()[1]
    return ""


def stage_hmm() -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    hmm = OUT / "hallmark.hmm"
    man_p = OUT / "hmm_manifest.tsv"
    rows = []
    chunks: list[str] = []

    for fn, cls, role in DF_FILES:
        src = DF / fn
        if not src.is_file():
            log(f"skip missing DF {fn}")
            continue
        txt = src.read_text()
        name = hmm_name_from_text(txt) or fn
        chunks.append(txt if txt.endswith("\n") else txt + "\n")
        rows.append({"hmm_name": name, "class": cls, "role": role, "source": f"DF:{fn}"})

    for fn, cls, role in CAS_FILES:
        src = CAS / fn
        if not src.is_file():
            log(f"skip missing Cas {fn}")
            continue
        txt = src.read_text()
        name = hmm_name_from_text(txt) or fn
        chunks.append(txt if txt.endswith("\n") else txt + "\n")
        rows.append({"hmm_name": name, "class": cls, "role": role, "source": f"CasFinder:{fn}"})

    for name, cls, role in PFAM_FETCH:
        r = subprocess.run([str(HMMFETCH), str(PFAM), name], capture_output=True, text=True)
        if r.returncode != 0 or "HMMER3" not in r.stdout:
            log(f"skip Pfam {name}: {r.stderr.strip()[:120]}")
            continue
        chunks.append(r.stdout if r.stdout.endswith("\n") else r.stdout + "\n")
        rows.append({"hmm_name": name, "class": cls, "role": role, "source": f"Pfam:{name}"})

    stat_p = OUT / "padloc_hmmstat.txt"
    if not stat_p.exists():
        log("hmmstat padlocdb…")
        run([str(HMMSTAT), str(PADLOC)], stdout=stat_p.open("w"))
    names = []
    for line in stat_p.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            names.append(parts[1])
    used_pad = set()
    for typ in PDC_TYPES:
        hit = next((n for n in names if n == typ or n.startswith(typ + "_")), None)
        if not hit or hit in used_pad:
            continue
        used_pad.add(hit)
        r = subprocess.run([str(HMMFETCH), str(PADLOC), hit], capture_output=True, text=True)
        if r.returncode != 0 or "HMMER3" not in r.stdout:
            log(f"skip PADLOC {typ}/{hit}")
            continue
        chunks.append(r.stdout if r.stdout.endswith("\n") else r.stdout + "\n")
        rows.append({"hmm_name": hit, "class": "PDC", "role": "core", "source": f"PADLOC:{typ}"})

    hmm.write_text("".join(chunks))
    pd.DataFrame(rows).to_csv(man_p, sep="\t", index=False)
    log(f"HMM n={len(rows)} -> {hmm}")
    run([str(HMMPRESS), "-f", str(hmm)])
    return hmm


def load_unknown_ids() -> set[str]:
    ie = pd.read_csv(IE_TSV, sep="\t")
    return set(ie.integrase_id.astype(str))


def stage_intervals(max_genomes: int | None) -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    unk = load_unknown_ids()
    idx = pd.read_csv(ATLAS, sep="\t")
    idx["gca"] = idx.rep_gbff.map(gca_of)
    idx["parsed"] = idx.element_id.map(parse_element)
    idx = idx[idx.parsed.notna() & idx.gca.notna()].copy()
    idx["contig"] = idx.parsed.map(lambda x: x[0])
    idx["lo"] = idx.parsed.map(lambda x: x[1])
    idx["hi"] = idx.parsed.map(lambda x: x[2])
    idx["is_unknown"] = idx.integrase_id.astype(str).isin(unk)

    gcas = sorted(idx.loc[idx.is_unknown, "gca"].unique())
    rows = []
    for g in gcas:
        faa = FINDER / f"{g}_genomic" / "orfs.faa"
        gff = FINDER / f"{g}_genomic" / "orfs.gff"
        rows.append({
            "gca": g,
            "has_orfs": faa.is_file() and gff.is_file(),
            "faa": str(faa) if faa.is_file() else "",
            "gff": str(gff) if gff.is_file() else "",
            "n_unk_ie": int(((idx.gca == g) & idx.is_unknown).sum()),
            "n_atlas_ie": int((idx.gca == g).sum()),
        })
    gen = pd.DataFrame(rows)
    keep = gen[gen.has_orfs].copy()
    if max_genomes:
        keep = keep.head(max_genomes)
    keep.to_csv(OUT / "genomes_used.tsv", sep="\t", index=False)
    used = set(keep.gca)
    iv = idx[idx.gca.isin(used)][
        ["integrase_id", "element_id", "gca", "contig", "lo", "hi", "is_unknown", "primary_label"]
    ]
    iv.to_csv(OUT / "ie_intervals.tsv", sep="\t", index=False)
    log(f"genomes with orfs {keep.has_orfs.sum()}/{len(gen)} used={len(keep)}; "
        f"unk IE on used genomes {iv.is_unknown.sum()}")
    return keep


def parse_gff_lengths(gff: Path) -> dict[str, int]:
    lens: dict[str, int] = {}
    rx = re.compile(r"seqlen=(\d+);seqhdr=\"(\S+)")
    with open(gff) as fh:
        for line in fh:
            if line.startswith("# Sequence Data:") or line.startswith("# Sequence Data"):
                m = rx.search(line)
                if m:
                    lens[m.group(2)] = int(m.group(1))
    return lens


def stage_lengths() -> pd.DataFrame:
    gen = pd.read_csv(OUT / "genomes_used.tsv", sep="\t")
    iv = pd.read_csv(OUT / "ie_intervals.tsv", sep="\t")
    iv["is_unknown"] = iv.is_unknown.astype(str).isin(("True", "true", "1"))
    by_g = {g: sub for g, sub in iv.groupby("gca")}
    recs = []
    for r in gen.itertuples(index=False):
        gff = Path(r.gff)
        contig_len = parse_gff_lengths(gff)
        genome_nt = sum(contig_len.values())
        n_cds = n_cds_unk = n_cds_other = 0
        cds_end_max: dict[str, int] = defaultdict(int)
        sub = by_g.get(r.gca, pd.DataFrame())
        unk_map: dict[str, list[tuple[int, int]]] = defaultdict(list)
        oth_map: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for row in sub.itertuples(index=False):
            (unk_map if row.is_unknown else oth_map)[row.contig].append((int(row.lo), int(row.hi)))
        with open(gff) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                p = line.split("\t")
                if len(p) < 5 or p[2] != "CDS":
                    continue
                try:
                    lo, hi = int(p[3]), int(p[4])
                except ValueError:
                    continue
                n_cds += 1
                ct = p[0]
                cds_end_max[ct] = max(cds_end_max[ct], hi)
                if overlaps(lo, hi, unk_map.get(ct, [])):
                    n_cds_unk += 1
                elif overlaps(lo, hi, oth_map.get(ct, [])):
                    n_cds_other += 1
        if genome_nt == 0:
            genome_nt = sum(cds_end_max.values())
            contig_len = dict(cds_end_max)
        L_unk = 0
        L_all_ie = 0
        for ct, clen in {**{c: 0 for c in set(unk_map) | set(oth_map)}, **contig_len}.items():
            cap = contig_len.get(ct, 10**12)
            u = [(max(1, a), min(b, cap)) for a, b in unk_map.get(ct, [])]
            o = [(max(1, a), min(b, cap)) for a, b in oth_map.get(ct, [])]
            L_unk += covered(u)
            L_all_ie += covered(u + o)
        recs.append({
            "gca": r.gca,
            "genome_nt": genome_nt,
            "n_contig": len(contig_len),
            "n_cds": n_cds,
            "n_cds_unk_ie": n_cds_unk,
            "n_cds_other_ie": n_cds_other,
            "n_cds_chrom": n_cds - n_cds_unk - n_cds_other,
            "L_unk_ie": L_unk,
            "L_all_ie": L_all_ie,
            "L_rest": max(0, genome_nt - L_unk),
            "L_chrom": max(0, genome_nt - L_all_ie),
            "n_unk_ie": int(r.n_unk_ie),
            "n_atlas_ie": int(r.n_atlas_ie),
        })
    df = pd.DataFrame(recs)
    df.to_csv(OUT / "genome_lengths.tsv", sep="\t", index=False)
    log(f"length table n={len(df)} genome_Mb={df.genome_nt.sum()/1e6:.1f} "
        f"unk_ie_Mb={df.L_unk_ie.sum()/1e6:.2f}")
    return df


def rewrite_faa(gca: str, src: Path, dest) -> int:
    n = 0
    with open(src) as fh:
        for line in fh:
            if line.startswith(">"):
                toks = line[1:].split()
                hid = toks[0]
                contig = hid.rsplit("_", 1)[0]
                try:
                    start, end = int(toks[2]), int(toks[4])
                except (IndexError, ValueError):
                    start, end = 0, 0
                dest.write(f">{gca}|{contig}|{start}|{end}\n")
                n += 1
            else:
                dest.write(line)
    return n


def stage_faa() -> list[Path]:
    gen = pd.read_csv(OUT / "genomes_used.tsv", sep="\t")
    d = OUT / "faa_batches"
    d.mkdir(parents=True, exist_ok=True)
    stamp = d / "stamp.txt"
    key = f"{len(gen)}\t{BATCH_N}\t{gen.gca.iloc[0]}\t{gen.gca.iloc[-1]}"
    if stamp.exists() and stamp.read_text().strip() != key:
        for p in d.glob("batch_*.faa"):
            p.unlink()
        hitdir = OUT / "hmm_tbl"
        if hitdir.is_dir():
            for p in hitdir.glob("batch_*.tbl"):
                p.unlink()
        log("cleared faa/hmm batches (genome list changed)")
    stamp.write_text(key + "\n")
    paths = []
    gcas = list(gen.gca)
    for i in range(0, len(gcas), BATCH_N):
        batch = gcas[i:i + BATCH_N]
        dest = d / f"batch_{i // BATCH_N:03d}.faa"
        if dest.exists() and dest.stat().st_size > 1000:
            paths.append(dest)
            continue
        n = 0
        with open(dest, "w") as out:
            for g in batch:
                src = FINDER / f"{g}_genomic" / "orfs.faa"
                n += rewrite_faa(g, src, out)
        log(f"faa {dest.name} proteins={n}")
        paths.append(dest)
    return paths


def _hmmsearch_one(args: tuple[str, str, str, int]) -> str:
    hmm, faa, tbl, cpu = args
    if Path(tbl).exists() and Path(tbl).stat().st_size > 100:
        return tbl
    tmp = tbl + ".tmp"
    cmd = [
        str(HMMSEARCH), "-E", str(E_CUT), "--tblout", tmp, "--noali",
        "--cpu", str(cpu), hmm, faa,
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    Path(tmp).replace(tbl)
    return tbl


def stage_search() -> Path:
    hmm = OUT / "hallmark.hmm"
    batches = sorted((OUT / "faa_batches").glob("batch_*.faa"))
    hitdir = OUT / "hmm_tbl"
    hitdir.mkdir(exist_ok=True)
    jobs = []
    for faa in batches:
        tbl = hitdir / (faa.stem + ".tbl")
        jobs.append((str(hmm), str(faa), str(tbl), CPU_PER_SEARCH))
    log(f"hmmsearch {len(jobs)} batches × {CPU_PER_SEARCH} cpu, {N_SEARCH_WORKERS} workers")
    with ProcessPoolExecutor(max_workers=N_SEARCH_WORKERS) as ex:
        futs = [ex.submit(_hmmsearch_one, j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            p = fut.result()
            log(f"  done {i}/{len(jobs)} {Path(p).name}")
    merged = OUT / "hmmsearch.tbl"
    with open(merged, "w") as out:
        for tbl in sorted(hitdir.glob("batch_*.tbl")):
            out.write(tbl.read_text())
    return merged


def parse_tbl(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 6:
                continue
            rows.append({
                "target": p[0],
                "hmm": p[2],
                "evalue": float(p[4]),
                "bits": float(p[5]),
            })
    return pd.DataFrame(rows)


def compartment_of(gca: str, contig: str, lo: int, hi: int,
                   unk: dict, oth: dict) -> str:
    key = (gca, contig)
    if overlaps(lo, hi, unk.get(key, [])):
        return "unk_ie"
    if overlaps(lo, hi, oth.get(key, [])):
        return "other_ie"
    return "chrom"


def stage_summarize() -> pd.DataFrame:
    man = pd.read_csv(OUT / "hmm_manifest.tsv", sep="\t")
    role = dict(zip(man.hmm_name, man.role))
    cls_of = dict(zip(man.hmm_name, man["class"]))
    hits = parse_tbl(OUT / "hmmsearch.tbl")
    if hits.empty:
        raise SystemExit("no hmmsearch hits")
    hits = hits[(hits.evalue <= E_CUT) & (hits.bits >= BIT_CUT)].copy()
    hits["class"] = hits.hmm.map(cls_of)
    hits["role"] = hits.hmm.map(role)
    hits = hits[hits["class"].notna() & hits.role.isin(("core", "rease", "control"))]
    # one protein × class
    hits = hits.sort_values("bits", ascending=False).drop_duplicates(["target", "class"])

    parts = hits.target.str.split("|", n=3, expand=True)
    hits["gca"] = parts[0]
    hits["contig"] = parts[1]
    hits["lo"] = pd.to_numeric(parts[2], errors="coerce")
    hits["hi"] = pd.to_numeric(parts[3], errors="coerce")
    hits = hits.dropna(subset=["lo", "hi"])

    iv = pd.read_csv(OUT / "ie_intervals.tsv", sep="\t")
    iv["is_unknown"] = iv.is_unknown.astype(str).isin(("True", "true", "1"))
    unk: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    oth: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    ie_by_key: dict[tuple[str, str], list[tuple[int, int, str]]] = defaultdict(list)
    for row in iv.itertuples(index=False):
        key = (row.gca, row.contig)
        lo, hi = int(row.lo), int(row.hi)
        if row.is_unknown:
            unk[key].append((lo, hi))
            ie_by_key[key].append((lo, hi, str(row.integrase_id)))
        else:
            oth[key].append((lo, hi))

    hits["comp"] = [
        compartment_of(g, c, int(a), int(b), unk, oth)
        for g, c, a, b in zip(hits.gca, hits.contig, hits.lo, hits.hi)
    ]

    def ie_id_for(g, c, a, b) -> str:
        mid = (int(a) + int(b)) // 2
        for lo, hi, iid in ie_by_key.get((g, c), []):
            if lo <= mid <= hi:
                return iid
        return ""

    hits["integrase_id"] = [
        ie_id_for(g, c, a, b) if comp == "unk_ie" else ""
        for g, c, a, b, comp in zip(hits.gca, hits.contig, hits.lo, hits.hi, hits.comp)
    ]
    hits.to_csv(OUT / "hits_assigned.tsv", sep="\t", index=False)

    lens = pd.read_csv(OUT / "genome_lengths.tsv", sep="\t")
    L_unk = int(lens.L_unk_ie.sum())
    L_rest = int(lens.L_rest.sum())
    L_chrom = int(lens.L_chrom.sum())
    L_all = int(lens.genome_nt.sum())
    n_cds_unk = int(lens.n_cds_unk_ie.sum())
    n_cds_rest = int(lens.n_cds.sum() - lens.n_cds_unk_ie.sum())
    n_cds_chrom = int(lens.n_cds_chrom.sum())
    n_gen = len(lens)
    n_unk_ie = int(lens.n_unk_ie.sum())

    recs = []
    for cls in ORDER:
        sub = hits[hits["class"] == cls]
        n_in = int((sub.comp == "unk_ie").sum())
        n_other = int((sub.comp == "other_ie").sum())
        n_chrom = int((sub.comp == "chrom").sum())
        n_out = n_other + n_chrom
        n_all = n_in + n_out
        rho_in = n_in / (L_unk / 1e6) if L_unk else math.nan
        rho_rest = n_out / (L_rest / 1e6) if L_rest else math.nan
        rho_chrom = n_chrom / (L_chrom / 1e6) if L_chrom else math.nan
        rho_all = n_all / (L_all / 1e6) if L_all else math.nan
        if n_out == 0 and n_in > 0:
            fold_rest = float("inf")
            fold_rest_pseudo = (n_in / (L_unk / 1e6)) / (0.5 / (L_rest / 1e6))
        elif n_in == 0 and n_out == 0:
            fold_rest = math.nan
            fold_rest_pseudo = math.nan
        elif rho_rest == 0:
            fold_rest = float("inf")
            fold_rest_pseudo = math.nan
        else:
            fold_rest = rho_in / rho_rest
            fold_rest_pseudo = fold_rest
        fold_chrom = (rho_in / rho_chrom) if rho_chrom else (float("inf") if n_in else math.nan)
        fold_genome = (rho_in / rho_all) if rho_all else math.nan

        table = [[n_in, max(n_cds_unk - n_in, 0)],
                 [n_out, max(n_cds_rest - n_out, 0)]]
        try:
            orat, p = fisher_exact(table, alternative="greater")
        except ValueError:
            orat, p = math.nan, math.nan

        # per-genome median fold: only genomes with ≥1 hit inside the unknown IE
        folds = []
        for gca, gs in sub.groupby("gca"):
            gl = lens.loc[lens.gca == gca]
            if gl.empty:
                continue
            li = int(gl.L_unk_ie.iloc[0])
            lo_ = int(gl.L_rest.iloc[0])
            if li <= 0 or lo_ <= 0:
                continue
            ni = int((gs.comp == "unk_ie").sum())
            no = int((gs.comp != "unk_ie").sum())
            if ni == 0:
                continue
            folds.append(((ni + 0.5) / (li / 1e6)) / ((no + 0.5) / (lo_ / 1e6)))
        med = float(pd.Series(folds).median()) if folds else math.nan
        n_ie_with = int(sub.loc[sub.comp == "unk_ie", "integrase_id"].replace("", pd.NA).nunique())

        recs.append({
            "class": cls,
            "layer": LAYER[cls],
            "n_hits_unk_ie": n_in,
            "n_hits_rest": n_out,
            "n_hits_other_ie": n_other,
            "n_hits_chrom": n_chrom,
            "n_unk_ie_with_hit": n_ie_with,
            "per_Mb_unk_ie": round(rho_in, 4) if rho_in == rho_in else math.nan,
            "per_Mb_rest": round(rho_rest, 4) if rho_rest == rho_rest else math.nan,
            "per_Mb_chrom": round(rho_chrom, 4) if rho_chrom == rho_chrom else math.nan,
            "per_Mb_genome": round(rho_all, 4) if rho_all == rho_all else math.nan,
            "fold_vs_rest": fold_rest,
            "fold_vs_rest_pseudo": fold_rest_pseudo,
            "fold_vs_chrom": fold_chrom,
            "fold_vs_genome": fold_genome,
            "median_per_genome_fold": med,
            "n_genomes_with_hit": int(sub.gca.nunique()),
            "fisher_or": orat,
            "fisher_p_greater": p,
        })

    out = pd.DataFrame(recs)
    mask = out.layer != "ctrl"
    pvals = out.loc[mask, "fisher_p_greater"].fillna(1.0).to_numpy()
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n)
    ranks[order] = np.arange(1, n + 1)
    raw = np.clip(pvals * n / ranks, 0, 1)
    mono = np.minimum.accumulate(raw[order[::-1]])[::-1]
    fdr = np.empty(n)
    fdr[order] = mono
    out.loc[mask, "fdr_bh"] = fdr
    out.loc[~mask, "fdr_bh"] = pd.NA
    meta = pd.DataFrame([{
        "n_genomes": n_gen,
        "n_unk_ie_on_genomes": n_unk_ie,
        "unknown_cohort": 4689,
        "genome_nt": L_all,
        "L_unk_ie": L_unk,
        "L_rest": L_rest,
        "L_chrom": L_chrom,
        "n_cds": int(lens.n_cds.sum()),
        "n_cds_unk_ie": n_cds_unk,
        "n_cds_rest": n_cds_rest,
        "n_cds_chrom": n_cds_chrom,
        "evalue": E_CUT,
        "bits": BIT_CUT,
        "detector": "compact hallmark HMM (DF + CasFinder + Pfam + 1 PADLOC/PDC type)",
        "contrast_rest": "genome minus union of unknown IE intervals (other MGEs stay in rest)",
        "contrast_chrom": "genome minus union of all atlas IE intervals",
    }])
    out.to_csv(OUT / "enrichment_by_class.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "enrichment_meta.tsv", sep="\t", index=False)
    log(out.to_string(index=False))
    log(meta.T.to_string())
    return out


def _fmt_fold(x) -> str:
    if x is None or (isinstance(x, float) and (math.isnan(x))):
        return "NA"
    if x == float("inf"):
        return "∞"
    if x >= 100:
        return f"{x:.0f}"
    if x >= 10:
        return f"{x:.1f}"
    return f"{x:.2f}"


def stage_plot() -> None:
    df = pd.read_csv(OUT / "enrichment_by_class.tsv", sep="\t")
    meta = pd.read_csv(OUT / "enrichment_meta.tsv", sep="\t").iloc[0]
    empty_lit = {"Gabija", "Septu", "Hachiman", "Zorya", "DISARM"}
    df = df[df.layer != "ctrl"].copy()
    df = df[~df["class"].isin(empty_lit)].copy()
    df = df[(df.n_hits_unk_ie + df.n_hits_rest) > 0].copy()
    df["fold_plot"] = df.fold_vs_rest_pseudo.map(lambda x: 200 if x == float("inf") else x)
    df = df.sort_values("fold_plot", ascending=True)
    colors = []
    for layer in df.layer:
        colors.append({"A": BLUE, "B": VERM, "C": GREEN}.get(layer, GREY))

    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    y = range(len(df))
    ax.barh(list(y), df.fold_plot, color=colors, height=0.72)
    ax.axvline(1, color=INK, lw=1.0, ls="--")
    ax.set_yticks(list(y))
    ax.set_yticklabels(df["class"])
    ax.set_xlabel("fold = (hits / Mb in unknown IE) / (hits / Mb in rest of genome)")
    ax.set_title(
        f"Unknown IE vs rest of host genome  ·  {int(meta.n_genomes)} genomes, "
        f"{int(meta.n_unk_ie_on_genomes)} IE"
    )
    xmax = max(8, float(df.fold_plot.max()) * 1.18)
    ax.set_xlim(0, min(xmax, 80) if float(df.fold_plot.max()) < 80 else xmax)
    for i, r in enumerate(df.itertuples(index=False)):
        lab = f"{_fmt_fold(r.fold_vs_rest)}   n={int(r.n_hits_unk_ie)}/{int(r.n_hits_rest)}"
        ax.text(min(r.fold_plot, ax.get_xlim()[1]) + 0.15, i, lab,
                va="center", fontsize=8, color=INK)
    ax.plot([], [], color=BLUE, lw=6, label="A maintenance")
    ax.plot([], [], color=VERM, lw=6, label="B defense")
    ax.plot([], [], color=GREEN, lw=6, label="C host-adaptive")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        dest = FIG / f"Figure_unknown_ie_genome_enrichment.{ext}"
        fig.savefig(dest, dpi=200 if ext == "png" else None, bbox_inches="tight")
        log(f"-> {dest}")
    plt.close(fig)

    # control panel
    ctrl = pd.read_csv(OUT / "enrichment_by_class.tsv", sep="\t")
    ctrl = ctrl[ctrl.layer == "ctrl"]
    if not ctrl.empty:
        fig, ax = plt.subplots(figsize=(6.4, 2.2))
        xs = list(ctrl["class"])
        ys = [200 if x == float("inf") else x for x in ctrl.fold_vs_rest_pseudo]
        ax.bar(xs, ys, color=GREY, width=0.55)
        ax.axhline(1, color=INK, ls="--", lw=1)
        ax.set_ylabel("fold vs rest")
        ax.set_title("Negative control (chromosomal PolA / DnaG)")
        for i, r in enumerate(ctrl.itertuples(index=False)):
            ax.text(i, ys[i], f" {_fmt_fold(r.fold_vs_rest)}", ha="center", va="bottom", fontsize=8)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            dest = FIG / f"Figure_unknown_ie_genome_enrichment_controls.{ext}"
            fig.savefig(dest, dpi=200 if ext == "png" else None, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["hmm", "intervals", "faa", "search", "summarize", "plot", "all"])
    ap.add_argument("--max-genomes", type=int, default=0)
    args = ap.parse_args()
    stages = (["hmm", "intervals", "faa", "search", "summarize", "plot"]
              if args.stage == "all" else [args.stage])
    mx = args.max_genomes or None
    if "hmm" in stages:
        stage_hmm()
    if "intervals" in stages:
        stage_intervals(mx)
        stage_lengths()
    if "faa" in stages:
        stage_faa()
    if "search" in stages:
        stage_search()
    if "summarize" in stages:
        stage_summarize()
    if "plot" in stages:
        stage_plot()


if __name__ == "__main__":
    main()
