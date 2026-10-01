"""Biological reading of Leiden-HD clusters for the unclassified cohort.

No new clusters are built here. The inputs are a ready feature matrix,
Leiden labels, and ``ie_rep_rows.tsv`` (rep_id, matrix row, ie_id, sample_id).
On top of that this module computes:

* per-cluster markers (binary X, p_in / p_out / log2 enrichment, marker_score),
* aggregated functional roles (hallmarks) from the Pfam clan plus rules,
* top ``product=`` names,
* host phyla, IE lengths, and IE atlas scores.

Pfam accessions missing from the marker table are not fetched from InterPro again.
``pfam_descriptions.tsv`` is written once and reloaded.
"""

from __future__ import annotations

import json
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, load_npz

# Paths. Pfam tables used by the current unknown clustering live under
# final_final/03_clustering/reference. The other names below belonged to the
# retired notebook and are kept so older callers still resolve a path.

REPO = Path("/home/lam34/MGE_finder").resolve()
NB3_BASE = REPO / "final_final" / "03_clustering"
NB3_OUT = NB3_BASE / "reference"
ROOT_ORG = REPO / "results_organized"
PADLOC_PIPELINE = REPO / "padloc_pipeline"

# Alternate padloc.csv locations. For gtdb the file is under results_organized.
# For thermaceae and deinococcales it sits in a separate padloc_pipeline run,
# named results_<group>_padloc.csv. Bins have no PADLOC run.
def _padloc_paths(group: str, sample_id: str) -> list[Path]:
    return [
        ROOT_ORG / group / sample_id / "padloc.csv",
        PADLOC_PIPELINE / f"results_{group}" / sample_id / f"results_{group}_padloc.csv",
    ]

PATH_FEATURE_SUBMAT = NB3_OUT / "filtered_ie_feature_subset.npz"
PATH_FEATURE_LABELS = ROOT_ORG / "genome_clustering" / "project" / "genome_matrix" / "genome_feature_labels.npz"
PATH_MAP = NB3_OUT / "ie_rep_rows.tsv"
PATH_LEIDEN = NB3_OUT / "filtered_ie_leiden_hd.tsv"
PATH_GTDB_TAX = ROOT_ORG / "misc" / "gtdb_taxonomy.tsv"
PATH_ATLAS = NB3_BASE / "filtered_ie_ie_atlas_precluster.tsv"
PATH_PROTEIN_LONG = ROOT_ORG / "misc" / "mge_protein_cluster_matrix_long.tsv"

PATH_PFAM_DESC = NB3_OUT / "pfam_descriptions.tsv"
PATH_PFAM_CATALOG = NB3_OUT / "pfam_catalog.tsv"
PATH_MARKERS_OUT = NB3_OUT / "cluster_markers_unclassified_leiden_hd.tsv"
PATH_HALLMARK_OUT = NB3_OUT / "cluster_hallmark_composition.tsv"
PATH_PRODUCTS_OUT = NB3_OUT / "cluster_top_products.tsv"
PATH_SUMMARY_OUT = NB3_OUT / "cluster_biology_summary.tsv"
PATH_PADLOC_PER_IE_OUT = NB3_OUT / "padloc_per_ie.tsv"
PATH_PADLOC_CLUSTER_LONG_OUT = NB3_OUT / "cluster_padloc_systems.tsv"
PATH_PADLOC_CLUSTER_SUMMARY_OUT = NB3_OUT / "cluster_padloc_summary.tsv"

# Pfam accession to functional role, for the accessions in this cohort.

# Role vocabulary:
#   - integrase_tyrosine: tyrosine recombinase (classic Phage_integrase)
#   - integrase_serine:    serine recombinase / resolvase
#   - regulator_HTH:       any HTH domain (CI/Cro, MarR, TetR, MerR, AraC...)
#   - replication:         replication helpers (Rep_trans, Helitron, RepSA, ...)
#   - methyltransferase:   N6_Mtase / Methyltransf_*
#   - restriction_T1:      Type I R-M (HsdM, HSDR, Methylase_S, MmeI*)
#   - defense_STAND:       NACHT/NB-ARC/nSTAND* (anti-phage defense)
#   - toxin_antitoxin:     AbiEi/MqsA/YdaS/Xre-MbcA/CopG (RHH_5)
#   - transposase:         HTH_Tnp_*, IstB, Zn_ribbon_IS1595
#   - aaa_atpase:          AAA domains (often helicase loaders / clamp)
#   - helicase:            named helicases (DEAD, RNA_helicase, SWI2_SNF2)
#   - tpr_repeat:          TPR repeats (scaffold)
#   - zn_ribbon:           Zn-ribbon domains (often Pol/primase)
#   - other:               everything else

PFAM_ROLE_RULES: dict[str, str] = {
    # tyrosine integrases / phage integrase family
    "PF00589": "integrase_tyrosine",
    "PF02899": "integrase_tyrosine",
    "PF12834": "integrase_tyrosine",
    "PF14659": "integrase_tyrosine",
    "PF13495": "integrase_tyrosine",
    "PF14657": "integrase_tyrosine",
    "PF12167": "integrase_tyrosine",
    "PF09003": "integrase_tyrosine",
    "PF13356": "integrase_tyrosine",
    "PF22022": "integrase_tyrosine",
    "PF24624": "integrase_tyrosine",
    "PF26003": "integrase_tyrosine",
    "PF13102": "integrase_tyrosine",
    "PF16795": "integrase_tyrosine",
    "PF12835": "integrase_tyrosine",
    "PF22823": "integrase_tyrosine",
    # serine resolvases
    "PF02796": "integrase_serine",
    # IntI/rve transposase / IS-related
    "PF00665": "transposase",
    "PF13683": "transposase",
    "PF01527": "transposase",
    "PF13022": "transposase",
    "PF13542": "transposase",
    "PF12760": "transposase",
    "PF01695": "transposase",  # IstB IS21 ATPase
    # replication-associated
    "PF00476": "replication",   # DNA pol A
    "PF01612": "replication",   # 3'-5' exonuclease
    "PF02486": "replication",   # Rep_trans
    "PF20199": "replication",   # RepSA
    "PF21860": "replication",   # Replitron_C
    "PF14214": "replication",   # Helitron_like_N
    "PF23343": "replication",   # REP_ORF2-G2P
    "PF26475": "replication",   # DUF8148_C / Rep clan
    "PF09250": "replication",   # Prim-Pol
    "PF01807": "replication",   # Zn_ribbon_DnaG (primase)
    "PF22548": "replication",   # AEP-TOTE primase
    "PF01751": "replication",   # Toprim
    "PF13155": "replication",   # Toprim_2
    "PF13362": "replication",
    "PF13662": "replication",
    "PF05840": "replication",   # Phage_GPA
    "PF03266": "replication",   # NTPase_1
    "PF19263": "replication",   # DUF5906 (P-loop, often Pol-helicase)
    # methyltransferases (DNA modification, often part of R-M)
    "PF02384": "methyltransferase",
    "PF01420": "restriction_T1",  # Methylase_S = TRD specificity
    "PF12161": "restriction_T1",  # HsdM_N
    "PF04313": "restriction_T1",  # HSDR_N
    "PF13588": "restriction_T1",
    "PF22679": "restriction_T1",
    "PF11867": "restriction_T1",
    "PF20464": "restriction_T1",  # MmeI_*
    "PF20465": "restriction_T1",
    "PF20466": "restriction_T1",
    "PF20467": "restriction_T1",
    "PF20473": "restriction_T1",
    "PF07669": "restriction_T1",  # Eco57I
    "PF04851": "restriction_T1",  # ResIII
    "PF25120": "restriction_T1",  # DUF7814 (MmeI_HD-like clan)
    "PF24441": "restriction_T1",  # DUF7560 (often near R-M)
    "PF11599": "methyltransferase",  # AviRa
    "PF08241": "methyltransferase",  # Methyltransf_11
    "PF08242": "methyltransferase",
    "PF13649": "methyltransferase",
    "PF13489": "methyltransferase",
    "PF13847": "methyltransferase",
    "PF05175": "methyltransferase",
    "PF01170": "methyltransferase",
    # STAND / NACHT / NB-ARC defense
    "PF20720": "defense_STAND",
    "PF25199": "defense_STAND",
    "PF05729": "defense_STAND",  # NACHT
    "PF00931": "defense_STAND",  # NB-ARC
    # toxin/antitoxin & abortive infection
    "PF13338": "toxin_antitoxin",  # AbiEi
    "PF15943": "toxin_antitoxin",  # YdaS toxin
    "PF15731": "toxin_antitoxin",  # MqsA antitoxin
    "PF20432": "toxin_antitoxin",  # Xre-like-HTH
    "PF23125": "toxin_antitoxin",  # Xre-MbcA-ParS
    "PF07878": "toxin_antitoxin",  # RHH_5 (CopG)
    # CRISPR / Cas
    "PF07282": "crispr_cas",  # Cas12f1-like_TNB
    # transcription factors / regulators (HTH)
    "PF01381": "regulator_HTH",  # HTH_3
    "PF12728": "regulator_HTH",  # HTH_17
    "PF12844": "regulator_HTH",  # HTH_19
    "PF13413": "regulator_HTH",  # HTH_25
    "PF13443": "regulator_HTH",  # HTH_26 (Cro/C1)
    "PF13412": "regulator_HTH",  # HTH_24 wHTH
    "PF13384": "regulator_HTH",  # HTH_23
    "PF13518": "regulator_HTH",  # HTH_28
    "PF13551": "regulator_HTH",  # HTH_29
    "PF13565": "regulator_HTH",  # HTH_32
    "PF13556": "regulator_HTH",
    "PF13744": "regulator_HTH",
    "PF13730": "regulator_HTH",
    "PF13693": "regulator_HTH",
    "PF13404": "regulator_HTH",  # HTH_AsnC
    "PF13309": "regulator_HTH",
    "PF22495": "regulator_HTH",  # RsaL
    "PF20038": "regulator_HTH",
    "PF13560": "regulator_HTH",  # HTH_31
    "PF13936": "regulator_HTH",
    "PF21531": "regulator_HTH",  # Rv2175c_wHTH
    "PF24250": "regulator_HTH",  # HVO_2718 wHTH
    "PF25213": "regulator_HTH",
    "PF14549": "regulator_HTH",  # P22_Cro
    "PF02954": "regulator_HTH",  # HTH_8 (Fis)
    "PF00440": "regulator_HTH",  # TetR_N
    "PF01047": "regulator_HTH",  # MarR
    "PF12802": "regulator_HTH",  # MarR_2
    "PF00376": "regulator_HTH",  # MerR
    "PF13411": "regulator_HTH",  # MerR_1
    "PF00165": "regulator_HTH",  # AraC
    "PF08220": "regulator_HTH",  # DeoR
    "PF08279": "regulator_HTH",  # HTH_11
    "PF08667": "regulator_HTH",  # BetR
    "PF08281": "regulator_HTH",  # Sigma70_r4
    "PF04545": "regulator_HTH",
    "PF07638": "regulator_HTH",  # ECF sigma
    "PF09339": "regulator_HTH",  # IclR
    "PF13545": "regulator_HTH",  # HTH_Crp
    "PF07022": "regulator_HTH",  # Phage_CI_repr
    "PF01978": "regulator_HTH",  # TrmB
    "PF03704": "regulator_HTH",  # BTAD
    "PF05930": "regulator_HTH",  # Phage_AlpA
    "PF17762": "regulator_HTH",  # ParB HTH
    "PF06785": "regulator_HTH",  # UPF0242
    "PF07374": "regulator_HTH",  # DUF1492
    "PF10825": "regulator_HTH",  # DUF2752
    "PF11917": "regulator_HTH",  # DUF3435 (DNA-mend clan)
    # AAA+ ATPases (often helicase loaders / clamp / unwinder)
    "PF00004": "aaa_atpase",
    "PF07728": "aaa_atpase",  # AAA_5
    "PF12775": "aaa_atpase",
    "PF13173": "aaa_atpase",
    "PF13175": "aaa_atpase",
    "PF13176": "aaa_atpase",
    "PF13191": "aaa_atpase",
    "PF13238": "aaa_atpase",
    "PF13245": "aaa_atpase",
    "PF13304": "aaa_atpase",
    "PF13401": "aaa_atpase",
    "PF13476": "aaa_atpase",
    "PF13521": "aaa_atpase",
    "PF13555": "aaa_atpase",
    "PF13604": "aaa_atpase",
    "PF13671": "aaa_atpase",
    "PF03193": "aaa_atpase",  # RsgA GTPase
    "PF18766": "aaa_atpase",  # SWI2/SNF2
    # transport / ABC
    "PF00005": "transport_abc",
    # helicases (RNA / DEAD)
    "PF00270": "helicase",
    "PF00910": "helicase",  # RNA_helicase
    # TPR repeats
    "PF00515": "tpr_repeat",
    "PF07719": "tpr_repeat",
    "PF07721": "tpr_repeat",
    "PF13174": "tpr_repeat",
    "PF13181": "tpr_repeat",
    "PF13371": "tpr_repeat",
    "PF13374": "tpr_repeat",
    "PF13414": "tpr_repeat",
    "PF13424": "tpr_repeat",
    "PF13428": "tpr_repeat",
    "PF13431": "tpr_repeat",
    "PF13432": "tpr_repeat",
    "PF14559": "tpr_repeat",
    "PF14561": "tpr_repeat",
    "PF20308": "tpr_repeat",
    "PF21033": "tpr_repeat",
    "PF24575": "tpr_repeat",
    "PF25063": "tpr_repeat",
    "PF17874": "tpr_repeat",
    "PF23914": "tpr_repeat",
    "PF12895": "tpr_repeat",
    # Zn ribbon (often DNA binding, including primase)
    "PF13240": "zn_ribbon",
    "PF13248": "zn_ribbon",
    "PF13717": "zn_ribbon",
    "PF09297": "zn_ribbon",
    "PF09723": "zn_ribbon",
    "PF17032": "zn_ribbon",
    "PF14803": "zn_ribbon",
    "PF12773": "zn_ribbon",
    "PF24463": "zn_ribbon",
    "PF24460": "zn_ribbon",
    "PF14353": "zn_ribbon",
    "PF21344": "zn_ribbon",
    "PF22747": "zn_ribbon",
    "PF23548": "zn_ribbon",
    "PF18912": "zn_ribbon",
    "PF26408": "zn_ribbon",
    "PF26788": "zn_ribbon",
    "PF18439": "zn_ribbon",
    "PF17915": "zn_ribbon",
    "PF17283": "zn_ribbon",
    "PF07754": "zn_ribbon",
    "PF08271": "zn_ribbon",
    "PF10571": "zn_ribbon",
    "PF01155": "zn_ribbon",  # HypA Zn beta ribbon
    # phage virion / packaging
    "PF06056": "phage_virion",  # Terminase_5 ATPase
    "PF16677": "phage_virion",  # GP3 packaging
    # SMC / chromosome maintenance
    "PF02463": "smc_chromosome",
    # exonuclease / DNA-mend rest
    # peptidase / metallopeptidase
    "PF01863": "peptidase",
    "PF08325": "peptidase",
    # iron-sulfur clusters
    "PF12837": "iron_sulfur",
    "PF12838": "iron_sulfur",
    # everything else => 'other' below
}


# Heuristics on full_name (case-insensitive regexes). Used for Pfam accessions
# that have no manual rule and no known clan. Order matters: the first match wins.
_NAME_RULES: list[tuple[re.Pattern, str]] = [
    (re.compile(r"\bphage\s*integrase|\btyr(osine)?\s*recombinase|\barm[-_ ]?dna|\bint(egrase)?_n", re.I), "integrase_tyrosine"),
    (re.compile(r"\bresolvase|\bserine\s*recombinase|\bgamma[-_ ]?delta\s*resolvase|\bxerd|\bxerc", re.I), "integrase_serine"),
    (re.compile(r"\btransposase|\bIS\d+|\bIstB|\bDDE_Tnp|\bRve\b|\binsertion\s*element|\bMu[-_ ]?A", re.I), "transposase"),
    (re.compile(r"\bDNA[-_ ]?pol|\bDNA[-_ ]?primase|\bRep_trans|\bhelitron|\breplication\s*init|\bRepSA|\bORF2[-_]?G2P|\bRepA|\bGPA", re.I), "replication"),
    (re.compile(r"\bToprim|\btopoisomerase|\bgyr(ase)?", re.I), "replication"),
    (re.compile(r"\bN6[-_ ]?Mtase|\bN[46][-_ ]?Mtase|\bMethyltransf|\bDNA[-_ ]?methyl|\bC5[-_ ]?MTase|\b[Mm]tase", re.I), "methyltransferase"),
    (re.compile(r"\bHsdM|\bHsdR|\bHsdS|\bMethylase_S|\bMmeI|\bResIII|\bEco57I|\brestriction[-_ ]?(endo|enzy|modi)", re.I), "restriction_T1"),
    (re.compile(r"\bNACHT|\bNB[-_ ]?ARC|\bSTAND|\bAvs|\bCBASS|\bThoeris|\bGabija|\bHachiman|\bDruantia|\bSeptu|\bWadjet", re.I), "defense_STAND"),
    (re.compile(r"\bantitoxin|\btoxin|\bAbi[A-Z]|\b(MqsA|YdaS|ParD|ParE|RelE|MazF|HipA|HipB)\b|\bRHH|\bCopG|\bXre[-_ ]?MbcA", re.I), "toxin_antitoxin"),
    (re.compile(r"\bCRISPR|\bCas\d+|\bCas[A-Z]\b|\bcas12|\bcas13|\bcas9", re.I), "crispr_cas"),
    (re.compile(r"\bHTH|\bhelix[-_ ]?turn[-_ ]?helix|\bMerR|\bMarR|\bTetR|\bAraC|\bDeoR|\bIclR|\bsigma[-_ ]?70|\bSigma[-_ ]?54|\bECF\s*sigma|\bCro[/_-]?C1|\bAlpA|\bBetR|\bphage[-_ ]?CI|\bRsaL|\bAsnC|\bFis\b|\bParB|\bTrmB", re.I), "regulator_HTH"),
    (re.compile(r"\bAAA\b|\bAAA[+_-]?\d|\bATPase\s*(family|associated)|\bP[-_ ]?loop|\bSWI2[/_]?SNF2|\bChlI|\bdynein|\bClpA|\bClpX|\bChsW", re.I), "aaa_atpase"),
    (re.compile(r"\bABC[-_ ]?tran|\bABC\s*transp", re.I), "transport_abc"),
    (re.compile(r"\bDEAD|\bDEAH|\b[Hh]elicase|\bSF1|\bSF2", re.I), "helicase"),
    (re.compile(r"\bTPR|\btetratric", re.I), "tpr_repeat"),
    (re.compile(r"\bZn[-_ ]?ribbon|\bzinc[-_ ]?ribbon|\bZn[-_ ]?finger|\bzinc[-_ ]?finger|\bC2H2", re.I), "zn_ribbon"),
    (re.compile(r"\bphage\s*(capsid|tail|portal|head|baseplate|major|minor)|\bterminase|\bGP[\dA-Z]|\bvirion|\blysin|\bholin|\btail\s*fiber", re.I), "phage_virion"),
    (re.compile(r"\bSMC\b|\bcondensin|\bMukB", re.I), "smc_chromosome"),
    (re.compile(r"\bpeptidase|\bprotease", re.I), "peptidase"),
    (re.compile(r"\b4Fe[-_ ]?4S|\biron[-_ ]?sulfur", re.I), "iron_sulfur"),
    (re.compile(r"\bcoiled[-_ ]?coil", re.I), "coiled_coil"),
]


def assign_role(pfam_acc: str, clan: str = "", full_name: str = "") -> str:
    """Functional role for a Pfam accession. Priority:
    1) the manual map PFAM_ROLE_RULES,
    2) a clan heuristic, when the clan is known,
    3) a full_name heuristic (regexes in _NAME_RULES),
    4) ``other``.
    """
    if pfam_acc in PFAM_ROLE_RULES:
        return PFAM_ROLE_RULES[pfam_acc]
    clan = (clan or "").strip()
    if clan:
        if clan == "HTH":
            return "regulator_HTH"
        if clan == "P-loop_NTPase":
            return "aaa_atpase"
        if clan == "TPR":
            return "tpr_repeat"
        if clan == "Zn_Beta_Ribbon":
            return "zn_ribbon"
        if clan == "NADP_Rossmann":
            return "methyltransferase"
        if clan in ("DNA-mend", "l-integrase_N", "MBD-like"):
            return "integrase_tyrosine"
        if clan == "RNase_H":
            return "transposase"
        if clan in ("Rep", "Toprim-like", "AEP"):
            return "replication"
        if clan in ("MmeI_HD-like", "TRD", "PDDEXK"):
            return "restriction_T1"
        if clan == "Met_repress":
            return "toxin_antitoxin"
        if clan == "FtsL":
            return "smc_chromosome"
        if clan == "Peptidase_MA":
            return "peptidase"
        if clan == "4Fe-4S":
            return "iron_sulfur"
    name = (full_name or "").strip()
    if name:
        for rgx, role in _NAME_RULES:
            if rgx.search(name):
                return role
    return "other"


# Short labels for the summary table and the figure
ROLE_SHORT_LABEL = {
    "integrase_tyrosine": "tyr-integrase",
    "integrase_serine": "ser-resolvase",
    "transposase": "transposase/IS",
    "replication": "replication",
    "methyltransferase": "methylase",
    "restriction_T1": "restriction (T-I/MmeI)",
    "defense_STAND": "STAND/NACHT defense",
    "toxin_antitoxin": "TA / Abi",
    "crispr_cas": "CRISPR/Cas",
    "regulator_HTH": "HTH regulator",
    "aaa_atpase": "AAA+ ATPase",
    "transport_abc": "ABC transport",
    "helicase": "helicase",
    "tpr_repeat": "TPR repeat",
    "zn_ribbon": "Zn ribbon",
    "phage_virion": "phage virion/pkg",
    "smc_chromosome": "SMC/condensin",
    "peptidase": "peptidase",
    "iron_sulfur": "[4Fe-4S] domain",
    "coiled_coil": "coiled-coil",
    "other": "other / unknown",
}


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def load_inputs():
    X = load_npz(PATH_FEATURE_SUBMAT).tocsr()
    lab = np.load(PATH_FEATURE_LABELS, allow_pickle=True)
    features = np.array([str(x) for x in lab["features"]])
    map_df = pd.read_csv(PATH_MAP, sep="\t").sort_values("matrix_row").reset_index(drop=True)
    lei = pd.read_csv(PATH_LEIDEN, sep="\t")
    lei["genome_id"] = lei["genome_id"].astype(str)
    map_df["rep_id"] = map_df["rep_id"].astype(str)
    meta = map_df.merge(lei[["genome_id", "leiden_hd_cluster"]], left_on="rep_id", right_on="genome_id", how="left")
    if meta["leiden_hd_cluster"].isna().any():
        raise ValueError("some IE have no leiden_hd_cluster; check the clustering output")
    meta["leiden_hd_cluster"] = meta["leiden_hd_cluster"].astype(int)
    if X.shape[0] != len(meta):
        raise ValueError(f"X rows ({X.shape[0]}) != meta rows ({len(meta)})")
    return X, features, meta


# ---------------------------------------------------------------------------
# 1. Per-cluster markers (binary X, log enrichment)
# ---------------------------------------------------------------------------


@dataclass
class MarkerParams:
    min_cluster_size: int = 20
    min_count_in: int = 5
    min_p_in_frac: float = 0.05
    top_n: int = 30


def compute_markers(X: csr_matrix, features: np.ndarray, meta: pd.DataFrame, params: MarkerParams | None = None) -> pd.DataFrame:
    params = params or MarkerParams()
    Xb = (X > 0).astype(np.int8)
    n_ie = Xb.shape[0]
    col_sums_global = np.asarray(Xb.sum(axis=0)).ravel()

    sizes = meta["leiden_hd_cluster"].value_counts()
    main = sorted([c for c, n in sizes.items() if n >= params.min_cluster_size])

    rows = []
    for c in main:
        idx = np.where(meta["leiden_hd_cluster"].values == c)[0]
        nin = len(idx)
        sub = Xb[idx]
        c_in = np.asarray(sub.sum(axis=0)).ravel()
        c_out = col_sums_global - c_in
        nout = n_ie - nin
        p_in = c_in / max(nin, 1)
        p_out = np.where(nout > 0, c_out / nout, 0.0)
        eps = 1e-3
        enrich = (p_in + eps) / (p_out + eps)
        score = np.log2(enrich) * np.sqrt(p_in)
        keep = (c_in >= max(params.min_count_in, int(params.min_p_in_frac * nin))) & (p_in > p_out)
        cand = np.where(keep)[0]
        if cand.size == 0:
            continue
        order = np.argsort(-score[cand])[: params.top_n]
        chosen = cand[order]
        for fi in chosen:
            rows.append({
                "cluster": int(c),
                "cluster_size": int(nin),
                "feature_id": features[fi],
                "count_in": int(c_in[fi]),
                "p_in": round(float(p_in[fi]), 3),
                "p_out": round(float(p_out[fi]), 3),
                "enrichment": round(float(enrich[fi]), 2),
                "marker_score": round(float(score[fi]), 3),
            })
    df = pd.DataFrame(rows)
    return df


# ---------------------------------------------------------------------------
# 2. Pfam descriptions and role labels on the markers
# ---------------------------------------------------------------------------


def _fetch_pfam(pf: str) -> tuple[str, str, str, str, str]:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}/"
    try:
        r = urllib.request.urlopen(url, timeout=15)
        d = json.loads(r.read())
        m = d.get("metadata", {})
        nm = m.get("name", {}) or {}
        clan = (m.get("set_info") or {}).get("name", "")
        return pf, nm.get("short", ""), nm.get("name", ""), clan, m.get("type", "")
    except Exception as e:
        return pf, "", "", "", f"ERR:{e}"


def ensure_pfam_descriptions(needed_pfam_ids: set[str]) -> pd.DataFrame:
    """Fetch descriptions for missing Pfam accessions from InterPro and update the TSV."""
    cur = pd.DataFrame(columns=["pfam_acc", "short_name", "full_name", "clan", "type"])
    if PATH_PFAM_DESC.is_file():
        cur = pd.read_csv(PATH_PFAM_DESC, sep="\t", dtype=str).fillna("")
    have = set(cur["pfam_acc"].astype(str))
    todo = sorted(needed_pfam_ids - have)
    if todo:
        out = []
        with ThreadPoolExecutor(max_workers=20) as ex:
            futs = [ex.submit(_fetch_pfam, pf) for pf in todo]
            for f in as_completed(futs):
                out.append(f.result())
        new = pd.DataFrame(out, columns=["pfam_acc", "short_name", "full_name", "clan", "type"])
        cur = pd.concat([cur, new], ignore_index=True).drop_duplicates("pfam_acc", keep="last")
        cur.to_csv(PATH_PFAM_DESC, sep="\t", index=False)
    return cur


def annotate_markers(markers: pd.DataFrame) -> pd.DataFrame:
    """Add short_name, full_name, and role columns to the markers."""
    pfam_ids = {f.replace("pfam_", "").split(".")[0] for f in markers["feature_id"] if f.startswith("pfam_")}
    desc = ensure_pfam_descriptions(pfam_ids)
    desc_idx = desc.set_index("pfam_acc")

    def _row(fid: str) -> tuple[str, str, str, str]:
        if not fid.startswith("pfam_"):
            return ("", "(protein family)", "", "other")
        pf = fid.replace("pfam_", "").split(".")[0]
        if pf not in desc_idx.index:
            return (pf, "", "", "other")
        r = desc_idx.loc[pf]
        clan = r.get("clan", "") or ""
        full = r.get("full_name", "") or ""
        return (r.get("short_name", ""), full, clan, assign_role(pf, clan=clan, full_name=full))

    annotated = markers.copy()
    annotated[["short_name", "full_name", "clan", "role"]] = annotated["feature_id"].apply(lambda x: pd.Series(_row(x)))
    return annotated


# ---------------------------------------------------------------------------
# 3. Hallmark composition per cluster (every Pfam feature in X, not only the markers)
# ---------------------------------------------------------------------------


def _load_pfam_lookup() -> pd.DataFrame:
    """Joined table pfam_acc → (full_name, clan). For the marker subset the clan
    comes from pfam_descriptions.tsv; the rest come from pfam_catalog.tsv
    (full_name only, clan='').
    """
    if not PATH_PFAM_CATALOG.is_file():
        raise FileNotFoundError(
            f"Missing {PATH_PFAM_CATALOG}. Download it once; "
            "see the Pfam catalog cell in the cluster-biology notebook."
        )
    cat = pd.read_csv(PATH_PFAM_CATALOG, sep="\t", dtype=str).fillna("")
    cat = cat[["pfam_acc", "full_name"]].rename(columns={"full_name": "full_name_cat"}).copy()
    if PATH_PFAM_DESC.is_file():
        desc = pd.read_csv(PATH_PFAM_DESC, sep="\t", dtype=str).fillna("")
        desc = desc[["pfam_acc", "short_name", "full_name", "clan"]]
        merged = cat.merge(desc, on="pfam_acc", how="left")
        merged["full_name"] = merged["full_name"].where(merged["full_name"].fillna("") != "", merged["full_name_cat"])
        merged["clan"] = merged["clan"].fillna("")
        merged["short_name"] = merged["short_name"].fillna("")
        return merged[["pfam_acc", "short_name", "full_name", "clan"]]
    cat["full_name"] = cat["full_name_cat"]
    cat["short_name"] = ""
    cat["clan"] = ""
    return cat[["pfam_acc", "short_name", "full_name", "clan"]]


def compute_hallmark_composition(X: csr_matrix, features: np.ndarray, meta: pd.DataFrame) -> pd.DataFrame:
    """Fraction of IE in each cluster that carry at least one Pfam of each role.

    Uses the local ``pfam_catalog.tsv`` (the full Pfam list). InterPro is not
    queried at runtime. Add an accession to ``pfam_descriptions.tsv`` when a
    more precise clan is needed.
    """
    pfam_mask = np.array([str(f).startswith("pfam_") for f in features])
    pfam_idx = np.where(pfam_mask)[0]
    Xb = (X > 0).astype(np.int8)
    col_sums = np.asarray(Xb.sum(axis=0)).ravel()
    used_pfam_idx = pfam_idx[col_sums[pfam_idx] > 0]
    lookup = _load_pfam_lookup().set_index("pfam_acc")

    role_for_col: dict[int, str] = {}
    for col in used_pfam_idx:
        pf = str(features[col]).replace("pfam_", "").split(".")[0]
        if pf in lookup.index:
            row = lookup.loc[pf]
            clan = row.get("clan", "") or ""
            full = row.get("full_name", "") or ""
        else:
            clan, full = "", ""
        role_for_col[int(col)] = assign_role(pf, clan=clan, full_name=full)

    by_role: dict[str, list[int]] = {}
    for col, role in role_for_col.items():
        by_role.setdefault(role, []).append(col)

    sizes = meta["leiden_hd_cluster"].value_counts()
    main = sorted([c for c, n in sizes.items() if n >= 20])

    rows = []
    for c in main:
        idx = np.where(meta["leiden_hd_cluster"].values == c)[0]
        nin = len(idx)
        sub = Xb[idx]
        for role, cols in by_role.items():
            if not cols:
                continue
            sums = np.asarray(sub[:, cols].sum(axis=1)).ravel()
            n_with = int(np.count_nonzero(sums))
            tot_hits = int(sums.sum())
            mean_hits = float(sums.mean())
            rows.append({
                "cluster": int(c),
                "cluster_size": nin,
                "role": role,
                "role_label": ROLE_SHORT_LABEL.get(role, role),
                "n_ie_with_role": n_with,
                "frac_ie_with_role": round(n_with / nin, 3),
                "mean_pfam_hits_per_ie": round(mean_hits, 2),
                "total_pfam_hits": tot_hits,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Top product= per cluster (hypothetical protein excluded)
# ---------------------------------------------------------------------------

_PRODUCT_RE = re.compile(r"product=([^|]+?)(?:$|\|)")


def _norm_product(p: str) -> str:
    p = (p or "").strip()
    pl = p.lower()
    if pl in {"hypothetical protein", "hypothetical", ""}:
        return ""
    p = re.sub(r"\bdomain-containing protein\b", "domain protein", pl)
    p = re.sub(r"\bcontaining protein\b", "containing", p)
    return p


def compute_top_products(meta: pd.DataFrame, top_n: int = 12) -> pd.DataFrame:
    long = pd.read_csv(PATH_PROTEIN_LONG, sep="\t")
    long["product"] = long["cluster_id"].astype(str).apply(lambda s: _PRODUCT_RE.search(s).group(1).strip() if _PRODUCT_RE.search(s) else "")
    long["product_norm"] = long["product"].apply(_norm_product)
    long = long[long["product_norm"] != ""]
    long["mge_id"] = long["mge_id"].astype(str)

    meta_loc = meta[["rep_id", "leiden_hd_cluster"]].copy()
    meta_loc["rep_id"] = meta_loc["rep_id"].astype(str)
    j = long.merge(meta_loc, left_on="mge_id", right_on="rep_id", how="inner")

    counts = j.groupby(["leiden_hd_cluster", "product_norm"]).size().reset_index(name="n_proteins")
    n_ie_per_cluster = j.groupby("leiden_hd_cluster")["mge_id"].nunique().to_dict()
    counts["n_ie_in_cluster"] = counts["leiden_hd_cluster"].map(n_ie_per_cluster)
    counts["per_ie"] = counts["n_proteins"] / counts["n_ie_in_cluster"]
    counts = counts.sort_values(["leiden_hd_cluster", "per_ie"], ascending=[True, False])
    top = counts.groupby("leiden_hd_cluster").head(top_n).reset_index(drop=True)
    return top


# ---------------------------------------------------------------------------
# 5. Host phylum, lengths, and IE atlas scores, joined into one summary
# ---------------------------------------------------------------------------


ATLAS_SCORE_COLS = ["IE_core_score", "phage_score", "inovirus_score", "conjugation_score", "replication_score"]


def attach_metadata(meta: pd.DataFrame) -> pd.DataFrame:
    df = meta.copy()
    df["sample_id"] = df["sample_id"].astype(str)
    tax = pd.read_csv(PATH_GTDB_TAX, sep="\t", dtype=str).fillna("")
    sample_col = "sample_id" if "sample_id" in tax.columns else tax.columns[0]
    cand_cols = [c for c in tax.columns if c == "host_gtdb_phylum"] or [
        c for c in tax.columns if c.lower() == "phylum"
    ]
    if cand_cols:
        tax = tax.rename(columns={cand_cols[0]: "host_gtdb_phylum"})
        tax_sub = tax[[sample_col, "host_gtdb_phylum"]].rename(columns={sample_col: "sample_id"})
        df = df.merge(tax_sub.drop_duplicates("sample_id"), on="sample_id", how="left")
    if "host_gtdb_phylum" not in df.columns:
        df["host_gtdb_phylum"] = ""

    atlas = pd.read_csv(PATH_ATLAS, sep="\t", dtype=str, low_memory=False)
    keep_cols = ["element_id"] + [c for c in ATLAS_SCORE_COLS if c in atlas.columns]
    atlas_num = atlas[keep_cols].copy()
    for c in ATLAS_SCORE_COLS:
        if c in atlas_num.columns:
            atlas_num[c] = pd.to_numeric(atlas_num[c], errors="coerce")
    atlas_num["element_id"] = atlas_num["element_id"].astype(str)
    df = df.merge(atlas_num, left_on="ie_id", right_on="element_id", how="left")
    return df


def cluster_summary(df: pd.DataFrame, hallmark: pd.DataFrame, products: pd.DataFrame) -> pd.DataFrame:
    """Final table: one row per cluster, with the key numbers."""
    sizes = df["leiden_hd_cluster"].value_counts()
    main = sorted([c for c, n in sizes.items() if n >= 20])
    rows = []
    for c in main:
        sub = df[df["leiden_hd_cluster"] == c]
        nin = len(sub)
        # lengths
        lens = pd.to_numeric(sub["ie_len_nt"], errors="coerce").dropna()
        len_med = float(lens.median()) if not lens.empty else float("nan")
        len_q25 = float(lens.quantile(0.25)) if not lens.empty else float("nan")
        len_q75 = float(lens.quantile(0.75)) if not lens.empty else float("nan")
        # host phylum
        gtdb_sub = sub[sub["group"].astype(str) == "gtdb"]
        phyl = gtdb_sub["host_gtdb_phylum"].dropna()
        phyl = phyl[phyl != ""]
        if not phyl.empty:
            vc = phyl.value_counts()
            top_phylum = vc.index[0]
            top_phylum_frac = round(float(vc.iloc[0] / len(phyl)), 2)
            distinct_phyla = int(phyl.nunique())
        else:
            top_phylum, top_phylum_frac, distinct_phyla = "", float("nan"), 0
        # group composition
        gcomp = sub["group"].value_counts(normalize=True).round(2).to_dict()
        # IE Atlas
        atlas_med = {}
        for col in ATLAS_SCORE_COLS:
            vals = pd.to_numeric(sub[col], errors="coerce") if col in sub.columns else pd.Series(dtype=float)
            atlas_med[col] = round(float(vals.median()), 3) if not vals.dropna().empty else 0.0
        # roles
        h = hallmark[hallmark["cluster"] == c].set_index("role")
        def fr(role):
            return float(h.loc[role, "frac_ie_with_role"]) if role in h.index else 0.0
        # products: top 3, hypothetical protein excluded
        p = products[products["leiden_hd_cluster"] == c].head(5)
        top_prod = "; ".join(f"{r['product_norm']} ({round(r['per_ie'],2)}/IE)" for _, r in p.iterrows())
        rows.append({
            "cluster": c,
            "n_ie": nin,
            "ie_len_median": int(len_med) if not np.isnan(len_med) else None,
            "ie_len_q25": int(len_q25) if not np.isnan(len_q25) else None,
            "ie_len_q75": int(len_q75) if not np.isnan(len_q75) else None,
            "group_composition": "; ".join(f"{k}:{v}" for k, v in gcomp.items()),
            "top_host_phylum": top_phylum,
            "top_host_phylum_frac": top_phylum_frac,
            "n_distinct_host_phyla": distinct_phyla,
            "median_IE_core_score": atlas_med.get("IE_core_score", 0.0),
            "median_replication_score": atlas_med.get("replication_score", 0.0),
            "frac_integrase_tyrosine": round(fr("integrase_tyrosine"), 2),
            "frac_integrase_serine": round(fr("integrase_serine"), 2),
            "frac_transposase": round(fr("transposase"), 2),
            "frac_replication": round(fr("replication"), 2),
            "frac_methyltransferase": round(fr("methyltransferase"), 2),
            "frac_restriction_T1": round(fr("restriction_T1"), 2),
            "frac_defense_STAND": round(fr("defense_STAND"), 2),
            "frac_toxin_antitoxin": round(fr("toxin_antitoxin"), 2),
            "frac_regulator_HTH": round(fr("regulator_HTH"), 2),
            "frac_aaa_atpase": round(fr("aaa_atpase"), 2),
            "frac_tpr_repeat": round(fr("tpr_repeat"), 2),
            "frac_helicase": round(fr("helicase"), 2),
            "frac_zn_ribbon": round(fr("zn_ribbon"), 2),
            "top_named_products": top_prod,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 6. PADLOC defense systems per cluster
# ---------------------------------------------------------------------------
#
# Where PADLOC output lives, by group:
#   * gtdb           — results_organized/gtdb/<sample>/padloc.csv
#                       (mirror of padloc_pipeline/results/<sample>/results_padloc.csv)
#   * thermaceae     — padloc_pipeline/results_thermaceae/<sample>/results_thermaceae_padloc.csv
#                       (there is no file under results_organized/thermaceae/)
#   * deinococcales  — both paths exist; the first one that is present is used
#   * bins           — PADLOC was not run; padloc_status='empty' is expected.
#
# seqid in padloc.csv matches ie_id (PADLOC was run on mge_region.fa,
# where each IE is its own record).


def _load_padloc_for_sample(group: str, sample_id: str) -> tuple[pd.DataFrame, str]:
    """Return (df, sample_status) for (group, sample_id).

    sample_status:
      * 'not_run'           — none of the expected padloc files exist; PADLOC was not run;
      * 'no_systems_sample' — a file exists but is 0 bytes (PADLOC ran and found no systems in the sample);
      * 'has_systems'       — a non-empty file; df holds the parsed hits.

    PADLOC wrote a 0-byte CSV when the sample had no defense systems, while the
    Prodigal and HMM steps had still run (orfs.faa and orfs.gff sit beside it).
    That has to stay distinct from "PADLOC was never started for this sample"
    (the bins cohort). Otherwise the denominator of "fraction of IE with a system"
    is inflated and cohorts cannot be compared.
    """
    files = _padloc_paths(group, sample_id)
    found_any = False
    for p in files:
        if not p.is_file():
            continue
        found_any = True
        if p.stat().st_size > 0:
            try:
                return pd.read_csv(p), "has_systems"
            except Exception:
                continue
    return pd.DataFrame(), ("no_systems_sample" if found_any else "not_run")


def _scan_padloc_per_ie(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """PADLOC systems inside each IE window.

    Returns (per_ie, hits_long):
      * per_ie: ie_id, padloc_status, n_padloc_rows, n_padloc_systems, padloc_systems
        where padloc_status is:
          - 'not_run'           — PADLOC was not run for this sample;
          - 'no_systems_sample' — it ran, and the whole sample has 0 systems;
          - 'no_rows_for_ie'    — the sample has systems, but not inside this IE window;
          - 'ok'                — at least one system falls inside this IE window.
      * hits_long: ie_id, group, sample_id, system — one row per IE system.
    """
    rows: list[dict] = []
    detail_rows: list[dict] = []
    cache: dict[tuple[str, str], tuple[pd.DataFrame, str]] = {}
    for _, r in meta.iterrows():
        ie_id = str(r["ie_id"])
        key = (str(r["group"]), str(r["sample_id"]))
        if key not in cache:
            cache[key] = _load_padloc_for_sample(*key)
        pad, sample_status = cache[key]
        if sample_status != "has_systems":
            rows.append({
                "ie_id": ie_id,
                "padloc_status": sample_status,  # not_run | no_systems_sample
                "n_padloc_rows": 0,
                "n_padloc_systems": 0,
                "padloc_systems": "",
            })
            continue
        if "seqid" not in pad.columns:
            rows.append({
                "ie_id": ie_id,
                "padloc_status": "no_systems_sample",
                "n_padloc_rows": 0,
                "n_padloc_systems": 0,
                "padloc_systems": "",
            })
            continue
        sub = pad[pad["seqid"].astype(str) == ie_id]
        if len(sub) == 0:
            rows.append({
                "ie_id": ie_id,
                "padloc_status": "no_rows_for_ie",
                "n_padloc_rows": 0,
                "n_padloc_systems": 0,
                "padloc_systems": "",
            })
            continue
        systems = sorted({str(x) for x in sub["system"].dropna().unique()})
        rows.append({
            "ie_id": ie_id,
            "padloc_status": "ok",
            "n_padloc_rows": int(len(sub)),
            "n_padloc_systems": int(len(systems)),
            "padloc_systems": "|".join(systems),
        })
        for sys_name in systems:
            detail_rows.append({
                "ie_id": ie_id,
                "group": key[0],
                "sample_id": key[1],
                "system": sys_name,
            })
    return pd.DataFrame(rows), pd.DataFrame(detail_rows)


def compute_padloc_summary(
    meta: pd.DataFrame,
    *,
    min_cluster_size: int = 20,
    top_n_systems: int = 5,
) -> dict[str, pd.DataFrame]:
    """Aggregate PADLOC systems by leiden_hd_cluster.

    Arguments:
      meta: DataFrame with ie_id, group, sample_id, leiden_hd_cluster
            (what load_inputs / attach_metadata return).
      min_cluster_size: keep only clusters of at least 20 IE, as in
                        compute_markers / compute_hallmark_composition.
      top_n_systems: how many top systems to list in the summary.

    Returns a dict: per_ie, cluster_long, cluster_summary.
    """
    needed = {"ie_id", "group", "sample_id", "leiden_hd_cluster"}
    miss = needed - set(meta.columns)
    if miss:
        raise ValueError(f"meta missing columns: {sorted(miss)}")

    per_ie, hits = _scan_padloc_per_ie(meta[["ie_id", "group", "sample_id"]])

    meta_loc = meta[["ie_id", "leiden_hd_cluster"]].copy()
    meta_loc["ie_id"] = meta_loc["ie_id"].astype(str)
    per_ie = per_ie.merge(meta_loc, on="ie_id", how="left")

    sizes = meta["leiden_hd_cluster"].value_counts()
    main = sorted([int(c) for c, n in sizes.items() if n >= min_cluster_size])

    # Long table: one row per (cluster, system).
    if hits.empty:
        cluster_long = pd.DataFrame(
            columns=[
                "cluster",
                "cluster_size",
                "system",
                "n_ie_with_system",
                "frac_ie_with_system",
                "n_total_hits",
            ]
        )
    else:
        hits = hits.merge(meta_loc, on="ie_id", how="left")
        hits = hits[hits["leiden_hd_cluster"].notna()]
        hits["leiden_hd_cluster"] = hits["leiden_hd_cluster"].astype(int)
        hits = hits[hits["leiden_hd_cluster"].isin(main)]

        n_ie_per_cluster = (
            meta[meta["leiden_hd_cluster"].isin(main)]
            .groupby("leiden_hd_cluster")["ie_id"]
            .nunique()
            .to_dict()
        )

        # number of PADLOC rows = number of gene hits assigned to a system
        agg = (
            hits.groupby(["leiden_hd_cluster", "system"])
            .agg(n_ie_with_system=("ie_id", "nunique"), n_total_hits=("ie_id", "size"))
            .reset_index()
            .rename(columns={"leiden_hd_cluster": "cluster"})
        )
        agg["cluster_size"] = agg["cluster"].map(n_ie_per_cluster).astype(int)
        agg["frac_ie_with_system"] = (agg["n_ie_with_system"] / agg["cluster_size"]).round(3)
        cluster_long = agg.sort_values(
            ["cluster", "n_ie_with_system"], ascending=[True, False]
        ).reset_index(drop=True)[
            [
                "cluster",
                "cluster_size",
                "system",
                "n_ie_with_system",
                "frac_ie_with_system",
                "n_total_hits",
            ]
        ]

    # Summary: one row per cluster.
    rows = []
    for c in main:
        sub = per_ie[per_ie["leiden_hd_cluster"] == c]
        n_ie = int(len(sub))
        n_with = int((sub["n_padloc_systems"].fillna(0) > 0).sum())
        frac_with = round(n_with / n_ie, 3) if n_ie else 0.0
        cl_systems = cluster_long[cluster_long["cluster"] == c]
        n_distinct = int(cl_systems["system"].nunique())
        top = cl_systems.head(top_n_systems)
        top_str = "; ".join(
            f"{r['system']} ({int(r['n_ie_with_system'])} IE, {r['frac_ie_with_system']:.2f})"
            for _, r in top.iterrows()
        )
        rows.append({
            "cluster": c,
            "n_ie": n_ie,
            "n_ie_with_padloc": n_with,
            "frac_ie_with_padloc": frac_with,
            "n_distinct_systems": n_distinct,
            "top_systems": top_str,
        })
    cluster_summary_df = pd.DataFrame(rows)

    return {
        "per_ie": per_ie,
        "cluster_long": cluster_long,
        "cluster_summary": cluster_summary_df,
    }


def run_all(save: bool = True):
    X, features, meta_raw = load_inputs()
    markers_raw = compute_markers(X, features, meta_raw)
    markers = annotate_markers(markers_raw)
    if save:
        markers.to_csv(PATH_MARKERS_OUT, sep="\t", index=False)
    hallmark = compute_hallmark_composition(X, features, meta_raw)
    if save:
        hallmark.to_csv(PATH_HALLMARK_OUT, sep="\t", index=False)
    products = compute_top_products(meta_raw)
    if save:
        products.to_csv(PATH_PRODUCTS_OUT, sep="\t", index=False)
    meta = attach_metadata(meta_raw)
    summary = cluster_summary(meta, hallmark, products)
    if save:
        summary.to_csv(PATH_SUMMARY_OUT, sep="\t", index=False)
    padloc = compute_padloc_summary(meta)
    if save:
        padloc["per_ie"].to_csv(PATH_PADLOC_PER_IE_OUT, sep="\t", index=False)
        padloc["cluster_long"].to_csv(PATH_PADLOC_CLUSTER_LONG_OUT, sep="\t", index=False)
        padloc["cluster_summary"].to_csv(PATH_PADLOC_CLUSTER_SUMMARY_OUT, sep="\t", index=False)
    return {
        "markers": markers,
        "hallmark": hallmark,
        "products": products,
        "summary": summary,
        "meta": meta,
        "padloc_per_ie": padloc["per_ie"],
        "padloc_cluster_long": padloc["cluster_long"],
        "padloc_cluster_summary": padloc["cluster_summary"],
    }


if __name__ == "__main__":
    out = run_all()
    print("markers:", out["markers"].shape)
    print("hallmark:", out["hallmark"].shape)
    print("products:", out["products"].shape)
    print("summary:", out["summary"].shape)
    print("padloc_per_ie:", out["padloc_per_ie"].shape)
    print("padloc_cluster_long:", out["padloc_cluster_long"].shape)
    print("padloc_cluster_summary:", out["padloc_cluster_summary"].shape)
