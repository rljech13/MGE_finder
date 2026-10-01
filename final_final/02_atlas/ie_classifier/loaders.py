from dataclasses import dataclass
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


@dataclass
class ParsedCDS:
    locus_tag: str
    start: int
    end: int
    strand: str
    product: str


def _parse_location(loc_str: str) -> Tuple[int, int, str]:
    loc_str = loc_str.strip()
    strand = "+"
    if loc_str.startswith("complement(") and loc_str.endswith(")"):
        loc_str = loc_str[11:-1].strip()
        strand = "-"
    nums = [int(x) for x in re.findall(r"\d+", loc_str)]
    if len(nums) >= 2:
        return min(nums), max(nums), strand
    if len(nums) == 1:
        return nums[0], nums[0], strand
    return 0, 0, strand


def parse_gbff_cds(gbff_path: Path) -> List[ParsedCDS]:
    """Minimal GenBank FEATURES parser for CDS + /locus_tag + /product."""
    text = gbff_path.read_text(encoding="utf-8", errors="replace")
    if "FEATURES" not in text or "ORIGIN" not in text:
        return []
    feat = text.split("FEATURES", 1)[1].split("ORIGIN", 1)[0]
    lines = feat.splitlines()
    out: List[ParsedCDS] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"^\s{5}(\S+)\s+(.+)$", line)
        if not m:
            i += 1
            continue
        ftype, rest = m.group(1), m.group(2).strip()
        if ftype != "CDS":
            i += 1
            continue
        loc_part = rest.split("/")[0].strip()
        start, end, strand = _parse_location(loc_part)
        quals: Dict[str, str] = {}
        i += 1
        while i < len(lines):
            q = lines[i]
            if re.match(r"^\s{5}\S+\s+", q) and not q.lstrip().startswith("/"):
                break
            mm = re.match(r'^\s+/([^=]+)="([^"]*)"', q)
            if mm:
                quals[mm.group(1)] = mm.group(2)
            else:
                mm2 = re.match(r"^\s+/([^=]+)=([^\s].*)$", q)
                if mm2:
                    quals[mm2.group(1)] = mm2.group(2).strip()
            i += 1
        lt = quals.get("locus_tag", "").strip()
        prod = quals.get("product", quals.get("note", "")).strip()
        if lt and start and end:
            out.append(ParsedCDS(locus_tag=lt, start=start, end=end, strand=strand, product=prod))
    out.sort(key=lambda x: (x.start, x.end))
    return out


def find_att_notes(gbff_path: Path) -> Dict[str, bool]:
    text = gbff_path.read_text(encoding="utf-8", errors="replace")
    low = text.lower()
    return {
        "attL_attR_present": "attl" in low and "attr" in low,
        "attL_or_attR": "attl" in low or "attr" in low,
    }


def gbff_path_from_rep_index_row(merged_dedup_dir: Path, file_col: str) -> Path:
    return merged_dedup_dir / file_col


def load_merge_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def merge_lookup_for_genome(df: pd.DataFrame, genome_id: str) -> Dict[str, Dict[str, str]]:
    sub = df[df["genome_id"].astype(str) == str(genome_id)]
    lut: Dict[str, Dict[str, str]] = {}
    for _, row in sub.iterrows():
        lt = str(row.get("locus_tag", "")).strip()
        if not lt:
            continue
        lut[lt] = {k: ("" if pd.isna(row[k]) else str(row[k])) for k in row.index}
    return lut


def hydrate_genes(
    cds_list: List[ParsedCDS],
    merge_lut: Dict[str, Dict[str, str]],
) -> List[Dict[str, Any]]:
    genes: List[Dict[str, Any]] = []
    for j, cds in enumerate(cds_list):
        row = merge_lut.get(cds.locus_tag, {})
        genes.append(
            {
                "idx": j,
                "locus_tag": cds.locus_tag,
                "start": cds.start,
                "end": cds.end,
                "strand": cds.strand,
                "product": cds.product or row.get("product", "") or row.get("predicted_protein_name", ""),
                "pfam_id": row.get("pfam_id", ""),
                "emapper_PFAMs": row.get("emapper_PFAMs", ""),
                "eggNOG_Description": row.get("eggNOG_Description", ""),
                "seed_eggNOG_ortholog": row.get("seed_eggNOG_ortholog", ""),
                "KEGG_ko": row.get("KEGG_ko", ""),
                "family_id": row.get("family_id", ""),
            }
        )
    return genes


def normalize_rep_cluster_key(cluster_id: str) -> str:
    """Strip |sample=... suffix so truth (locus-only) matches rep_index cluster_id."""
    s = str(cluster_id).strip()
    if "|sample=" in s:
        return s.split("|sample=", 1)[0].strip()
    return s


def build_rep_locus_file_lookup(rep_df: pd.DataFrame) -> Dict[str, str]:
    """Map normalized locus/member key -> rep gbff filename (one row per key)."""
    out: Dict[str, str] = {}
    if rep_df is None or rep_df.empty:
        return out
    for _, row in rep_df.iterrows():
        rid = normalize_rep_cluster_key(str(row.get("cluster_id", "") or ""))
        if not rid:
            continue
        fn = str(row.get("file", "") or "").strip()
        if fn and rid not in out:
            out[rid] = fn
    return out


def truth_row_meta(row: pd.Series) -> Dict[str, Any]:
    ie_len = ""
    for col in ("ie_len_nt", "ie_length_nt", "mge_len", "ie_len"):
        v = row.get(col, "")
        if v not in (None, "", "nan"):
            ie_len = str(v).strip()
            break
    return {
        "element_id": str(row.get("mge_canonical_id", "")).strip(),
        "genome_id": str(row.get("matrix_row_genome_id", "")).strip(),
        "leiden_hd_cluster": str(row.get("leiden_hd_cluster", "")).strip(),
        "sample_id": str(row.get("sample_id", "")).strip(),
        "tRNA_type_assigned": str(row.get("tRNA_type_assigned", "")).strip(),
        "integrase_id_lookup": str(row.get("integrase_id_lookup", "")).strip(),
        "finder_data_missing": str(row.get("finder_data_missing", "")).strip(),
        "ie_len_nt": ie_len,
    }


def truth_cluster_key(row: pd.Series) -> str:
    for k in ("mge_canonical_id", "cluster_id", "member"):
        v = str(row.get(k, "") or "").strip()
        if v:
            return v
    return ""


def rep_filename_for_lookup(rep_file_by_locus: Dict[str, str], cluster_key: str) -> str:
    if not cluster_key or not rep_file_by_locus:
        return ""
    k = normalize_rep_cluster_key(cluster_key)
    return str(rep_file_by_locus.get(k, "") or "").strip()


def rep_filename_for_cluster(rep_df: pd.DataFrame, cluster_key: str) -> str:
    """Deprecated path: builds lookup each call — prefer build_rep_locus_file_lookup + rep_filename_for_lookup."""
    return rep_filename_for_lookup(build_rep_locus_file_lookup(rep_df), cluster_key)
