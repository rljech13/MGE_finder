#!/usr/bin/env python3
"""Classify IE via ≤40 kb attL BLAST search window (Bakta slice on contig)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mge_paths import FF, REPO

SCRIPTS = FF / "02_atlas"
NB7 = REPO / "results_organized/genome_clustering/project/notebook7_extended"

WINDOW_BP = 40_000
OUT = FF / "artifacts/blast_window_40kb_atlas_classification.tsv"

sys.path.insert(0, str(SCRIPTS))
from ie_classifier.classify import classify_from_genes  # noqa: E402
from ie_classifier.config import load_config  # noqa: E402
from ie_classifier.loaders import hydrate_genes, merge_lookup_for_genome, parse_gbff_cds  # noqa: E402
from ie_classifier.schema import TSV_COLUMNS, flatten_record  # noqa: E402
from ie_classifier.rules import coarse_label_from_primary  # noqa: E402


def find_merged_path(sample_id: str, source: str) -> Path | None:
    ro = REPO / "results_organized"
    bakta = REPO / "bakta_pipeline"
    for base in (ro / source / sample_id, bakta / f"results_bakta_{source}" / sample_id):
        p = base / "annotation.gbff.merged"
        if p.is_file():
            return p.resolve()
    return None


def blast_window(trna_start: int, trna_end: int, trna_strand: str) -> tuple[int, int]:
    if str(trna_strand).strip() == "+":
        return int(trna_end), int(trna_end) + WINDOW_BP
    return max(1, int(trna_start) - WINDOW_BP), int(trna_start)


def pick_contig_record(records: list, contig: str) -> SeqRecord | None:
    for r in records:
        if (r.id or "").split()[0] == contig or contig in (r.id or ""):
            return r
    for r in records:
        if contig in (r.name or "") or contig in (r.description or ""):
            return r
    return None


def features_in_window(rec: SeqRecord, lo: int, hi: int) -> list:
    out = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        try:
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
        except Exception:
            continue
        if end >= lo and start <= hi:
            out.append(feat)
    return out


def cds_from_features(feats: list) -> list:
    """Minimal CDS list compatible with parse_gbff_cds / hydrate."""
    from ie_classifier.loaders import ParsedCDS

    rows = []
    for f in feats:
        quals = f.qualifiers
        lt = (quals.get("locus_tag") or quals.get("protein_id") or [""])[0]
        prod = (quals.get("product") or quals.get("Note") or [""])[0]
        start = int(f.location.start) + 1
        end = int(f.location.end)
        strand = "+" if f.location.strand != -1 else "-"
        if lt:
            rows.append(ParsedCDS(locus_tag=lt, start=start, end=end, strand=strand, product=prod))
    rows.sort(key=lambda x: (x.start, x.end))
    return rows


def main() -> None:
    replay = pd.read_csv(
        REPO / "snapshots/v3pS_pre_rebuild/replay_summary.tsv", sep="\t", low_memory=False,
    )
    replay = replay.drop_duplicates("integrase_id", keep="first")

    rep = pd.read_csv(REPO / "results_organized/merged_deduplicated/rep_index.tsv", sep="\t", dtype=str)
    rep["integrase_id"] = (
        rep["cluster_id"].astype(str).str.split("|sample=", regex=False).str[0]
        .str.split(":", n=1, regex=False).str[0]
    )
    rep_lut = rep.drop_duplicates("integrase_id", keep="first").set_index("integrase_id")

    merge_df = pd.read_csv(NB7 / "locus_eggnog_merged.tsv", sep="\t", dtype=str, low_memory=False)
    cfg = load_config(FF / "02_atlas" / "ie_classifier" / "ie_atlas_labels.yaml")

    merged_cache: dict[str, list] = {}
    rows_out = []
    total = len(replay)
    for i, (_, row) in enumerate(replay.iterrows(), start=1):
        iid = str(row["integrase_id"])
        if iid not in rep_lut.index:
            continue
        r = rep_lut.loc[iid]
        rep_id = str(r["rep_id"])
        sample_id = str(r["sample_id"])
        source = str(r["source"])
        contig = str(row["contig"])
        try:
            ts, te = int(row["trna_start"]), int(row["trna_end"])
        except Exception:
            continue
        lo, hi = blast_window(ts, te, str(row["trna_strand"]))

        cache_key = f"{sample_id}|{source}"
        if cache_key not in merged_cache:
            mp = find_merged_path(sample_id, source)
            if mp is None:
                merged_cache[cache_key] = []
            else:
                try:
                    merged_cache[cache_key] = list(SeqIO.parse(mp, "genbank"))
                except Exception:
                    merged_cache[cache_key] = []
        recs = merged_cache[cache_key]
        contig_rec = pick_contig_record(recs, contig)
        if contig_rec is None:
            continue

        seq_len = len(contig_rec.seq)
        lo = max(1, lo)
        hi = min(seq_len, hi)
        if hi <= lo:
            continue

        feats = features_in_window(contig_rec, lo, hi)
        cds = cds_from_features(feats)
        if not cds:
            continue

        lut = merge_lookup_for_genome(merge_df, rep_id)
        genes = hydrate_genes(cds, lut)
        element_id = f"{iid}:{contig}:{lo}-{hi}"
        meta = {
            "element_id": element_id,
            "genome_id": rep_id,
            "sample_id": sample_id,
            "leiden_hd_cluster": "",
            "tRNA_type_assigned": "",
            "finder_data_missing": "",
            "ie_len_nt": str(hi - lo + 1),
        }
        rec = classify_from_genes(genes, meta, cfg, att={"attL_attR_present": False, "attL_or_attR": False})
        flat = flatten_record(rec)
        flat["integrase_id"] = iid
        flat["coarse_label"] = coarse_label_from_primary(flat.get("primary_label", ""))
        flat["blast_window_lo"] = lo
        flat["blast_window_hi"] = hi
        rows_out.append(flat)

        if i == 1 or i % 500 == 0 or i == total:
            print(f"  progress {i}/{total}  classified={len(rows_out)}", flush=True)

    out_df = pd.DataFrame(rows_out)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT, sep="\t", index=False)
    print(f"→ {OUT} ({len(out_df):,} rows)")


if __name__ == "__main__":
    main()
