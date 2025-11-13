#!/usr/bin/env python3
"""
Merge attL/attR misc_feature annotations from an annotated .gbk into a Bakta .gbff.
"""

import argparse
import os
import sys
import warnings
from typing import Dict, List

from Bio import SeqIO, BiopythonWarning
from Bio.SeqFeature import SeqFeature

# allow import of top‑level scripts/logger.py
SCRIPT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

sys.path.append("../finder_pipeline/scripts")
from logger import Logger as LoggerWrapper

warnings.filterwarnings("ignore", category=BiopythonWarning)

_wrapper = LoggerWrapper(
    name="merge_att_features",
    log_to_console=True,
    log_to_file=True
)
logger = _wrapper.get_logger()
trace = _wrapper.trace_call


@trace
def extract_att_features(
    annotated_gbk: str
) -> Dict[str, List[SeqFeature]]:
    """
    Extract attL/attR misc_features from an annotated GBK (.gbk/.gbff).

    Args:
        annotated_gbk (str): Path to the GenBank file with existing
            misc_feature /note="attL"/"attR".

    Returns:
        Dict[str, List[SeqFeature]]: record.id → list of SeqFeature
    """
    logger.info(f"Loading annotated GBK: {annotated_gbk!r}")
    recs = list(SeqIO.parse(annotated_gbk, "genbank"))
    logger.info(f"Found {len(recs)} records in annotated file")

    feats_by_id: Dict[str, List[SeqFeature]] = {}
    total = 0
    for rec in recs:
        logger.debug(f"Record {rec.id!r}, {len(rec.features)} total features")
        att = [
            f for f in rec.features
            if f.type == "misc_feature"
            and "note" in f.qualifiers
            and any(n in ("attL", "attR") for n in f.qualifiers["note"])
        ]
        if att:
            feats_by_id[rec.id] = att
            logger.info(f"  → {len(att)} att-features from {rec.id!r}")
            total += len(att)

    logger.info(f"Total extracted: {total} att-features from {len(feats_by_id)} records")
    return feats_by_id


@trace
def merge_features(
    bakta_gbff: str,
    feats: Dict[str, List[SeqFeature]],
    out_path: str
) -> None:
    """
    Merge previously extracted att-features into the Bakta-generated .gbff.

    Args:
        bakta_gbff (str): Path to Bakta .gbff.
        feats (dict): Output of extract_att_features().
        out_path (str): Where to write the merged GBFF.
    """
    logger.info(f"Loading primary GBFF: {bakta_gbff!r}")
    main = list(SeqIO.parse(bakta_gbff, "genbank"))
    logger.info(f"Found {len(main)} records in primary GBFF")

    merged_count = 0
    recs_updated = 0
    for rec in main:
        if rec.id in feats:
            cnt = len(feats[rec.id])
            rec.features.extend(feats[rec.id])
            merged_count += cnt
            recs_updated += 1
            logger.info(f"Merged {cnt} features into record {rec.id!r}")

    missing = set(feats) - {r.id for r in main}
    if missing:
        logger.warning(f"Features for records not in GBFF: {sorted(missing)}")

    logger.info(f"Writing {merged_count} features into {recs_updated} records → {out_path!r}")
    SeqIO.write(main, out_path, "genbank")
    logger.info("Merge finished")


def main():
    parser = argparse.ArgumentParser(
        description="Merge attL/attR from .gbk into Bakta .gbff"
    )
    parser.add_argument(
        "--att", required=True,
        help="Annotated .gbk (with misc_feature attL/attR)"
    )
    parser.add_argument(
        "--gbff", required=True,
        help="Bakta-generated .gbff to receive the features"
    )
    parser.add_argument(
        "--out", required=True,
        help="Output path for merged GBFF"
    )
    args = parser.parse_args()

    logger.info(f"Args: att={args.att!r}, gbff={args.gbff!r}, out={args.out!r}")
    for fn in (args.att, args.gbff):
        if not os.path.isfile(fn):
            logger.error(f"File not found: {fn!r}")
            sys.exit(1)

    feats = extract_att_features(args.att)
    if not feats:
        logger.warning("No att-features extracted; copying GBFF unchanged")
        recs = list(SeqIO.parse(args.gbff, "genbank"))
        SeqIO.write(recs, args.out, "genbank")
        logger.info(f"Copied {args.gbff!r} → {args.out!r}")
    else:
        merge_features(args.gbff, feats, args.out)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Fatal: {e}", exc_info=True)
        sys.exit(1)