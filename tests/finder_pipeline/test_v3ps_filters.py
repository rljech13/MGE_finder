"""Unit tests for V_3pS strict selection and the manuscript confidence cascade."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "finder_pipeline" / "scripts"))

from v3ps_filters import (
    FilterThresholds,
    classify_prodigal_overlap,
    evaluate_ie_candidate,
    exact_anchored_run,
    select_v3ps_strict_hit,
)

# attL and tRNA region sharing a 20 nt perfect run from the anchored end (no repeats)
_ATTL = "AACCGGTTAACCGGTTACGT"
_TRNA = "NNNN" + _ATTL


def _raw(length=20, mismatch=0, gapped=False):
    """Anchored at the tRNA 3' end. Ungapped unless ``gapped`` stretches the query span."""
    qend = 70
    qstart = 1 if gapped else qend - length + 1
    sstart, send = 100, 100 + length - 1
    return pd.DataFrame([
        {"qseqid": "i1:c1:1-70:+", "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
         "length": length, "mismatch": mismatch, "gapopen": 0, "pident": 100, "wstart": 1000},
    ])


def _audit(raw, **kw):
    args = dict(
        integrase_id="i1", contig="c1", trna_start=1, trna_end=70, trna_strand="+", trna_len=70,
        integrase_start=1000, integrase_end=1902, ie_id="i1:c1:100-3500", ie_len_nt=2500,
        ie_has_n=False, raw_hits=raw, cds_by_contig={"c1": [(2000, 3000, "+", "cds1")]},
        thresholds=FilterThresholds(), attl_seq=_ATTL, trna_seq=_TRNA,
    )
    args.update(kw)
    return evaluate_ie_candidate(**args)


class TestV3psStrict(unittest.TestCase):
    def test_select_longest_strict_hit(self):
        raw = pd.DataFrame([
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 68, "sstart": 100, "send": 200,
             "length": 68, "wstart": 1000},
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 69, "sstart": 100, "send": 220,
             "length": 85, "wstart": 1000},
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 50, "sstart": 100, "send": 150,
             "length": 50, "wstart": 1000},
        ])
        hit = select_v3ps_strict_hit(raw, trna_len=70, trna_strand="+", shift=3)
        self.assertIsNotNone(hit)
        self.assertEqual(hit["best_attl_len_bp"], 85)

    def test_rejects_short_3p_anchor(self):
        raw = pd.DataFrame([
            {"qseqid": "i1:c1:1-70:+", "qstart": 1, "qend": 60, "sstart": 100, "send": 200,
             "length": 60, "wstart": 1000},
        ])
        self.assertIsNone(select_v3ps_strict_hit(raw, trna_len=70, trna_strand="+", shift=3))

    def test_prodigal_intergenic(self):
        info = classify_prodigal_overlap(500, 520, [(600, 700, "+", "cds1")])
        self.assertEqual(info["prodigal_hit_class"], "intergenic")

    def test_prodigal_partial_overlap(self):
        info = classify_prodigal_overlap(500, 520, [(510, 700, "+", "cds1")])
        self.assertEqual(info["prodigal_hit_class"], "partial_CDS")

    def test_exact_run_from_anchor(self):
        self.assertEqual(exact_anchored_run(_ATTL, _TRNA, "+"), 20)
        broken = "NNNN" + "C" * len(_ATTL)
        self.assertLess(exact_anchored_run(_ATTL, broken, "+"), 8)

    def test_confident_passes_all_filters(self):
        audit = _audit(_raw())
        self.assertTrue(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "")
        self.assertEqual(audit["attL_tier"], "exact")

    def test_confident_rejects_short_integrase(self):
        audit = _audit(_raw(), integrase_start=1000, integrase_end=1500)
        self.assertFalse(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "integrase_too_short")

    def test_exact_run_below_candidate_minimum(self):
        audit = _audit(_raw(length=20), attl_seq="AAAA", trna_seq="TTTT")
        self.assertFalse(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "exact_run_too_short")

    def test_relaxed_length_keeps_short_exact_run(self):
        # exact run 12 is below 14, BLAST length 20 is at least 17; mismatches are allowed
        attl, trna = "AAAAAAAA" + _ATTL[8:], "TTTTTTTT" + _TRNA[8:]
        self.assertEqual(exact_anchored_run(attl, trna, "+"), 12)
        audit = _audit(_raw(length=20, mismatch=5), attl_seq=attl, trna_seq=trna)
        self.assertTrue(audit["passed_confident"])
        self.assertEqual(audit["attL_tier"], "relaxed")

    def test_length_gate_rejects_when_both_tiers_fail(self):
        attl, trna = "AAAAAAAA" + _ATTL[8:], "TTTTTTTT" + _TRNA[8:]
        audit = _audit(_raw(length=16), attl_seq=attl, trna_seq=trna)
        self.assertFalse(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "attl_too_short")

    def test_gapped_alignment_rejected(self):
        audit = _audit(_raw(gapped=True))
        self.assertFalse(audit["passed_confident"])
        self.assertEqual(audit["reject_reason"], "attl_gapped")


if __name__ == "__main__":
    unittest.main()
