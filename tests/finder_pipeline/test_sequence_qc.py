"""Unit tests for IE_finder sequence QC helpers."""
import tempfile
import unittest
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "finder_pipeline" / "scripts"))

from sequence_qc import fasta_contains_ambiguous_n, parse_fasta_lengths_and_n_flags


class TestSequenceQC(unittest.TestCase):
    def _write_fasta(self, path: Path, records: list[tuple[str, str]]) -> None:
        lines = []
        for header, seq in records:
            lines.append(f">{header}")
            lines.append(seq)
        path.write_text("\n".join(lines) + "\n")

    def test_no_ambiguous_n(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "clean.fa"
            self._write_fasta(fa, [("contig1", "ACGTACGT"), ("contig2", "TTTT")])
            self.assertFalse(fasta_contains_ambiguous_n(fa))

    def test_has_ambiguous_n(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "dirty.fa"
            self._write_fasta(fa, [("contig1", "ACGTNACGT")])
            self.assertTrue(fasta_contains_ambiguous_n(fa))

    def test_parse_lengths_and_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "mixed.fa"
            self._write_fasta(
                fa,
                [("a", "ACGT"), ("b", "NNNN"), ("c", "ATGC")],
            )
            lengths, has_n = parse_fasta_lengths_and_n_flags(fa)
            self.assertEqual(lengths, {"a": 4, "b": 4, "c": 4})
            self.assertEqual(has_n, {"a": False, "b": True, "c": False})


if __name__ == "__main__":
    unittest.main()
