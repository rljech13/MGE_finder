"""Tests for prepare_fastas N-filter and bin ingestion."""
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "finder_pipeline" / "scripts"))

from prepare_fastas import process_input_sources


class TestPrepareFastas(unittest.TestCase):
    def _write(self, path: Path, content: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    def test_reject_ambiguous_n_skips_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src" / "GCA_test" / "GCA_test_genomic.fna"
            out = Path(tmp) / "out"
            self._write(src, ">contig\nACGTN\n")

            processed = process_input_sources(
                [str(src.parent.parent)],
                str(out),
                reject_ambiguous_n=True,
            )
            self.assertEqual(processed, [])
            self.assertFalse((out / "GCA_test.fna").exists())

    def test_accepts_clean_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src" / "GCA_clean" / "GCA_clean_genomic.fna"
            out = Path(tmp) / "out"
            self._write(src, ">contig\nACGTACGT\n")

            processed = process_input_sources(
                [str(src.parent.parent)],
                str(out),
                reject_ambiguous_n=True,
            )
            self.assertEqual(processed, ["GCA_clean"])
            self.assertTrue((out / "GCA_clean.fna").exists())

    def test_bin_pattern(self):
        with tempfile.TemporaryDirectory() as tmp:
            src_dir = Path(tmp) / "bins"
            out = Path(tmp) / "out"
            self._write(src_dir / "bin_1.fa", ">bin\nATGC\n")

            processed = process_input_sources(
                [str(src_dir)],
                str(out),
                bins_pattern="*.fa",
                reject_ambiguous_n=False,
            )
            self.assertEqual(processed, ["bin_1"])
            self.assertTrue((out / "bin_1.fna").exists())


if __name__ == "__main__":
    unittest.main()
