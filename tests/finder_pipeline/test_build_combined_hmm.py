"""Smoke tests for build_combined_hmm.py."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "finder_pipeline" / "scripts" / "build_combined_hmm.py"
PFAM_DIR = REPO_ROOT / "pfam"


def _tool_env() -> dict[str, str]:
    env = os.environ.copy()
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        bin_dir = str(Path(prefix) / "bin")
        env["PATH"] = f"{bin_dir}:{env.get('PATH', '')}"
    return env


class TestBuildCombinedHmm(unittest.TestCase):
    @unittest.skipUnless(PFAM_DIR.is_dir(), "pfam/ directory not found")
    @unittest.skipUnless(shutil.which("hmmpress", path=_tool_env().get("PATH", "")), "hmmpress not in PATH")
    def test_combine_and_hmmpress(self):
        pfams = sorted(PFAM_DIR.glob("*.hmm"))
        self.assertGreaterEqual(len(pfams), 1, "Need at least one Pfam HMM")

        with tempfile.TemporaryDirectory() as tmp:
            out_hmm = Path(tmp) / "combined.hmm"
            log_file = Path(tmp) / "build.log"
            cmd = [
                sys.executable,
                str(SCRIPT),
                "--output", str(out_hmm),
                "--log", str(log_file),
            ]
            for pfam in pfams[:2]:
                cmd.extend(["--pfam", str(pfam)])

            result = subprocess.run(cmd, capture_output=True, text=True, env=_tool_env())
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(out_hmm.is_file())
            self.assertGreater(out_hmm.stat().st_size, 0)
            for suffix in (".h3f", ".h3p", ".h3i", ".h3m"):
                self.assertTrue(Path(str(out_hmm) + suffix).is_file(), f"Missing hmmpress index {suffix}")


if __name__ == "__main__":
    unittest.main()
