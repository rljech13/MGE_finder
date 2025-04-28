# padloc_wrapper.py
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import shutil
import glob

# --- parse arguments ---
parser = argparse.ArgumentParser(description="Padloc wrapper with logging")
parser.add_argument("--logger-dir", required=True)
parser.add_argument("--faa",        required=True)
parser.add_argument("--gff",        required=True)
parser.add_argument("--outdir",     required=True)
parser.add_argument("--log",        required=True)
parser.add_argument("--force",      action="store_true")
parser.add_argument("--cpu",        type=int, default=1)
args = parser.parse_args()

# --- import Logger ---
sys.path.insert(0, args.logger_dir)
from logger import Logger
logger = Logger(name="padloc_wrapper", log_to_console=True, log_to_file=False).get_logger()

# ensure log exists
os.makedirs(os.path.dirname(args.log), exist_ok=True)
open(args.log, "a").close()

# ensure outdir exists
os.makedirs(args.outdir, exist_ok=True)
logger.info(f"Running Padloc in: {args.outdir}")

# run padloc
cmd = [
    "padloc", "--fix-prodigal",
    "--faa",  args.faa,
    "--gff",  args.gff,
    "--cpu",  str(args.cpu),
    "--outdir", args.outdir
]
if args.force:
    cmd.append("--force")

logger.info("CMD: " + " ".join(cmd))
with open(args.log, "ab") as lf:
    proc = subprocess.run(cmd, stdout=lf, stderr=lf)
if proc.returncode != 0:
    logger.error(f"Padloc failed (exit {proc.returncode}); see {args.log}")
    sys.exit(proc.returncode)
logger.info("✅ Padloc completed successfully.")

# compute prefix from the parent folder name of outdir (e.g. 'results')
prefix = os.path.basename(os.path.dirname(args.outdir))

# pick up the CSV Padloc generated (<anything>_padloc.csv) and rename to {prefix}_padloc.csv
wild = glob.glob(os.path.join(args.outdir, "*_padloc.csv"))
dst_csv = os.path.join(args.outdir, f"{prefix}_padloc.csv")
if wild:
    shutil.move(wild[0], dst_csv)
    logger.info(f"Renamed {os.path.basename(wild[0])} → {os.path.basename(dst_csv)}")
else:
    open(dst_csv, "w").close()
    logger.warning(f"No CSV produced; created empty stub {os.path.basename(dst_csv)}")