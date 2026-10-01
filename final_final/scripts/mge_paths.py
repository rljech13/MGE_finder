"""Locations for the module experiments.

The repository root is the directory that contains ``finder_pipeline/Snakefile``,
found by walking up from this file. Scripts keep working when they are started
through a link in ``final_final/02_atlas``, ``03_clustering``, ``04_cargo`` or ``05_abundance``.

External databases and HMMER binaries can be overridden:

  MGE_HMMSEARCH, MGE_HMMFETCH, MGE_HMMPRESS, MGE_HMMSTAT
  MGE_PFAM, MGE_PADLOC_HMM, MGE_DEFENSE_FINDER, MGE_CASFINDER
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path


def repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "finder_pipeline" / "Snakefile").is_file():
            return parent
    raise SystemExit(f"repository root not found above {here}")


def tool(name: str, *fallbacks: str) -> Path:
    env = os.environ.get("MGE_" + name.upper())
    if env:
        return Path(env)
    found = shutil.which(name)
    if found:
        return Path(found)
    for fallback in fallbacks:
        if Path(fallback).is_file():
            return Path(fallback)
    return Path(name)


def data_file(env: str, *candidates: str) -> Path:
    chosen = os.environ.get(env)
    if chosen:
        return Path(chosen)
    for candidate in candidates:
        if Path(candidate).exists():
            return Path(candidate)
    if not candidates:
        raise SystemExit(f"{env} is not set and no default path was given")
    return Path(candidates[0])


REPO = repo_root()
FF = REPO / "final_final"

HMMSEARCH = tool("hmmsearch", "/home/lam34/miniforge3/envs/MGE_finder/bin/hmmsearch")
HMMFETCH = tool("hmmfetch", "/home/lam34/miniforge3/envs/MGE_finder/bin/hmmfetch")
HMMPRESS = tool("hmmpress", "/home/lam34/miniforge3/envs/MGE_finder/bin/hmmpress")
HMMSTAT = tool("hmmstat", "/home/lam34/miniforge3/envs/MGE_finder/bin/hmmstat")
PFAM = data_file("MGE_PFAM", "/home/lam34/bakta_db/db/pfam.h3m")
PADLOC = data_file("MGE_PADLOC_HMM", "/home/lam34/miniforge3/envs/padloc/data/hmm/padlocdb.hmm")
DF = data_file("MGE_DEFENSE_FINDER", "/home/lam34/.macsyfinder/models/defense-finder-models/profiles")
CAS = data_file("MGE_CASFINDER", "/home/lam34/.macsyfinder/models/CasFinder/profiles")
