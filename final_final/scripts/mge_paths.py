"""Locations for the module experiments.

The repository root is the directory that contains ``finder_pipeline/Snakefile``,
found by walking up from this file. Scripts keep working when they are started
through a link in ``final_final/02_atlas``, ``03_clustering``, ``04_cargo`` or ``05_abundance``.

HMMER binaries are taken from ``PATH``, or from MGE_HMMSEARCH, MGE_HMMFETCH,
MGE_HMMPRESS, MGE_HMMSTAT. Databases are only read from the environment:

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


def data_file(env: str) -> Path | None:
    """Database path from the environment. Nothing is assumed about this machine."""
    chosen = os.environ.get(env)
    return Path(chosen) if chosen else None


REPO = repo_root()
FF = REPO / "final_final"

HMMSEARCH = tool("hmmsearch")
HMMFETCH = tool("hmmfetch")
HMMPRESS = tool("hmmpress")
HMMSTAT = tool("hmmstat")
PFAM = data_file("MGE_PFAM")
PADLOC = data_file("MGE_PADLOC_HMM")
DF = data_file("MGE_DEFENSE_FINDER")
CAS = data_file("MGE_CASFINDER")
