#!/usr/bin/env python3
"""Render top PFAM cluster markers as a readable table figure.

Input:
  artifacts/nb04_cluster_unknown_6960/cluster_pfam_markers_c6960_unk.tsv

Outputs:
  figures/04_cluster_unknown_6960/09_cluster_pfam_markers_table.{png,svg}
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams

FF = Path(__file__).resolve().parents[1]
ART = FF / "artifacts/nb04_cluster_unknown_6960"
FIG = FF / "figures/04_cluster_unknown_6960"
PATH_MARKERS = ART / "cluster_pfam_markers_c6960_unk.tsv"
OUT_STEM = FIG / "09_cluster_pfam_markers_table"

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def _clip(text: str, n: int) -> str:
    text = str(text or "").strip()
    return text if len(text) <= n else text[: n - 1] + "…"


def build_table_rows(df: pd.DataFrame) -> tuple[list[list[str]], list[str]]:
    df = df.sort_values(["leiden_cluster", "marker_score"], ascending=[True, False])
    rows: list[list[str]] = []
    for cl, sub in df.groupby("leiden_cluster", sort=True):
        n_ie = int(sub["cluster_size"].iloc[0])
        for rank, r in enumerate(sub.itertuples(index=False), start=1):
            short = r.short_name if pd.notna(r.short_name) and str(r.short_name).strip() else r.pfam_acc
            rows.append([
                f"C{int(cl)}" if rank == 1 else "",
                f"{n_ie:,}" if rank == 1 else "",
                str(rank),
                short,
                _clip(r.full_name, 78),
                f"{100 * r.p_in:.0f}%",
                f"{r.enrichment:.1f}×",
            ])
    cols = ["Cluster", "n IE", "#", "PFAM", "Description", "p_in", "Enrichment"]
    return rows, cols


def render_markers_table() -> Path:
    FIG.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PATH_MARKERS, sep="\t")
    rows, cols = build_table_rows(df)

    n_rows = len(rows)
    fig_h = max(6.5, 0.30 * n_rows + 2.2)
    fig, ax = plt.subplots(figsize=(16.0, fig_h))
    ax.axis("off")

    tbl = ax.table(
        cellText=rows,
        colLabels=cols,
        loc="upper center",
        cellLoc="left",
        colLoc="left",
        bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1.0, 1.40)

    col_widths = [0.07, 0.06, 0.04, 0.12, 0.50, 0.07, 0.10]
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d9d9d9")
        if row == 0:
            cell.set_facecolor("#eef2f7")
            cell.set_text_props(weight="bold", color="#1a1a1a")
            cell.set_height(0.045)
            continue
        data_row = row - 1
        cluster_start = rows[data_row][0] != ""
        if cluster_start:
            cell.set_facecolor("#f7fbff")
        elif data_row % 2 == 0:
            cell.set_facecolor("#fafafa")
        else:
            cell.set_facecolor("#ffffff")
        if col < len(col_widths):
            cell.set_width(col_widths[col])

    ax.set_title(
        "6960 Unknown · top PFAM markers per Leiden cluster (Fisher enrichment)",
        fontsize=11, weight="bold", pad=14, loc="left", x=0,
    )
    fig.text(
        0.0, 0.008,
        "p_in = fraction of IE in cluster carrying PFAM; enrichment = p_in / p_out",
        fontsize=8, color="#555",
    )
    plt.subplots_adjust(left=0.02, right=0.98, top=0.93, bottom=0.05)

    for ext in ("png", "svg"):
        fig.savefig(f"{OUT_STEM}.{ext}", dpi=160, bbox_inches="tight", pad_inches=0.30)
    plt.close(fig)

    # Drop legacy heatmap if present.
    for legacy in (FIG / "09_heatmap_top_markers.png", FIG / "09_heatmap_top_markers.svg"):
        if legacy.exists():
            legacy.unlink()

    print(f"→ {OUT_STEM}.png")
    return OUT_STEM.with_suffix(".png")


if __name__ == "__main__":
    render_markers_table()
