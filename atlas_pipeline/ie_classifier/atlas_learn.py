
from pathlib import Path

import pandas as pd


def learn_atlas_markers_stub(out_path: Path) -> None:
    """Write placeholder TSV for future enrichment from high-confidence seeds."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        [
            {
                "marker_key": "stub",
                "family_id": "",
                "supporting_elements": 0,
                "note": "replace with atlas-derived associations when seeds are available",
            }
        ]
    )
    df.to_csv(out_path, sep="\t", index=False)
