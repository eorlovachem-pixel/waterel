#!/usr/bin/env python3
"""Create local Google Drive payload trees. Does not invent QC files."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "drive_staging"
PARENT_ID = "14ltUb63LmmihKGgzMf-geoXxPUmNLSp2"
SUBDIRS = (
    "final_structures", "tables", "figures", "validation", "reports", "archive",
)


def write_tree(ion_folder: str, ion: str) -> Path:
    dest = STAGE / ion_folder
    for name in SUBDIRS:
        (dest / name).mkdir(parents=True, exist_ok=True)
        readme = dest / name / "README.md"
        if not readme.is_file():
            readme.write_text(
                f"# {ion} / {name}\n\n"
                "No frequency-confirmed geometries yet. "
                f"`CASE_STATUS` for {ion}_04–{ion}_14 is `NOT_ASSIGNED`.\n"
                "This placeholder exists so the Drive folder layout is ready.\n",
                encoding="utf-8",
            )
    man = {
        "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
        "google_drive_parent_id": PARENT_ID,
        "folder": ion_folder,
        "ion": ion,
        "qc_uploaded": False,
        "reason": "NOT_ASSIGNED; no GWH/xTB results to upload",
        "local_packings": [
            f"inputs/{ion}_{n:02d}/packing.vecfle" for n in range(4, 15)
        ],
        "subdirs": list(SUBDIRS),
        "stoichiometry_dirs": [f"{ion}_{n:02d}" for n in range(4, 15)],
    }
    (dest / "MANIFEST.json").write_text(
        json.dumps(man, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for n in range(4, 15):
        (dest / f"{ion}_{n:02d}").mkdir(parents=True, exist_ok=True)
        (dest / f"{ion}_{n:02d}" / "README.md").write_text(
            f"# {ion}_{n:02d}\n\nNo production GWH/xTB for this stoichiometry yet.\n",
            encoding="utf-8",
        )
    (dest / "UPLOAD_STATUS.md").write_text(
        f"# {ion_folder}\n\n"
        f"Parent Drive folder: https://drive.google.com/drive/folders/{PARENT_ID}\n\n"
        "OAuth upload is attempted by `tools/drive_upload.py`. "
        "Until credentials exist, only this local staging tree is present.\n",
        encoding="utf-8",
    )
    return dest


def main() -> None:
    STAGE.mkdir(parents=True, exist_ok=True)
    paths = [
        write_tree("K+_eorlovachem-pixel", "K"),
        write_tree("Cl-_eorlovachem-pixel", "Cl"),
        write_tree("Br-_eorlovachem-pixel", "Br"),
    ]
    for p in paths:
        print("staged", p)


if __name__ == "__main__":
    main()
