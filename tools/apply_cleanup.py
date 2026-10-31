#!/usr/bin/env python3
"""Archive unused session helpers and the pipeline-check copy. Never deletes unique QC data.

Canonical calculations stay in results/K_06/pilot_2/, inputs/, pilots/, cases/.
CSV CASE_STATUS.csv remains the assignment source of truth.
"""
from __future__ import annotations
import csv, json, shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKUP = ROOT / "archive" / "cleanup_backup_2026-10-31"
HELPERS = [
    "_verify_originals.py",
    "_redraw_table.py",
    "_copy_postprocess.py",
    "_verify_pilot.py",
    "_run_verify_and_copy.py",
    "_copy_pilot_report.py",
    "_pilot_status.py",
    "_test_mode_follow.py",
    "_install_pilot_files.py",
    "_peek_g98.py",
    "_pilot_precheck.py",
    "_run_hold_tests.sh",
    "_run_pipeline_tests.sh",
    "_inventory_cleanup.py",
    "_pilot_inspect.sh",
]
KEEP_FOREVER = (
    ROOT / "results" / "K_06" / "pilot_2",
    ROOT / "inputs",
)


def csv_rows():
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_case_status_json():
    rows = csv_rows()
    dest = ROOT / "CASE_STATUS.json"
    dest.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return dest, rows


def move_to_backup(src: Path):
    if not src.exists():
        return None
    rel = src.relative_to(ROOT)
    dest = BACKUP / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return dest
    shutil.move(str(src), str(dest))
    return dest


def write_pipeline_check_stub():
    stub_dir = ROOT / "results" / "K_06" / "pilot_2_pipeline_check"
    stub_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "archived": True,
        "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
        "note": (
            "Full pipeline-check tree (duplicate of pilot_2 xtb_runs + CoCoS copies) "
            "was moved to archive/cleanup_backup_2026-10-31. Canonical geometries remain "
            "in results/K_06/pilot_2/."
        ),
        "canonical_pilot": "results/K_06/pilot_2",
        "confirmed_ids": ["frame_001", "frame_000_minus", "frame_000_plus"],
        "unique_representative": "frame_001",
        "unique_energy_Eh": -30.517135487228,
        "cocos_stable_grouping": True,
        "n_unique": 1,
        "backup": "archive/cleanup_backup_2026-10-31/results/K_06/pilot_2_pipeline_check",
    }
    (stub_dir / "PIPELINE_CHECK_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main():
    BACKUP.mkdir(parents=True, exist_ok=True)
    moved = []
    old_json = ROOT / "CASE_STATUS.json"
    if old_json.is_file() and "assignment_status" not in old_json.read_text(encoding="utf-8"):
        dest = BACKUP / "CASE_STATUS.json.pre_assignment_column"
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(old_json, dest)
            moved.append("CASE_STATUS.json (stale copy to backup)")
    for name in HELPERS:
        dest = move_to_backup(ROOT / name)
        if dest:
            moved.append(str(Path(name)))
    check = ROOT / "results" / "K_06" / "pilot_2_pipeline_check"
    if check.is_dir() and not (check / "PIPELINE_CHECK_SUMMARY.json").is_file():
        dest = move_to_backup(check)
        if dest:
            moved.append("results/K_06/pilot_2_pipeline_check")
    summary = write_pipeline_check_stub()
    json_path, rows = write_case_status_json()
    assert all(r["assignment_status"] in ("ON_HOLD", "NOT_ASSIGNED") for r in rows)
    assert next(r for r in rows if r["case"] == "K_06")["assignment_status"] == "ON_HOLD"
    assert not any(r["ion"] in ("Na", "Rb") for r in rows)
    for keep in KEEP_FOREVER:
        if not keep.exists():
            print("warning missing canonical", keep)
    report = {
        "backup": str(BACKUP),
        "moved": moved,
        "case_status_json": str(json_path),
        "pipeline_check_stub": summary,
        "untouched": [
            "results/K_06/pilot_2",
            "inputs/",
            "CASE_STATUS.csv assignment_status values",
        ],
    }
    (BACKUP / "CLEANUP_MANIFEST.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
