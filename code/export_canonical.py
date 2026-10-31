#!/usr/bin/env python3
"""Copy unique K6 pilot products into results/<ion>/<n>/ without rewriting originals.

The canonical original tree remains results/K_06/pilot_2/. This layout is a
publication-facing mirror with SHA-256 provenance.
"""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    import pipeline as P
    import qc_checks
    import shell_descriptors as SD
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import pipeline as P
    import qc_checks
    import shell_descriptors as SD

ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "results" / "K_06" / "pilot_2"
HARTREE_TO_KCAL = 627.5094740631

CONFIRMED = [
    {
        "id": "frame_001",
        "origin": "free optimization of GWH structure 0001",
        "xyz": PILOT / "xtb_runs" / "frame_001" / "xtbopt.xyz",
        "opt_log": PILOT / "xtb_runs" / "frame_001" / "xtb.out",
        "hess_log": PILOT / "xtb_runs" / "frame_001" / "hess" / "xtb_hess.out",
        "parent": None,
        "energy_Eh": -30.517135487228,
        "delta_kcal_mol": 0.0,
        "lowest_nu_cm1": 31.12,
        "cocos_group": 0,
    },
    {
        "id": "frame_000_minus",
        "origin": "frame_000 imaginary-mode displacement -0.15 A, reoptimized",
        "xyz": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "opt" / "xtbopt.xyz",
        "opt_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "opt" / "xtb.out",
        "hess_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "hess" / "xtb_hess.out",
        "parent": "frame_000",
        "energy_Eh": -30.517100856916,
        "delta_kcal_mol": 0.021730848870468323,
        "lowest_nu_cm1": 23.62,
        "cocos_group": 0,
    },
    {
        "id": "frame_000_plus",
        "origin": "frame_000 imaginary-mode displacement +0.15 A, reoptimized",
        "xyz": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "opt" / "xtbopt.xyz",
        "opt_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "opt" / "xtb.out",
        "hess_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "hess" / "xtb_hess.out",
        "parent": "frame_000",
        "energy_Eh": -30.517067438374,
        "delta_kcal_mol": 0.042701300584500446,
        "lowest_nu_cm1": 22.12,
        "cocos_group": 0,
    },
]


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def copy_if(src: Path, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_file():
        shutil.copy2(src, dest)
        return True
    return False


def verify_energy(xyz: Path, expected: float):
    fr = P.read_xyz(xyz)[0]
    e = qc_checks.energy_from_comment(fr["comment"])
    if e is None:
        e = P.energy_in_comment(fr["comment"])
    if e is None:
        raise RuntimeError(f"{xyz}: no energy in comment")
    if abs(e - expected) > 1e-9:
        raise RuntimeError(f"{xyz}: energy {e} != {expected}")
    return fr, e


def export_k06():
    if not PILOT.is_dir():
        raise FileNotFoundError(PILOT)
    dest = ROOT / "results" / "K" / "6"
    for sub in ("gwh", "xtb", "hessian", "mode_follow", "cocos", "environment", "final"):
        (dest / sub).mkdir(parents=True, exist_ok=True)

    # Pointers to unmodified originals.
    (dest / "SOURCE_PILOT.txt").write_text(
        "Original unmodified tree: results/K_06/pilot_2/\n"
        "This directory is a publication-facing copy of unique confirmed products.\n"
        "Do not treat it as a second independent calculation.\n",
        encoding="utf-8",
    )
    copy_if(PILOT / "gwh" / "gwh_opt_str0000" / "mol_str0000_fin.xyz", dest / "gwh" / "gwh_0000_fin.xyz")
    copy_if(PILOT / "gwh" / "gwh_opt_str0001" / "mol_str0001_fin.xyz", dest / "gwh" / "gwh_0001_fin.xyz")
    copy_if(PILOT / "source" / "constrained_multiframe.xyz", dest / "gwh" / "constrained_multiframe.xyz")

    rows = []
    frames_final = []
    mdatom_note = dest / "environment" / "MDAtomLocalEnvironment_STATUS.txt"
    mdatom_log = PILOT / "postprocess_confirmed_minima" / "MDAtomLocalEnvironment.log"
    if mdatom_log.is_file():
        shutil.copy2(mdatom_log, dest / "environment" / "MDAtomLocalEnvironment.log")
        mdatom_note.write_text("Copied original MDAtomLocalEnvironment stdout. Python descriptors are separate.\n", encoding="utf-8")
    else:
        mdatom_note.write_text(
            "MDAtomLocalEnvironment was not executed on this pilot (optional --mdatom). "
            "Python first-shell descriptors in final/descriptors.csv are separate from author stdout.\n",
            encoding="utf-8",
        )

    copy_if(
        PILOT / "postprocess_confirmed_minima" / "cocos_primary" / "H12KO6_quick_cleanup.trj",
        dest / "cocos" / "H12KO6_quick_cleanup.trj",
    )
    copy_if(
        PILOT / "postprocess_confirmed_minima" / "POSTPROCESS_REPORT.md",
        dest / "cocos" / "POSTPROCESS_REPORT.md",
    )

    for item in CONFIRMED:
        xyz = Path(item["xyz"])
        fr, e = verify_energy(xyz, item["energy_Eh"])
        qc_checks.assert_cluster_composition(fr, "K", 6, charge=1)
        digest = sha256(xyz)
        copy_if(xyz, dest / "xtb" / f"{item['id']}.xyz")
        copy_if(xyz, dest / "final" / f"{item['id']}.xyz")
        copy_if(Path(item["opt_log"]), dest / "xtb" / f"{item['id']}_xtb.out")
        copy_if(Path(item["hess_log"]), dest / "hessian" / f"{item['id']}_hess.out")
        if item["parent"]:
            copy_if(xyz, dest / "mode_follow" / f"{item['id']}.xyz")
        desc = SD.cluster_descriptors(fr, "K", n_water=6, xyz_path=xyz)
        row = {
            "id": item["id"],
            "ion": "K",
            "charge": 1,
            "n_water": 6,
            "origin": item["origin"],
            "energy_Eh": f"{e:.14f}",
            "delta_kcal_mol": item["delta_kcal_mol"],
            "lowest_nu_cm1": item["lowest_nu_cm1"],
            "n_significant_imaginary": 0,
            "freq_confirmed": True,
            "cocos_group": item["cocos_group"],
            "xyz_sha256": digest,
            "source_xyz": str(xyz.relative_to(ROOT)),
            "visualization": f"figures/burette/scenes/{item['id']}/original.xyz",
            **SD.csv_row(desc),
        }
        rows.append(row)
        frames_final.append({
            "symbols": fr["symbols"],
            "xyz": fr["xyz"],
            "comment": (
                f"id={item['id']} origin={item['origin']} E = {e:.14f} Eh "
                f"cocos_group={item['cocos_group']} CN={desc['CN']} sha256={digest}"
            ),
        })

    P.write_csv(dest / "final" / "descriptors.csv", rows)
    P.write_xyz(dest / "final" / "unique_representatives.xyz", frames_final)
    (dest / "final" / "manifest.json").write_text(json.dumps({
        "ion": "K",
        "n_water": 6,
        "charge": 1,
        "method": "GFN2-xTB / ALPB(water)",
        "n_gwh_starts": 2,
        "n_freq_confirmed": 3,
        "n_cocos_groups": 1,
        "cocos_note": "All three confirmed minima belong to CoCoS group 0; representative is frame_001.",
        "cn_cut_status": SD.CN_CUT_STATUS,
        "cn_cut_A": SD.CN_CUT_PROVISIONAL_A["K"],
        "pilot_original": "results/K_06/pilot_2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "rows": rows,
    }, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    return dest, rows
