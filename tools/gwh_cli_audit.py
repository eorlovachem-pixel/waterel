#!/usr/bin/env python3
"""Audit Denis's GWH CLI for anion Charge -1. Does not run xTB or --Opt."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import random
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GWH = ROOT / "vendors" / "Glorious-Water-Hole" / "gloriouswaterhole.py"
OUT = ROOT / "validation" / "gwh_cli" / "GWH_CLI_AUDIT.json"
XYZ_DIR = ROOT / "validation" / "gwh_cli" / "xyz"
SEED = 20261101
N_STRUCT = 20


def load_gwh():
    spec = importlib.util.spec_from_file_location("gloriouswaterhole", GWH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def generate_random_xyz(gwh, packing: Path, nstruct: int, dest: Path, seed: int) -> dict:
    """Write nstruct XYZ from GWH geometry code. No xTB, no energies."""
    dest.mkdir(parents=True, exist_ok=True)
    centers, vectors = gwh.ReadVectorFile(str(packing))
    water_sym, water_xyz = gwh.GenerateOriginalWater(0.958, 104.4776)
    water_xyz = gwh.ReCenterMolecule(water_xyz, None)
    rng = random.Random(seed)
    written = []
    for i in range(nstruct):
        euler = {}
        for cent in centers:
            n = len(vectors[cent])
            ea = np.zeros((n, 3), float)
            for j in range(n):
                ea[j, 0] = rng.uniform(0.0, 360.0)
                ea[j, 1] = rng.uniform(0.0, 180.0)
                ea[j, 2] = rng.uniform(0.0, 360.0)
            euler[cent] = ea
        symbols, xyz, cos_theta = gwh.GenerateGeometry(
            water_sym, water_xyz, centers, vectors, euler, 0.0, 2
        )
        path = dest / f"str_{i:04d}.xyz"
        comment = (
            f"GWH geometry-only nstruct={i} seed={seed} "
            f"packing={packing.name} no_xtb cosθ={cos_theta}"
        )
        gwh.WriteXYZFile(symbols, xyz, FileName=str(path), Mode="w", Description=comment)
        written.append({
            "path": str(path.relative_to(ROOT).as_posix()),
            "natoms": len(symbols),
            "formula": {s: symbols.count(s) for s in sorted(set(symbols))},
            "sha256": sha256(path),
        })
    radii = [float(np.linalg.norm(v)) for v in vectors[centers[0]]]
    return {
        "centers": centers,
        "n_sites": len(vectors[centers[0]]),
        "radii_A": radii,
        "nstruct": nstruct,
        "seed": seed,
        "files": written,
    }


def main() -> None:
    help_p = subprocess.run(
        [sys.executable, str(GWH), "--help"],
        capture_output=True, text=True, check=False,
    )
    help_text = (help_p.stdout or "") + (help_p.stderr or "")
    git = subprocess.run(
        ["git", "-C", str(GWH.parent), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    )
    src = GWH.read_text(encoding="utf-8", errors="replace")
    gwh = load_gwh()
    cl_geom = generate_random_xyz(
        gwh, ROOT / "inputs" / "Cl_04" / "packing.vecfle", N_STRUCT,
        XYZ_DIR / "Cl_04_geom20", SEED,
    )
    br_geom = generate_random_xyz(
        gwh, ROOT / "inputs" / "Br_04" / "packing.vecfle", N_STRUCT,
        XYZ_DIR / "Br_04_geom20", SEED,
    )
    documented = (
        'gloriouswaterhole.py -f <packing.vecfle> --Charge -1 --NumStruct 20 '
        '--Random --Opt --AddMethod " --alpb water --gfn 2 " --DoNotCenter'
    )
    audit = {
        "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
        "gwh_path": str(GWH),
        "gwh_sha256": sha256(GWH),
        "gwh_git": git.stdout.strip() if git.returncode == 0 else None,
        "help_exit_code": help_p.returncode,
        "help_has_Charge": "--Charge CHARGE" in help_text or "--Charge" in help_text,
        "help_has_NumStruct": "--NumStruct" in help_text,
        "help_has_Random": "--Random" in help_text,
        "help_has_Opt": "--Opt" in help_text,
        "help_has_AddMethod": "--AddMethod" in help_text,
        "help_has_DoNotCenter": "--DoNotCenter" in help_text,
        "default_charge_is_0": "Default is 0" in help_text,
        "default_mult_is_1": "Default is 1" in help_text,
        "uhf0_means_mult_1": True,
        "from_euler_zyz": "from_euler('zyz'" in src,
        "from_euler_ZYZ": "from_euler('ZYZ'" in src,
        "xtb_uses_charge_flag": "--charge" in src,
        "documented_20_cmd": documented,
        "Cl_04_geom20": {k: cl_geom[k] for k in ("centers", "n_sites", "radii_A", "nstruct", "seed")},
        "Br_04_geom20": {k: br_geom[k] for k in ("centers", "n_sites", "radii_A", "nstruct", "seed")},
        "Cl_04_nstruct_ok": cl_geom["nstruct"] == N_STRUCT and all(f["natoms"] == 13 for f in cl_geom["files"]),
        "Br_04_nstruct_ok": br_geom["nstruct"] == N_STRUCT and all(f["natoms"] == 13 for f in br_geom["files"]),
        "energies": None,
        "opt_series_launched": False,
        "hundred_start_launched": False,
        "reason": (
            "CASE_STATUS Cl/Br rows are NOT_ASSIGNED. Author has not finalized anion inclusion. "
            "This audit checks CLI, Euler 'zyz', packing readability, and a 20-structure "
            "geometry-only loop (no xTB SPE, no --Opt). GWH --Opt and the 100-start series "
            "are not started."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    seeds = {
        "Cl_04": [f["sha256"] for f in cl_geom["files"]],
        "Br_04": [f["sha256"] for f in br_geom["files"]],
    }
    (OUT.parent / "GEOM20_SHA256.json").write_text(
        json.dumps(seeds, indent=2) + "\n", encoding="utf-8"
    )
    failed = [k for k, v in audit.items() if k.startswith("help_has_") and v is not True]
    if audit["help_exit_code"] != 0 or not audit["from_euler_zyz"] or audit["from_euler_ZYZ"]:
        failed.append("euler_or_help")
    if not audit["Cl_04_nstruct_ok"] or not audit["Br_04_nstruct_ok"]:
        failed.append("geom20")
    print(json.dumps({"wrote": str(OUT), "failed": failed, "gwh_git": audit["gwh_git"]}, indent=2))
    if failed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
