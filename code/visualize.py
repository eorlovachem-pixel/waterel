#!/usr/bin/env python3
"""Prepare Burette scenes and export molecular illustrations.

Burette (https://github.com/SergeiNikolenko/Burette) is a macOS desktop app
plus a browser demo. xyzrender is the documented XYZ fallback renderer used by
Burette. This module:

* never modifies input coordinates of the source XYZ (originals are copied);
* writes Burette-openable scene packages (XYZ + camera JSON);
* renders publication images with xyzrender when the CLI is available;
* never labels a matplotlib or other non-Burette image as a Burette render.

If Burette.app / Mol* export is unavailable on this host, the blocker is
recorded in figures/burette/scenes/BURETTE_STATUS.json.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    import qc_checks
    import pipeline as P
    import shell_descriptors as SD
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import qc_checks
    import pipeline as P
    import shell_descriptors as SD

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
SCENE_DIR = FIG / "burette" / "scenes"
STRUCT_DIR = FIG / "burette" / "structures"
COMP_DIR = FIG / "burette" / "comparisons"
COORD_DIR = FIG / "burette" / "coordination"
PUB_DIR = FIG / "publication"

# CPK-ish legend used in documentation and xyzrender defaults.
ATOM_LEGEND = {
    "Li": "violet",
    "Na": "blue",
    "K": "purple / violet (central ion, enlarged)",
    "Rb": "dark red / violet",
    "F": "green",
    "Cl": "green",
    "Br": "dark red",
    "O": "red",
    "H": "white",
}

CANVAS_PX = 2400
XYZRENDER_CONFIG = "tube"
CAMERA = {
    "align": "ion_at_origin_then_Kabsch_on_O_for_comparison_copy_only",
    "original_xyz_untouched": True,
    "xyzrender_flags": {
        "config": XYZRENDER_CONFIG,
        "canvas_px": CANVAS_PX,
        "hy": True,
        "nci": True,
        "orient": False,
        "background": "white",
        "atom_scale": 1.15,
        "hl_central": True,
    },
}


def _which_xyzrender():
    return shutil.which("xyzrender") or shutil.which("xyzrender.exe")


def _burette_status():
    return {
        "burette_repo": "https://github.com/SergeiNikolenko/Burette",
        "burette_docs": "https://burette-landing.vercel.app/docs",
        "burette_demo": "https://burette-landing.vercel.app/demo",
        "host": sys.platform,
        "burette_app_available": False,
        "burette_app_reason": (
            "Burette.app is a macOS 12+ Apple Silicon desktop application "
            "(Homebrew cask SergeiNikolenko/burette). This host is Windows/WSL. "
            "Finder Quick Look, the Codex plugin, and the iPhone target do not "
            "run in Cursor on Windows. The public browser demo is a preview "
            "surface, not a headless exporter; screenshots would include UI chrome, "
            "which is forbidden for publication figures."
        ),
        "xyzrender_role": (
            "Burette documents xyzrender as the XYZ fallback renderer "
            "(https://github.com/SergeiNikolenko/Burette README, External xyzrender). "
            "CLI exports from xyzrender are labeled renderer=xyzrender, not Burette.app."
        ),
        "xyzrender_path": _which_xyzrender(),
        "checked_utc": datetime.now(timezone.utc).isoformat(),
    }


def copy_xyz_verbatim(src: Path, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    data = Path(src).read_bytes()
    dest.write_bytes(data)
    if dest.read_bytes() != data:
        raise RuntimeError(f"XYZ copy changed bytes: {src} -> {dest}")
    return hashlib.sha256(data).hexdigest()


def align_ion_and_oxygens(frame, ion):
    """Return a viewing copy: ion at origin, Kabsch-align O to a reference octahedron.

    Original coordinates are not written back to the source file.
    """
    symbols = np.array(frame["symbols"])
    xyz = np.asarray(frame["xyz"], float).copy()
    ci = int(np.where(symbols == ion)[0][0])
    xyz -= xyz[ci]
    o_ids = np.where(symbols == "O")[0]
    o = xyz[o_ids]
    # Principal axes of the O cloud: reproducible, no random seed.
    cov = o.T @ o
    w, v = np.linalg.eigh(cov)
    order = np.argsort(w)[::-1]
    rot = v[:, order]
    if np.linalg.det(rot) < 0:
        rot[:, -1] *= -1
    xyz = xyz @ rot
    # Flip so the farthest O has +z, then +x.
    o = xyz[o_ids]
    if o[np.argmax(np.abs(o[:, 2])), 2] < 0:
        xyz[:, 2] *= -1
        xyz[:, 1] *= -1
    o = xyz[o_ids]
    if o[np.argmax(np.abs(o[:, 0])), 0] < 0:
        xyz[:, 0] *= -1
        xyz[:, 1] *= -1
    return {"symbols": list(frame["symbols"]), "xyz": xyz, "comment": frame.get("comment", "")}


def write_scene(entry, dest_dir: Path):
    dest_dir.mkdir(parents=True, exist_ok=True)
    src = Path(entry["xyz_path"])
    original = dest_dir / "original.xyz"
    digest = copy_xyz_verbatim(src, original)
    frames = P.read_xyz(original)
    if len(frames) != 1:
        raise ValueError(f"{src}: expected a single frame, got {len(frames)}")
    aligned = align_ion_and_oxygens(frames[0], entry["ion"])
    P.write_xyz(dest_dir / "aligned_for_view.xyz", [aligned])
    nci_pairs = []
    desc = entry.get("descriptors") or {}
    hb = json.loads(desc.get("hbonds_json") or "[]")
    xh = json.loads(desc.get("anion_xh_json") or "[]")
    for bond in hb:
        nci_pairs.append(f"{bond['h_index']+1}-{bond['acceptor_O']+1}")
    for bond in xh:
        nci_pairs.append(f"{bond['h_index']+1}-1")  # central ion is atom 1 after typical K-first XYZ; overwritten below
    # Use actual central index (1-based) from descriptors.
    central_1b = int(desc.get("central_index", 0)) + 1
    nci_pairs = []
    for bond in hb:
        nci_pairs.append(f"{int(bond['h_index'])+1}-{int(bond['acceptor_O'])+1}")
    for bond in xh:
        nci_pairs.append(f"{int(bond['h_index'])+1}-{central_1b}")
    scene = {
        "id": entry["id"],
        "ion": entry["ion"],
        "n_water": entry.get("n_water"),
        "charge": qc_checks.protocol_charge(entry["ion"]),
        "energy_Eh": entry.get("energy_Eh"),
        "cocos_group": entry.get("cocos_group"),
        "origin": entry.get("origin"),
        "source_xyz": str(src),
        "original_xyz": "original.xyz",
        "aligned_xyz": "aligned_for_view.xyz",
        "original_sha256": digest,
        "method": "GFN2-xTB / ALPB(water), frequencies confirmed ν ≥ −5 cm⁻¹",
        "camera": CAMERA,
        "atom_legend": ATOM_LEGEND,
        "nci_pairs_1based": nci_pairs,
        "burette_open": {
            "macos": "Open original.xyz in Burette.app (Mol* or xyzrender renderer).",
            "windows_wsl": "Burette.app not available; use xyzrender CLI on aligned_for_view.xyz.",
        },
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (dest_dir / "scene.json").write_text(json.dumps(scene, indent=2, ensure_ascii=False), encoding="utf-8")
    return scene, dest_dir / "aligned_for_view.xyz", nci_pairs


def xyzrender_cmd(inp: Path, out: Path, nci_pairs, highlight_ion=True, canvas=CANVAS_PX):
    cmd = [
        "xyzrender", str(inp),
        "-o", str(out),
        "--config", XYZRENDER_CONFIG,
        "-S", str(canvas),
        "--hy",
        "--no-orient",
        "-B", "white",
        "-a", "1.15",
        "--nci",
    ]
    if highlight_ion:
        cmd.extend(["--hl", "M", "#5b2c6f"])
        cmd.extend(["--radius-scale", "M", "1.35"])
    for pair in nci_pairs:
        cmd.extend(["--nci-bond", pair])
    return cmd


def run_xyzrender(inp: Path, out: Path, nci_pairs):
    exe = _which_xyzrender()
    if not exe:
        return {"ok": False, "reason": "xyzrender CLI not on PATH", "path": None}
    cmd = xyzrender_cmd(inp, out, nci_pairs)
    cmd[0] = exe
    out.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    info = {
        "ok": proc.returncode == 0 and out.is_file() and out.stat().st_size > 0,
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": (proc.stdout or "")[-2000:],
        "stderr": (proc.stderr or "")[-2000:],
        "path": str(out) if out.exists() else None,
        "bytes": out.stat().st_size if out.exists() else 0,
        "renderer": "xyzrender",
        "not_burette_app": True,
    }
    (out.parent / (out.stem + "_xyzrender.log")).write_text(
        json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return info


def render_entry(entry, figure_id):
    scene_dir = SCENE_DIR / entry["id"]
    scene, aligned, nci_pairs = write_scene(entry, scene_dir)
    products = []
    for suffix, dest_root in (("overview", STRUCT_DIR), ("shell", COORD_DIR)):
        for ext in (".svg", ".png"):
            dest = dest_root / f"{figure_id}_{entry['id']}_{suffix}{ext}"
            info = run_xyzrender(aligned, dest, nci_pairs)
            products.append({"file": str(dest), "kind": suffix, "ext": ext, **info})
    return {"scene": scene, "products": products}


def comparison_overlay(entries, figure_id):
    """Overlay is only valid for same stoichiometry; xyzrender --ensemble on a multi-XYZ."""
    if len(entries) < 2:
        return None
    COMP_DIR.mkdir(parents=True, exist_ok=True)
    SCENE_DIR.mkdir(parents=True, exist_ok=True)
    frames = []
    for e in entries:
        src = SCENE_DIR / e["id"] / "aligned_for_view.xyz"
        if not src.exists():
            write_scene(e, SCENE_DIR / e["id"])
        frames.extend(P.read_xyz(SCENE_DIR / e["id"] / "aligned_for_view.xyz"))
    multi = SCENE_DIR / f"{figure_id}_ensemble.xyz"
    P.write_xyz(multi, frames)
    products = []
    exe = _which_xyzrender()
    if not exe:
        return {"ok": False, "reason": "xyzrender CLI not on PATH", "xyz": str(multi)}
    for ext in (".svg", ".png"):
        dest = COMP_DIR / f"{figure_id}_ensemble{ext}"
        cmd = [
            exe, str(multi), "-o", str(dest),
            "--config", XYZRENDER_CONFIG, "-S", str(CANVAS_PX),
            "--hy", "--no-orient", "-B", "white", "--ensemble",
            "--hl", "M", "#5b2c6f",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        products.append({
            "file": str(dest),
            "ok": proc.returncode == 0 and dest.is_file(),
            "returncode": proc.returncode,
            "renderer": "xyzrender",
            "cmd": cmd,
        })
    return {"xyz": str(multi), "products": products}


def write_status():
    SCENE_DIR.mkdir(parents=True, exist_ok=True)
    status = _burette_status()
    (SCENE_DIR / "BURETTE_STATUS.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return status


def visualize_confirmed(entries, figure_prefix="Fig_K06_Pilot"):
    """entries: list of dicts with id, ion, xyz_path, energy_Eh, cocos_group, origin, n_water, descriptors?."""
    status = write_status()
    FIG.mkdir(parents=True, exist_ok=True)
    for d in (SCENE_DIR, STRUCT_DIR, COMP_DIR, COORD_DIR, PUB_DIR):
        d.mkdir(parents=True, exist_ok=True)
    reports = []
    for e in entries:
        if "descriptors" not in e:
            fr = P.read_xyz(e["xyz_path"])[0]
            e["descriptors"] = SD.cluster_descriptors(
                fr, e["ion"], n_water=e.get("n_water"), xyz_path=e["xyz_path"]
            )
        reports.append(render_entry(e, figure_prefix))
    overlay = comparison_overlay(entries, figure_prefix)
    index_rows = []
    n_ok = 0
    for e, rep in zip(entries, reports):
        for prod in rep["products"]:
            if prod.get("ok"):
                n_ok += 1
            index_rows.append({
                "figure_id": f"{figure_prefix}_{e['id']}_{prod.get('kind')}",
                "geometry_id": e["id"],
                "ion": e["ion"],
                "xyz_source": e["xyz_path"],
                "xyz_scene": str(SCENE_DIR / e["id"] / "original.xyz"),
                "method": "GFN2-xTB/ALPB(water), Hessian-confirmed",
                "energy_Eh": e.get("energy_Eh"),
                "cocos_group": e.get("cocos_group"),
                "CN": (e.get("descriptors") or {}).get("CN"),
                "renderer": "xyzrender" if prod.get("ok") else "NONE",
                "burette_app": False,
                "file": prod.get("path") or prod.get("file"),
                "ok": prod.get("ok"),
                "export": f"canvas={CANVAS_PX}px config={XYZRENDER_CONFIG} PNG/SVG",
            })
    summary = {
        "burette": status,
        "n_geometries": len(entries),
        "n_successful_exports": n_ok,
        "overlay": overlay,
        "index": index_rows,
    }
    (SCENE_DIR / f"{figure_prefix}_render_report.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    return summary
