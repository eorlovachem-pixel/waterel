#!/usr/bin/env python3
"""Write FIGURE_INDEX.md, RESULTS_INDEX.md and related publication docs from on-disk files."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
PUB = FIG / "publication"


def rel(p: Path) -> str:
    try:
        return p.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(p)


def existing(paths):
    return [p for p in paths if p.is_file()]


def write_figure_index():
    FIG.mkdir(parents=True, exist_ok=True)
    burette_status = {}
    st = FIG / "burette" / "scenes" / "BURETTE_STATUS.json"
    if st.is_file():
        burette_status = json.loads(st.read_text(encoding="utf-8"))
    k6_csv = ROOT / "results" / "K" / "6" / "final" / "descriptors.csv"
    man = ROOT / "results" / "K" / "6" / "final" / "manifest.json"
    manifest = json.loads(man.read_text(encoding="utf-8")) if man.is_file() else {}
    rows = manifest.get("rows") or []
    xyzrender_png = list((FIG / "burette").rglob("*.png"))
    xyzrender_svg = list((FIG / "burette").rglob("*.svg"))
    mpl = sorted(PUB.glob("Fig_*")) if PUB.is_dir() else []
    lines = [
        "# Figure index",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Renderer policy",
        "",
        "- **Burette.app** (macOS Mol* workspace) is **not available** on this Windows/WSL host.",
        f"- Reason: {burette_status.get('burette_app_reason', 'see figures/burette/scenes/BURETTE_STATUS.json')}",
        "- **xyzrender** is Burette's documented XYZ fallback renderer. Files under `figures/burette/` produced by the xyzrender CLI are labeled `renderer=xyzrender`, not Burette.app.",
        "- **Matplotlib** ball-and-stick and numerical plots live under `figures/publication/` and are labeled as Matplotlib.",
        "- No dummy Li/Na/Rb/F/Cl/Br structures were drawn.",
        "",
        "## Atom colour legend (ball-and-stick)",
        "",
        "| Element | Colour |",
        "|---|---|",
        "| K (central) | purple `#8e44ad` / xyzrender highlight `#5b2c6f` |",
        "| O | red `#e74c3c` |",
        "| H | white |",
        "| covalent bonds | dark grey |",
        "| O–H···O (if drawn) | dashed blue |",
        "| anion X···H (if drawn) | dashed teal |",
        "",
        f"xyzrender exports found: {len(xyzrender_png)} PNG, {len(xyzrender_svg)} SVG.",
        f"xyzrender on PATH at generation: `{burette_status.get('xyzrender_path')}`",
        "",
        "## Geometry-backed figures",
        "",
        "| Figure ID | Geometry ID | Source XYZ | Method | E / Eh | CoCoS | CN | Renderer | File |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    geom_files = {
        "Fig_K06_Pilot": PUB / "Fig_K06_Pilot.png",
        "Fig_Conformers": PUB / "Fig_Conformers.png",
        "Fig_Coordination_Shells": PUB / "Fig_Coordination_Shells.png",
        "Fig_Structures_Cations": PUB / "Fig_Structures_Cations.png",
        "Fig_K06_DeltaE": PUB / "Fig_K06_DeltaE.png",
        "Fig_K06_R_ion_O": PUB / "Fig_K06_R_ion_O.png",
        "Fig_K06_orientation": PUB / "Fig_K06_orientation.png",
    }
    for fig_id, path in geom_files.items():
        status = rel(path) if path.is_file() else "NOT GENERATED"
        renderer = "matplotlib"
        if fig_id.startswith("Fig_K06_") and "DeltaE" not in fig_id and "R_ion" not in fig_id and "orientation" not in fig_id:
            renderer = "matplotlib (structure); xyzrender if present under figures/burette/"
        e = rows[0]["energy_Eh"] if rows else ""
        gid = "frame_001 + frame_000_minus + frame_000_plus" if "Pilot" in fig_id or "Conformers" in fig_id else (rows[0]["id"] if rows else "")
        xyz = "results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz"
        lines.append(
            f"| {fig_id} | {gid} | {xyz} | GFN2-xTB/ALPB(water), Hessian-confirmed | {e} | 0 | 6 | {renderer} | {status} |"
        )
    lines += [
        "",
        "## Not generated (no calculated structures)",
        "",
        "- `Fig_Structures_Anions` — no F−/Cl−/Br− XYZ in this kit.",
        "- `R(n)`, `CN(n)`, `<cosθ>(n)` series — only n = 6 exists.",
        "",
        "## Per-geometry scenes",
        "",
    ]
    for r in rows:
        scene = FIG / "burette" / "scenes" / r["id"] / "original.xyz"
        lines.append(
            f"- `{r['id']}` SHA-256 `{r.get('xyz_sha256','')}` E = {r.get('energy_Eh')} Eh "
            f"CoCoS {r.get('cocos_group')} CN {r.get('CN')} "
            f"scene `{rel(scene) if scene.is_file() else 'missing'}`"
        )
    (FIG / "FIGURE_INDEX.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return FIG / "FIGURE_INDEX.md"


def write_results_index():
    dest = ROOT / "RESULTS_INDEX.md"
    man_p = ROOT / "results" / "K" / "6" / "final" / "manifest.json"
    man = json.loads(man_p.read_text(encoding="utf-8")) if man_p.is_file() else {}
    text = f"""# Results index

Generated: {datetime.now(timezone.utc).isoformat()}

## What was actually calculated

| System | n | Starts | Hessian-confirmed minima | CoCoS groups | Status |
|---|---|---|---|---|---|
| K+ (H2O)6 pilot | 6 | 2 | 3 | 1 | complete pilot; production 100-start ON_HOLD |
| Li+, Na+, Rb+ production | — | — | — | — | not run (NOT_ASSIGNED / claimed elsewhere) |
| F−, Cl−, Br− | — | — | — | — | code supports charge −1; no packing/assignment; no QC |

Known K6 pilot regression: 2 GWH starts, 3 frequency-confirmed geometries, 1 CoCoS group.
Original files remain in `results/K_06/pilot_2/` (not rewritten).
Canonical mirror: `results/K/6/`.

## Final unique representatives

- Multi-XYZ: `results/K/6/final/unique_representatives.xyz`
- Descriptors: `results/K/6/final/descriptors.csv`
- Manifest: `results/K/6/final/manifest.json`

CN cutoff for K is 3.20 Å, status **BLOCKING_UNCONFIRMED** until Denis agrees the selection rule.
All three confirmed K6 geometries have CN = 6 at that cutoff (observed K–O ≤ 2.706 Å).

## Figures

See `figures/FIGURE_INDEX.md`.

n_freq_confirmed = {man.get('n_freq_confirmed')}
n_cocos_groups = {man.get('n_cocos_groups')}
"""
    dest.write_text(text, encoding="utf-8")
    return dest
