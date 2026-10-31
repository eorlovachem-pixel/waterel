#!/usr/bin/env python3
"""Post-process the three K+(H2O)6 pilot geometries that have no significant imaginary modes.

Does not invent energies or frequencies. Re-reads completed xTB files, runs the
author CoCoS snapshot with seeded_cocos.py, and writes XYZ / CSV / SVG / PNG.
Does not overwrite the original GWH, free-opt, or Hessian directories.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
import freq_parser  # noqa: E402
import mode_follow  # noqa: E402
import pipeline as P  # noqa: E402
import qc_checks  # noqa: E402

PILOT = ROOT / "results" / "K_06" / "pilot_2"
OUT = PILOT / "postprocess_confirmed_minima"
ION = "K"
N_WATER = 6
CHARGE = 1
IMAG_THRESH = 5.0
HARTREE_TO_KCAL = 627.5094740631
FIXED_SEED = 20261030
FIXED_SAMPLES = 50
ROBUST_SEEDS = (20261030, 101, 202, 303, 7)
ROBUST_SAMPLES = (50, 100)
COCOS = ROOT / "code" / "cocos_uploaded_snapshot.py"
WRAPPER = ROOT / "code" / "seeded_cocos.py"

CANDIDATES = [
    {
        "id": "frame_001",
        "label": "frame_001",
        "origin": "free optimization of GWH structure 0001",
        "xyz": PILOT / "xtb_runs" / "frame_001" / "xtbopt.xyz",
        "opt_log": PILOT / "xtb_runs" / "frame_001" / "xtb.out",
        "hess_log": PILOT / "xtb_runs" / "frame_001" / "hess" / "xtb_hess.out",
        "vibspectrum": PILOT / "xtb_runs" / "frame_001" / "hess" / "vibspectrum",
        "g98": PILOT / "xtb_runs" / "frame_001" / "hess" / "g98.out",
    },
    {
        "id": "frame_000_plus",
        "label": "frame_000/mode_follow/plus",
        "origin": "frame_000 imaginary-mode displacement +0.15 A, reoptimized",
        "xyz": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "opt" / "xtbopt.xyz",
        "opt_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "opt" / "xtb.out",
        "hess_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "hess" / "xtb_hess.out",
        "vibspectrum": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "hess" / "vibspectrum",
        "g98": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "plus_0.15A" / "hess" / "g98.out",
    },
    {
        "id": "frame_000_minus",
        "label": "frame_000/mode_follow/minus",
        "origin": "frame_000 imaginary-mode displacement -0.15 A, reoptimized",
        "xyz": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "opt" / "xtbopt.xyz",
        "opt_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "opt" / "xtb.out",
        "hess_log": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "hess" / "xtb_hess.out",
        "vibspectrum": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "hess" / "vibspectrum",
        "g98": PILOT / "xtb_runs" / "frame_000" / "mode_follow" / "minus_0.15A" / "hess" / "g98.out",
    },
]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_one(spec):
    for key in ("xyz", "opt_log", "hess_log", "vibspectrum", "g98"):
        p = spec[key]
        if not Path(p).is_file():
            raise RuntimeError(f"missing {key}: {p}")
    frames = qc_checks.read_xyz(spec["xyz"])
    if len(frames) != 1:
        raise RuntimeError(f"{spec['xyz']}: expected 1 XYZ frame, got {len(frames)}")
    fr = frames[0]
    qc_checks.assert_cluster_composition(fr, ION, N_WATER, CHARGE)
    log = Path(spec["opt_log"]).read_text(errors="replace")
    energy, opt_info = qc_checks.require_completed_opt(log, fr, fr["symbols"])
    hess = freq_parser.require_completed_hessian(spec["hess_log"], imag_thresh_cm1=IMAG_THRESH)
    vib = mode_follow.parse_vibspectrum(spec["vibspectrum"])
    vib_cls = mode_follow.classify_modes(vib, IMAG_THRESH)
    g98_modes = mode_follow.parse_g98_modes(spec["g98"])
    n_expected = 3 * len(fr["symbols"])
    if vib_cls["n_modes"] != n_expected:
        raise RuntimeError(f"{spec['id']}: vibspectrum has {vib_cls['n_modes']} modes, expected {n_expected}")
    if hess["n_frequencies_parsed"] != n_expected:
        raise RuntimeError(f"{spec['id']}: Hessian log parsed {hess['n_frequencies_parsed']} frequencies, expected {n_expected}")
    if vib_cls["n_significant_imaginary"] != 0 or hess["imaginary_cm1"]:
        raise RuntimeError(f"{spec['id']} still has significant imaginary modes; not a confirmed minimum")
    geom = P.geometry_features(fr, ION)
    return {
        "id": spec["id"],
        "label": spec["label"],
        "origin": spec["origin"],
        "frame": fr,
        "energy_Eh": float(energy),
        "opt_converged": bool(opt_info["opt_converged"]),
        "normal_termination": bool(opt_info["normal_termination"]),
        "n_atoms": len(fr["symbols"]),
        "formula": qc_checks.formula_of(fr["symbols"]),
        "xyz_sha256": sha256(spec["xyz"]),
        "xyz_path": str(spec["xyz"]),
        "opt_log": str(spec["opt_log"]),
        "hess_log": str(spec["hess_log"]),
        "n_frequencies": hess["n_frequencies_parsed"],
        "n_significant_imaginary": vib_cls["n_significant_imaginary"],
        "n_near_zero": vib_cls["n_near_zero"],
        "lowest_vibrational_cm1": min(
            (f for f in vib_cls["frequencies_cm1"] if f > IMAG_THRESH), default=None
        ),
        "near_zero_cm1": vib_cls["near_zero_cm1"],
        "g98_n_modes": len(g98_modes),
        "geometry": geom,
        "xtbopt_comment": fr["comment"],
    }


def run_cocos(xyz_path, workdir, seed, samples):
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    local_xyz = workdir / "free_sorted_for_cocos.xyz"
    shutil.copy2(xyz_path, local_xyz)
    cmd = [
        sys.executable,
        str(WRAPPER),
        str(seed),
        str(COCOS),
        "-f",
        "free_sorted_for_cocos.xyz",
        "--DoNotUseGlobalOpt",
        "--NumSamples",
        str(samples),
        "--DoNotSaveDoubles",
    ]
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    log = proc.stdout + "\n" + proc.stderr
    (workdir / "cocos_log.txt").write_text(log, encoding="utf-8")
    (workdir / "cocos_command.json").write_text(
        json.dumps({"cmd": cmd, "returncode": proc.returncode, "seed": seed, "samples": samples}, indent=2),
        encoding="utf-8",
    )
    if proc.returncode != 0:
        raise RuntimeError(f"CoCoS failed seed={seed} samples={samples}: {log[-2000:]}")
    added = [int(x) for x in re.findall(r"added structure (\d+)\b", proc.stdout)]
    removed = re.findall(r"Removed structure (\d+) as a copy of (\d+)", proc.stdout)
    trjs = list(workdir.glob("*_quick_cleanup.trj"))
    if len(trjs) != 1:
        raise RuntimeError(f"{workdir}: expected one *_quick_cleanup.trj, found {trjs}")
    kept = P.read_xyz(trjs[0])
    return {
        "seed": seed,
        "samples": samples,
        "returncode": proc.returncode,
        "added_indices": added,
        "n_unique": len(kept),
        "n_added_log": len(added),
        "removed": [{"structure": int(a), "copy_of": int(b)} for a, b in removed],
        "cleanup_trj": str(trjs[0]),
        "log": str(workdir / "cocos_log.txt"),
        "cmd": cmd,
    }


def write_multiframe_xyz(path, records, emin):
    frames = []
    for r in records:
        delta = (r["energy_Eh"] - emin) * HARTREE_TO_KCAL
        comment = (
            f"id={r['id']} origin={r['label']} "
            f"E = {r['energy_Eh']:.14f} Eh "
            f"delta_E = {delta:.8f} kcal/mol "
            f"n_imag = {r['n_significant_imaginary']} "
            f"method=GFN2-xTB/ALPB(water) chrg=+1"
        )
        frames.append({"symbols": r["frame"]["symbols"], "xyz": r["frame"]["xyz"], "comment": comment})
    P.write_xyz(path, frames)


def align_for_view(symbols, xyz):
    pts = np.array(xyz, float)
    k = np.where(np.array(symbols) == "K")[0][0]
    pts = pts - pts[k]
    o = pts[np.array(symbols) == "O"]
    cov = np.cov(o.T)
    evals, evecs = np.linalg.eigh(cov)
    axes = evecs[:, np.argsort(evals)[::-1]]
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1
    return pts @ axes


def draw_cluster(ax, symbols, xyz, title, subtitle):
    pts = align_for_view(symbols, xyz)
    order = np.argsort(pts[:, 1])  # draw back-to-front
    display = {"K": ("#5B3A91", 220), "O": ("#C44E3A", 95), "H": ("#E6E1D8", 28)}
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            d = float(np.linalg.norm(pts[i] - pts[j]))
            si, sj = symbols[i], symbols[j]
            if {si, sj} == {"O", "H"} and d < 1.25:
                ax.plot([pts[i, 0], pts[j, 0]], [pts[i, 2], pts[j, 2]], color="#444444", lw=1.2, zorder=1)
            elif {si, sj} == {"K", "O"} and d < 3.6:
                ax.plot([pts[i, 0], pts[j, 0]], [pts[i, 2], pts[j, 2]], color="#8A8A8A", lw=0.8, ls="--", zorder=0)
    for idx in order:
        s = symbols[idx]
        color, size = display[s]
        ax.scatter(pts[idx, 0], pts[idx, 2], s=size, c=color, edgecolors="#222222", linewidths=0.4, zorder=2 + idx * 0.01)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10, pad=4)
    ax.text(0.5, -0.08, subtitle, transform=ax.transAxes, ha="center", va="top", fontsize=8)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def figure_geometries(records, emin, png, svg):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.8))
    for ax, r in zip(axes, records):
        delta = (r["energy_Eh"] - emin) * HARTREE_TO_KCAL
        subtitle = f"E = {r['energy_Eh']:.8f} Eh\nΔE = {delta:.4f} kcal/mol · n_imag = {r['n_significant_imaginary']}"
        draw_cluster(ax, r["frame"]["symbols"], r["frame"]["xyz"], r["label"], subtitle)
    fig.suptitle("[K(H2O)6]+ confirmed minima, GFN2-xTB/ALPB(water), charge +1", fontsize=11, y=0.98)
    fig.tight_layout(rect=(0, 0.04, 1, 0.92))
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)


def figure_table(rows, png, svg):
    fig, ax = plt.subplots(figsize=(11.2, 2.6))
    ax.axis("off")
    col_labels = [
        "id",
        "origin",
        "E / Eh",
        "ΔE / kcal mol⁻¹",
        "n_imag (|ν|>5 cm⁻¹)",
        "n_near-zero",
        "lowest vib / cm⁻¹",
        "CoCoS group",
    ]
    cell = [
        [
            r["id"],
            r["label"],
            f"{r['energy_Eh']:.12f}",
            f"{r['delta_kcal_mol']:.6f}",
            str(r["n_significant_imaginary"]),
            str(r["n_near_zero"]),
            f"{r['lowest_vibrational_cm1']:.2f}",
            str(r["cocos_group"]),
        ]
        for r in rows
    ]
    table = ax.table(cellText=cell, colLabels=col_labels, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1, 1.6)
    for (row, col), cell_obj in table.get_celld().items():
        cell_obj.set_edgecolor("#CCCCCC")
        if row == 0:
            cell_obj.set_facecolor("#EEEAE3")
            cell_obj.set_text_props(weight="medium")
        elif rows[row - 1]["delta_kcal_mol"] == 0.0:
            cell_obj.set_facecolor("#F4F1EA")
    ax.set_title("[K(H2O)6]+ pilot confirmed minima · GFN2-xTB/ALPB(water), +1", fontsize=10, pad=12)
    fig.tight_layout()
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def write_csv(path, rows):
    fields = [
        "id",
        "label",
        "origin",
        "energy_Eh",
        "delta_kcal_mol",
        "n_significant_imaginary",
        "n_near_zero",
        "lowest_vibrational_cm1",
        "n_frequencies",
        "formula",
        "opt_converged",
        "normal_termination",
        "distance_min_A",
        "distance_mean_A",
        "distance_max_A",
        "cocos_group",
        "xyz_path",
        "xyz_sha256",
    ]
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_report(path, records, emin, primary, robust, group_note):
    lines = [
        "# Confirmed minima post-processing — [K(H2O)6]+ 2-structure pilot",
        "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        "Numbers below are re-read from completed xTB files. No fabricated QC values.",
        "",
        "## Geometries kept",
        "",
        "Only structures with verified free optimization, full 3N frequencies, and **zero** significant imaginary modes (|ν| > 5 cm⁻¹).",
        "Translational/rotational near-zero modes (|ν| ≤ 5 cm⁻¹) are not counted as imaginary.",
        "",
    ]
    for r in records:
        lines += [
            f"### {r['id']} (`{r['label']}`)",
            "",
            f"- origin: {r['origin']}",
            f"- formula: {r['formula']} ({r['n_atoms']} atoms)",
            f"- E = {r['energy_Eh']:.14f} Eh",
            f"- ΔE = {r['delta_kcal_mol']:.8f} kcal/mol vs lowest of this set",
            f"- opt_converged / normal_termination: {r['opt_converged']} / {r['normal_termination']}",
            f"- n_frequencies = {r['n_frequencies']} (expected 57)",
            f"- n_significant_imaginary = {r['n_significant_imaginary']}",
            f"- n_near_zero = {r['n_near_zero']}",
            f"- lowest vibrational = {r['lowest_vibrational_cm1']} cm⁻¹",
            f"- K–O min/mean/max = {r['geometry']['distance_min_A']:.4f} / {r['geometry']['distance_mean_A']:.4f} / {r['geometry']['distance_max_A']:.4f} Å",
            f"- XYZ: `{r['xyz_path']}`",
            "",
        ]
    lines += [
        f"Lowest energy in this set: **{emin:.14f} Eh** (`{records[0]['id']}`).",
        "",
        "## CoCoS (author snapshot, fixed protocol)",
        "",
        f"- code: `{COCOS}` SHA256 `{sha256(COCOS)}`",
        f"- wrapper: `{WRAPPER}` with `random.seed` and `np.random.seed`",
        f"- production settings: seed **{FIXED_SEED}**, NumSamples **{FIXED_SAMPLES}**, `--DoNotUseGlobalOpt`, `--DoNotSaveDoubles`",
        f"- CritFreq = 83.4 cm⁻¹, NLRMSDFreq = 0.035 cm⁻¹ (CoCoS defaults)",
        f"- input order: energy-sorted, as in `pipeline.dedup`",
        f"- unique representatives at production settings: **{primary['n_unique']}**",
        f"- added indices: {primary['added_indices']}",
        f"- removed: {primary['removed']}",
        f"- {group_note}",
        "",
        "## Seed robustness (same XYZ, same flags)",
        "",
        "| samples | seed | n_unique | added indices | removed |",
        "|---|---|---|---|---|",
    ]
    for r in robust:
        rem = "; ".join(f"{x['structure']}→{x['copy_of']}" for x in r["removed"]) or "—"
        lines.append(
            f"| {r['samples']} | {r['seed']} | {r['n_unique']} | {r['added_indices']} | {rem} |"
        )
    nuniq = {r["n_unique"] for r in robust}
    added_sets = {tuple(r["added_indices"]) for r in robust}
    if len(nuniq) == 1 and len(added_sets) == 1:
        lines.append("")
        lines.append(
            f"For the {len(robust)} runs above, the number of groups and the retained indices did not change. "
            "This is not a proof of algorithmic completeness."
        )
    else:
        lines.append("")
        lines.append(
            "WARNING: CoCoS grouping changed with seed and/or NumSamples. "
            f"n_unique values {sorted(nuniq)}; added-index sets {sorted(added_sets)}."
        )
    lines += [
        "",
        "## Files",
        "",
        f"- `{OUT / 'K06_pilot_confirmed_minima.xyz'}`",
        f"- `{OUT / 'confirmed_minima_table.csv'}`",
        f"- `{OUT / 'confirmed_minima_geometries.png'}` / `.svg`",
        f"- `{OUT / 'confirmed_minima_table.png'}` / `.svg`",
        f"- `{OUT / 'cocos_primary'}`",
        f"- `{OUT / 'cocos_sensitivity'}`",
        "",
        "These three geometries are a 2-structure pilot follow-up, not a 100-structure search.",
        "",
    ]
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise SystemExit(f"{OUT} is not empty; refusing overwrite")
    OUT.mkdir(parents=True)
    verified = [verify_one(spec) for spec in CANDIDATES]
    emin = min(r["energy_Eh"] for r in verified)
    verified.sort(key=lambda r: (r["energy_Eh"], r["id"]))
    for r in verified:
        r["delta_kcal_mol"] = (r["energy_Eh"] - emin) * HARTREE_TO_KCAL

    xyz_sorted = OUT / "energy_sorted_for_cocos.xyz"
    write_multiframe_xyz(xyz_sorted, verified, emin)

    primary = run_cocos(xyz_sorted, OUT / "cocos_primary", FIXED_SEED, FIXED_SAMPLES)
    # "added structure i" uses the original energy-sorted index.
    # "Removed structure i as a copy of j" uses j as the unique-list index (NewTrj).
    survivors = list(dict.fromkeys(primary["added_indices"]))
    groups = {}
    for g, idx in enumerate(survivors):
        groups[verified[idx]["id"]] = g
    for rem in primary["removed"]:
        orig_i = int(rem["structure"])
        uniq_j = int(rem["copy_of"])
        rep_orig = survivors[uniq_j]
        groups[verified[orig_i]["id"]] = groups[verified[rep_orig]["id"]]
    missing = [rec["id"] for rec in verified if rec["id"] not in groups]
    if missing:
        raise RuntimeError(f"CoCoS did not assign a group for {missing}")
    for rec in verified:
        rec["cocos_group"] = groups[rec["id"]]

    robust = [primary]
    for ns in ROBUST_SAMPLES:
        for seed in ROBUST_SEEDS:
            if ns == FIXED_SAMPLES and seed == FIXED_SEED:
                continue
            job = OUT / "cocos_sensitivity" / f"N{ns}_seed{seed}"
            robust.append(run_cocos(xyz_sorted, job, seed, ns))

    rows = []
    for r in verified:
        row = {
            **{k: r[k] for k in (
                "id", "label", "origin", "energy_Eh", "delta_kcal_mol",
                "n_significant_imaginary", "n_near_zero", "lowest_vibrational_cm1",
                "n_frequencies", "opt_converged", "normal_termination",
                "xyz_path", "xyz_sha256", "cocos_group",
            )},
            "formula": "K1O6H12",
            "distance_min_A": r["geometry"]["distance_min_A"],
            "distance_mean_A": r["geometry"]["distance_mean_A"],
            "distance_max_A": r["geometry"]["distance_max_A"],
        }
        rows.append(row)

    n_groups = len(set(r["cocos_group"] for r in verified))
    group_note = (
        f"CoCoS placed the three geometries into {n_groups} equivalence class(es) "
        f"(seed {FIXED_SEED}, NumSamples {FIXED_SAMPLES})."
    )

    write_multiframe_xyz(OUT / "K06_pilot_confirmed_minima.xyz", verified, emin)
    write_csv(OUT / "confirmed_minima_table.csv", rows)
    figure_geometries(
        verified, emin,
        OUT / "confirmed_minima_geometries.png",
        OUT / "confirmed_minima_geometries.svg",
    )
    figure_table(
        rows,
        OUT / "confirmed_minima_table.png",
        OUT / "confirmed_minima_table.svg",
    )
    write_report(OUT / "POSTPROCESS_REPORT.md", verified, emin, primary, robust, group_note)

    serializable = []
    for r in verified:
        d = {k: v for k, v in r.items() if k != "frame"}
        serializable.append(d)
    (OUT / "postprocess_summary.json").write_text(
        json.dumps(
            {
                "emin_Eh": emin,
                "records": serializable,
                "cocos_primary": primary,
                "cocos_robustness": robust,
                "cocos_sha256": sha256(COCOS),
                "fixed_seed": FIXED_SEED,
                "fixed_samples": FIXED_SAMPLES,
            },
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )
    # Also place the annotated XYZ at the pilot root, without touching gwh/xtb_runs.
    dest_xyz = PILOT / "K06_pilot_confirmed_minima.xyz"
    if dest_xyz.exists():
        raise SystemExit(f"refusing to overwrite {dest_xyz}")
    shutil.copy2(OUT / "K06_pilot_confirmed_minima.xyz", dest_xyz)
    print("WROTE", OUT)
    print("XYZ", dest_xyz)
    print("groups", n_groups, "primary unique", primary["n_unique"])


if __name__ == "__main__":
    main()
