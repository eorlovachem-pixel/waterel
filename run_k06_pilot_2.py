#!/usr/bin/env python3
"""End-to-end 2-structure K+(H2O)6 pilot. Never fabricates QC numbers.

Writes under results/K_06/pilot_2/ via a symlink from pilots/K_06 so the
existing workflow.py is used unchanged. Does not touch Li_07 or cases/.
Does not run 100 structures.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
import freq_parser  # noqa: E402
import mode_follow  # noqa: E402
import qc_checks  # noqa: E402

XTB = Path.home() / "xtb-standalone" / "xtb-dist" / "bin" / "xtb"
SYSTEM_XTB = Path("/usr/bin/xtb")
NSTRUCT = 2
IMAG_THRESH = 5.0
STEP_A = 0.15
ION = "K"
N_WATER = 6
CHARGE = 1
VEC_SHA = "06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5"
GWH_SHA = "7411811374b6422a15f23d6d86a47f2145f8ac47"


def run(cmd, cwd, log_path, env=None, timeout=None):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    with log_path.open("w", encoding="utf-8") as fh:
        fh.write(f"# cmd: {cmd}\n# cwd: {cwd}\n# start_UTC: {started}\n")
        fh.flush()
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            env=env or os.environ.copy(),
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=timeout,
        )
    ended = datetime.now(timezone.utc).isoformat()
    rec = {
        "cmd": cmd,
        "cwd": str(cwd),
        "log": str(log_path),
        "returncode": proc.returncode,
        "start_UTC": started,
        "end_UTC": ended,
    }
    print(json.dumps(rec), flush=True)
    return rec


def imag_score(frame):
    vib = frame.get("vibspectrum_classification") or {}
    sig = vib.get("significant_imaginary") or []
    if sig:
        return min(x["frequency_cm1"] for x in sig)
    im = frame.get("hess", {}).get("imaginary_cm1") or []
    if im:
        return min(im)
    return 0.0


def write_report(dest, vec_info, euler, gwh_head, frames, follow, commands):
    lines = []
    lines.append("# K+(H2O)6 2-structure pilot report")
    lines.append("")
    lines.append(f"Generated UTC: {datetime.now(timezone.utc).isoformat()}")
    lines.append("This file contains only values read from completed xTB/GWH files. No fabricated energies.")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("- Case: `K_06`")
    lines.append("- Structures requested: **2**")
    lines.append("- 100-structure series: **not run**")
    lines.append("- Previous `pilots/Li_07` and `cases/Li_07`: not modified")
    lines.append(f"- Work directory: `{dest}`")
    lines.append(f"- Pipeline symlink: `{ROOT / 'pilots' / 'K_06'} -> {dest}`")
    lines.append("")
    lines.append("## Inputs and programs")
    lines.append("")
    lines.append(f"- packing.vecfle SHA256: `{vec_info['sha256']}`")
    lines.append(f"- ion / n_water / radius: {vec_info['ion']} / {vec_info['n_water']} / {vec_info['radius_A']} Å")
    lines.append(f"- GWH git: `{gwh_head}`")
    lines.append(
        f"- GWH Euler: author `from_euler('zyz')` present={euler.get('author_method_intact')}; "
        f"SciPy interprets lowercase as {euler.get('scipy_interprets_author_zyz_as')} "
        "(source not rewritten)"
    )
    lines.append(f"- xTB binary: `{XTB}` SHA256 `{qc_checks.sha256_file(XTB)}`")
    lines.append("- Method: GFN2-xTB, ALPB(water), charge +1")
    lines.append("- Free opt: unconstrained `--opt tight`")
    lines.append("- Hessian threshold for **significant** imaginary modes: |ν| > 5.0 cm⁻¹")
    lines.append("- Near-zero trans/rot: |ν| ≤ 5.0 cm⁻¹, not counted as imaginary")
    lines.append("")
    lines.append("## Commands")
    lines.append("")
    for i, rec in enumerate(commands, 1):
        lines.append(f"{i}. returncode={rec['returncode']} `{' '.join(str(x) for x in rec['cmd'])}`")
        lines.append(f"   log: `{rec['log']}`")
    lines.append("")
    lines.append("## Constrained GWH stage")
    lines.append("")
    gwh_log = dest / "gwh_stdout_stderr.log"
    gwh_exit = dest / "gwh_exit_code.txt"
    lines.append(f"- GWH exit file: `{gwh_exit}` contents `{gwh_exit.read_text().strip() if gwh_exit.exists() else 'MISSING'}`")
    lines.append(f"- GWH combined log: `{gwh_log}` exists={gwh_log.is_file()}")
    src_xyz = dest / "source" / "constrained_multiframe.xyz"
    src_en = dest / "source" / "energies.csv"
    lines.append(f"- Constrained multi-XYZ: `{src_xyz}` exists={src_xyz.is_file()}")
    lines.append(f"- Constrained energies CSV: `{src_en}` exists={src_en.is_file()}")
    if src_xyz.is_file():
        constr = qc_checks.read_xyz(src_xyz)
        lines.append(f"- Constrained frames parsed: {len(constr)}")
        for i, fr in enumerate(constr):
            qc_checks.assert_cluster_composition(fr, ION, N_WATER, CHARGE)
            e = qc_checks.energy_from_comment(fr["comment"])
            lines.append(f"  - structure {i:04d}: {len(fr['symbols'])} atoms, formula {qc_checks.formula_of(fr['symbols'])}, E_constrained = {e} Eh")
            lines.append(f"    file: `{dest / 'gwh' / f'gwh_opt_str{i:04d}' / f'mol_str{i:04d}_fin.xyz'}`")
    lines.append("")
    lines.append("## Free optimization (GFN2-xTB / ALPB water, charge +1)")
    lines.append("")
    lines.append("Success requires: process exit 0 **and** `GEOMETRY OPTIMIZATION CONVERGED` **and** `normal termination of xtb` **and** matching TOTAL ENERGY vs `xtbopt.xyz`.")
    lines.append("")
    for fr in frames:
        h = fr["hess"]
        vib = fr.get("vibspectrum_classification") or {}
        lines.append(f"### {fr['id']}")
        lines.append("")
        lines.append(f"- E_free = **{fr['opt_energy_Eh']:.14f} Eh**")
        lines.append(f"- opt_converged: {fr['opt_info'].get('opt_converged')}")
        lines.append(f"- normal_termination: {fr['opt_info'].get('normal_termination')}")
        lines.append(f"- geometry: `{fr['xtbopt_xyz']}`")
        lines.append(f"- opt log: `{fr['opt_log']}`")
        lines.append(f"- hess log: `{fr['hess_log']}`")
        lines.append(f"- Hessian classification (parser): `{h.get('classification')}`")
        lines.append(f"- n_frequencies parsed from xtb log: {h.get('n_frequencies_parsed')}")
        lines.append(f"- significant imaginary from log parser (ν < -5 cm⁻¹): {h.get('imaginary_cm1')}")
        lines.append(f"- near-zero from log parser (|ν| ≤ 5 cm⁻¹): {h.get('near_zero_cm1')}")
        if vib:
            lines.append(f"- vibspectrum n_modes: {vib.get('n_modes')}")
            lines.append(f"- vibspectrum near-zero count: {vib.get('n_near_zero')} values {vib.get('near_zero_cm1')}")
            lines.append(f"- vibspectrum significant imaginary count: {vib.get('n_significant_imaginary')}")
            lines.append(f"- vibspectrum significant imaginary: {vib.get('significant_imaginary')}")
            lines.append(f"- vibspectrum lowest frequency: {vib.get('lowest_frequency_cm1')} cm⁻¹")
            lines.append(f"- vibspectrum all frequencies (cm⁻¹): {vib.get('frequencies_cm1')}")
        lines.append(f"- g98.out vibrational frequencies (cm⁻¹): {fr.get('g98_frequencies_cm1')}")
        lines.append("")
    lines.append("## Mode following (±0.15 Å)")
    lines.append("")
    if follow is None:
        lines.append("No structure had a significant imaginary frequency (ν < -5 cm⁻¹). Displacement was not required and was not run.")
    else:
        lines.append(f"Selected structure: `{follow['id']}` (most negative significant imaginary mode).")
        disp = follow.get("displacement") or {}
        lines.append(f"- Chosen g98 frequency: {disp.get('chosen_g98_frequency_cm1')} cm⁻¹")
        lines.append(f"- Scale convention: {(disp.get('scale') or {}).get('convention')}")
        lines.append(f"- step: {(disp.get('scale') or {}).get('step_A')} Å")
        for br in follow.get("branches") or []:
            lines.append(f"### branch `{br['name']}`")
            lines.append(f"- opt returncode: {(br.get('opt') or {}).get('returncode')}")
            ver = br.get("verified") or {}
            lines.append(f"- verified: {json.dumps(ver, default=str)}")
            hs = br.get("hess") or {}
            if hs:
                lines.append(f"- re-hess parser classification: {(hs.get('parser') or {}).get('classification')}")
                lines.append(f"- re-hess significant imag: {(hs.get('parser') or {}).get('imaginary_cm1')}")
                lines.append(f"- re-hess vibspectrum: {json.dumps(hs.get('vibspectrum'), default=str)}")
            lines.append("")
    lines.append("## File map")
    lines.append("")
    lines.append(f"- `{dest / 'gwh'}` constrained GWH workdirs")
    lines.append(f"- `{dest / 'xtb_runs'}` free opt + Hessian")
    lines.append(f"- `{dest / 'provenance'}` copied packing.vecfle, commands, environment")
    lines.append(f"- `{dest / 'pilot_summary.json'}` machine-readable copy of this report's numbers")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (dest / "PILOT_REPORT.md").write_text(text, encoding="utf-8")
    print("WROTE", dest / "PILOT_REPORT.md")


def main():
    dest = ROOT / "results" / "K_06" / "pilot_2"
    pilots = ROOT / "pilots"
    pilots_k06 = pilots / "K_06"
    cases_k06 = ROOT / "cases" / "K_06"
    li07 = ROOT / "pilots" / "Li_07"
    vec = ROOT / "inputs" / "K_06" / "packing.vecfle"
    gwh = ROOT / "vendors" / "Glorious-Water-Hole" / "gloriouswaterhole.py"

    if cases_k06.exists():
        raise SystemExit(f"Refusing to run: {cases_k06} already exists (would mix with production).")
    if not li07.is_dir():
        print("WARNING: pilots/Li_07 missing; continuing without touching it.")
    if not XTB.is_file():
        raise SystemExit(f"Working xtb missing: {XTB}")
    if XTB.resolve() == SYSTEM_XTB.resolve():
        raise SystemExit("Refusing /usr/bin/xtb for free opt")
    if pilots_k06.exists() and not pilots_k06.is_symlink():
        raise SystemExit(f"{pilots_k06} already exists as a real directory; refusing overwrite.")
    if dest.exists() and any(dest.iterdir()):
        raise SystemExit(f"{dest} is not empty; refusing overwrite of a previous pilot.")

    dest.mkdir(parents=True, exist_ok=True)
    if not pilots_k06.exists():
        rel = os.path.relpath(dest, pilots)
        pilots_k06.symlink_to(rel, target_is_directory=True)

    provenance = dest / "provenance"
    provenance.mkdir(exist_ok=True)
    shutil.copy2(vec, provenance / "packing.vecfle")
    (provenance / "packing.vecfle.sha256").write_text(
        qc_checks.sha256_file(vec) + "\n", encoding="utf-8"
    )
    vec_info = qc_checks.validate_vecfle_bytes(vec, ION, N_WATER, 3.04, VEC_SHA)
    euler = qc_checks.gwh_uses_intrinsic_zyz(gwh)
    if not euler["intrinsic_zyz"]:
        raise SystemExit(f"GWH is not intrinsic ZYZ: {euler}")
    gwh_head = subprocess.run(
        ["git", "-C", str(gwh.parent), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    if gwh_head != GWH_SHA:
        raise SystemExit(f"GWH SHA {gwh_head} != pinned {GWH_SHA}")

    env = os.environ.copy()
    env["PATH"] = str(XTB.parent) + os.pathsep + env.get("PATH", "")
    env["XTB_BIN"] = str(XTB)
    env["XTBPATH"] = str(XTB.parent.parent / "share" / "xtb")
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    # Do not seed GWH: author --Random uses unseeded random.uniform.

    py = sys.executable
    wf = str(ROOT / "workflow.py")
    xtb = str(XTB)
    commands = []

    commands.append(
        run(
            [py, wf, "gwh", "--case", "K_06", "--pilot", "--nstruct", str(NSTRUCT), "--xtb", xtb],
            ROOT,
            dest / "logs" / "01_gwh.log",
            env=env,
            timeout=7200,
        )
    )
    if commands[-1]["returncode"] != 0:
        raise SystemExit("GWH failed; see logs/01_gwh.log")
    commands.append(run([py, wf, "collect", "--case", "K_06", "--pilot"], ROOT, dest / "logs" / "02_collect.log", env=env))
    if commands[-1]["returncode"] != 0:
        raise SystemExit("collect failed")
    commands.append(run([py, wf, "prepare", "--case", "K_06", "--pilot"], ROOT, dest / "logs" / "03_prepare.log", env=env))
    if commands[-1]["returncode"] != 0:
        raise SystemExit("prepare failed")
    commands.append(
        run(
            [py, wf, "relax", "--case", "K_06", "--pilot", "--xtb", xtb, "--workers", "1"],
            ROOT,
            dest / "logs" / "04_relax.log",
            env=env,
            timeout=7200,
        )
    )
    if commands[-1]["returncode"] != 0:
        raise SystemExit("relax failed")
    commands.append(
        run(
            [py, wf, "hess", "--case", "K_06", "--pilot", "--xtb", xtb, "--workers", "1", "--imag-thresh", str(IMAG_THRESH)],
            ROOT,
            dest / "logs" / "05_hess.log",
            env=env,
            timeout=7200,
        )
    )
    if commands[-1]["returncode"] != 0:
        raise SystemExit("hess failed")

    (provenance / "commands.json").write_text(json.dumps(commands, indent=2), encoding="utf-8")
    (provenance / "environment.json").write_text(
        json.dumps(
            {
                "xtb": xtb,
                "xtb_sha256": qc_checks.sha256_file(XTB),
                "gwh_git": gwh_head,
                "gwh_euler": euler,
                "vecfle": vec_info,
                "python": sys.version,
                "which_xtb_during_run": shutil.which("xtb", path=env["PATH"]),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    # Verify each frame and classify Hessian.
    frames = []
    expected = None
    xtb_runs = dest / "xtb_runs"
    for i in range(NSTRUCT):
        fid = f"frame_{i:03d}"
        folder = xtb_runs / fid
        opt = folder / "xtbopt.xyz"
        log = folder / "xtb.out"
        hdir = folder / "hess"
        hlog = hdir / "xtb_hess.out"
        xyz = qc_checks.read_xyz(opt)[0]
        if expected is None:
            expected = xyz["symbols"]
        qc_checks.assert_cluster_composition(xyz, ION, N_WATER, CHARGE)
        energy, opt_info = qc_checks.require_completed_opt(log.read_text(errors="replace"), xyz, expected)
        hess = freq_parser.require_completed_hessian(hlog, imag_thresh_cm1=IMAG_THRESH)
        vib_path = hdir / "vibspectrum"
        vib_rows = mode_follow.parse_vibspectrum(vib_path) if vib_path.is_file() else []
        vib_cls = mode_follow.classify_modes(vib_rows, IMAG_THRESH) if vib_rows else None
        g98 = hdir / "g98.out"
        g98_modes = mode_follow.parse_g98_modes(g98) if g98.is_file() else []
        rec = {
            "id": fid,
            "opt_energy_Eh": energy,
            "opt_info": {k: opt_info[k] for k in opt_info if k != "version_line"} | {"version_line": opt_info.get("version_line")},
            "xtbopt_xyz": str(opt),
            "opt_log": str(log),
            "hess_log": str(hlog),
            "hess": hess,
            "vibspectrum_classification": vib_cls,
            "g98_frequencies_cm1": [m["frequency_cm1"] for m in g98_modes],
            "n_g98_modes": len(g98_modes),
        }
        frames.append(rec)

    follow = None
    imag_frames = [
        f
        for f in frames
        if (f["vibspectrum_classification"] or {}).get("n_significant_imaginary", 0) > 0
        or f["hess"]["imaginary_cm1"]
    ]
    if imag_frames:
        target = min(imag_frames, key=imag_score)
        folder = xtb_runs / target["id"]
        hdir = folder / "hess"
        mfdir = folder / "mode_follow"
        mfdir.mkdir(exist_ok=True)
        mf_cmd = [
            py,
            str(ROOT / "code" / "mode_follow.py"),
            "--xyz", str(folder / "xtbopt.xyz"),
            "--g98", str(hdir / "g98.out"),
            "--vibspectrum", str(hdir / "vibspectrum"),
            "--outdir", str(mfdir),
            "--step-A", str(STEP_A),
            "--imag-thresh", str(IMAG_THRESH),
        ]
        commands.append(run(mf_cmd, ROOT, dest / "logs" / "06_mode_follow_displace.log", env=env))
        follow = {"id": target["id"], "displacement": json.loads((mfdir / "mode_follow.json").read_text()), "branches": []}
        for name in ("plus", "minus"):
            xyz_in = mfdir / f"displaced_{name}_{STEP_A:.2f}A.xyz"
            if not xyz_in.is_file():
                continue
            bdir = mfdir / f"{name}_{STEP_A:.2f}A"
            (bdir / "opt").mkdir(parents=True, exist_ok=True)
            shutil.copy2(xyz_in, bdir / "opt" / "input.xyz")
            opt_cmd = [xtb, "input.xyz", "--gfn", "2", "--chrg", str(CHARGE), "--alpb", "water", "--opt", "tight"]
            rec_opt = run(opt_cmd, bdir / "opt", bdir / "opt" / "xtb.out", env=env, timeout=3600)
            branch = {"name": name, "opt": rec_opt, "hess": None, "verified": None}
            opt_xyz = bdir / "opt" / "xtbopt.xyz"
            opt_log = bdir / "opt" / "xtb.out"
            if opt_xyz.is_file():
                fr = qc_checks.read_xyz(opt_xyz)[0]
                try:
                    e, info = qc_checks.require_completed_opt(opt_log.read_text(errors="replace"), fr, expected)
                    branch["verified"] = {"energy_Eh": e, "opt_converged": info["opt_converged"], "normal_termination": info["normal_termination"]}
                    h2 = bdir / "hess"
                    h2.mkdir(exist_ok=True)
                    shutil.copy2(opt_xyz, h2 / "xtbopt.xyz")
                    hess_cmd = [xtb, "xtbopt.xyz", "--gfn", "2", "--chrg", str(CHARGE), "--alpb", "water", "--hess"]
                    rec_h = run(hess_cmd, h2, h2 / "xtb_hess.out", env=env, timeout=3600)
                    parsed = freq_parser.require_completed_hessian(h2 / "xtb_hess.out", imag_thresh_cm1=IMAG_THRESH)
                    vib2 = mode_follow.parse_vibspectrum(h2 / "vibspectrum") if (h2 / "vibspectrum").is_file() else []
                    branch["hess"] = {
                        "returncode": rec_h["returncode"],
                        "parser": parsed,
                        "vibspectrum": mode_follow.classify_modes(vib2, IMAG_THRESH) if vib2 else None,
                    }
                except Exception as exc:
                    branch["verified"] = {"error": str(exc)}
            follow["branches"].append(branch)
        (mfdir / "follow_summary.json").write_text(json.dumps(follow, indent=2, default=str), encoding="utf-8")

    (provenance / "commands.json").write_text(json.dumps(commands, indent=2, default=str), encoding="utf-8")
    summary = {
        "case": "K_06",
        "nstruct": NSTRUCT,
        "charge": CHARGE,
        "method": "GFN2-xTB / ALPB(water)",
        "imag_thresh_cm1": IMAG_THRESH,
        "frames": frames,
        "mode_follow": follow,
        "commands": commands,
    }
    (dest / "pilot_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    write_report(dest, vec_info, euler, gwh_head, frames, follow, commands)
    print("PILOT_RAW_SUMMARY_WRITTEN", dest / "pilot_summary.json")


if __name__ == "__main__":
    main()
