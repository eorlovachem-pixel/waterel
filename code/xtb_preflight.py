#!/usr/bin/env python3
"""Real xTB geometry-optimization and Hessian smoke tests.

xtb --version is not accepted as proof. This script runs ANCOPT and a Hessian
on tiny systems and records the actual logs. It never fabricates energies.
It does not replace /usr/bin/xtb.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import freq_parser
import qc_checks

STANDALONE_DEFAULT = Path.home() / "xtb-standalone" / "xtb-dist" / "bin" / "xtb"
SYSTEM_XTB = Path("/usr/bin/xtb")

WATER_XYZ = """3
preflight water
O   0.0000000000   0.0000000000   0.1173000000
H   0.0000000000   0.7572000000  -0.4692000000
H   0.0000000000  -0.7572000000  -0.4692000000
"""

K_WATER_XYZ = """4
preflight K+(H2O) charge +1
K   0.0000000000   0.0000000000   0.0000000000
O   0.0000000000   0.0000000000   2.7000000000
H   0.7572000000   0.0000000000   3.2865000000
H  -0.7572000000   0.0000000000   3.2865000000
"""


def run_xtb(xtb_bin, workdir: Path, args, timeout):
    workdir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    # Do not hide /usr/bin/xtb globally; only the explicit binary is invoked.
    cmd = [str(xtb_bin), *args]
    log = workdir / "xtb.out"
    with log.open("w", encoding="utf-8") as fh:
        proc = subprocess.run(
            cmd,
            cwd=workdir,
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    text = log.read_text(errors="replace")
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "log": str(log),
        "log_text": text,
        "workdir": str(workdir),
        "parsed": qc_checks.parse_xtb_log(text),
    }


def summarize_opt(result, expected_symbols):
    text = result["log_text"]
    out = {
        "cmd": result["cmd"],
        "returncode": result["returncode"],
        "workdir": result["workdir"],
        "format_crash": result["parsed"]["format_crash"],
        "normal_termination": result["parsed"]["normal_termination"],
        "opt_converged": result["parsed"]["opt_converged"],
        "total_energy_Eh": result["parsed"]["total_energy_Eh"],
        "version_line": result["parsed"]["version_line"],
        "xtbopt_exists": (Path(result["workdir"]) / "xtbopt.xyz").is_file(),
        "ok": False,
        "error": None,
    }
    try:
        if result["parsed"]["format_crash"]:
            raise RuntimeError("Fortran format crash in ANCOPT")
        opt = Path(result["workdir"]) / "xtbopt.xyz"
        if not opt.is_file():
            raise RuntimeError("xtbopt.xyz missing")
        frame = qc_checks.read_xyz(opt)[0]
        energy, _info = qc_checks.require_completed_opt(text, frame, expected_symbols)
        out["total_energy_Eh"] = energy
        out["ok"] = True
    except Exception as exc:
        out["error"] = str(exc)
    return out


def summarize_hess(result):
    out = {
        "cmd": result["cmd"],
        "returncode": result["returncode"],
        "workdir": result["workdir"],
        "format_crash": result["parsed"]["format_crash"],
        "normal_termination": result["parsed"]["normal_termination"],
        "ok": False,
        "error": None,
        "hessian": None,
    }
    try:
        hess = freq_parser.require_completed_hessian(result["log"])
        out["hessian"] = {
            "n_frequencies_parsed": hess["n_frequencies_parsed"],
            "imaginary_cm1": hess["imaginary_cm1"],
            "classification": hess["classification"],
            "is_genuine_minimum_candidate": hess["is_genuine_minimum_candidate"],
            "frequencies_cm1": hess["frequencies_cm1"],
        }
        out["ok"] = True
    except Exception as exc:
        out["error"] = str(exc)
        try:
            out["hessian"] = freq_parser.classify_hessian(result["log_text"])
        except Exception:
            pass
    return out


def write_xyz(path: Path, body: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--xtb",
        default=None,
        help="working xTB binary; default is env_xtb.sh / ~/xtb-standalone, never /usr/bin/xtb",
    )
    p.add_argument("--also-test-system-xtb", action="store_true", default=True)
    p.add_argument("--skip-system-xtb", action="store_true")
    p.add_argument("--timeout", type=int, default=180)
    args = p.parse_args()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    root = ROOT / "preflight_runs" / stamp
    root.mkdir(parents=True, exist_ok=True)
    report = {
        "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
        "system_xtb": str(SYSTEM_XTB) if SYSTEM_XTB.exists() else None,
        "working_xtb": None,
        "tests": {},
        "usable_for_free_optimization": False,
        "note": "Values below are parsed from real xTB logs. Missing jobs are failures, not invented numbers.",
    }

    working = Path(qc_checks.resolve_working_xtb(args.xtb, kit=ROOT))
    report["working_xtb"] = str(working)
    if SYSTEM_XTB.exists() and working == SYSTEM_XTB.resolve():
        report["warning"] = (
            "Requested working xtb is /usr/bin/xtb. Ubuntu 26.04 libgfortran 16 "
            "is known to crash ANCOPT in xtb 6.7.1."
        )

    if args.also_test_system_xtb and not args.skip_system_xtb and SYSTEM_XTB.exists():
        d = root / "system_xtb_water_opt"
        write_xyz(d / "input.xyz", WATER_XYZ)
        raw = run_xtb(SYSTEM_XTB, d, ["input.xyz", "--gfn", "2", "--opt", "tight"], args.timeout)
        (d / "excerpt.txt").write_text("\n".join(raw["log_text"].splitlines()[-40:]), encoding="utf-8")
        report["tests"]["system_xtb_water_opt"] = summarize_opt(raw, ["O", "H", "H"])

    d = root / "standalone_water_opt"
    write_xyz(d / "input.xyz", WATER_XYZ)
    raw = run_xtb(working, d, ["input.xyz", "--gfn", "2", "--opt", "tight"], args.timeout)
    water_opt = summarize_opt(raw, ["O", "H", "H"])
    report["tests"]["standalone_water_opt"] = water_opt

    if water_opt["ok"]:
        hd = root / "standalone_water_hess"
        hd.mkdir(parents=True, exist_ok=True)
        shutil.copy2(d / "xtbopt.xyz", hd / "xtbopt.xyz")
        rawh = run_xtb(working, hd, ["xtbopt.xyz", "--gfn", "2", "--hess"], args.timeout)
        report["tests"]["standalone_water_hess"] = summarize_hess(rawh)

    kd = root / "standalone_k_water_opt"
    write_xyz(kd / "input.xyz", K_WATER_XYZ)
    rawk = run_xtb(
        working,
        kd,
        ["input.xyz", "--gfn", "2", "--chrg", "1", "--alpb", "water", "--opt", "tight"],
        args.timeout,
    )
    kopt = summarize_opt(rawk, ["K", "O", "H", "H"])
    report["tests"]["standalone_k_water_alpb_opt"] = kopt
    if kopt["ok"]:
        kh = root / "standalone_k_water_hess"
        kh.mkdir(parents=True, exist_ok=True)
        shutil.copy2(kd / "xtbopt.xyz", kh / "xtbopt.xyz")
        rawkh = run_xtb(
            working,
            kh,
            ["xtbopt.xyz", "--gfn", "2", "--chrg", "1", "--alpb", "water", "--hess"],
            args.timeout,
        )
        report["tests"]["standalone_k_water_alpb_hess"] = summarize_hess(rawkh)

    needed = [
        "standalone_water_opt",
        "standalone_water_hess",
        "standalone_k_water_alpb_opt",
        "standalone_k_water_alpb_hess",
    ]
    report["usable_for_free_optimization"] = all(
        report["tests"].get(name, {}).get("ok") for name in needed
    )
    report["root"] = str(root)
    out = root / "preflight_report.json"
    # Do not dump full logs into JSON; they remain as xtb.out files.
    slim = json.loads(json.dumps(report))
    out.write_text(json.dumps(slim, indent=2, ensure_ascii=False), encoding="utf-8")
    latest = ROOT / "preflight_runs" / "latest.json"
    latest.write_text(out.read_text(encoding="utf-8"), encoding="utf-8")
    print(json.dumps(slim, indent=2, ensure_ascii=False))
    if not report["usable_for_free_optimization"]:
        raise SystemExit("xTB preflight failed; do not start free optimization or a 100-structure series.")


if __name__ == "__main__":
    main()
