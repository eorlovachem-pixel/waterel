#!/usr/bin/env python3
"""Displace a geometry along an xTB Cartesian normal mode.

Does not invent frequencies. Reads g98.out written by xtb --hess.
Displacement convention (documented, not a change of the GWH/xTB method):
the chosen mode vector is scaled so the largest atomic Cartesian
displacement is exactly `step_A` angstrom, then applied as ±step_A.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

FREQ_LINE_RE = re.compile(r"^\s*Frequencies --(.+)$", re.I)
ATOM_LINE_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*$"
)
VIB_ROW_RE = re.compile(
    r"^\s*(\d+)\s+\S*\s+([-+]?\d+\.\d+)\s+"
)


def read_xyz(path):
    p = Path(path)
    lines = p.read_text(encoding="utf-8").splitlines()
    n = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    symbols = []
    xyz = []
    for line in lines[2 : 2 + n]:
        toks = line.split()
        symbols.append(toks[0])
        xyz.append([float(toks[1]), float(toks[2]), float(toks[3])])
    return {"symbols": symbols, "xyz": np.array(xyz, float), "comment": comment}


def write_xyz(path, frame):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(frame["symbols"])), frame.get("comment", "")]
    for s, (x, y, z) in zip(frame["symbols"], frame["xyz"]):
        lines.append(f"{s:2s} {x:.10f} {y:.10f} {z:.10f}")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_vibspectrum(path):
    p = Path(path)
    rows = []
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.strip().startswith("$") or line.strip().startswith("#"):
            continue
        toks = line.split()
        if len(toks) < 3:
            continue
        try:
            idx = int(toks[0])
        except ValueError:
            continue
        # xtb: mode [sym] wavenumber intensity ...
        freq = None
        for tok in toks[1:]:
            try:
                freq = float(tok)
                break
            except ValueError:
                continue
        if freq is None:
            continue
        rows.append({"mode_index": idx, "frequency_cm1": freq})
    return rows


def parse_g98_modes(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        m = FREQ_LINE_RE.match(lines[i])
        if not m:
            i += 1
            continue
        freqs = [float(x) for x in m.group(1).split()]
        n_cols = len(freqs)
        # seek Atom AN header
        j = i + 1
        while j < len(lines) and "Atom AN" not in lines[j]:
            j += 1
        j += 1
        disp = [[] for _ in range(n_cols)]
        while j < len(lines):
            am = ATOM_LINE_RE.match(lines[j])
            if not am:
                break
            vals = [float(am.group(k)) for k in range(3, 12)]
            for col in range(n_cols):
                disp[col].append(vals[3 * col : 3 * col + 3])
            j += 1
        for col, freq in enumerate(freqs):
            arr = np.array(disp[col], float)
            blocks.append({"frequency_cm1": freq, "cartesian": arr})
        i = j
    if not blocks:
        raise RuntimeError(f"{path}: no g98.out frequency/mode blocks parsed")
    return blocks


def classify_modes(vib_rows, imag_thresh_cm1=5.0):
    thresh = abs(float(imag_thresh_cm1))
    freqs = [r["frequency_cm1"] for r in vib_rows]
    near_zero = [r for r in vib_rows if abs(r["frequency_cm1"]) <= thresh]
    significant_imag = [r for r in vib_rows if r["frequency_cm1"] < -thresh]
    real_vib = [r for r in vib_rows if r["frequency_cm1"] > thresh]
    return {
        "n_modes": len(vib_rows),
        "imag_thresh_cm1": thresh,
        "n_near_zero": len(near_zero),
        "near_zero_cm1": [r["frequency_cm1"] for r in near_zero],
        "n_significant_imaginary": len(significant_imag),
        "significant_imaginary": significant_imag,
        "n_real_vibrational": len(real_vib),
        "lowest_frequency_cm1": min(freqs) if freqs else None,
        "frequencies_cm1": freqs,
    }


def choose_imaginary_g98_mode(g98_modes, imag_thresh_cm1=5.0):
    imag = [m for m in g98_modes if m["frequency_cm1"] < -abs(float(imag_thresh_cm1))]
    if not imag:
        return None
    return min(imag, key=lambda m: m["frequency_cm1"])


def scale_mode_max_atom(cartesian, step_A):
    atom_norms = np.linalg.norm(cartesian, axis=1)
    max_norm = float(atom_norms.max())
    if max_norm <= 0.0 or not math.isfinite(max_norm):
        raise RuntimeError("mode vector has zero/non-finite max atomic displacement")
    scaled = cartesian * (float(step_A) / max_norm)
    return scaled, {
        "convention": "largest_atomic_cartesian_displacement_equals_step_A",
        "step_A": float(step_A),
        "unscaled_max_atom_A": max_norm,
        "unscaled_rms_atom_A": float(np.sqrt(np.mean(atom_norms**2))),
        "scaled_max_atom_A": float(np.linalg.norm(scaled, axis=1).max()),
        "scaled_rms_atom_A": float(np.sqrt(np.mean(np.linalg.norm(scaled, axis=1) ** 2))),
    }


def displace(frame, scaled_mode, sign, comment):
    xyz = frame["xyz"] + sign * scaled_mode
    return {
        "symbols": list(frame["symbols"]),
        "xyz": xyz,
        "comment": comment,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--g98", required=True)
    ap.add_argument("--vibspectrum")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--step-A", type=float, default=0.15)
    ap.add_argument("--imag-thresh", type=float, default=5.0)
    args = ap.parse_args()
    frame = read_xyz(args.xyz)
    modes = parse_g98_modes(args.g98)
    vib = parse_vibspectrum(args.vibspectrum) if args.vibspectrum else []
    cls = classify_modes(vib, args.imag_thresh) if vib else None
    chosen = choose_imaginary_g98_mode(modes, args.imag_thresh)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    report = {
        "xyz": str(Path(args.xyz).resolve()),
        "g98": str(Path(args.g98).resolve()),
        "vibspectrum": str(Path(args.vibspectrum).resolve()) if args.vibspectrum else None,
        "g98_n_modes": len(modes),
        "g98_frequencies_cm1": [m["frequency_cm1"] for m in modes],
        "vibspectrum_classification": cls,
        "chosen_g98_frequency_cm1": None,
        "displaced": False,
    }
    if chosen is None:
        (out / "mode_follow.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print("No significant imaginary mode in g98.out; no displacement written.")
        return
    scaled, meta = scale_mode_max_atom(chosen["cartesian"], args.step_A)
    report["chosen_g98_frequency_cm1"] = chosen["frequency_cm1"]
    report["scale"] = meta
    report["displaced"] = True
    for sign, name in ((+1.0, "plus"), (-1.0, "minus")):
        tagged = displace(
            frame,
            scaled,
            sign,
            comment=(
                f"displaced {name} along imag mode nu={chosen['frequency_cm1']:.4f} cm-1 "
                f"max-atom {args.step_A} A"
            ),
        )
        xyz_path = out / f"displaced_{name}_{args.step_A:.2f}A.xyz"
        write_xyz(xyz_path, tagged)
        report.setdefault("files", {})[name] = str(xyz_path)
    (out / "mode_follow.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
