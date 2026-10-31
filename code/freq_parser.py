#!/usr/bin/env python3
"""Parse xTB Hessian / frequency logs. Does not invent frequencies."""
from __future__ import annotations

import math
import re
from pathlib import Path

try:
    from qc_checks import FORMAT_CRASH_RE, parse_xtb_log
except ImportError:  # package-style import from repo root
    from code.qc_checks import FORMAT_CRASH_RE, parse_xtb_log

IMAG_WORD_RE = re.compile(r"\bimaginary\b", re.I)
FREQ_BLOCK_RE = re.compile(
    r"(?:projected\s+)?(?:vibrational\s+)?frequencies(?:\s*\(cm(?:\*\*-1|\u207b\u00b9)?\))?",
    re.I,
)
FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")


def _floats_in_line(line: str):
    return [float(x) for x in FLOAT_RE.findall(line)]


def parse_frequencies_cm1(log_text: str) -> list[float]:
    """Collect vibrational frequencies from an xTB Hessian log.

    xTB prints a table after a frequencies header. Index-only lines are ignored.
    Values are returned in cm^-1 in the printed order. Empty list means the
    parser could not find a table — that is a failure, not zero imaginary modes.
    """
    lines = log_text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if FREQ_BLOCK_RE.search(line) and "reduced" not in line.lower():
            start = i + 1
            break
    if start is None:
        return []
    freqs = []
    for line in lines[start: start + 80]:
        low = line.lower()
        if any(k in low for k in ("reduced mass", "ir intensity", "raman", "character", "-----")):
            if freqs:
                break
            continue
        values = _floats_in_line(line)
        if not values:
            if freqs and line.strip() == "":
                continue
            if freqs:
                break
            continue
        # Header index rows are small positive integers 1..3N-6
        if all(v.is_integer() and 1 <= v <= 500 for v in values) and max(values) <= len(values) + 12:
            continue
        freqs.extend(values)
    return freqs


def classify_hessian(log_text: str, imag_thresh_cm1: float = 5.0) -> dict:
    if FORMAT_CRASH_RE.search(log_text):
        raise RuntimeError("Hessian log contains the ANCOPT Fortran format crash")
    info = parse_xtb_log(log_text)
    freqs = parse_frequencies_cm1(log_text)
    imag = [f for f in freqs if f < -abs(float(imag_thresh_cm1))]
    near_zero = [f for f in freqs if abs(f) <= abs(float(imag_thresh_cm1))]
    n_imag_words = len(IMAG_WORD_RE.findall(log_text))
    status = "NOT_A_MINIMUM" if imag else "NO_IMAGINARY_ABOVE_THRESHOLD"
    if not freqs:
        status = "FREQUENCIES_NOT_PARSED"
    return {
        "normal_termination": info["normal_termination"],
        "format_crash": info["format_crash"],
        "has_frequency_section": info["has_frequency_section"] or bool(freqs),
        "n_frequencies_parsed": len(freqs),
        "frequencies_cm1": freqs,
        "imaginary_cm1": imag,
        "near_zero_cm1": near_zero,
        "imaginary_word_hits": n_imag_words,
        "imag_thresh_cm1": float(imag_thresh_cm1),
        "classification": status,
        "nu_min_cm1": -abs(float(imag_thresh_cm1)),
        "all_nu_ge_minus_thresh": bool(freqs) and all(
            f >= -abs(float(imag_thresh_cm1)) for f in freqs
        ),
        "is_genuine_minimum_candidate": (
            info["normal_termination"]
            and bool(freqs)
            and not imag
            and not info["format_crash"]
        ),
    }


def is_confirmed_minimum(result, nu_min_cm1: float = -5.0) -> bool:
    """Confirmed minimum for CoCoS: all parsed ν ≥ nu_min_cm1 (default −5 cm⁻¹)."""
    freqs = result.get("frequencies_cm1") or []
    if not freqs:
        return False
    if result.get("format_crash"):
        return False
    if not result.get("normal_termination"):
        return False
    return all(f >= float(nu_min_cm1) for f in freqs)


def require_completed_hessian(log_path, imag_thresh_cm1: float = 5.0) -> dict:
    p = Path(log_path)
    text = p.read_text(errors="replace")
    result = classify_hessian(text, imag_thresh_cm1=imag_thresh_cm1)
    result["log"] = str(p)
    if not result["normal_termination"]:
        raise RuntimeError(f"{p}: Hessian job did not terminate normally")
    if not result["has_frequency_section"] and result["n_frequencies_parsed"] == 0:
        raise RuntimeError(f"{p}: no frequency table in Hessian log")
    if result["n_frequencies_parsed"] == 0:
        raise RuntimeError(f"{p}: frequency table header found but no numbers parsed")
    if any(not math.isfinite(x) for x in result["frequencies_cm1"]):
        raise RuntimeError(f"{p}: non-finite frequency")
    return result
