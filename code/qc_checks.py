#!/usr/bin/env python3
"""Strict checks for XYZ, composition, charge, geometry and xTB logs.

This module never invents quantum-chemical numbers. It only parses files
that already exist and raises if they are incomplete or inconsistent.
"""
from __future__ import annotations

import hashlib
import math
import os
import re
import shutil
from collections import Counter
from pathlib import Path

import numpy as np

ENERGY_RE = re.compile(
    r"\bE\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*Eh"
)
ENERGY_COLON_RE = re.compile(
    r"energy\s*:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)", re.I
)
TOTAL_RE = re.compile(
    r"TOTAL ENERGY\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)", re.I
)
OPT_RE = re.compile(r"GEOMETRY OPTIMIZATION CONVERGED", re.I)
NORMAL_TERM_RE = re.compile(r"normal termination of xtb", re.I)
FORMAT_CRASH_RE = re.compile(r"Missing comma between descriptors", re.I)
HESS_HEADER_RE = re.compile(r"FREQUENCY CALCULATION|vibrational frequencies", re.I)

SUPPORTED_CATIONS = ("Li", "Na", "K", "Rb")
SUPPORTED_ANIONS = ("F", "Cl", "Br")
SUPPORTED_IONS = SUPPORTED_CATIONS + SUPPORTED_ANIONS
# Default argparse/cation charge. Use protocol_charge(ion) for the ion-specific value.
PROTOCOL_CHARGE = 1
PROTOCOL_CHARGE_CATION = 1
PROTOCOL_CHARGE_ANION = -1
# Radii taken only from packing files in this kit. Na/Rb/anion radii are not invented.
KNOWN_INITIAL_M_O_A = {"Li": 2.42, "K": 3.04}
SYSTEM_XTB = Path("/usr/bin/xtb")
SCIPY_FROM_EULER_DOCS = (
    "SciPy Rotation.from_euler: seq in {'X','Y','Z'} is intrinsic (body-fixed); "
    "seq in {'x','y','z'} is extrinsic (space-fixed). "
    "from_euler('ZYZ') is SciPy-intrinsic ZYZ; from_euler('zyz') is SciPy-extrinsic ZYZ. "
    "They are not the same rotation for the same angles. "
    "Docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_euler.html"
)


def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kit_root() -> Path:
    return Path(__file__).resolve().parents[1]


def standalone_xtb_candidates(root=None):
    root = Path(root) if root else kit_root()
    return [
        root / "vendors" / "xtb-dist" / "bin" / "xtb",
        root / "bin" / "xtb",
        Path.home() / "xtb-standalone" / "xtb-dist" / "bin" / "xtb",
    ]


def is_system_xtb(path) -> bool:
    try:
        p = Path(path).expanduser().resolve()
    except OSError:
        return False
    if SYSTEM_XTB.exists():
        return p == SYSTEM_XTB.resolve()
    return str(p) == "/usr/bin/xtb"


def resolve_working_xtb(explicit=None, kit=None, allow_system=False):
    """Return the Intel-static xTB path. Do not silently use /usr/bin/xtb for ANCOPT.

    Search order: explicit path (unless it is just the name 'xtb'), $XTB_BIN,
    kit vendors/bin, ~/xtb-standalone, then PATH. /usr/bin/xtb is skipped unless
    allow_system=True. Missing binary is a hard error; no energies are invented.
    """
    unspecified = explicit in (None, "", "xtb")
    ordered = []
    if not unspecified:
        ordered.append(Path(str(explicit)).expanduser())
    env_bin = os.environ.get("XTB_BIN")
    if env_bin:
        ordered.append(Path(env_bin).expanduser())
    ordered.extend(standalone_xtb_candidates(kit))
    which = shutil.which("xtb")
    if which:
        ordered.append(Path(which))

    seen = set()
    found_system = None
    for cand in ordered:
        p = cand
        if not p.exists():
            w = shutil.which(str(p))
            p = Path(w) if w else p
        if not p.exists() or not os.access(p, os.X_OK):
            continue
        p = p.resolve()
        if p in seen:
            continue
        seen.add(p)
        if is_system_xtb(p):
            found_system = p
            if not unspecified and not allow_system:
                raise RuntimeError(
                    "Refusing /usr/bin/xtb (Ubuntu libgfortran 16 crashes ANCOPT in 6.7.1). "
                    "Use the Intel-static binary via env_xtb.sh / ~/xtb-standalone. "
                    "Do not replace /usr/bin/xtb."
                )
            if allow_system:
                return p
            continue
        return p
    if found_system and not allow_system:
        raise RuntimeError(
            "Only /usr/bin/xtb is available. Ubuntu libgfortran 16 crashes ANCOPT. "
            "Install Intel-static 6.7.1 at ~/xtb-standalone or source env_xtb.sh. "
            "Do not replace /usr/bin/xtb."
        )
    raise RuntimeError(
        "xTB executable unavailable; install xTB first. No simulated energy will be generated."
    )


def read_xyz(path):
    p = Path(path)
    lines = p.read_text(encoding="utf-8").splitlines()
    i = 0
    frames = []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].strip())
        except ValueError as ex:
            raise ValueError(f"{p}: bad atom count at line {i+1}") from ex
        if n <= 0 or i + 2 + n > len(lines):
            raise ValueError(f"{p}: incomplete XYZ at {i+1}")
        comment = lines[i + 1]
        symbols = []
        xyz = []
        for line in lines[i + 2 : i + 2 + n]:
            toks = line.split()
            if len(toks) < 4:
                raise ValueError(f"{p}: invalid atom row: {line}")
            symbols.append(toks[0])
            xyz.append([float(k) for k in toks[1:4]])
        coords = np.array(xyz, float)
        if not np.isfinite(coords).all():
            raise ValueError(f"{p}: non-finite coordinates in a frame")
        frames.append({"symbols": symbols, "xyz": coords, "comment": comment})
        i += n + 2
    if not frames:
        raise ValueError(f"No structures in {p}")
    return frames


def formula_of(symbols):
    return dict(Counter(symbols))


def expected_formula(ion: str, n_water: int):
    return {ion: 1, "O": int(n_water), "H": 2 * int(n_water)}


def expected_natoms(n_water: int) -> int:
    return 1 + 3 * int(n_water)


def assert_supported_ion(ion: str):
    if ion not in SUPPORTED_IONS:
        raise ValueError(
            f"unsupported ion {ion!r}; this hydrate protocol supports "
            f"cations {SUPPORTED_CATIONS} and anions {SUPPORTED_ANIONS}"
        )
    return ion


def is_cation(ion: str) -> bool:
    return assert_supported_ion(ion) in SUPPORTED_CATIONS


def is_anion(ion: str) -> bool:
    return assert_supported_ion(ion) in SUPPORTED_ANIONS


def protocol_charge(ion: str) -> int:
    """Formal cluster charge: +1 for Li/Na/K/Rb, −1 for F/Cl/Br."""
    if is_cation(ion):
        return PROTOCOL_CHARGE_CATION
    return PROTOCOL_CHARGE_ANION


def xtb_chrg_flag(ion: str) -> str:
    return str(protocol_charge(ion))


def assert_cluster_composition(frame, ion: str, n_water: int, charge: int):
    assert_supported_ion(ion)
    wanted = expected_formula(ion, n_water)
    found = formula_of(frame["symbols"])
    if found != wanted:
        raise ValueError(f"wrong formula {found}, wanted {wanted}")
    n = len(frame["symbols"])
    if n != expected_natoms(n_water):
        raise ValueError(f"atom count {n}, wanted {expected_natoms(n_water)}")
    expected_charge = protocol_charge(ion)
    if int(charge) != expected_charge:
        raise ValueError(
            f"this hydrate protocol requires charge {expected_charge:+d} for {ion}, got {charge}"
        )
    if frame["symbols"].count(ion) != 1:
        raise ValueError(f"expected exactly one {ion}")
    return wanted


def energy_from_comment(comment):
    m = ENERGY_RE.search(comment)
    if m:
        return float(m.group(1))
    m = ENERGY_COLON_RE.search(comment)
    return float(m.group(1)) if m else None


def parse_xtb_log(text: str) -> dict:
    """Parse one xTB stdout/stderr log. Returns a dict; never fills energies if absent."""
    info = {
        "normal_termination": bool(NORMAL_TERM_RE.search(text)),
        "opt_converged": bool(OPT_RE.search(text)),
        "format_crash": bool(FORMAT_CRASH_RE.search(text)),
        "total_energy_Eh": None,
        "n_total_energy_prints": 0,
        "has_frequency_section": bool(HESS_HEADER_RE.search(text)),
        "version_line": None,
    }
    totals = [float(x) for x in TOTAL_RE.findall(text)]
    info["n_total_energy_prints"] = len(totals)
    if totals:
        info["total_energy_Eh"] = totals[-1]
    for line in text.splitlines():
        if "xtb version" in line.lower() or "* xtb version" in line:
            info["version_line"] = line.strip()
            break
    return info


def require_completed_opt(log_text: str, xyz_frame, expected_symbols, energy_match_atol=1e-6):
    """Geometry from xtbopt.xyz and energy from the same completed log."""
    if FORMAT_CRASH_RE.search(log_text):
        raise RuntimeError(
            "xTB ANCOPT crashed with Fortran format error "
            "(Missing comma between descriptors). This binary is not usable."
        )
    info = parse_xtb_log(log_text)
    if not info["normal_termination"]:
        raise RuntimeError("xTB log lacks 'normal termination of xtb'")
    if not info["opt_converged"]:
        raise RuntimeError("xTB log lacks 'GEOMETRY OPTIMIZATION CONVERGED'")
    if info["total_energy_Eh"] is None or not math.isfinite(info["total_energy_Eh"]):
        raise RuntimeError("no finite TOTAL ENERGY in the completed xTB log")
    if xyz_frame["symbols"] != expected_symbols:
        raise RuntimeError(
            f"optimized atom order {xyz_frame['symbols']} != {expected_symbols}"
        )
    comment_e = energy_from_comment(xyz_frame["comment"])
    if comment_e is not None and not math.isclose(
        comment_e, info["total_energy_Eh"], rel_tol=0.0, abs_tol=energy_match_atol
    ):
        raise RuntimeError(
            f"energy mismatch between xtbopt.xyz comment ({comment_e}) "
            f"and log TOTAL ENERGY ({info['total_energy_Eh']})"
        )
    return info["total_energy_Eh"], info


def validate_vecfle_bytes(path, ion, n_water, radius, expected_sha256=None):
    p = Path(path)
    raw = p.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if expected_sha256 and digest != expected_sha256:
        raise RuntimeError(f"{p}: SHA256 {digest} != {expected_sha256}")
    lines = [x.split() for x in raw.decode("utf-8").splitlines() if x.strip()]
    if len(lines) != int(n_water):
        raise RuntimeError(f"{p}: {len(lines)} sites, expected {n_water}")
    if any(row[0] != ion or len(row) != 5 for row in lines):
        raise RuntimeError(f"{p}: invalid ion label or column count")
    vec = np.array([[float(s) for s in row[1:4]] for row in lines], float)
    rads = np.array([float(row[4]) for row in lines], float)
    norms = np.linalg.norm(vec, axis=1)
    if np.max(np.abs(norms - 1.0)) > 1e-7:
        raise RuntimeError(f"{p}: non-unit direction vectors")
    if np.max(np.abs(rads - float(radius))) > 1e-10:
        raise RuntimeError(f"{p}: radius mismatch vs {radius}")
    return {
        "path": str(p),
        "sha256": digest,
        "n_water": int(n_water),
        "ion": ion,
        "radius_A": float(radius),
        "unit_vectors": True,
    }


def gwh_uses_intrinsic_zyz(gwh_source: Path) -> dict:
    """Audit GWH Euler strings. Does not rewrite the author source.

    Author commit 7411811 uses from_euler('zyz', ...). SciPy documentation
    (Rotation.from_euler) treats lowercase sequences as extrinsic and uppercase
    as intrinsic, so author 'zyz' is SciPy-extrinsic ZYZ, not SciPy-intrinsic
    'ZYZ'. That naming mismatch is reported; the GWH file is left unchanged
    until the author confirms which convention was intended.
    """
    text = Path(gwh_source).read_text(encoding="utf-8", errors="replace")
    has_zyz = "from_euler('zyz'" in text or 'from_euler("zyz"' in text
    has_ZYZ = "from_euler('ZYZ'" in text or 'from_euler("ZYZ"' in text
    has_zyx = "from_euler('zyx'" in text or 'from_euler("zyx"' in text
    has_ZYX = "from_euler('ZYX'" in text or 'from_euler("ZYX"' in text
    author_intact = bool(has_zyz and not has_zyx and not has_ZYZ)
    return {
        "path": str(gwh_source),
        "sha256": sha256_file(gwh_source),
        "has_zyz": has_zyz,
        "has_ZYZ_uppercase": has_ZYZ,
        "has_zyx": has_zyx,
        "has_ZYX_uppercase": has_ZYX,
        "zyz_ZYZ_discrepancy": bool(has_zyz and has_ZYZ),
        "author_method_intact": author_intact,
        # Back-compat alias: True iff the author lowercase 'zyz' string is present
        # and unmixed. This is NOT "SciPy intrinsic". See scipy_interprets_author_zyz_as.
        "intrinsic_zyz": author_intact,
        "scipy_interprets_author_zyz_as": "extrinsic" if has_zyz else None,
        "scipy_docs_uppercase_ZYZ_is": "intrinsic",
        "scipy_docs_lowercase_zyz_is": "extrinsic",
        "scipy_note": SCIPY_FROM_EULER_DOCS,
        "method_policy": (
            "Leave the author source unchanged. Do not rewrite 'zyz' to 'ZYZ' "
            "or the reverse without an explicit author decision. "
            "SciPy docs: lowercase = extrinsic, uppercase = intrinsic."
        ),
        "open_question": (
            "Author commit 7411811 writes from_euler('zyz'). SciPy treats that as "
            "extrinsic ZYZ. Confirm with Denis whether the intended rigid-water "
            "rotation is SciPy-extrinsic 'zyz' or SciPy-intrinsic 'ZYZ'."
        ),
    }
