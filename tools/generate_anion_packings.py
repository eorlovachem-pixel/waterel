#!/usr/bin/env python3
"""Build Cl/Br n=4–14 and K_04/K_05 packing.vecfle from Sloane dim-3 codes + mendeleev radii.

Does not launch GWH or xTB. Does not mark cases ASSIGNED.
Does not overwrite existing author K_06–K_14 or Li packings.
"""
from __future__ import annotations

import csv
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from mendeleev import element
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "inputs"
SLOANE_DIR = INPUTS / "sloane_dim3"
CSV_PATH = ROOT / "CASE_STATUS.csv"
JSON_PATH = ROOT / "CASE_STATUS.json"
PROV_PATH = INPUTS / "ANION_PACKING_PROVENANCE.json"
N_RANGE = range(4, 15)
PM_TO_A = 0.01
SLOANE_URL = "https://neilsloane.com/packings/dim3/pack.3.{n}.txt"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def ionic_radius_table(symbol: str) -> list[dict]:
    el = element(symbol)
    rows = []
    for ir in el.ionic_radii or []:
        rows.append({
            "symbol": symbol,
            "name": el.name,
            "charge": int(ir.charge),
            "coordination": str(ir.coordination),
            "spin": ir.spin or "",
            "ionic_radius_pm": float(ir.ionic_radius) if ir.ionic_radius is not None else None,
            "crystal_radius_pm": float(ir.crystal_radius) if ir.crystal_radius is not None else None,
            "ionic_radius_A": (
                None if ir.ionic_radius is None else round(float(ir.ionic_radius) * PM_TO_A, 6)
            ),
            "econf": ir.econf,
            "origin": (ir.origin or "").strip().rstrip(","),
            "most_reliable": bool(ir.most_reliable),
        })
    return rows


def pick_radius(rows: list[dict], charge: int, coordination: str) -> dict:
    hits = [
        r for r in rows
        if r["charge"] == charge and r["coordination"] == coordination
        and r["ionic_radius_pm"] is not None
    ]
    if len(hits) != 1:
        raise RuntimeError(f"expected one radius for charge={charge} CN={coordination}, got {hits}")
    return hits[0]


def download_sloane(n: int) -> Path:
    SLOANE_DIR.mkdir(parents=True, exist_ok=True)
    dest = SLOANE_DIR / f"pack.3.{n}.txt"
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    url = SLOANE_URL.format(n=n)
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    dest.write_bytes(data)
    return dest


def parse_sloane(path: Path, n: int) -> np.ndarray:
    nums = [float(x) for x in path.read_text(encoding="utf-8").split() if x.strip()]
    if len(nums) != 3 * n:
        raise RuntimeError(f"{path}: expected {3*n} numbers, got {len(nums)}")
    vec = np.array(nums, float).reshape(n, 3)
    norms = np.linalg.norm(vec, axis=1)
    if np.min(norms) < 1e-8:
        raise RuntimeError(f"{path}: zero-length vector")
    unit = vec / norms[:, None]
    if np.max(np.abs(np.linalg.norm(unit, axis=1) - 1.0)) > 1e-8:
        raise RuntimeError(f"{path}: failed to unit-normalize")
    return unit


def read_author_unit_vectors(ion: str, n: int) -> np.ndarray | None:
    path = INPUTS / f"{ion}_{n:02d}" / "packing.vecfle"
    if not path.is_file():
        return None
    lines = [ln.split() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    vec = np.array([[float(x) for x in row[1:4]] for row in lines], float)
    norms = np.linalg.norm(vec, axis=1)
    return vec / norms[:, None]


def rmsd_after_assignment(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD of two unit-vector sets after Hungarian matching (no rotation search)."""
    cost = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(cost)
    paired = b[ci]
    return float(np.sqrt(np.mean(np.sum((a - paired) ** 2, axis=1))))


def format_vecfle(ion: str, unit: np.ndarray, radius: float) -> str:
    lines = []
    for ux, uy, uz in unit:
        lines.append(f"{ion:<5s}{ux:16.10f}{uy:16.10f}{uz:16.10f}{radius:8.2f}")
    return "\n".join(lines) + "\n"


def load_csv_rows() -> list[dict]:
    with CSV_PATH.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def write_case_status(rows: list[dict]) -> None:
    fieldnames = [
        "case", "ion", "charge", "n_water", "initial_M_O_A", "n_start_target",
        "input_sha256", "start_status", "claim_status", "assignment_status",
    ]
    with CSV_PATH.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in fieldnames})
    JSON_PATH.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def new_row(case: str, ion: str, charge: int, n: int, radius: float, digest: str) -> dict:
    return {
        "case": case,
        "ion": ion,
        "charge": str(charge),
        "n_water": str(n),
        "initial_M_O_A": f"{radius:.2f}",
        "n_start_target": "100",
        "input_sha256": digest,
        "start_status": "NO_GWH_RESULTS__ASSIGNMENT_PENDING",
        "claim_status": "TECHNICAL_PREP_ONLY__CONFIRM_WITH_DENIS",
        "assignment_status": "NOT_ASSIGNED",
    }


def write_packing(ion: str, n: int, unit: np.ndarray, radius: float) -> dict:
    case = f"{ion}_{n:02d}"
    dest = INPUTS / case / "packing.vecfle"
    dest.parent.mkdir(parents=True, exist_ok=True)
    text = format_vecfle(ion, unit, radius)
    data = text.encode("utf-8")
    dest.write_bytes(data)
    digest = sha256_bytes(data)
    return {
        "case": case,
        "ion": ion,
        "n_water": n,
        "initial_M_O_A": radius,
        "path": str(dest.relative_to(ROOT).as_posix()),
        "sha256": digest,
    }


def main() -> None:
    cl_tab = ionic_radius_table("Cl")
    br_tab = ionic_radius_table("Br")
    o_tab = ionic_radius_table("O")
    cl_vi = pick_radius(cl_tab, charge=-1, coordination="VI")
    br_vi = pick_radius(br_tab, charge=-1, coordination="VI")
    o_vi = pick_radius(o_tab, charge=-2, coordination="VI")
    r_cl = round(cl_vi["ionic_radius_A"] + o_vi["ionic_radius_A"], 2)
    r_br = round(br_vi["ionic_radius_A"] + o_vi["ionic_radius_A"], 2)
    r_k = 3.04
    if r_cl != 3.21 or r_br != 3.36:
        raise RuntimeError(f"unexpected sums: Cl {r_cl} Br {r_br}")

    sloane = {}
    sloane_meta = {}
    for n in N_RANGE:
        path = download_sloane(n)
        unit = parse_sloane(path, n)
        sloane[n] = unit
        kvec = read_author_unit_vectors("K", n)
        meta = {
            "file": str(path.relative_to(ROOT).as_posix()),
            "url": SLOANE_URL.format(n=n),
            "sha256": sha256_bytes(path.read_bytes()),
            "n": n,
        }
        if kvec is not None:
            meta["rmsd_vs_author_K_hungarian"] = rmsd_after_assignment(unit, kvec)
        sloane_meta[str(n)] = meta

    packing_rows = []
    for n in N_RANGE:
        packing_rows.append(write_packing("Cl", n, sloane[n], r_cl))
        packing_rows.append(write_packing("Br", n, sloane[n], r_br))
    packing_rows.append(write_packing("K", 4, sloane[4], r_k))
    packing_rows.append(write_packing("K", 5, sloane[5], r_k))

    existing = load_csv_rows()
    keep = [
        r for r in existing
        if r["case"] not in {p["case"] for p in packing_rows}
    ]
    charge_of = {"Cl": -1, "Br": -1, "K": 1}
    added = [
        new_row(p["case"], p["ion"], charge_of[p["ion"]], p["n_water"], p["initial_M_O_A"], p["sha256"])
        for p in packing_rows
    ]
    # Keep Li then K_06–14 then new K_04/05 then Cl then Br for readability.
    added_sorted = sorted(
        added,
        key=lambda r: ({"K": 0, "Cl": 1, "Br": 2}[r["ion"]], int(r["n_water"])),
    )
    write_case_status(keep + added_sorted)

    provenance = {
        "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
        "mendeleev": "1.3.0",
        "docs": "https://mendeleev.readthedocs.io/en/stable/",
        "sloane": "https://neilsloane.com/packings/dim3/",
        "pm_to_A": "ionic_radius_A = ionic_radius_pm / 100",
        "selection_rule": (
            "X- uses the only tabulated anion entry (charge -1, coordination VI). "
            "O2- VI is the oxygen radius. Packing R = r(X-) + r(O2-). "
            "K_04/K_05 use the author K series radius 3.04 A, not a Shannon sum. "
            "Existing K_06–K_14 and Li packings were not rewritten."
        ),
        "not_shannon_label_caveat": (
            "mendeleev origin string for Cl-/Br- VI ionic_radius is "
            "Pauling's (1960) crystal radius. Numeric pairs 181/167 pm (Cl) and "
            "196/182 pm (Br) match Shannon 1976 CN=6 ionic/crystal radii. "
            "most_reliable is False for both halide anions."
        ),
        "Cl_minus_VI": cl_vi,
        "Br_minus_VI": br_vi,
        "O2_minus_VI": o_vi,
        "packing_R_A": {"Cl": r_cl, "Br": r_br, "K_04_K_05": r_k},
        "full_ionic_radius_tables": {"Cl": cl_tab, "Br": br_tab, "O": o_tab},
        "sloane_by_n": sloane_meta,
        "cases": packing_rows,
        "assignment": (
            "New CASE_STATUS rows are NOT_ASSIGNED. DENIS_ASSIGNED was not set. "
            "Production 100-start GWH/xTB is refused until Denis assigns."
        ),
    }
    PROV_PATH.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Cl R={r_cl} A  Br R={r_br} A  K_04/K_05 R={r_k} A")
    print(f"wrote {len(packing_rows)} packing files")
    print(f"provenance: {PROV_PATH}")


if __name__ == "__main__":
    main()
