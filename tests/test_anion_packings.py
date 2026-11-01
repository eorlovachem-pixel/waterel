#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import qc_checks


def _rows():
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        return {r["case"]: r for r in csv.DictReader(fh)}


def test_cl_br_known_radii_match_mendeleev_sums():
    assert qc_checks.KNOWN_INITIAL_M_O_A["Cl"] == pytest.approx(3.21)
    assert qc_checks.KNOWN_INITIAL_M_O_A["Br"] == pytest.approx(3.36)
    prov = json.loads((ROOT / "inputs" / "ANION_PACKING_PROVENANCE.json").read_text(encoding="utf-8"))
    assert prov["pm_to_A"].startswith("ionic_radius_A")
    assert prov["Cl_minus_VI"]["coordination"] == "VI"
    assert prov["Cl_minus_VI"]["charge"] == -1
    assert prov["Br_minus_VI"]["coordination"] == "VI"
    assert prov["O2_minus_VI"]["coordination"] == "VI"
    assert prov["packing_R_A"]["Cl"] == 3.21
    assert prov["packing_R_A"]["Br"] == 3.36
    assert prov["packing_R_A"]["K_04_K_05"] == 3.04
    cl_cns = {r["coordination"] for r in prov["full_ionic_radius_tables"]["Cl"] if r["charge"] == -1}
    br_cns = {r["coordination"] for r in prov["full_ionic_radius_tables"]["Br"] if r["charge"] == -1}
    assert cl_cns == {"VI"}
    assert br_cns == {"VI"}


def test_cl_br_packings_match_case_status_and_unit_vectors():
    rows = _rows()
    for ion, radius in (("Cl", 3.21), ("Br", 3.36)):
        for n in range(4, 15):
            case = f"{ion}_{n:02d}"
            row = rows[case]
            info = qc_checks.validate_vecfle_bytes(
                ROOT / "inputs" / case / "packing.vecfle",
                ion=ion,
                n_water=n,
                radius=radius,
                expected_sha256=row["input_sha256"],
            )
            assert info["n_water"] == n
            assert int(row["charge"]) == qc_checks.protocol_charge(ion) == -1


def test_cl_and_br_share_sloane_directions():
    cl = np.array(
        [[float(x) for x in ln.split()[1:4]]
         for ln in (ROOT / "inputs" / "Cl_06" / "packing.vecfle").read_text().splitlines() if ln.strip()]
    )
    br = np.array(
        [[float(x) for x in ln.split()[1:4]]
         for ln in (ROOT / "inputs" / "Br_06" / "packing.vecfle").read_text().splitlines() if ln.strip()]
    )
    k = np.array(
        [[float(x) for x in ln.split()[1:4]]
         for ln in (ROOT / "inputs" / "K_06" / "packing.vecfle").read_text().splitlines() if ln.strip()]
    )
    assert np.allclose(cl, br)
    prov = json.loads((ROOT / "inputs" / "ANION_PACKING_PROVENANCE.json").read_text(encoding="utf-8"))
    rmsd6 = prov["sloane_by_n"]["6"]["rmsd_vs_author_K_hungarian"]
    # Author K_06 matches the Sloane n=6 spherical code (octahedral-like).
    assert rmsd6 == pytest.approx(0.0, abs=1e-8)
    assert np.allclose(cl, k)


def test_n4_is_tetrahedron_n5_is_tbp():
    cl4 = np.array(
        [[float(x) for x in ln.split()[1:4]]
         for ln in (ROOT / "inputs" / "Cl_04" / "packing.vecfle").read_text().splitlines() if ln.strip()]
    )
    assert cl4.shape == (4, 3)
    dots = cl4 @ cl4.T
    off = dots[np.triu_indices(4, 1)]
    assert np.allclose(off, -1.0 / 3.0, atol=1e-6)
    cl5 = np.array(
        [[float(x) for x in ln.split()[1:4]]
         for ln in (ROOT / "inputs" / "Cl_05" / "packing.vecfle").read_text().splitlines() if ln.strip()]
    )
    assert cl5.shape == (5, 3)
    assert np.max(np.abs(np.linalg.norm(cl5, axis=1) - 1.0)) < 1e-8
    off5 = (cl5 @ cl5.T)[np.triu_indices(5, 1)]
    # Sloane n=5 is a numerical trigonal bipyramid (one nearly antipodal pair).
    assert np.min(off5) < -0.95


def test_k04_k05_use_author_k_radius_and_are_not_assigned():
    rows = _rows()
    for case, n in (("K_04", 4), ("K_05", 5)):
        row = rows[case]
        assert row["assignment_status"] == "NOT_ASSIGNED"
        assert int(row["charge"]) == 1
        qc_checks.validate_vecfle_bytes(
            ROOT / "inputs" / case / "packing.vecfle",
            ion="K",
            n_water=n,
            radius=3.04,
            expected_sha256=row["input_sha256"],
        )


def test_k_li_packings_untouched():
    import hashlib
    k06 = ROOT / "inputs" / "K_06" / "packing.vecfle"
    digest = hashlib.sha256(k06.read_bytes()).hexdigest()
    assert digest == "06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5"
    rows = _rows()
    assert rows["K_06"]["assignment_status"] == "ON_HOLD"
    assert rows["Li_07"]["assignment_status"] == "NOT_ASSIGNED"
    assert "Na_06" not in rows
    assert "F_04" not in rows
