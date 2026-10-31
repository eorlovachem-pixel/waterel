#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import qc_checks

K06_SHA = "06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5"


def test_k06_vecfle_matches_author_bytes_and_geometry():
    info = qc_checks.validate_vecfle_bytes(
        ROOT / "inputs" / "K_06" / "packing.vecfle",
        ion="K",
        n_water=6,
        radius=3.04,
        expected_sha256=K06_SHA,
    )
    assert info["n_water"] == 6
    assert qc_checks.expected_natoms(6) == 19
    assert qc_checks.expected_formula("K", 6) == {"K": 1, "O": 6, "H": 12}


def test_case_status_k06_charge_and_target():
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        rows = {r["case"]: r for r in csv.DictReader(fh)}
    k = rows["K_06"]
    assert k["ion"] == "K"
    assert int(k["charge"]) == 1
    assert int(k["n_water"]) == 6
    assert float(k["initial_M_O_A"]) == pytest.approx(3.04)
    assert int(k["n_start_target"]) == 100
    assert k["input_sha256"] == K06_SHA
    assert k["assignment_status"] == "ON_HOLD"
    ions = {r["ion"] for r in rows.values()}
    assert "Na" not in ions
    assert "Rb" not in ions
    assert ions <= {"Li", "K"}
    for name, row in rows.items():
        if name == "K_06":
            continue
        assert row["assignment_status"] == "NOT_ASSIGNED", name


def test_li_packings_match_case_status():
    import csv
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["ion"] == "Li"]
    assert rows
    for row in rows:
        qc_checks.validate_vecfle_bytes(
            ROOT / "inputs" / row["case"] / "packing.vecfle",
            ion="Li",
            n_water=int(row["n_water"]),
            radius=float(row["initial_M_O_A"]),
            expected_sha256=row["input_sha256"],
        )
        assert float(row["initial_M_O_A"]) == pytest.approx(qc_checks.KNOWN_INITIAL_M_O_A["Li"])


def test_case_status_json_matches_csv_and_does_not_invent_assignments():
    import json
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        csv_rows = list(csv.DictReader(fh))
    json_rows = json.loads((ROOT / "CASE_STATUS.json").read_text(encoding="utf-8"))
    assert [r["case"] for r in csv_rows] == [r["case"] for r in json_rows]
    for c, j in zip(csv_rows, json_rows):
        assert c["ion"] == j["ion"]
        assert c["input_sha256"] == j["input_sha256"]
        assert c["assignment_status"] == j["assignment_status"]
        assert j["assignment_status"] in ("ON_HOLD", "NOT_ASSIGNED")
    assert next(r for r in csv_rows if r["case"] == "K_06")["assignment_status"] == "ON_HOLD"


def test_na_rb_supported_in_code_but_not_added_without_assignment():
    for ion in ("Na", "Rb"):
        qc_checks.assert_supported_ion(ion)
        assert qc_checks.expected_formula(ion, 6) == {ion: 1, "O": 6, "H": 12}
    assert "Na" not in qc_checks.KNOWN_INITIAL_M_O_A
    assert "Rb" not in qc_checks.KNOWN_INITIAL_M_O_A


def test_k06_directions_are_octahedral_like_unit_vectors():
    p = ROOT / "inputs" / "K_06" / "packing.vecfle"
    rows = [line.split() for line in p.read_text().splitlines() if line.strip()]
    vec = np.array([[float(x) for x in r[1:4]] for r in rows])
    assert vec.shape == (6, 3)
    assert np.max(np.abs(np.linalg.norm(vec, axis=1) - 1)) < 1e-8
    # Opposite poles along z are present in the author packing.
    z = vec[:, 2]
    assert np.any(np.isclose(z, 1.0))
    assert np.any(np.isclose(z, -1.0))
