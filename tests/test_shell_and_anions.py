#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import pipeline as P
import qc_checks
import shell_descriptors as SD


def _xyz(symbols, coords, comment="E = -1.0 Eh"):
    return {"symbols": list(symbols), "xyz": np.array(coords, float), "comment": comment}


def test_anion_charge_accepted_cation_charge_rejected_for_f():
    # F- + 1 water: F, O, H, H — not a production packing, only a composition unit test.
    fr = _xyz(["F", "O", "H", "H"], [[0, 0, 0], [2.5, 0, 0], [3.0, 0.8, 0], [3.0, -0.8, 0]])
    qc_checks.assert_cluster_composition(fr, "F", 1, charge=-1)
    with pytest.raises(ValueError, match="charge"):
        qc_checks.assert_cluster_composition(fr, "F", 1, charge=1)


def test_prepare_accepts_anion_charge(tmp_path):
    fr = _xyz(
        ["F", "O", "H", "H"],
        [[0, 0, 0], [2.5, 0, 0], [3.0, 0.8, 0], [3.0, -0.8, 0]],
        "E = -1.000000000000 Eh",
    )
    xyz = tmp_path / "in.xyz"
    P.write_xyz(xyz, [fr])
    csvp = tmp_path / "e.csv"
    csvp.write_text("id,batch,structure,E_final_Eh,nfev,success\nx,b,0,-1.000000000000,,True\n", encoding="utf-8")
    P.prepare(type("A", (), {
        "input": str(xyz), "energies_csv": str(csvp), "output": str(tmp_path / "o"),
        "central": "F", "charge": -1, "allow_existing": False,
    })())
    man = (tmp_path / "o" / "manifest.json").read_text(encoding="utf-8")
    assert '"charge": -1' in man


def test_cn_cut_is_marked_blocking():
    assert SD.CN_CUT_STATUS == "BLOCKING_UNCONFIRMED"
    assert "K" in SD.CN_CUT_PROVISIONAL_A


def test_k6_pilot_descriptors_full_shell():
    xyz = ROOT / "results" / "K_06" / "pilot_2" / "xtb_runs" / "frame_001" / "xtbopt.xyz"
    if not xyz.is_file():
        pytest.skip("K6 pilot XYZ not present")
    fr = P.read_xyz(xyz)[0]
    desc = SD.cluster_descriptors(fr, "K", n_water=6, xyz_path=xyz)
    assert desc["water_connectivity"] == "nominal_2H_per_O"
    assert desc["CN"] == 6
    assert desc["full_first_shell"] is True
    assert desc["proton_transfer_candidate"] is False
    assert desc["R_min_A"] < desc["R_mean_A"] < desc["R_max_A"]
    assert desc["orientation_cos_mean"] is not None
    assert desc["xyz_sha256"]
    assert len(desc["xyz_sha256"]) == 64


def test_anion_xh_contact_and_no_false_proton_transfer():
    # F at origin, water with H pointing toward F but still bound to O.
    fr = _xyz(
        ["F", "O", "H", "H"],
        [
            [0.0, 0.0, 0.0],
            [2.60, 0.0, 0.0],
            [1.70, 0.0, 0.0],
            [3.10, 0.85, 0.0],
        ],
    )
    desc = SD.cluster_descriptors(fr, "F", n_water=1)
    assert desc["n_XH_O_contacts"] >= 1
    assert desc["proton_transfer_candidate"] is False


def test_matplotlib_does_not_draw_ion_oxygen_as_covalent():
    import mol_matplotlib as MM
    symbols = ["K", "O", "H"]
    xyz = np.array([[0, 0, 0], [2.64, 0, 0], [3.2, 0.8, 0]], float)
    bonds = MM.covalent_bonds(symbols, xyz)
    assert all({symbols[i], symbols[j]} != {"K", "O"} for i, j, _ in bonds)
    assert any({symbols[i], symbols[j]} == {"O", "H"} for i, j, _ in bonds)
    coords = MM.coordination_contacts(symbols, xyz, "K", 3.20)
    assert len(coords) == 1


def test_case_status_cl_br_not_assigned_no_na_rb_f():
    import csv
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    ions = {r["ion"] for r in rows}
    assert {"Cl", "Br"} <= ions
    assert not {"F", "Na", "Rb"} & ions
    anion = [r for r in rows if r["ion"] in ("Cl", "Br")]
    assert {r["case"] for r in anion} == {
        f"{ion}_{n:02d}" for ion in ("Cl", "Br") for n in range(4, 15)
    }
    for row in anion:
        assert row["assignment_status"] == "NOT_ASSIGNED"
        assert int(row["charge"]) == -1
        assert int(row["n_start_target"]) == 100
