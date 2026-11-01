#!/usr/bin/env python3
"""Pipeline confirmed-minima filter and overwrite guards. No fabricated QC numbers."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import freq_parser
import mode_follow
import pipeline as P
import qc_checks


def _xyz(symbols, coords, comment="E = -1.0 Eh"):
    return {"symbols": list(symbols), "xyz": np.array(coords, float), "comment": comment}


def test_prepare_rejects_wrong_charge_and_ion(tmp_path):
    fr = _xyz(["K"] + ["O"] * 6 + ["H"] * 12, [[0, 0, 0]] * 19, "E = -1.000000000000 Eh")
    xyz = tmp_path / "in.xyz"
    P.write_xyz(xyz, [fr])
    csvp = tmp_path / "e.csv"
    csvp.write_text("id,batch,structure,E_final_Eh,nfev,success\nx,b,0,-1.000000000000,,True\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported ion"):
        P.prepare(type("A", (), {"input": str(xyz), "energies_csv": str(csvp), "output": str(tmp_path / "o"),
                                 "central": "Cs", "charge": 1, "allow_existing": False})())
    with pytest.raises(ValueError, match="charge"):
        P.prepare(type("A", (), {"input": str(xyz), "energies_csv": str(csvp), "output": str(tmp_path / "o2"),
                                 "central": "K", "charge": 0, "allow_existing": False})())


def test_mode_follow_scale_and_nu_cutoff():
    cart = np.array([[0.2, 0, 0], [0, 0, 0], [0, 0, 0]], float)
    scaled, meta = mode_follow.scale_mode_max_atom(cart, 0.15)
    assert meta["scaled_max_atom_A"] == pytest.approx(0.15)
    assert scaled[0, 0] == pytest.approx(0.15)
    rows = [{"mode_index": i + 1, "frequency_cm1": f} for i, f in enumerate([-0.0, 0.0, -66.06, 31.12])]
    cls = mode_follow.classify_modes(rows, 5.0)
    assert cls["n_significant_imaginary"] == 1
    assert cls["significant_imaginary"][0]["frequency_cm1"] == pytest.approx(-66.06)


def test_collect_confirmed_excludes_saddle_keeps_refine_branches(tmp_path):
    out = tmp_path / "xtb_runs"
    (out / "frame_000").mkdir(parents=True)
    (out / "frame_001").mkdir(parents=True)
    k6 = ["K"] + ["O"] * 6 + ["H"] * 12
    coords = [[float(i), 0.0, 0.0] for i in range(19)]
    e_saddle = -30.517082221236
    e_min = -30.517135487228
    e_plus = -30.517067438374
    for folder, energy, comment in (
        (out / "frame_000", e_saddle, f" energy: {e_saddle:.12f} gnorm: 0.0005 xtb: 6.7.1 (edcfbbe)"),
        (out / "frame_001", e_min, f" energy: {e_min:.12f} gnorm: 0.0002 xtb: 6.7.1 (edcfbbe)"),
    ):
        P.write_xyz(folder / "input.xyz", [_xyz(k6, coords, f"E = {energy:.14f} Eh")])
        P.write_xyz(folder / "xtbopt.xyz", [_xyz(k6, coords, comment)])
        (folder / "xtb.out").write_text(
            f"GEOMETRY OPTIMIZATION CONVERGED\nTOTAL ENERGY              {energy:.12f} Eh\nnormal termination of xtb\n",
            encoding="utf-8",
        )
    saddle_log = """FREQUENCY CALCULATION
frequencies (cm**-1)
eigval :       -0.00    -0.00    -0.00     0.00     0.00     0.00
eigval :      -66.06    28.17    33.34
normal termination of xtb
"""
    min_log = """FREQUENCY CALCULATION
frequencies (cm**-1)
eigval :       -0.00    -0.00    -0.00     0.00     0.00     0.00
eigval :       31.12    41.18    42.98
normal termination of xtb
"""
    (out / "frame_000" / "hess").mkdir()
    (out / "frame_000" / "hess" / "xtb_hess.out").write_text(saddle_log, encoding="utf-8")
    (out / "frame_001" / "hess").mkdir()
    (out / "frame_001" / "hess" / "xtb_hess.out").write_text(min_log, encoding="utf-8")
    plus = out / "frame_000" / "mode_follow" / "plus_0.15A" / "opt"
    plus.mkdir(parents=True)
    P.write_xyz(plus / "xtbopt.xyz", [_xyz(k6, coords, f" energy: {e_plus:.12f} gnorm: 0.0002 xtb: 6.7.1 (edcfbbe)")])
    (plus / "xtb.out").write_text(
        f"GEOMETRY OPTIMIZATION CONVERGED\nTOTAL ENERGY              {e_plus:.12f} Eh\nnormal termination of xtb\n",
        encoding="utf-8",
    )
    (plus.parent / "hess").mkdir()
    (plus.parent / "hess" / "xtb_hess.out").write_text(min_log, encoding="utf-8")
    (out / "manifest.json").write_text(json.dumps({
        "charge": 1, "central": "K", "method": "GFN2-xTB", "solvent": "ALPB(water)",
        "entries": [
            {"id": "frame_000", "constrained_energy_Eh": -30.51, "source_run_id": "gwh_opt_str0000"},
            {"id": "frame_001", "constrained_energy_Eh": -30.51, "source_run_id": "gwh_opt_str0001"},
        ],
    }), encoding="utf-8")
    confirmed, notes, expected, manifest = P.collect_confirmed_minima(out, imag_thresh=5.0)
    ids = [c["id"] for c in confirmed]
    assert "frame_000" not in ids
    assert "frame_001" in ids
    assert "frame_000_plus" in ids
    assert "frame_000_minus" not in ids  # incomplete branch excluded
    assert all(c["n_imag"] == 0 for c in confirmed)
    statuses = {n["id"]: n["status"] for n in notes}
    assert statuses["frame_000"] == "SADDLE"
    assert statuses["frame_001"] == "CONFIRMED_FREE_OPT"


def test_collect_nested_refine_round_two(tmp_path):
    out = tmp_path / "xtb_runs"
    (out / "frame_000").mkdir(parents=True)
    k6 = ["K"] + ["O"] * 6 + ["H"] * 12
    coords = [[float(i), 0.0, 0.0] for i in range(19)]
    e0 = -30.517082221236
    e1 = -30.517090000000
    e2 = -30.517135487228
    saddle_log = """FREQUENCY CALCULATION
frequencies (cm**-1)
eigval :       -0.00    -0.00    -0.00     0.00     0.00     0.00
eigval :      -66.06    28.17    33.34
normal termination of xtb
"""
    min_log = """FREQUENCY CALCULATION
frequencies (cm**-1)
eigval :       -0.00    -0.00    -0.00     0.00     0.00     0.00
eigval :       31.12    41.18    42.98
normal termination of xtb
"""

    def _write_opt(folder, energy, log):
        folder.mkdir(parents=True, exist_ok=True)
        P.write_xyz(folder / "xtbopt.xyz", [_xyz(k6, coords, f" energy: {energy:.12f} gnorm: 0.0002 xtb: 6.7.1 (edcfbbe)")])
        (folder / "xtb.out").write_text(
            f"GEOMETRY OPTIMIZATION CONVERGED\nTOTAL ENERGY              {energy:.12f} Eh\nnormal termination of xtb\n",
            encoding="utf-8",
        )
        (folder.parent / "hess" if folder.name == "opt" else folder / "hess").mkdir(exist_ok=True)
        hdir = folder.parent / "hess" if folder.name == "opt" else folder / "hess"
        hdir.mkdir(exist_ok=True)
        (hdir / "xtb_hess.out").write_text(log, encoding="utf-8")

    P.write_xyz(out / "frame_000" / "input.xyz", [_xyz(k6, coords, f"E = {e0:.14f} Eh")])
    _write_opt(out / "frame_000", e0, saddle_log)
    plus = out / "frame_000" / "mode_follow" / "plus_0.15A"
    _write_opt(plus / "opt", e1, saddle_log)
    plus2 = plus / "mode_follow" / "plus_0.15A"
    _write_opt(plus2 / "opt", e2, min_log)
    (out / "manifest.json").write_text(json.dumps({
        "charge": 1, "central": "K",
        "entries": [{"id": "frame_000", "constrained_energy_Eh": -30.51, "source_run_id": "gwh_opt_str0000"}],
    }), encoding="utf-8")
    confirmed, notes, expected, manifest = P.collect_confirmed_minima(out, imag_thresh=5.0, max_iters=3)
    ids = [c["id"] for c in confirmed]
    assert "frame_000" not in ids
    assert "frame_000_plus" not in ids
    assert "frame_000_plus_plus" in ids
    limited, notes1, _, _ = P.collect_confirmed_minima(out, imag_thresh=5.0, max_iters=1)
    assert limited == []
    assert any(n.get("status") == "MAX_ITERS_STILL_SADDLE" for n in notes1)


def test_resolve_xtb_uses_explicit_file_not_system(tmp_path, monkeypatch):
    fake = tmp_path / "xtb"
    fake.write_text("#!/bin/sh\n", encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setattr(
        qc_checks, "standalone_xtb_candidates",
        lambda root=None: [fake],
    )
    got = qc_checks.resolve_working_xtb(str(fake), kit=tmp_path)
    assert got == fake.resolve()
    monkeypatch.delenv("XTB_BIN", raising=False)
    monkeypatch.setenv("PATH", "/usr/bin")
    resolved = qc_checks.resolve_working_xtb(None, kit=tmp_path)
    assert resolved == fake.resolve()
    assert not qc_checks.is_system_xtb(resolved)
    monkeypatch.setattr(
        qc_checks, "standalone_xtb_candidates",
        lambda root=None: [tmp_path / "nope" / "xtb"],
    )
    monkeypatch.setenv("PATH", str(tmp_path / "empty_path"))
    with pytest.raises(RuntimeError, match="unavailable"):
        qc_checks.resolve_working_xtb(None, kit=tmp_path)


def test_resolve_xtb_raises_when_nothing_usable(tmp_path, monkeypatch):
    monkeypatch.delenv("XTB_BIN", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path / "empty_path"))
    monkeypatch.setattr(
        qc_checks, "standalone_xtb_candidates",
        lambda root=None: [tmp_path / "nope" / "xtb"],
    )
    with pytest.raises(RuntimeError, match="unavailable"):
        qc_checks.resolve_working_xtb(None, kit=tmp_path)


def test_k06_on_hold_is_not_bypassed_by_allow_full_batch(monkeypatch):
    sys.path.insert(0, str(ROOT))
    import workflow as W
    monkeypatch.delenv("DENIS_ASSIGNED", raising=False)
    args = type("A", (), {"case": "K_06"})()
    row = W.select(args)
    assert row["assignment_status"] == "ON_HOLD"
    with pytest.raises(RuntimeError, match="ON_HOLD"):
        W.require_launch_allowed(args, row)
    monkeypatch.setenv("DENIS_ASSIGNED", "yes")
    W.require_launch_allowed(args, row)


def test_refuse_overwrite_helper(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("a", encoding="utf-8")
    with pytest.raises(RuntimeError, match="Refusing overwrite"):
        P.refuse_overwrite(p, allow=False)
    assert P.refuse_overwrite(p, allow=True) == p


def test_k06_pilot_copy_does_not_include_saddle(tmp_path):
    src = ROOT / "results" / "K_06" / "pilot_2" / "xtb_runs"
    if not (src / "manifest.json").is_file():
        pytest.skip("K6 pilot files not present in this tree")
    dest = tmp_path / "xtb_runs"
    shutil.copytree(src, dest)
    original_opt = (src / "frame_000" / "xtbopt.xyz").read_bytes()
    confirmed, notes, expected, manifest = P.collect_confirmed_minima(dest, imag_thresh=5.0)
    ids = [c["id"] for c in confirmed]
    assert "frame_001" in ids
    assert "frame_000" not in ids
    assert "frame_000_plus" in ids
    assert "frame_000_minus" in ids
    assert all(c["n_imag"] == 0 for c in confirmed)
    assert (src / "frame_000" / "xtbopt.xyz").read_bytes() == original_opt
    for c in confirmed:
        qc_checks.assert_cluster_composition(c["frame"], "K", 6, 1)
        assert freq_parser.is_confirmed_minimum(
            freq_parser.require_completed_hessian(c["hess_log"], 5.0), nu_min_cm1=-5.0
        )
