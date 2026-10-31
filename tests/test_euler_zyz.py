#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import qc_checks


def test_gwh_source_is_author_intrinsic_zyz():
    src = ROOT / "vendors" / "Glorious-Water-Hole" / "gloriouswaterhole.py"
    if not src.is_file():
        import pytest
        pytest.skip("GWH vendor tree not cloned in this copy")
    info = qc_checks.gwh_uses_intrinsic_zyz(src)
    assert info["has_zyz"]
    assert not info["has_zyx"]
    assert not info["has_ZYZ_uppercase"]
    assert not info["zyz_ZYZ_discrepancy"]
    assert info["author_method_intact"]
    # Alias kept for callers; it means author 'zyz' is present, not SciPy-intrinsic.
    assert info["intrinsic_zyz"]
    assert info["scipy_interprets_author_zyz_as"] == "extrinsic"
    assert info["scipy_docs_uppercase_ZYZ_is"] == "intrinsic"
    assert "Leave the author source unchanged" in info["method_policy"]


def test_gwh_constrained_stage_is_rigid_water_euler_not_xtb_ancopt():
    src = ROOT / "vendors" / "Glorious-Water-Hole" / "gloriouswaterhole.py"
    if not src.is_file():
        import pytest
        pytest.skip("GWH vendor tree not cloned in this copy")
    text = src.read_text(encoding="utf-8", errors="replace")
    assert "spo.minimize" in text or "scipy.optimize" in text
    assert "rOH=0.958" in text
    assert "aHOH=104.4776" in text
    assert "def RunXTB" in text
    assert "--charge" in text
    # GWH constrained stage calls xtb as a single-point engine; ANCOPT --opt is later.
    cmd_lines = [ln for ln in text.splitlines() if "xtb " in ln and "Command" in ln]
    assert cmd_lines
    assert all("--opt" not in ln for ln in cmd_lines)


def test_zyz_vs_ZYZ_are_not_the_same_scipy_convention(tmp_path):
    fake = tmp_path / "gwh.py"
    fake.write_text("M = spst.Rotation.from_euler('ZYZ', [phi, theta, chi], degrees=True)\n", encoding="utf-8")
    info = qc_checks.gwh_uses_intrinsic_zyz(fake)
    assert info["has_ZYZ_uppercase"]
    assert not info["has_zyz"]
    assert not info["author_method_intact"]
    assert not info["intrinsic_zyz"]
    both = tmp_path / "both.py"
    both.write_text(
        "from_euler('zyz', a)\nfrom_euler('ZYZ', b)\n",
        encoding="utf-8",
    )
    disc = qc_checks.gwh_uses_intrinsic_zyz(both)
    assert disc["zyz_ZYZ_discrepancy"]
    assert not disc["author_method_intact"]


def test_scipy_from_euler_lowercase_zyz_differs_from_uppercase_ZYZ():
    import numpy as np
    from scipy.spatial.transform import Rotation as R
    angles = [12.0, 34.0, 56.0]
    m_lower = R.from_euler("zyz", angles, degrees=True).as_matrix()
    m_upper = R.from_euler("ZYZ", angles, degrees=True).as_matrix()
    assert not np.allclose(m_lower, m_upper)
    # Documented relation: extrinsic ABC(a,b,c) == intrinsic CBA(c,b,a)
    m_int_rev = R.from_euler("ZYZ", list(reversed(angles)), degrees=True).as_matrix()
    assert np.allclose(m_lower, m_int_rev)


def test_workflow_does_not_dedup_before_free_opt():
    wf = (ROOT / "workflow.py").read_text(encoding="utf-8")
    assert "No CoCoS filtering or energy-window preselection allowed before free optimization" in wf
    assert "--alpb water --gfn 2" in wf
    assert "--DoNotCenter" in wf
    assert "--Charge" in wf
