#!/usr/bin/env python3
"""Parser tests on labelled fixtures. These are not quantum-chemical results."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import freq_parser
import qc_checks

FIXTURE_XYZ = """19
fixture only E = -12.345678901234 Eh
K   0.0 0.0 0.0
O   2.0 0.0 0.0
H   2.5 0.7 0.0
H   2.5 -0.7 0.0
O  -2.0 0.0 0.0
H  -2.5 0.7 0.0
H  -2.5 -0.7 0.0
O   0.0 2.0 0.0
H   0.7 2.5 0.0
H  -0.7 2.5 0.0
O   0.0 -2.0 0.0
H   0.7 -2.5 0.0
H  -0.7 -2.5 0.0
O   0.0 0.0 2.0
H   0.7 0.0 2.5
H  -0.7 0.0 2.5
O   0.0 0.0 -2.0
H   0.7 0.0 -2.5
H  -0.7 0.0 -2.5
"""

OPT_LOG_OK = """
           * xtb version 6.7.1 (edcfbbe)
           GEOMETRY OPTIMIZATION CONVERGED
           TOTAL ENERGY              -12.345678901234 Eh
           TOTAL ENERGY              -12.345678901234 Eh
           normal termination of xtb
"""

OPT_LOG_CRASH = """
           * total energy  :     -7.3318435 Eh
At line 852 of file ../src/optimizer.f90 (unit = 6, file = 'stdout')
Fortran runtime error: Missing comma between descriptors
"""

HESS_LOG = """
     FREQUENCY CALCULATION
     frequencies (cm**-1)
          1          2          3          4
      12.10      45.20     1601.3     3650.1
     reduced masses (amu)
           normal termination of xtb
"""

HESS_IMAG = """
     FREQUENCY CALCULATION
     frequencies (cm**-1)
          1          2          3
     -85.10      40.00     1600.0
           normal termination of xtb
"""


def test_supported_ions_cations_and_anions():
    assert qc_checks.SUPPORTED_CATIONS == ("Li", "Na", "K", "Rb")
    assert qc_checks.SUPPORTED_ANIONS == ("F", "Cl", "Br")
    assert qc_checks.SUPPORTED_IONS == ("Li", "Na", "K", "Rb", "F", "Cl", "Br")
    for ion in qc_checks.SUPPORTED_CATIONS:
        assert qc_checks.protocol_charge(ion) == 1
        wanted = qc_checks.expected_formula(ion, 6)
        assert wanted == {ion: 1, "O": 6, "H": 12}
        qc_checks.assert_supported_ion(ion)
    for ion in qc_checks.SUPPORTED_ANIONS:
        assert qc_checks.protocol_charge(ion) == -1
        assert qc_checks.xtb_chrg_flag(ion) == "-1"
        qc_checks.assert_supported_ion(ion)
    with pytest.raises(ValueError, match="unsupported ion"):
        qc_checks.assert_supported_ion("Cs")


def test_k6_formula_and_charge_rule():
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "k6.xyz"
        p.write_text(FIXTURE_XYZ, encoding="utf-8")
        fr = qc_checks.read_xyz(p)[0]
        got = qc_checks.assert_cluster_composition(fr, "K", 6, charge=1)
        assert got == {"K": 1, "O": 6, "H": 12}
        assert qc_checks.energy_from_comment(fr["comment"]) == pytest.approx(-12.345678901234)
        with pytest.raises(ValueError):
            qc_checks.assert_cluster_composition(fr, "K", 6, charge=0)


def test_opt_log_and_geometry_must_come_from_same_run():
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "xtbopt.xyz"
        p.write_text(FIXTURE_XYZ, encoding="utf-8")
        fr = qc_checks.read_xyz(p)[0]
        e, info = qc_checks.require_completed_opt(OPT_LOG_OK, fr, fr["symbols"])
        assert e == pytest.approx(-12.345678901234)
        assert info["opt_converged"]
        with pytest.raises(RuntimeError, match="Fortran format"):
            qc_checks.require_completed_opt(OPT_LOG_CRASH, fr, fr["symbols"])


def test_frequency_parser_minimum_vs_saddle():
    ok = freq_parser.classify_hessian(HESS_LOG)
    assert ok["is_genuine_minimum_candidate"]
    assert ok["n_frequencies_parsed"] == 4
    bad = freq_parser.classify_hessian(HESS_IMAG)
    assert not bad["is_genuine_minimum_candidate"]
    assert bad["classification"] == "NOT_A_MINIMUM"
    assert bad["imaginary_cm1"][0] == pytest.approx(-85.10)
    assert freq_parser.is_confirmed_minimum(ok, nu_min_cm1=-5.0)
    assert not freq_parser.is_confirmed_minimum(bad, nu_min_cm1=-5.0)
    near = freq_parser.classify_hessian(
        """
     FREQUENCY CALCULATION
     frequencies (cm**-1)
         -0.00      12.0     1600.0
           normal termination of xtb
"""
    )
    assert freq_parser.is_confirmed_minimum(near, nu_min_cm1=-5.0)


def test_na_reference_is_19_atoms_but_not_used_as_k():
    ref = ROOT.parents[1] / "reference" / "Na6_12_minima_reference.xyz"
    if not ref.is_file():
        ref = ROOT.parent.parent / "reference" / "Na6_12_minima_reference.xyz"
    alt = Path("/home/antiz/hydration-research/Na6_12_minima_reference.xyz")
    for candidate in (
        ROOT / "all-12-minima-na-centered.xyz",
        ROOT.parents[1] / "reference" / "Na6_12_minima_reference.xyz",
        Path("/mnt/c/Users/antiz/OneDrive/Документы/Cursor_K6_research_handoff/Cursor_research/reference/Na6_12_minima_reference.xyz"),
    ):
        if candidate.is_file():
            frames = qc_checks.read_xyz(candidate)
            assert len(frames) == 12
            for fr in frames:
                assert qc_checks.formula_of(fr["symbols"]) == {"Na": 1, "O": 6, "H": 12}
                assert "K" not in fr["symbols"]
            return
    pytest.skip("Na reference XYZ not mounted in this tree")
