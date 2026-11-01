#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_gwh_cli_audit_json_if_present():
    p = ROOT / "validation" / "gwh_cli" / "GWH_CLI_AUDIT.json"
    if not p.is_file():
        pytest.skip("GWH CLI audit not generated in this tree")
    audit = json.loads(p.read_text(encoding="utf-8"))
    assert audit["from_euler_zyz"] is True
    assert audit["from_euler_ZYZ"] is False
    assert audit["help_has_Charge"] is True
    assert audit["opt_series_launched"] is False
    assert audit["hundred_start_launched"] is False
    assert audit["Cl_04_nstruct_ok"] is True
    assert audit["Br_04_nstruct_ok"] is True
