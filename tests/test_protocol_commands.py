#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_run_one_script_does_not_overwrite_and_has_preflight_gate():
    text = (ROOT / "RUN_ONE_APPROVED_CASE.sh").read_text(encoding="utf-8")
    assert "nstruct 100" in text
    assert "environment_check.py" in text or "xtb_preflight.py" in text
    assert "confirm" in text.lower() or "ASSIGNED" in text or "Denis" in text
    assert "workflow.py hess" in text
    assert "workflow.py refine" in text
    assert "workflow.py dedup" in text
    assert "ALLOW_FULL_BATCH" in text
    assert "ON_HOLD" in text
    assert "DENIS_ASSIGNED" in text
    assert "source ./env_xtb.sh" in text
    assert "XTB_BIN" in text


def test_denis_assigned_is_not_hardcoded_yes():
    sh = (ROOT / "RUN_ONE_APPROVED_CASE.sh").read_text(encoding="utf-8")
    wf = (ROOT / "workflow.py").read_text(encoding="utf-8")
    assert "DENIS_ASSIGNED" in sh
    assert "require_launch_allowed" in wf
    assert "os.environ.get('DENIS_ASSIGNED')" in wf


def test_original_pipeline_refuses_empty_xtb():
    text = (ROOT / "code" / "pipeline.py").read_text(encoding="utf-8")
    assert "No simulated energy will be generated" in text
    assert "--alpb" in text and "water" in text
    assert "--gfn" in text
    assert "allow_incomplete" in text
    assert "def refine" in text
    assert "confirmed_minima_for_cocos" in text
    assert "ν ≥" in text or "nu_min" in text or "is_confirmed_minimum" in text
    assert "DEFAULT_MAX_REFINE_ITERS" in text
    assert "DEFAULT_COCOS_SEED = 20261030" in text
    assert "resolve_working_xtb" in text or "require_xtb" in text
    wf = (ROOT / "workflow.py").read_text(encoding="utf-8")
    assert "process-existing" in wf
    assert "visualize" in wf
    assert "protocol_charge" in wf
