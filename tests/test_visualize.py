#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import pipeline as P
import visualize as V


def test_copy_xyz_verbatim_preserves_bytes(tmp_path):
    src = ROOT / "results" / "K_06" / "pilot_2" / "xtb_runs" / "frame_001" / "xtbopt.xyz"
    if not src.is_file():
        pytest.skip("K6 pilot XYZ not present")
    dest = tmp_path / "original.xyz"
    digest = V.copy_xyz_verbatim(src, dest)
    assert dest.read_bytes() == src.read_bytes()
    assert digest == hashlib.sha256(src.read_bytes()).hexdigest()
    fr0 = P.read_xyz(src)[0]
    fr1 = P.read_xyz(dest)[0]
    import numpy as np
    assert np.allclose(fr0["xyz"], fr1["xyz"])
    assert fr0["symbols"] == fr1["symbols"]


def test_aligned_copy_does_not_rewrite_source(tmp_path):
    src = ROOT / "results" / "K_06" / "pilot_2" / "xtb_runs" / "frame_001" / "xtbopt.xyz"
    if not src.is_file():
        pytest.skip("K6 pilot XYZ not present")
    before = src.read_bytes()
    entry = {"id": "frame_001", "ion": "K", "n_water": 6, "xyz_path": str(src),
             "energy_Eh": -30.517135487228, "cocos_group": 0, "origin": "test"}
    scene, aligned, _ = V.write_scene(entry, tmp_path / "scene")
    assert src.read_bytes() == before
    assert (tmp_path / "scene" / "original.xyz").read_bytes() == before
    assert aligned.is_file()
    assert scene["original_sha256"] == hashlib.sha256(before).hexdigest()
    assert scene["camera"]["original_xyz_untouched"] is True


def test_burette_status_does_not_claim_macos_app():
    st = V._burette_status()
    assert st["burette_app_available"] is False
    assert "macOS" in st["burette_app_reason"]


def test_workflow_visualize_does_not_call_launch_gate():
    text = (ROOT / "workflow.py").read_text(encoding="utf-8")
    # process-existing must not require DENIS_ASSIGNED / require_launch_allowed
    assert "def process_existing" in text
    block = text.split("def process_existing")[1].split("def visualize_cmd")[0]
    assert "require_launch_allowed" not in block
    assert "GWH" in block or "xTB" in block or "pilot" in block.lower() or "export" in block
