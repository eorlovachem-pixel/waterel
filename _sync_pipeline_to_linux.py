#!/usr/bin/env python3
"""Backward-compatible entry point. Implementation: tools/sync_to_linux.py."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).resolve().parent / "tools" / "sync_to_linux.py"), run_name="__main__")
