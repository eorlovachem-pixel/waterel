#!/usr/bin/env python3
from pathlib import Path
import os
import subprocess
import sys

KIT = Path(__file__).resolve().parents[1]
os.chdir(KIT)
print("KIT", KIT, flush=True)
r = subprocess.run([sys.executable, "-m", "pytest", "tests", "-q", "--tb=short"], cwd=KIT)
print("PYTEST", r.returncode, flush=True)
if r.returncode != 0:
    sys.exit(r.returncode)
r = subprocess.run([sys.executable, "workflow.py", "process-existing"], cwd=KIT)
print("PROCESS", r.returncode, flush=True)
sys.exit(r.returncode)
