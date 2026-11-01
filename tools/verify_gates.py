#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT))
import pipeline as P
import workflow as W


def main() -> None:
    u = P.read_xyz(ROOT / "results/K/6/final/unique_representatives.xyz")
    c = P.read_xyz(ROOT / "results/K/6/final/confirmed_minima.xyz")
    print("unique_frames", len(u))
    print("confirmed_frames", len(c))
    print("unique_comment", u[0]["comment"][:120])

    report = {"unique_n": len(u), "confirmed_n": len(c), "gates": {}}
    for case in ("Cl_04", "K_06", "Br_04"):
        p = subprocess.run(
            [sys.executable, str(ROOT / "workflow.py"), "gwh", "--case", case, "--nstruct", "20"],
            cwd=ROOT, capture_output=True, text=True,
        )
        text = (p.stdout or "") + "\n" + (p.stderr or "")
        report["gates"][case] = {
            "returncode": p.returncode,
            "refused": p.returncode != 0 and ("NOT_ASSIGNED" in text or "ON_HOLD" in text),
            "snippet": text.strip()[-400:],
        }
        print(case, "rc", p.returncode, "refused", report["gates"][case]["refused"])

    oauth_cands = [
        Path.home() / "hydration-research" / "google_oauth_client.json",
        ROOT / "google_oauth_client.json",
    ]
    report["oauth_clients"] = [str(p) for p in oauth_cands if p.is_file()]
    up = subprocess.run(
        [sys.executable, str(ROOT / "workflow.py"), "upload-drive"],
        cwd=ROOT, capture_output=True, text=True,
    )
    report["upload_drive"] = {
        "returncode": up.returncode,
        "out": ((up.stdout or "") + (up.stderr or ""))[-800:],
    }
    dest = ROOT / "reports" / "PRODUCTION_GATE.json"
    dest.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("wrote", dest)


if __name__ == "__main__":
    main()
