#!/usr/bin/env python3
"""Honest per-stoichiometry counts. Zeros stay zeros."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    with (ROOT / "CASE_STATUS.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        if r["ion"] not in ("K", "Cl", "Br"):
            continue
        case = r["case"]
        n = int(r["n_water"])
        gwh = 2 if case == "K_06" else 0
        confirmed = 3 if case == "K_06" else 0
        cocos = 1 if case == "K_06" else 0
        unique = 1 if case == "K_06" else 0
        xtb_ok = 2 if case == "K_06" else 0  # two free-opt parents; plus two mode-follow branches
        out.append({
            "case": case,
            "ion": r["ion"],
            "n_water": n,
            "charge": r["charge"],
            "assignment_status": r["assignment_status"],
            "n_gwh_starts": gwh,
            "n_xtb_ok": xtb_ok,
            "n_confirmed": confirmed,
            "n_cocos_groups": cocos,
            "n_unique": unique,
            "note": "K_06 numbers are the 2-start pilot, not a 100-start series" if case == "K_06" else "NOT_ASSIGNED; no production QC",
        })
    dest = ROOT / "reports" / "SERIES_COUNTS.csv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("wrote", dest, "rows", len(out))


if __name__ == "__main__":
    main()
