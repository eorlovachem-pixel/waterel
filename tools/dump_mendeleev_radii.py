#!/usr/bin/env python3
"""Print Shannon-style ionic radii from mendeleev for Cl, Br, O."""
from __future__ import annotations

import json
from mendeleev import element

rows = []
for sym in ("Cl", "Br", "O"):
    el = element(sym)
    print(f"==== {sym} {el.name} mendeleev={el} ====")
    print("atomic_radius_pm", el.atomic_radius)
    for ir in el.ionic_radii or []:
        rec = {
            "symbol": sym,
            "charge": ir.charge,
            "coordination": str(ir.coordination),
            "spin": ir.spin,
            "ionic_radius_pm": ir.ionic_radius,
            "crystal_radius_pm": ir.crystal_radius,
            "econf": ir.econf,
            "origin": ir.origin,
            "most_reliable": ir.most_reliable,
        }
        rows.append(rec)
        print(rec)
    print()

print(json.dumps(rows, indent=2, ensure_ascii=False))
