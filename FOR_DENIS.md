# For Denis

This file is the public-repo counterpart of `MESSAGE_TO_DENIS.md`. Nothing in it sets `DENIS_ASSIGNED`.

## Current gate

- `K_06` remains `ON_HOLD`. The 100-start series has **not** been launched.
- K_04, K_05, K_07–K_14, Li_07–Li_14, Cl_04–Cl_14, Br_04–Br_14 are `NOT_ASSIGNED`.
- Na⁺ and Rb⁺ are **not** in this case table.
- Cl⁻/Br⁻ packings exist (Sloane directions, mendeleev VI radii). **No** GWH `--Opt`, **no** xTB free opt.
- Launch refuses `ON_HOLD` and `NOT_ASSIGNED` unless `DENIS_ASSIGNED=yes`.

## Existing data that *was* run (method development, 2-start pilot)

K⁺(H₂O)₆, GFN2-xTB/ALPB(water):

1. `frame_001` — E = −30.517135487228 Eh, genuine minimum, CoCoS group 0 (**unique representative**).
2. `frame_000` minus branch — E = −30.517100856916 Eh, confirmed minimum, same CoCoS group.
3. `frame_000` plus branch — E = −30.517067438374 Eh, confirmed minimum, same CoCoS group.

`unique_representatives.xyz` now holds only `frame_001`. All three geometries remain in `confirmed_minima.xyz`. `frame_000` free-opt was a saddle (ν = −66.06 cm⁻¹) and is excluded from grouping.

## Decisions still needed

See `MESSAGE_TO_DENIS.md` (six questions: K_06 100-start; remaining K; Drive GWH reuse; Cl/Br assignment; Euler `'zyz'`; CN/publication).

Public `waterel` contains the K₆ **pilot** plus code, tests, and Cl/Br/K_04–K_05 packings. Vendor trees and raw `cases/` stay local.
