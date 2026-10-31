# For Denis

This file is the public-repo counterpart of `MESSAGE_TO_DENIS.md`. Nothing in it changes `CASE_STATUS.csv` or sets `DENIS_ASSIGNED`.

## Current gate

- `K_06` remains `ON_HOLD`. The 100-start series has **not** been launched.
- Li_07–Li_14 and K_07–K_14 remain `NOT_ASSIGNED`.
- Na⁺ and Rb⁺ are **not** in this case table (Na verbally claimed; Rb already running elsewhere).
- F⁻/Cl⁻/Br⁻ are supported in code (`--chrg -1`) but have **no** packing files and **no** calculations.

## Existing data that *was* run (method development, 2-start pilot)

K⁺(H₂O)₆, GFN2-xTB/ALPB(water):

1. `frame_001` — E = −30.517135487228 Eh, genuine minimum, CoCoS group 0 (representative).
2. `frame_000` minus branch — E = −30.517100856916 Eh, confirmed minimum, same CoCoS group.
3. `frame_000` plus branch — E = −30.517067438374 Eh, confirmed minimum, same CoCoS group.

`frame_000` free-opt was a saddle (ν = −66.06 cm⁻¹) and is excluded from grouping.

## Decisions still needed (same as MESSAGE_TO_DENIS.md)

1. May this workstation take K⁺(H₂O)₆, 100 starts, packing SHA256 `06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5`?
2. Keep GWH `from_euler('zyz')` (SciPy-extrinsic) or switch to `'ZYZ'` (SciPy-intrinsic)? Source is unchanged until you say.
3. Which of Li_07–14 / K_07–14, if any, are assigned here?
4. CN cutoff for publication (provisional K cut 3.20 Å is marked `BLOCKING_UNCONFIRMED`).
5. Confirmation that the K₆ **pilot** geometries and figures may stay in the public `waterel` repository.

Until (1) and (5) are explicit, production chemistry stays blocked and unpublished co-author series are not pushed.
