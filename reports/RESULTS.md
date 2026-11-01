# Results (draft)

Do not treat this file as a completed 33-stoichiometry survey. Only the K₆ pilot exists.

## K⁺(H₂O)₆ pilot (2 GWH starts)

Three Hessian-confirmed minima, one CoCoS group. Unique representative: `frame_001`.

| id | E / Eh | ΔE / kcal mol⁻¹ | ν_min / cm⁻¹ | CN (3.20 Å, unconfirmed) |
|---|---|---|---|---|
| frame_001 | −30.517135487228 | 0.000 | 31.12 | 6 |
| frame_000_minus | −30.517100856916 | 0.022 | 23.62 | 6 |
| frame_000_plus | −30.517067438374 | 0.043 | 22.12 | 6 |

`unique_representatives.xyz` contains **only** `frame_001`. All three geometries are in `confirmed_minima.xyz` / `all_group_members.xyz`.

MDAtomLocalEnvironment was **not** run on this pilot. Python descriptors are in `results/K/6/final/descriptors.csv`.

## Cl⁻ and Br⁻

Packings exist. No GWH `--Opt`, no free xTB, no frequencies. Series plots are not generated.
