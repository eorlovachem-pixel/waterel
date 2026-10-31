# Methodology audit — K+(H2O)6 protocol

Recorded 2026-10-30 against the author sources cloned on this machine. This is an audit, not a claim that K6 production results exist.

## Charge, formula, method

- `CASE_STATUS.csv` sets `K_06` charge **+1**, ion **K**, **n_water=6**.
- Expected composition is **19 atoms**: 1 K + 6 O + 12 H.
- `workflow.py gwh` passes `--Charge` from that table and `--AddMethod " --alpb water --gfn 2 "`.
- Free relaxation in `pipeline.relax` is `xtb input.xyz --gfn 2 --chrg <manifest charge> --alpb water --opt tight`.
- `inputs/K_06/packing.vecfle` SHA256 `06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5` matches `CASE_STATUS.csv`. Six unit vectors, radius 3.04 Å. Not reconstructed.

## GWH constrained stage (rigid water, fixed oxygens)

Author `RunOptimizationOfTheCluster` does **not** freeze atoms inside xTB ANCOPT. It:

1. Builds each water from experimental geometry (`rOH=0.958` Å, `aHOH=104.4776°`, NIST CCCBDB defaults in that function).
2. Places oxygen along the vecfle direction at the listed M–O distance (fixed for this stage).
3. Optimizes **only Euler angles** with SciPy `minimize` (Powell).
4. Calls xTB as a **single-point** energy engine (`xtb … --charge N <AddMethod>`), not `--opt`.

That is the rigid-water / fixed-oxygen first stage. Because ANCOPT is not used here, the Ubuntu `/usr/bin/xtb` format crash does not appear during GWH. It appears at **free** optimization.

`--DoNotCenter` is passed by `workflow.py`, so the ion stays the geometric centre rather than the SPC/E charge centre of water.

## Euler convention (author `zyz`; SciPy lowercase = extrinsic)

Local clone was `bfc474c0ea56c1d0c9cfae6285deac9aff065dff` (`from_euler('zyx', …)`).

Author commit **`7411811374b6422a15f23d6d86a47f2145f8ac47`** (2026-10-30, Denis Tikhonov) changes that one line to `from_euler('zyz', …)`.

SciPy `Rotation.from_euler` documentation (current): **uppercase** `{'X','Y','Z'}` is **intrinsic** (body-fixed); **lowercase** `{'x','y','z'}` is **extrinsic** (space-fixed). Therefore author `'zyz'` is SciPy-**extrinsic** ZYZ, and `'ZYZ'` would be SciPy-**intrinsic** ZYZ. They are not the same rotation for the same angles.

The kit **does not rewrite** `'zyz'` → `'ZYZ'`. Confirm with Denis which convention was intended. An earlier audit note that “lowercase is intrinsic” was incorrect relative to the SciPy docs.

The vendor tree was fast-forwarded to that commit. Existing `cases/Li_07` (100 structures) and `pilots/Li_07` were produced with **old `zyx`** and must not be mixed with later `zyz` results or overwritten.

## Free optimization: all constraints off, no pre-filter

- `collect` keeps every GWH structure that has `normal termination of xtb`, the right formula, and `E = … Eh`.
- There is no CoCoS call and no energy window before `prepare`/`relax`.
- `relax` runs unconstrained `--opt tight` on every prepared frame in its own `xtb_runs/frame_XXX/` directory.
- CoCoS runs only in `dedup`, after free optimization.

## Energy and geometry from the same completed job

`pipeline.get_result` reads `xtbopt.xyz` and `xtb.out` from the **same frame folder**. Tightened checks reject the Fortran format crash, require `GEOMETRY OPTIMIZATION CONVERGED` plus `normal termination of xtb`, and cross-check the XYZ comment energy against the last `TOTAL ENERGY` in that log.

## Hessian / imaginary frequencies

The original orchestrator did **not** run `--hess`. START_HERE requires that imaginary-frequency structures are not counted as minima. `code/freq_parser.py` and `workflow.py hess` implement that as a **separate** stage. It is not auto-run on 100 structures. Treat this as a protocol addition that still needs coordinator confirmation before production use.

## Separate directories, no overwrite

- GWH refuses if `gwh/` already has files.
- `source/` collect refuses if outputs exist.
- `prepare` refuses a non-empty `xtb_runs/` unless `--allow-existing` and byte-identical inputs.
- Pilots live under `pilots/`, production under `cases/`.
- Saved Li_07 calculations were left in place.

## xTB binaries

- `/usr/bin/xtb` 6.7.1 (unknown), Ubuntu package, linked to **libgfortran5 16**. ANCOPT writes `(1x,"("f7.2"%)")` at `optimizer.f90:852` of tag v6.7.1 (`26b28010e805f7d1aeeef39813feb473e69cc4be`). That format is rejected at runtime. Do not replace this file.
- Working binary: Intel-static `~/xtb-standalone/xtb-dist/bin/xtb`, reports **6.7.1 (edcfbbe)**, SHA256 `debf27a9e0fa4bfb5ca75aafe4b90d8211f08ec2f4a482f375a4987212eaa12a`. Same version string as the Na reference XYZ comments. Enable with `source env_xtb.sh`.

## CoCoS / MDAtomLocalEnvironment

- GitLab CoCoS `47f6680da2eee90d9b18536fd94d2bc8f6d05aad` file SHA256 equals the uploaded snapshot `code/cocos_uploaded_snapshot.py` (`4c538fb83812cef33a81db25de9b973d258e781c75703d959e7e2966ad4e7cbc`). Use the GitLab clone; the snapshot is a bit-identical backup.
- MDAtom `3c43de70abef1f64fb31a2fe2693e4e73b88083c`. Pipeline correctly calls `-e` (element), not `-n` (index). `--DoNotSaveDoubles` is `store_false` in CoCoS (passing the flag turns saving doubles off).

## What was already computed on this Linux tree (do not delete)

- `pilots/Li_07`: 2 GWH structures, free opt, CoCoS, MDAtom. GWH SHA `bfc474c` (old Euler), system xtb recorded in provenance.
- `cases/Li_07`: 100 GWH constrained structures, same old Euler, exit 0. Not a K6 result. Not overwritten.

## Not done here

- No 100-structure K_06 series.
- No fabricated GFN2 energies.
- Assignment of K6 still needs Denis, as in `ON_HOLD_AND_SCOPE.md`.
