# Confirmed minima post-processing — [K(H2O)6]+ 2-structure pilot

Generated UTC: 2026-10-30T13:46:04.730059+00:00
Numbers below are re-read from completed xTB files. No fabricated QC values.

## Geometries kept

Only structures with verified free optimization, full 3N frequencies, and **zero** significant imaginary modes (|ν| > 5 cm⁻¹).
Translational/rotational near-zero modes (|ν| ≤ 5 cm⁻¹) are not counted as imaginary.

### frame_001 (`frame_001`)

- origin: free optimization of GWH structure 0001
- formula: {'K': 1, 'O': 6, 'H': 12} (19 atoms)
- E = -30.51713548722800 Eh
- ΔE = 0.00000000 kcal/mol vs lowest of this set
- opt_converged / normal_termination: True / True
- n_frequencies = 57 (expected 57)
- n_significant_imaginary = 0
- n_near_zero = 6
- lowest vibrational = 31.12 cm⁻¹
- K–O min/mean/max = 2.6146 / 2.6378 / 2.6752 Å
- XYZ: `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/xtb_runs/frame_001/xtbopt.xyz`

### frame_000_minus (`frame_000/mode_follow/minus`)

- origin: frame_000 imaginary-mode displacement -0.15 A, reoptimized
- formula: {'K': 1, 'O': 6, 'H': 12} (19 atoms)
- E = -30.51710085691600 Eh
- ΔE = 0.02173085 kcal/mol vs lowest of this set
- opt_converged / normal_termination: True / True
- n_frequencies = 57 (expected 57)
- n_significant_imaginary = 0
- n_near_zero = 6
- lowest vibrational = 23.62 cm⁻¹
- K–O min/mean/max = 2.6295 / 2.6471 / 2.7028 Å
- XYZ: `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/xtb_runs/frame_000/mode_follow/minus_0.15A/opt/xtbopt.xyz`

### frame_000_plus (`frame_000/mode_follow/plus`)

- origin: frame_000 imaginary-mode displacement +0.15 A, reoptimized
- formula: {'K': 1, 'O': 6, 'H': 12} (19 atoms)
- E = -30.51706743837400 Eh
- ΔE = 0.04270130 kcal/mol vs lowest of this set
- opt_converged / normal_termination: True / True
- n_frequencies = 57 (expected 57)
- n_significant_imaginary = 0
- n_near_zero = 6
- lowest vibrational = 22.12 cm⁻¹
- K–O min/mean/max = 2.6465 / 2.6603 / 2.7059 Å
- XYZ: `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/xtb_runs/frame_000/mode_follow/plus_0.15A/opt/xtbopt.xyz`

Lowest energy in this set: **-30.51713548722800 Eh** (`frame_001`).

## CoCoS (author snapshot, fixed protocol)

- code: `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/code/cocos_uploaded_snapshot.py` SHA256 `4c538fb83812cef33a81db25de9b973d258e781c75703d959e7e2966ad4e7cbc`
- wrapper: `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/code/seeded_cocos.py` with `random.seed` and `np.random.seed`
- production settings: seed **20261030**, NumSamples **50**, `--DoNotUseGlobalOpt`, `--DoNotSaveDoubles`
- CritFreq = 83.4 cm⁻¹, NLRMSDFreq = 0.035 cm⁻¹ (CoCoS defaults)
- input order: energy-sorted, as in `pipeline.dedup`
- unique representatives at production settings: **1**
- added indices: [0]
- removed: [{'structure': 1, 'copy_of': 0}, {'structure': 2, 'copy_of': 0}]
- CoCoS placed the three geometries into 1 equivalence class(es) (seed 20261030, NumSamples 50).

## Seed robustness (same XYZ, same flags)

| samples | seed | n_unique | added indices | removed |
|---|---|---|---|---|
| 50 | 20261030 | 1 | [0] | 1→0; 2→0 |
| 50 | 101 | 1 | [0] | 1→0; 2→0 |
| 50 | 202 | 1 | [0] | 1→0; 2→0 |
| 50 | 303 | 1 | [0] | 1→0; 2→0 |
| 50 | 7 | 1 | [0] | 1→0; 2→0 |
| 100 | 20261030 | 1 | [0] | 1→0; 2→0 |
| 100 | 101 | 1 | [0] | 1→0; 2→0 |
| 100 | 202 | 1 | [0] | 1→0; 2→0 |
| 100 | 303 | 1 | [0] | 1→0; 2→0 |
| 100 | 7 | 1 | [0] | 1→0; 2→0 |

For the 10 runs above, the number of groups and the retained indices did not change. This is not a proof of algorithmic completeness.

## Files

- `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/postprocess_confirmed_minima/K06_pilot_confirmed_minima.xyz`
- `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/postprocess_confirmed_minima/confirmed_minima_table.csv`
- `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/postprocess_confirmed_minima/confirmed_minima_geometries.png` / `.svg`
- `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/postprocess_confirmed_minima/confirmed_minima_table.png` / `.svg`
- `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/postprocess_confirmed_minima/cocos_primary`
- `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30/results/K_06/pilot_2/postprocess_confirmed_minima/cocos_sensitivity`

These three geometries are a 2-structure pilot follow-up, not a 100-structure search.
