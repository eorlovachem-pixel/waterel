# Cleanup report (2026-11-01)

Final audit of `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30` before publishing `waterel`.
No quantum-chemical jobs were launched.

## Backup (outside the published Git tree)

- Timestamp UTC: `2026-10-31T23:18:31Z`
- Full tree copy: `/home/antiz/hydration-research/_backups/Hydration_unassigned_2026-10-30_2026-11-01`
- Local archive of unused helpers: `/home/antiz/hydration-research/_local_archive/cleanup_2026-11-01`

## Method

- Recursively hashed files (SHA-256).
- Searched imports and tests for each candidate.
- Byte-identical copies inside `results/`, `figures/`, and `inputs/` are **kept**: they are GWH/xTB/Hessian/CoCoS/canonical-mirror/Burette provenance, not unused drafts.
- Unique XYZ, `packing.vecfle`, Hessians, logs, checkpoints, and original kit files were not deleted.

## Archived (moved, not destroyed)

| Path | SHA-256 | Reason |
|---|---|---|
| `tools/finish_git.py` | `6a515c4d2561b28dab505916e148419f593f26ba1e812cfb1635c6d5beea11a0` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |
| `tools/git_waterel.py` | `6e71896cf64e796f7efda837bbd63839d92a34ff0544ca8983b2d0d56ef90cd4` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |
| `tools/wsl_go.py` | `b29a0a820e9e83946e4ec23b9a65cf30d707484e5b44dd9378c4c145067b125e` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |
| `tools/run_all.py` | `c52e96d3e98b2f30d9d3650369f5ab6afe48891caf5de970e20c5b4dfb922cea` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |
| `tools/apply_cleanup.py` | `cb788e6148585701306d3758343f6b4472058fd2817ede7f28057402ec553858` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |
| `tools/audit_cleanup.py` | `870d066aff5cc7e25f41cea439798ebba060eaaf75ab75b5b4899a2d611d0b09` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |
| `_sync_pipeline_to_linux.py` | `5406b6f2020a8ea30344ae173669a1c06423af46431896129a402c014be48833` | session helper; not imported by `workflow.py` / `code/pipeline.py` / `tests/` |

## Deleted (regenerable caches only)

- `__pycache__/` and `.pytest_cache/` trees under the kit.
- Temporary inventory scripts live under `/home/antiz/hydration-research/` and are not part of the published repo.

## Kept (not drafts)

- Canonical K6 pilot: `results/K_06/pilot_2/` (XYZ, xtb logs, Hessians, g98, provenance, CoCoS, mode follow).
- Canonical mirror: `results/K/6/{gwh,xtb,hessian,mode_follow,cocos,environment,final}/`.
- Kit packing files: `inputs/*/packing.vecfle`.
- Original snapshots: `originals/`.
- Pilot launch and postprocess scripts: `run_k06_pilot_2.py`, `postprocess_k06_confirmed_minima.py`.
- Packing analysis: `analyze_packings.py` (documented in the Russian README).
- CoCoS robustness helper: `check_cocos_robustness.py`.
- Windows→Linux sync: `tools/sync_to_linux.py`.
- Figure snapshot `figures/publication/figures_matplotlib.py` (byte-identical to `code/figures_matplotlib.py`; listed in `MATPLOTLIB_MANIFEST.json`).
- Local, unpublished trees remain gitignored and were **not** copied to GitHub: `vendors/`, `cases/`, `pilots/`, previous `archive/`.

## Duplicate SHA-256 groups that were retained

CoCoS `free_sorted_for_cocos.xyz` / `H12KO6_quick_cleanup.trj` copies across seeds, `hess/xtbopt.xyz` copies of free-opt geometries, and `results/K/6/` mirrors of `pilot_2` are job provenance. Identical xyzrender PNGs in `figures/burette/coordination/` and `figures/burette/structures/` are kept as the two documented figure slots.

SHA-256 manifest rewritten with 335 published-path files (`SHA256SUMS.txt`).

## Publication scope

Public `waterel` does not include `pilots/Li_07` (old local Euler trial), vendor checkouts, or the previous cleanup backup. Those stay in the local WSL tree / `_local_archive`.

