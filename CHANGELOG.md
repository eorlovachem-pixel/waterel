# Changelog (2026-10-31)

## Final audit / cleanup (2026-11-01)

- Backed up the Linux kit to `/home/antiz/hydration-research/_backups/` (outside this Git tree).
- Archived unused session helpers (`tools/finish_git.py`, `tools/git_waterel.py`, `tools/wsl_go.py`, `tools/run_all.py`, `tools/apply_cleanup.py`, `_sync_pipeline_to_linux.py`) to a local archive; they are not imported by `workflow.py` / `pipeline.py` / tests.
- Regenerable caches (`__pycache__`, `.pytest_cache`) removed. Unique XYZ, packing files, Hessians, logs, CoCoS, MDAtom, figures, and original kit files were kept.
- `SHA256SUMS.txt` now covers published paths only (not `vendors/`, `cases/`, `pilots/`, `archive/`).
- Details: `CLEANUP_REPORT.md`.

## Pipeline extension (same modules, no alternate pipeline)

- Code now accepts cations Li/Na/K/Rb (`--chrg 1`) and anions F/Cl/Br (`--chrg -1`) via `qc_checks.protocol_charge`. **No anion or Na/Rb rows added to `CASE_STATUS.csv`.**
- `workflow.py process-existing` / `visualize` post-process the existing K₆ pilot only (no GWH/xTB launch, no `ON_HOLD` bypass).
- First-shell descriptors: CN, ion–O distances, radial deformation, cosθ, P₂, O–ion–O, O–H···O, anion X···H–O, proton-transfer flag, XYZ SHA-256. CN cut is `BLOCKING_UNCONFIRMED`.
- Canonical mirror `results/K/6/{gwh,xtb,hessian,mode_follow,cocos,environment,final}/` without rewriting `results/K_06/pilot_2/`.
- Burette: scene packages under `figures/burette/scenes/`. Burette.app is macOS-only on this host; xyzrender CLI used when present; Matplotlib structures explicitly labeled as Matplotlib.
- Matplotlib numerical plots in `figures/publication/`. Missing n-series and anion figures are documented, not invented.
- Tests: anion charge, descriptors, visualizer byte-preserving XYZ copy. GitHub Actions `tests.yml`.
- Docs: `README.md`, `FOR_DENIS.md`, `RESULTS_INDEX.md`, `figures/FIGURE_INDEX.md`.

## Assignments / methodology (no calculation launch)

## Assignments / methodology (no calculation launch)

- Audited `CASE_STATUS.csv`, kit docs, and GWH vendor git. No confirmed Denis assignment of K_06. **Did not** change `ON_HOLD` / `NOT_ASSIGNED`. **Did not** set `DENIS_ASSIGNED=yes`.
- Na⁺ / Rb⁺ not added to the case table.
- Euler: documented SciPy `zyz` (extrinsic) vs `ZYZ` (intrinsic) vs author `'zyz'`. GWH source not rewritten.
- Added `MESSAGE_TO_DENIS.md` and `ASSIGNMENT_AUDIT.md`.

## Structure cleanup

- SHA-256 inventory: duplicate groups are mostly CoCoS working copies of the same XYZ and `hess/xtbopt.xyz` copies of free-opt geometries. Those job directories were **kept** (provenance).
- Session-only `_*.py` / `_*.sh` helpers moved to `archive/cleanup_backup_2026-10-31/` after backup. Not imported by `workflow.py` / `pipeline.py`.
- `results/K_06/pilot_2_pipeline_check/` (pipeline replay copy + duplicated CoCoS inputs) moved to the same backup. Stub `PIPELINE_CHECK_SUMMARY.json` left in place.
- Canonical data **untouched**: `results/K_06/pilot_2/` (XYZ, xtb logs, Hessian, provenance), `inputs/*.vecfle`, `pilots/`, `cases/`.
- `CASE_STATUS.json` regenerated from `CASE_STATUS.csv` (CSV is source of truth).
- Sync helper lives in `tools/sync_to_linux.py`.
- Added `.gitignore` for `__pycache__`, pytest cache, generated environment report, backup trees.

## Compatibility

- `workflow.py` still reads `CASE_STATUS.csv`.
- Production paths (`code/`, `inputs/`, `results/K_06/pilot_2/`) unchanged.
- CoCoS seed remains 20261030. Overwrite guards unchanged.
