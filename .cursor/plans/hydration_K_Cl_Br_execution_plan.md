# K⁺, Cl⁻, Br⁻ hydrate execution plan

Canonical workdir: `/home/antiz/hydration-research/Hydration_unassigned_2026-10-30`. Public repo: `https://github.com/eorlovachem-pixel/waterel` `main` at `226af299` (CI green). Do not touch the xTB resolver.

**Hard gate:** production GWH `--Opt`, free xTB, Hessian, refine, and CoCoS stay refused while `assignment_status` is `ON_HOLD` or `NOT_ASSIGNED` unless Denis explicitly assigns that case. `DENIS_ASSIGNED=yes` is the only documented bypass; `ALLOW_FULL_BATCH` is not. Independent stages below complete without assignment. The 33×100 series is **not** marked done if the gate is still closed.

See the approved Cursor plan `k_cl_br_hydrates_9120c7c2` for the full specification. This file is the in-repo checkpoint copy for resume after interruptions.
