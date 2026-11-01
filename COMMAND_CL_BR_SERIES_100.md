# Exact commands for Cl⁻(H₂O)ₙ and Br⁻(H₂O)ₙ, n = 4–14

**Do not run this until Denis assigns the anion cases.**
`assignment_status=NOT_ASSIGNED` for `Cl_04`–`Cl_14` and `Br_04`–`Br_14`.
`ALLOW_FULL_BATCH=yes` does **not** bypass that status.

Charge −1, multiplicity 1 (UHF=0). GWH Euler convention remains `from_euler('zyz')`.
Do not CoCoS-filter before free optimization.

GWH CLI (verified against `gloriouswaterhole.py --help`, commit `7411811374b6422a15f23d6d86a47f2145f8ac47`):

```bash
gloriouswaterhole.py -f <packing.vecfle> --Charge -1 --NumStruct 20 --Random --Opt --AddMethod " --alpb water --gfn 2 " --DoNotCenter
```

`--Charge` is `int`, default 0; −1 is accepted. `--Mult` default 1. `--Opt` runs SciPy Powell on Euler angles with an xTB SPE at every step. Even without `--Opt`, GWH still calls one xTB SPE for the initial random Euler geometry. That is why this recipe is gated.

After Denis assigns, from WSL Ubuntu:

```bash
cd ~/hydration-research/Hydration_unassigned_2026-10-30
source env_xtb.sh
python3 environment_check.py
python3 -m pytest tests -q
python3 code/xtb_preflight.py
python3 workflow.py check-inputs
python3 workflow.py status --case Cl_04
```

Twenty-structure CLI check (still blocked today):

```bash
# DENIS_ASSIGNED=yes python3 workflow.py gwh --case Cl_04 --nstruct 20 --pilot
```

One hundred starts per stoichiometry (still blocked today):

```bash
# DENIS_ASSIGNED=yes ALLOW_FULL_BATCH=yes bash RUN_ONE_APPROVED_CASE.sh Cl_04
# … repeat for Cl_05–Cl_14 and Br_04–Br_14
```

Free-opt protocol after GWH, unchanged: `normal → Hessian → mode following ±0.15 Å → tight → final Hessian`. CoCoS and MDAtomLocalEnvironment only on frequency-confirmed minima (`ν ≥ −5 cm⁻¹`).

Do not modify K, Li, Na, Rb, or other participants’ results.
