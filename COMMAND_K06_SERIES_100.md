# Exact command for the 100-structure K+(H2O)6 series

**Do not run this until Denis assigns K_06.** `assignment_status=ON_HOLD`.
`ALLOW_FULL_BATCH=yes` does **not** bypass that hold.

The 2-structure pilot under `results/K_06/pilot_2/` and `pilots/Li_07` / `cases/Li_07` must stay untouched. Production writes to `cases/K_06/`, which currently does not exist.

From a clean PowerShell / cmd:

```bat
wsl -d Ubuntu
```

Then:

```bash
cd ~/hydration-research/Hydration_unassigned_2026-10-30
source env_xtb.sh
echo "$XTB_BIN"
which xtb
# must NOT print /usr/bin/xtb; XTB_BIN must be the Intel-static 6.7.1 (edcfbbe)

python3 environment_check.py
python3 -m pytest tests -q
python3 code/xtb_preflight.py
python3 workflow.py check-inputs
python3 workflow.py status --case K_06
```

`workflow.py gwh --case K_06` currently raises `ON_HOLD`. The launch recipe after Denis confirms:

```bash
DENIS_ASSIGNED=yes ALLOW_FULL_BATCH=yes bash RUN_ONE_APPROVED_CASE.sh K_06
```

Equivalent explicit steps, still **not** to be launched now:

```bash
source env_xtb.sh
# Python resolves the standalone binary automatically; --xtb is optional.
DENIS_ASSIGNED=yes python3 workflow.py gwh --case K_06 --nstruct 100
DENIS_ASSIGNED=yes python3 workflow.py collect --case K_06
DENIS_ASSIGNED=yes python3 workflow.py prepare --case K_06
DENIS_ASSIGNED=yes python3 workflow.py relax --case K_06 --workers 2
DENIS_ASSIGNED=yes python3 workflow.py hess --case K_06
DENIS_ASSIGNED=yes python3 workflow.py refine --case K_06
DENIS_ASSIGNED=yes python3 workflow.py dedup --case K_06 \
  --mdatom vendors/MDAtomLocalEnvironment/atom_local_environment_computer.py
python3 workflow.py aggregate
```

`refine` repeats ±0.15 Å displacement / reopt / Hessian up to `--max-iters` (default 3) until ν ≥ −5 cm⁻¹. CoCoS seed remains **20261030**. CoCoS receives confirmed minima only. Do not pass `--pilot` for the 100-structure job.
