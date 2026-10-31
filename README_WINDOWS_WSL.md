See the copy at the project root `README_WINDOWS_WSL.md`, or follow the commands below from Ubuntu.

```bash
wsl -d Ubuntu
cd ~/hydration-research/Hydration_unassigned_2026-10-30
source env_xtb.sh
# exports XTB_BIN to the Intel-static binary; Python also auto-resolves it
python3 environment_check.py
python3 -m pytest tests -q
python3 code/xtb_preflight.py --xtb "$HOME/xtb-standalone/xtb-dist/bin/xtb"
python3 workflow.py check-inputs
```

Do not replace `/usr/bin/xtb`. Do not run the 100-structure series until preflight passes and K6 is assigned.
