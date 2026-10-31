#!/usr/bin/env bash
# Runs the author's GWH -> full xTB -> CoCoS -> MDAtomLocalEnvironment protocol.
# Run ONLY after explicit assignment of this case by Denis.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
CASE=${1:-}
if [[ -z "$CASE" ]]; then
  echo "Usage: bash RUN_ONE_APPROVED_CASE.sh Li_07";exit 2
fi
python workflow.py status --case "$CASE"
if [[ ! -f vendors/Glorious-Water-Hole/gloriouswaterhole.py ]]; then
  echo "Missing GWH. Run: bash setup_vendors.sh";exit 1
fi
if ! command -v xtb >/dev/null 2>&1; then
  echo "Missing xtb executable. Install xtb before actual chemistry calculations.";exit 1
fi
python workflow.py gwh --case "$CASE" --nstruct 100
python workflow.py collect --case "$CASE"
python workflow.py prepare --case "$CASE"
python workflow.py relax --case "$CASE" --workers 2
python workflow.py dedup --case "$CASE" \
   --mdatom vendors/MDAtomLocalEnvironment/atom_local_environment_computer.py
python check_cocos_robustness.py --case "$CASE"
python workflow.py aggregate
printf '\nFinished case: %s. Independently inspect final structures and logs before publication.\n' "$CASE"
