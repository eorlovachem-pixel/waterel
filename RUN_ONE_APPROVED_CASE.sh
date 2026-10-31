#!/usr/bin/env bash
# Runs the author's GWH -> full xTB -> CoCoS -> MDAtomLocalEnvironment protocol.
# Run ONLY after explicit assignment of this case by Denis.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
CASE=${1:-}
if [[ -z "$CASE" ]]; then
  echo "Usage: ALLOW_FULL_BATCH=yes bash RUN_ONE_APPROVED_CASE.sh K_06";exit 2
fi
# shellcheck disable=SC1091
source ./env_xtb.sh
ASSIGN=$(python3 -c "import csv; from pathlib import Path; c='$CASE';
rows={r['case']:r for r in csv.DictReader(Path('CASE_STATUS.csv').open(encoding='utf-8'))};
print((rows.get(c) or {}).get('assignment_status','UNKNOWN'))")
if [[ "$ASSIGN" == "ON_HOLD" && "${DENIS_ASSIGNED:-}" != "yes" ]]; then
  echo "Refusing $CASE: assignment_status=ON_HOLD. Denis has not assigned this case."
  echo "ALLOW_FULL_BATCH does not bypass ON_HOLD. After confirmation: DENIS_ASSIGNED=yes ALLOW_FULL_BATCH=yes."
  exit 1
fi
if [[ -z "${XTB_BIN:-}" ]] || [[ ! -x "${XTB_BIN}" ]]; then
  echo "Standalone xtb not set (XTB_BIN). env_xtb.sh must locate ~/xtb-standalone or vendors/xtb-dist.";exit 1
fi
if [[ "$(command -v xtb)" == "/usr/bin/xtb" ]] || [[ "$XTB_BIN" == "/usr/bin/xtb" ]]; then
  echo "Refusing /usr/bin/xtb for ANCOPT. source env_xtb.sh first.";exit 1
fi
python3 environment_check.py
python3 code/xtb_preflight.py --xtb "$XTB_BIN"
python3 workflow.py check-inputs
python3 workflow.py status --case "$CASE"
if [[ ! -f vendors/Glorious-Water-Hole/gloriouswaterhole.py ]]; then
  echo "Missing GWH. Run: bash setup_vendors.sh";exit 1
fi
if [[ "${ALLOW_FULL_BATCH:-}" != "yes" ]]; then
  echo "Refusing the 100-structure series. Pass preflight, get assignment, then ALLOW_FULL_BATCH=yes."
  exit 1
fi
python3 workflow.py gwh --case "$CASE" --nstruct 100
python3 workflow.py collect --case "$CASE"
python3 workflow.py prepare --case "$CASE"
python3 workflow.py relax --case "$CASE" --workers 2
python3 workflow.py hess --case "$CASE"
python3 workflow.py refine --case "$CASE"
python3 workflow.py dedup --case "$CASE" \
   --mdatom vendors/MDAtomLocalEnvironment/atom_local_environment_computer.py
python3 workflow.py aggregate
printf '\nFinished case: %s. Independently inspect final structures and logs before publication.\n' "$CASE"
