# Source from bash:  source env_xtb.sh
# Prepends the Intel-static xTB 6.7.1 (edcfbbe) without touching /usr/bin/xtb.
# Also exports XTB_BIN so Python (workflow.py / pipeline.py) uses the same binary
# even if a later PATH lookup would hit /usr/bin/xtb.
_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_STANDALONE="${HOME}/xtb-standalone/xtb-dist"
_LOCAL="${_ROOT}/vendors/xtb-dist"
_CHOSEN=""
if [ -x "${_ROOT}/bin/xtb" ]; then
  _CHOSEN="${_ROOT}/bin"
  export XTBPATH="${_ROOT}/share/xtb"
elif [ -x "${_LOCAL}/bin/xtb" ]; then
  _CHOSEN="${_LOCAL}/bin"
  export XTBPATH="${_LOCAL}/share/xtb"
elif [ -x "${_STANDALONE}/bin/xtb" ]; then
  _CHOSEN="${_STANDALONE}/bin"
  export XTBPATH="${_STANDALONE}/share/xtb"
fi
if [ -n "${_CHOSEN}" ]; then
  export PATH="${_CHOSEN}:${PATH}"
  export XTB_BIN="${_CHOSEN}/xtb"
else
  echo "env_xtb.sh: no static xtb found. Keep /usr/bin/xtb unused for ANCOPT." >&2
  unset XTB_BIN || true
fi
hash -r 2>/dev/null || true
command -v xtb >/dev/null && xtb --version | sed -n '1,20p'
echo "which xtb: $(command -v xtb)"
echo "XTB_BIN: ${XTB_BIN:-unset}"
