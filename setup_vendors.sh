#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
mkdir -p vendors
# Pinned author revisions (2026-10-30). GWH includes the intrinsic ZYZ commit.
GWH_SHA=7411811374b6422a15f23d6d86a47f2145f8ac47
COCOS_SHA=47f6680da2eee90d9b18536fd94d2bc8f6d05aad
MDATOM_SHA=3c43de70abef1f64fb31a2fe2693e4e73b88083c
XTB_TAG=v6.7.1
if [ ! -d vendors/Glorious-Water-Hole ]; then git clone https://gitlab.com/madschumacher/glorious-water-hole.git vendors/Glorious-Water-Hole; fi
if [ ! -d vendors/CoCoS ]; then git clone https://gitlab.com/madschumacher/cocoso.git vendors/CoCoS; fi
if [ ! -d vendors/MDAtomLocalEnvironment ]; then git clone https://gitlab.com/madschumacher/mdatomlocalenvironment.git vendors/MDAtomLocalEnvironment; fi
if [ ! -d vendors/xtb ]; then git clone --depth 1 --branch "$XTB_TAG" https://github.com/grimme-lab/xtb.git vendors/xtb; fi
git -C vendors/Glorious-Water-Hole fetch --all --tags
git -C vendors/Glorious-Water-Hole merge --ff-only "$GWH_SHA" || git -C vendors/Glorious-Water-Hole checkout "$GWH_SHA"
git -C vendors/CoCoS fetch --all --tags
git -C vendors/CoCoS merge --ff-only "$COCOS_SHA" || git -C vendors/CoCoS checkout "$COCOS_SHA"
git -C vendors/MDAtomLocalEnvironment fetch --all --tags
git -C vendors/MDAtomLocalEnvironment merge --ff-only "$MDATOM_SHA" || git -C vendors/MDAtomLocalEnvironment checkout "$MDATOM_SHA"
printf '\nGWH revision: '; git -C vendors/Glorious-Water-Hole rev-parse HEAD
printf 'CoCoS revision: '; git -C vendors/CoCoS rev-parse HEAD
printf 'MDAtomLocalEnvironment revision: '; git -C vendors/MDAtomLocalEnvironment rev-parse HEAD
printf 'xtb source: '; git -C vendors/xtb describe --tags --always
printf '\nDo not replace /usr/bin/xtb. Use the Intel-static 6.7.1 (edcfbbe) binary via source env_xtb.sh.\n'
