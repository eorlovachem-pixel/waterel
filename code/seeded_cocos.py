#!/usr/bin/env python3
"""Run the ORIGINAL author CoCoS using a deterministic random seed."""
import sys, random, runpy
import numpy as np
if len(sys.argv)<4: raise SystemExit('usage: seeded_cocos.py SEED /path/to/cocos.py [-f ...]')
seed=int(sys.argv[1]); script=sys.argv[2]
np.random.seed(seed);random.seed(seed)
sys.argv=[script,*sys.argv[3:]]
runpy.run_path(script,run_name='__main__')
