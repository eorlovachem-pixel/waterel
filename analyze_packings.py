#!/usr/bin/env python3
"""Analyze the supplied initial packing vectors WITHOUT computing electronic energies."""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
with (root/'CASE_STATUS.csv').open() as f: cases=list(csv.DictReader(f))
results=[]
for case in cases:
    p=root/'inputs'/case['case']/'packing.vecfle'
    records=[line.split() for line in p.read_text().splitlines() if line.strip()]
    u=np.array([[float(x) for x in rec[1:4]] for rec in records]);rad=float(case['initial_M_O_A']);
    dots=np.clip(u@u.T,-1,1)
    i,j=np.triu_indices(len(u),1)
    theta=np.degrees(np.arccos(dots[i,j])); d=rad*np.linalg.norm(u[i]-u[j],axis=1)
    results.append({'case':case['case'],'ion':case['ion'],'n_water':case['n_water'],
        'initial_M_O_A':rad,'min_pair_angle_deg':float(theta.min()),
        'mean_pair_angle_deg':float(theta.mean()),
        'min_initial_O_O_A':float(d.min()),'mean_initial_O_O_A':float(d.mean()),
        'max_initial_O_O_A':float(d.max()),
        'note':'Sphere-packing input geometry only; NOT optimized xTB geometry'})
with (root/'input_geometry_summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=results[0].keys());w.writeheader();w.writerows(results)
fig,ax=plt.subplots(figsize=(8.4,4.8))
for ion in ['Li','K']:
    rows=[r for r in results if r['ion']==ion]
    ax.plot([int(r['n_water']) for r in rows],[r['min_pair_angle_deg'] for r in rows],marker='o',label=ion)
ax.set(xlabel='Number of water molecules in input, n',ylabel='Minimum pairwise angle between input directions (deg)',
       title='Validated sphere-packing inputs (not quantum-chemical results)')
ax.legend();ax.grid(alpha=.18)
fig.tight_layout();fig.savefig(root/'figures/input_min_angular_spacing.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(8.4,4.8))
for ion in ['Li','K']:
    rows=[r for r in results if r['ion']==ion]
    ax.plot([int(r['n_water']) for r in rows],[r['min_initial_O_O_A'] for r in rows],marker='o',label=ion)
ax.set(xlabel='Number of water molecules, n',ylabel='Minimum initial O-O distance (angstrom)',
       title='Geometric input O-O separations before optimization')
ax.legend();ax.grid(alpha=.18)
fig.tight_layout();fig.savefig(root/'figures/input_min_OO_distance.png',dpi=170);plt.close(fig)
print('Analyzed',len(results),'author-supplied input packings; did not compute energies or optimize.')
print('Li7:',next(r for r in results if r['case']=='Li_07'))
print('K6:',next(r for r in results if r['case']=='K_06'))
