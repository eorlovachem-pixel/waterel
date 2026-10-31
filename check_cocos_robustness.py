#!/usr/bin/env python3
"""Sensitivity of CoCoS fast-mode results to seeds and rotational starting trials.
Run only on freely optimized structures; no scientific interpretation of constrained data.
"""
import argparse,csv,json,re,shutil,subprocess,sys,hashlib
from pathlib import Path

root=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--case',required=True)
    p.add_argument('--seeds',nargs='+',type=int,default=[101,202,303])
    p.add_argument('--samples',nargs='+',type=int,default=[50,100])
    p.add_argument('--cocos',default=str(root/'code/cocos_uploaded_snapshot.py'))
    a=p.parse_args()
    path=root/'cases'/a.case/'xtb_runs/postprocess'
    src=path/'confirmed_minima_for_cocos.xyz'
    if not src.exists():
        src=path/'free_sorted_for_cocos.xyz'
    ok=path/'COMPLETENESS.txt'
    if not src.exists() or not ok.exists() or 'COMPLETE: all requested inputs checked.' not in ok.read_text():
        raise RuntimeError('Only completed free-optimization datasets may be validated.')
    code=Path(a.cocos).resolve();wrapper=root/'code/seeded_cocos.py'
    data=[]
    for ns in a.samples:
        for seed in a.seeds:
            job=path/'cocos_sensitivity'/f'N{ns}_seed{seed}';job.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,job/'free_sorted_for_cocos.xyz')
            cmd=[sys.executable,str(wrapper),str(seed),str(code),'-f','free_sorted_for_cocos.xyz',
                 '--DoNotUseGlobalOpt','--NumSamples',str(ns),'--DoNotSaveDoubles']
            proc=subprocess.run(cmd,cwd=job,capture_output=True,text=True)
            (job/'cocos_log.txt').write_text(proc.stdout+'\n'+proc.stderr)
            if proc.returncode:raise RuntimeError(f'CoCoS failed ({ns}, {seed}), log {job}')
            survivors=[int(x) for x in re.findall(r'added structure (\d+)\b',proc.stdout)]
            vals={'samples':ns,'seed':seed,'num_unique':len(survivors),'survivor_sorted_indices':';'.join(map(str,survivors)),
                  'cocos_sha256':hashlib.sha256(code.read_bytes()).hexdigest()}
            data.append(vals)
            print(a.case,vals['samples'],vals['seed'],vals['num_unique'])
    f=path/'cocos_sensitivity.csv'
    with f.open('w',newline='') as g:
        writer=csv.DictWriter(g,fieldnames=data[0].keys());writer.writeheader();writer.writerows(data)
    print('Saved',f)
    if len(set(r['num_unique'] for r in data))>1:print('WARNING: CoCoS number of unique structures depends on parameters!')
    elif len(set(r['survivor_sorted_indices'] for r in data))>1:print('WARNING: same number but different identities retained!')
    else:print('No variation detected for this finite set of runs; still not a proof of algorithmic completeness.')

if __name__=='__main__':main()
