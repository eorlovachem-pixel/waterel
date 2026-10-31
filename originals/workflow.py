#!/usr/bin/env python3
"""Orchestrate *real* GWH -> free xTB -> CoCoS -> MDAtomLocalEnvironment runs.
This file never fabricates quantum-chemical results.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, shutil, subprocess, sys, datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
import pipeline as P


def cases():
    with (ROOT/'CASE_STATUS.csv').open(encoding='utf-8') as f:
        return {r['case']:r for r in csv.DictReader(f)}


def select(args):
    db=cases()
    if args.case not in db:
        raise ValueError('Unknown case '+repr(args.case)+'; available: '+', '.join(sorted(db)))
    return db[args.case]


def git_revision(d):
    p=subprocess.run(['git','-C',str(d),'rev-parse','HEAD'],capture_output=True,text=True)
    return p.stdout.strip() if p.returncode==0 else 'UNAVAILABLE (non-git or uncommitted snapshot)'


def case_root(args):
    return ROOT/('pilots' if getattr(args,'pilot',False) else 'cases')/args.case


def validate_vecfle(args):
    allcases=cases()
    for cid,row in allcases.items():
        p=ROOT/'inputs'/cid/'packing.vecfle'
        if not p.exists():raise RuntimeError(f'Missing {p}')
        import numpy as np
        lines=[x.split() for x in p.read_text().splitlines() if x.strip()]
        n=int(row['n_water'])
        if len(lines)!=n:raise RuntimeError(f'{cid}: {len(lines)} sites, expected {n}')
        v=np.array([[float(s) for s in r[1:4]] for r in lines])
        rads=[float(s[4]) for s in lines]
        if any(l[0]!=row['ion'] or len(l)!=5 for l in lines): raise RuntimeError(f'{cid}: invalid elements or columns')
        if max(abs(np.linalg.norm(v,axis=1)-1))>1e-7:raise RuntimeError(f'{cid}: non-unit direction')
        if max(abs(x-float(row['initial_M_O_A'])) for x in rads)>1e-10:raise RuntimeError(f'{cid}: radius mismatch')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row['input_sha256']:raise RuntimeError(f'{cid}: input SHA256 mismatch')
        print(f'{cid}: {n} normalized directions, initial {row["initial_M_O_A"]} A, SHA256 OK')


def gwh(args):
    row=select(args);cr=case_root(args);cr.mkdir(parents=True,exist_ok=True)
    out=cr/'gwh';gwhpath=Path(args.gwh).expanduser().resolve()
    if not gwhpath.is_file():raise FileNotFoundError(f'GWH source not found: {gwhpath}. Run setup_vendors.sh first or use --gwh.')
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f'{out} already has files; do not restart a GWH batch on top of an earlier batch. Archive it first.')
    out.mkdir(exist_ok=True)
    inp=(ROOT/'inputs'/args.case/'packing.vecfle').resolve()
    cmd=[sys.executable,str(gwhpath),'-f',str(inp),'--Charge',str(row['charge']),
         '--NumStruct',str(args.nstruct),'--Random','--Opt','--AddMethod',' --alpb water --gfn 2 ','--DoNotCenter']
    info={'case':args.case,'cmd':cmd,'timestamp_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'source_vecfle_sha256':row['input_sha256'],'GWH_git_revision':git_revision(gwhpath.parent),
          'xtb_version':None,'nstruct_target':args.nstruct,'stage':'CONSTRAINED',
          'note':'No CoCoS filtering or energy-window preselection allowed before free optimization.'}
    ver=subprocess.run(['xtb','--version'],capture_output=True,text=True) if shutil.which('xtb') else None
    if ver:info['xtb_version']=(ver.stdout+'\n'+ver.stderr).strip()[:2000]
    (cr/'gwh_command.json').write_text(json.dumps(info,indent=2,ensure_ascii=False))
    print('STARTING REAL GWH:',args.case,'number of structures:',args.nstruct,flush=True)
    with (cr/'gwh_stdout_stderr.log').open('w') as log:
        p=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=False)
    (cr/'gwh_exit_code.txt').write_text(str(p.returncode)+'\n')
    print('GWH exit code:',p.returncode,'; inspect',cr/'gwh_stdout_stderr.log')
    if p.returncode!=0:raise RuntimeError('GWH exited nonzero; do not collect without checking logs')


def collect(args):
    row=select(args);cr=case_root(args);out=cr/'gwh'
    if not out.is_dir():raise RuntimeError('GWH result folder missing: '+str(out))
    if not (cr/'gwh_command.json').exists():raise RuntimeError('No GWH provenance gwh_command.json: do not treat arbitrary files as this run')
    provenance=json.loads((cr/'gwh_command.json').read_text())
    if provenance.get('case')!=args.case:raise RuntimeError('Case ID mismatch')
    target=provenance['nstruct_target'];results=[];errors=[]
    for i in range(target):
        d=out/f'gwh_opt_str{i:04d}'
        xyz=d/f'mol_str{i:04d}_fin.xyz'
        err=d/'qc_tmp.err'
        if not xyz.is_file():errors.append(f'{i:04d}: missing {xyz.name}');continue
        if not err.is_file() or 'normal termination of xtb' not in err.read_text(errors='replace').lower():
            errors.append(f'{i:04d}: xTB normal termination not confirmed in qc_tmp.err');continue
        try:
            frames=P.read_xyz(xyz)
            if len(frames)!=1:raise ValueError('expected one frame')
            fr=frames[0]
            wanted={row['ion']:1,'O':int(row['n_water']),'H':2*int(row['n_water'])}
            found={s:fr['symbols'].count(s) for s in set(fr['symbols'])}
            if found!=wanted:raise ValueError(f'wrong formula {found}, wanted {wanted}')
            en=P.energy_in_comment(fr['comment'])
            if en is None:raise ValueError('missing E = ... Eh in xyz comment')
            results.append((i,fr,float(en)))
        except Exception as e:
            errors.append(f'{i:04d}: parse or validation error: {e}')
    cr.mkdir(exist_ok=True)
    (cr/'collect_report.json').write_text(json.dumps({'case':args.case,'target':target,'valid':len(results),
        'invalid':errors,'STAGE':'CONSTRAINED_GWH_ONLY'},ensure_ascii=False,indent=2))
    print(args.case,'verified GWH outputs',len(results),'/',target)
    if errors:
        print('\n'.join(errors[:20]));raise RuntimeError('Incomplete GWH run. No publication-grade data constructed. See collect_report.json.')
    assert len(results)==target
    source=cr/'source';source.mkdir(exist_ok=True)
    xyzout=source/'constrained_multiframe.xyz'
    enout=source/'energies.csv'
    if xyzout.exists() or enout.exists():raise RuntimeError('Existing source data: refusing overwrite')
    P.write_xyz(xyzout,[x[1] for x in results]);min_e=min(x[2] for x in results)
    with enout.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['id','batch','structure','E_initial_Eh','E_final_Eh','nfev','success','delta_E_kcal_mol']);w.writeheader()
        for i,fr,e in results:w.writerow({'id':f'gwh_opt_str{i:04d}','batch':'single_GWH_batch','structure':i,
            'E_initial_Eh':'','E_final_Eh':format(e,'.14f'),'nfev':'','success':'True',
            'delta_E_kcal_mol':format((e-min_e)*P.HARTREE_TO_KCAL_MOL,'.8f')})
    print('Wrote',xyzout,'and',enout)


def prep(args):
    select(args);cr=case_root(args)
    src=cr/'source';inp=src/'constrained_multiframe.xyz';en=src/'energies.csv'
    if not inp.is_file() or not en.is_file(): raise RuntimeError('First run gwh and collect')
    row=cases()[args.case]
    out=cr/'xtb_runs'
    if not out.is_dir():
        P.prepare(argparse.Namespace(input=str(inp),energies_csv=str(en),output=str(out),
            central=row['ion'],charge=1,allow_existing=False))
    else:
        print('xtb_runs exists. Will not overwrite any prepared input; confirm it is complete.')
    P.audit(argparse.Namespace(input=str(inp),energies_csv=str(en),output=str(cr/'constrained_audit'),
              central=row['ion'],shell_cut=None))


def relax(args):
    select(args);cr=case_root(args);out=cr/'xtb_runs'
    if not (out/'manifest.json').is_file():raise RuntimeError('Run `prepare` first')
    P.relax(argparse.Namespace(output=str(out),xtb=args.xtb,workers=args.workers,limit=args.limit))


def dedup(args):
    select(args);cr=case_root(args);out=cr/'xtb_runs'
    if not (out/'manifest.json').is_file():raise RuntimeError('Run prepare and relax first')
    coc=ROOT/'code/cocos_uploaded_snapshot.py'
    if args.cocos:coc=Path(args.cocos).expanduser().resolve()
    md=None
    if args.mdatom:md=Path(args.mdatom).expanduser().resolve()
    P.dedup(argparse.Namespace(output=str(out),cocos=str(coc),seed=args.seed,
        samples=args.samples,shell_cut=None,mdatom=str(md) if md else None,
        allow_incomplete=False))


def status(args):
    db=cases(); selected={args.case:db[args.case]} if args.case else db
    for name,row in selected.items():
        cr=ROOT/'cases'/name
        gwh=cr/'gwh';s=cr/'source';xtb=cr/'xtb_runs'
        gwhfiles=len(list(gwh.glob('gwh_opt_str*/mol_str*_fin.xyz'))) if gwh.is_dir() else 0
        collect=('yes' if (s/'constrained_multiframe.xyz').exists() else 'no')
        planned= int(row['n_start_target'])
        prepared=len(list(xtb.glob('frame_*/input.xyz'))) if xtb.exists() else 0
        free=len(list(xtb.glob('frame_*/xtbopt.xyz'))) if xtb.exists() else 0
        post=('yes' if (xtb/'postprocess/unique_representatives.csv').exists() else 'no')
        print(f'{name}: GWH={gwhfiles}/{planned}; collected={collect}; xTB inputs={prepared}; xTB XYZ present={free}; postprocessed={post}')


def aggregate(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    db=cases(); records=[]
    for name,row in db.items():
        post=ROOT/'cases'/name/'xtb_runs'/'postprocess'
        comple=post/'COMPLETENESS.txt';rep=post/'unique_representatives.csv';allfree=post/'all_freely_optimized.csv'
        if not (comple.exists() and rep.exists() and allfree.exists()):continue
        if 'COMPLETE: all requested inputs checked.' not in comple.read_text():continue
        with rep.open() as f:unique=list(csv.DictReader(f))
        with allfree.open() as f:allrows=list(csv.DictReader(f))
        if len(allrows)!=int(row['n_start_target']):continue
        records.append({'case':name,'ion':row['ion'],'n_water':int(row['n_water']),'n_initial':len(allrows),
             'n_unique_after_free_relaxation':len(unique),
             'lowest_free_energy_Eh':min(float(x['energy_free_Eh']) for x in allrows),
             'max_relative_energy_unique_kcal_mol':max(float(x['delta_free_kcal_mol']) for x in unique)})
    if not records:raise RuntimeError('No completed case is available; do not generate made-up cross-case summaries.')
    dest=ROOT/'aggregate_results';dest.mkdir(exist_ok=True)
    P.write_csv(dest/'completed_cases_summary.csv',records)
    for ion in sorted(set(r['ion'] for r in records)):
        items=sorted([r for r in records if r['ion']==ion],key=lambda r:r['n_water'])
        fig,ax=plt.subplots(figsize=(7,4))
        ax.bar([str(r['n_water']) for r in items],[r['n_unique_after_free_relaxation'] for r in items])
        ax.set(xlabel='Number of water molecules, n',ylabel='Unique structures after free relaxation',
               title=f'{ion}+ (H2O)n: unique conformers found, GFN2-xTB/ALPB')
        fig.tight_layout();fig.savefig(dest/f'{ion}_unique_conformers.png',dpi=170);plt.close(fig)
    print('Wrote summary for',len(records),'completed cases:',dest)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('status',help='read-only status of each case');a.add_argument('--case');a.set_defaults(func=status)
    a=sub.add_parser('check-inputs');a.set_defaults(func=validate_vecfle)
    a=sub.add_parser('aggregate',help='only fully completed cases');a.set_defaults(func=aggregate)
    for cmd,fn in [('gwh',gwh),('collect',collect),('prepare',prep),('relax',relax),('dedup',dedup)]:
        a=sub.add_parser(cmd);a.add_argument('--case',required=True,choices=sorted(cases()));a.add_argument('--pilot',action='store_true',help='keep small pilot in separate pilots/ folder')
        if cmd=='gwh':
            a.add_argument('--gwh',default=str(ROOT/'vendors/Glorious-Water-Hole/gloriouswaterhole.py'))
            a.add_argument('--nstruct',type=int,default=100)
        elif cmd=='relax':
            a.add_argument('--xtb',default='xtb');a.add_argument('--workers',type=int,default=2);a.add_argument('--limit',type=int)
        elif cmd=='dedup':
            a.add_argument('--cocos',help='optional current source cocos.py; default bundled submitted snapshot')
            a.add_argument('--mdatom',help='path to installed MDAtomLocalEnvironment atom_local_environment_computer.py')
            a.add_argument('--seed',type=int,default=20261030);a.add_argument('--samples',type=int,default=50)
        a.set_defaults(func=fn)
    x=p.parse_args();x.func(x)

if __name__=='__main__':main()
