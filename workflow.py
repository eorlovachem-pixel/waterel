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
    row=db[args.case]
    import qc_checks
    qc_checks.assert_supported_ion(row['ion'])
    expected=qc_checks.protocol_charge(row['ion'])
    if int(row['charge']) != expected:
        raise ValueError(
            f"{args.case}: protocol requires charge {expected:+d} for {row['ion']}, got {row['charge']}"
        )
    row=dict(row)
    row['packing_path']=str((ROOT/'inputs'/args.case/'packing.vecfle').resolve())
    row['seed']=int(getattr(args,'seed', 20261030) or 20261030)
    row['n_starts']=int(row['n_start_target'])
    row['case_status']=row.get('assignment_status')
    row['species']=row['ion']
    return row


def require_launch_allowed(args, row):
    """Refuse ON_HOLD and NOT_ASSIGNED production chemistry. ALLOW_FULL_BATCH is not a bypass."""
    status=(row.get('assignment_status') or '').strip().upper()
    if status in ('ON_HOLD', 'NOT_ASSIGNED') and os.environ.get('DENIS_ASSIGNED')!='yes':
        raise RuntimeError(
            f"{args.case} assignment_status={status}. "
            "Do not start GWH/opt/hess/refine/CoCoS until Denis assigns this case. "
            "ALLOW_FULL_BATCH is not sufficient. After confirmation set DENIS_ASSIGNED=yes."
        )


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


def resolve_xtb(explicit=None):
    import qc_checks
    return qc_checks.resolve_working_xtb(explicit, kit=ROOT)


def prepend_xtb_path(xtb_bin):
    if xtb_bin is None:
        return
    bin_dir=str(Path(xtb_bin).resolve().parent)
    os.environ['PATH']=bin_dir+os.pathsep+os.environ.get('PATH','')


def gwh(args):
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args);cr.mkdir(parents=True,exist_ok=True)
    out=cr/'gwh';gwhpath=Path(args.gwh).expanduser().resolve()
    if not gwhpath.is_file():raise FileNotFoundError(f'GWH source not found: {gwhpath}. Run setup_vendors.sh first or use --gwh.')
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f'{out} already has files; do not restart a GWH batch on top of an earlier batch. Archive it first.')
    out.mkdir(exist_ok=True)
    xtb_bin=resolve_xtb(getattr(args,'xtb',None)); prepend_xtb_path(xtb_bin)
    gwh_euler=None
    try:
        sys.path.insert(0,str(ROOT/'code'))
        import qc_checks
        gwh_euler=qc_checks.gwh_uses_intrinsic_zyz(gwhpath)
        if gwh_euler.get('zyz_ZYZ_discrepancy'):
            raise RuntimeError(
                "GWH source contains both from_euler('zyz') and from_euler('ZYZ'). "
                "SciPy lowercase is extrinsic, uppercase is intrinsic. Do not silently pick one."
            )
        if not gwh_euler.get('author_method_intact', gwh_euler.get('intrinsic_zyz')):
            raise RuntimeError(
                "GWH source does not contain the author from_euler('zyz') string "
                "(commit 7411811374b6422a15f23d6d86a47f2145f8ac47). "
                "SciPy docs: 'zyz' is extrinsic, 'ZYZ' is intrinsic. "
                "Do not rewrite zyz↔ZYZ without an author decision."
            )
    except ImportError:
        pass
    inp=(ROOT/'inputs'/args.case/'packing.vecfle').resolve()
    cmd=[sys.executable,str(gwhpath),'-f',str(inp),'--Charge',str(row['charge']),
         '--NumStruct',str(args.nstruct),'--Random','--Opt','--AddMethod',' --alpb water --gfn 2 ','--DoNotCenter']
    info={'case':args.case,'cmd':cmd,'timestamp_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'source_vecfle_sha256':row['input_sha256'],'GWH_git_revision':git_revision(gwhpath.parent),
          'gwh_euler':gwh_euler,
          'xtb_binary':str(xtb_bin) if xtb_bin else None,
          'xtb_version':None,'nstruct_target':args.nstruct,'stage':'CONSTRAINED',
          'note':'No CoCoS filtering or energy-window preselection allowed before free optimization.'}
    xtb_exe=str(xtb_bin) if xtb_bin else 'xtb'
    ver=subprocess.run([xtb_exe,'--version'],capture_output=True,text=True) if (xtb_bin or shutil.which('xtb')) else None
    if ver:info['xtb_version']=(ver.stdout+'\n'+ver.stderr).strip()[:2000]
    (cr/'gwh_command.json').write_text(json.dumps(info,indent=2,ensure_ascii=False))
    print('STARTING REAL GWH:',args.case,'number of structures:',args.nstruct,flush=True)
    with (cr/'gwh_stdout_stderr.log').open('w') as log:
        p=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=False)
    (cr/'gwh_exit_code.txt').write_text(str(p.returncode)+'\n')
    print('GWH exit code:',p.returncode,'; inspect',cr/'gwh_stdout_stderr.log')
    if p.returncode!=0:raise RuntimeError('GWH exited nonzero; do not collect without checking logs')


def collect(args):
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args);out=cr/'gwh'
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
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args)
    src=cr/'source';inp=src/'constrained_multiframe.xyz';en=src/'energies.csv'
    if not inp.is_file() or not en.is_file(): raise RuntimeError('First run gwh and collect')
    row=cases()[args.case]
    out=cr/'xtb_runs'
    if not out.is_dir():
        P.prepare(argparse.Namespace(input=str(inp),energies_csv=str(en),output=str(out),
            central=row['ion'],charge=int(row['charge']),allow_existing=False))
    else:
        print('xtb_runs exists. Will not overwrite any prepared input; confirm it is complete.')
    P.audit(argparse.Namespace(input=str(inp),energies_csv=str(en),output=str(cr/'constrained_audit'),
              central=row['ion'],shell_cut=None))


def relax(args):
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args);out=cr/'xtb_runs'
    if not (out/'manifest.json').is_file():raise RuntimeError('Run `prepare` first')
    xtb_bin=resolve_xtb(getattr(args,'xtb',None)); prepend_xtb_path(xtb_bin)
    P.relax(argparse.Namespace(output=str(out),xtb=str(xtb_bin),workers=args.workers,limit=args.limit))


def hess(args):
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args);out=cr/'xtb_runs'
    if not (out/'manifest.json').is_file():raise RuntimeError('Run prepare and relax first')
    xtb_bin=resolve_xtb(getattr(args,'xtb',None)); prepend_xtb_path(xtb_bin)
    P.hess(argparse.Namespace(output=str(out),xtb=str(xtb_bin),
                              workers=args.workers,limit=args.limit,imag_thresh=args.imag_thresh))


def refine(args):
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args);out=cr/'xtb_runs'
    if not (out/'manifest.json').is_file():raise RuntimeError('Run prepare, relax and hess first')
    xtb_bin=resolve_xtb(getattr(args,'xtb',None)); prepend_xtb_path(xtb_bin)
    P.refine(argparse.Namespace(output=str(out),xtb=str(xtb_bin),
                                limit=args.limit,imag_thresh=args.imag_thresh,step_A=args.step_A,
                                max_iters=getattr(args,'max_iters', P.DEFAULT_MAX_REFINE_ITERS)))


def preflight(args):
    env=subprocess.run([sys.executable,str(ROOT/'environment_check.py')],check=False)
    cmd=[sys.executable,str(ROOT/'code'/'xtb_preflight.py')]
    xtb_bin=None
    try:
        xtb_bin=resolve_xtb(getattr(args,'xtb',None))
    except RuntimeError as exc:
        raise RuntimeError(str(exc)) from exc
    if xtb_bin:
        cmd += ['--xtb', str(xtb_bin)]
        prepend_xtb_path(xtb_bin)
    xtb=subprocess.run(cmd,check=False)
    if env.returncode!=0 or xtb.returncode!=0:
        raise RuntimeError('Preflight failed. Do not start a production series. See environment_check_report.json and preflight_runs/.')
    print('Preflight passed: real opt+hess completed with the working xTB binary.')


def dedup(args):
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args);out=cr/'xtb_runs'
    if not (out/'manifest.json').is_file():raise RuntimeError('Run prepare and relax first')
    coc=ROOT/'code/cocos_uploaded_snapshot.py'
    if args.cocos:coc=Path(args.cocos).expanduser().resolve()
    md=None
    if args.mdatom:md=Path(args.mdatom).expanduser().resolve()
    P.dedup(argparse.Namespace(output=str(out),cocos=str(coc),seed=args.seed,
        samples=args.samples,shell_cut=None,mdatom=str(md) if md else None,
        allow_incomplete=False,imag_thresh=getattr(args,'imag_thresh',5.0),
        step_A=getattr(args,'step_A',0.15),
        max_iters=getattr(args,'max_iters', P.DEFAULT_MAX_REFINE_ITERS),
        skip_robustness=getattr(args,'skip_robustness', False)))


def status(args):
    db=cases(); selected={args.case:db[args.case]} if args.case else db
    for name,row in selected.items():
        for kind in ('cases','pilots'):
            cr=ROOT/kind/name
            if not cr.exists() and kind=='pilots':
                continue
            if kind=='cases' and not cr.exists() and (ROOT/'pilots'/name).exists():
                continue
            gwh=cr/'gwh';s=cr/'source';xtb=cr/'xtb_runs'
            gwhfiles=len(list(gwh.glob('gwh_opt_str*/mol_str*_fin.xyz'))) if gwh.is_dir() else 0
            collect=('yes' if (s/'constrained_multiframe.xyz').exists() else 'no')
            planned= int(row['n_start_target'])
            prepared=len(list(xtb.glob('frame_*/input.xyz'))) if xtb.exists() else 0
            free=len(list(xtb.glob('frame_*/xtbopt.xyz'))) if xtb.exists() else 0
            hess=len(list(xtb.glob('frame_*/hess/xtb_hess.out'))) if xtb.exists() else 0
            conf=xtb/'postprocess'/'confirmed_minima.csv'
            nconf=0
            if conf.exists():
                with conf.open(encoding='utf-8') as fh:
                    nconf=max(0, sum(1 for _ in fh)-1)
            post=('yes' if (xtb/'postprocess/unique_representatives.csv').exists() else 'no')
            print(f'{kind}/{name}: assignment={row.get("assignment_status","?")}; GWH={gwhfiles}/{planned}; collected={collect}; xTB inputs={prepared}; xTB XYZ present={free}; hess={hess}; confirmed_csv={nconf}; postprocessed={post}')


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
               title=f'{ion}(H2O)n: unique conformers found, GFN2-xTB/ALPB')
        fig.tight_layout();fig.savefig(dest/f'{ion}_unique_conformers.png',dpi=170);plt.close(fig)
    print('Wrote summary for',len(records),'completed cases:',dest)


def process_existing(args):
    """Post-process already computed results. Does not launch GWH/xTB and does not bypass ON_HOLD."""
    import export_canonical
    import visualize
    import figures_matplotlib as FM
    dest, rows = export_canonical.export_k06()
    entries = []
    for r in rows:
        entries.append({
            'id': r['id'],
            'ion': r['ion'],
            'n_water': int(r['n_water']),
            'xyz_path': str(ROOT / r['source_xyz']),
            'energy_Eh': float(r['energy_Eh']),
            'cocos_group': r['cocos_group'],
            'origin': r['origin'],
            'descriptors': None,
        })
    vis = visualize.visualize_confirmed(entries, figure_prefix='Fig_K06_Pilot')
    table = dest / 'final' / 'descriptors.csv'
    plots = FM.make_k06_publication_plots(table)
    import mol_matplotlib
    import write_indexes
    mol_files = mol_matplotlib.k06_matplotlib_structures(entries)
    write_indexes.write_figure_index()
    write_indexes.write_results_index()
    print('Canonical results:', dest)
    print('Visualize successful exports:', vis.get('n_successful_exports'))
    print('Matplotlib numeric:', plots.get('files'))
    print('Matplotlib structures (not Burette):', [str(p) for p in mol_files])


def visualize_cmd(args):
    process_existing(args)


def _count_glob(path, pattern):
    p=Path(path)
    return len(list(p.glob(pattern))) if p.exists() else 0


def write_status_json(args, extra=None):
    row=select(args)
    cr=case_root(args)
    xtb=cr/'xtb_runs'
    post=xtb/'postprocess'
    info={
        'case':args.case,
        'ion':row['ion'],
        'n_water':int(row['n_water']),
        'charge':int(row['charge']),
        'assignment_status':row.get('assignment_status'),
        'timestamp_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'gwh_ok':_count_glob(cr/'gwh','gwh_opt_str*/mol_str*_fin.xyz'),
        'n_start_target':int(row['n_start_target']),
        'relax_ok':_count_glob(xtb,'frame_*/xtbopt.xyz'),
        'hess_ok':_count_glob(xtb,'frame_*/hess/xtb_hess.out'),
        'saddles':None,
        'confirmed':None,
        'cocos_groups':None,
        'packing_sha256':row.get('input_sha256'),
        'note':'Counts are file presence, not invented energies.',
    }
    conf=post/'confirmed_minima.csv'
    if conf.is_file():
        with conf.open(encoding='utf-8') as fh:
            info['confirmed']=max(0,sum(1 for _ in fh)-1)
    uniq=post/'unique_representatives.csv'
    if uniq.is_file():
        with uniq.open(encoding='utf-8') as fh:
            info['cocos_groups']=max(0,sum(1 for _ in fh)-1)
    if extra:
        info.update(extra)
    dest=cr/'STATUS.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(info,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('Wrote',dest)
    return dest


def resume(args):
    """Continue unfinished relax/hess/refine. Does not restart a non-empty GWH batch."""
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args)
    if not (cr/'xtb_runs'/'manifest.json').is_file():
        raise RuntimeError('Nothing to resume: run prepare first. GWH is not restarted into a non-empty folder.')
    relax(args); hess(args); refine(args)
    write_status_json(args, {'stage':'resume'})


def run_series(args):
    """GWH -> collect -> prepare -> relax -> hess -> refine. Still gated."""
    row=select(args); require_launch_allowed(args, row)
    cr=case_root(args)
    gwh_dir=cr/'gwh'
    nstruct=int(getattr(args,'nstruct',100) or 100)
    if not gwh_dir.exists() or not any(gwh_dir.iterdir()):
        args.nstruct=nstruct
        if not hasattr(args,'gwh'):
            args.gwh=str(ROOT/'vendors/Glorious-Water-Hole/gloriouswaterhole.py')
        if not hasattr(args,'xtb'):
            args.xtb=None
        gwh(args)
    else:
        print(f'{gwh_dir} already has files; skipping GWH (resume).')
    collect(args); prep(args); relax(args); hess(args); refine(args)
    write_status_json(args, {'stage':'run-series','nstruct':nstruct})


def export_results(args):
    import export_canonical
    dest, rows = export_canonical.export_case(args.case)
    write_status_json(args, {'stage':'export-results','exported':str(dest),'n_rows':len(rows)})
    print('Exported', dest)


def upload_drive(args):
    script=ROOT/'tools'/'drive_upload.py'
    p=subprocess.run([sys.executable,str(script)],check=False)
    if p.returncode!=0:
        raise RuntimeError('Drive upload did not complete. Staging remains in drive_staging/.')


def verify_delivery(args):
    stage=ROOT/'drive_staging'
    report=stage/'DRIVE_UPLOAD_REPORT.json'
    if not report.is_file():
        raise RuntimeError('No DRIVE_UPLOAD_REPORT.json; upload has not succeeded.')
    data=json.loads(report.read_text(encoding='utf-8'))
    print(json.dumps(data,indent=2)[:4000])
    missing=[]
    for folder, info in (data.get('folders') or {}).items():
        if not info.get('id'):
            missing.append(folder)
    if missing:
        raise RuntimeError('Drive folders missing ids: '+', '.join(missing))
    print('Delivery report present for', list((data.get('folders') or {}).keys()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('status',help='read-only status of each case');a.add_argument('--case');a.set_defaults(func=status)
    a=sub.add_parser('check-inputs');a.set_defaults(func=validate_vecfle)
    a=sub.add_parser('aggregate',help='only fully completed cases');a.set_defaults(func=aggregate)
    a=sub.add_parser('preflight');a.add_argument('--xtb');a.set_defaults(func=preflight)
    a=sub.add_parser('process-existing', help='descriptors + Burette scenes + Matplotlib from existing K6 pilot; no GWH/xTB')
    a.set_defaults(func=process_existing)
    a=sub.add_parser('visualize', help='same as process-existing (postprocess only)')
    a.set_defaults(func=visualize_cmd)
    a=sub.add_parser('resume', help='continue unfinished relax/hess/refine; does not restart GWH')
    a.add_argument('--case',required=True);a.add_argument('--pilot',action='store_true')
    a.add_argument('--xtb',default=None);a.add_argument('--workers',type=int,default=2)
    a.add_argument('--limit',type=int);a.add_argument('--imag-thresh',type=float,default=5.0)
    a.add_argument('--step-A',type=float,default=0.15)
    a.add_argument('--max-iters',type=int,default=P.DEFAULT_MAX_REFINE_ITERS)
    a.set_defaults(func=resume)
    a=sub.add_parser('run-series', help='gated GWH+opt+hess+refine for one case')
    a.add_argument('--case',required=True);a.add_argument('--pilot',action='store_true')
    a.add_argument('--nstruct',type=int,default=100);a.add_argument('--xtb',default=None)
    a.add_argument('--gwh',default=str(ROOT/'vendors/Glorious-Water-Hole/gloriouswaterhole.py'))
    a.add_argument('--workers',type=int,default=2);a.add_argument('--limit',type=int)
    a.add_argument('--imag-thresh',type=float,default=5.0);a.add_argument('--step-A',type=float,default=0.15)
    a.add_argument('--max-iters',type=int,default=P.DEFAULT_MAX_REFINE_ITERS)
    a.set_defaults(func=run_series)
    a=sub.add_parser('export-results');a.add_argument('--case',required=True);a.add_argument('--pilot',action='store_true')
    a.set_defaults(func=export_results)
    a=sub.add_parser('upload-drive');a.set_defaults(func=upload_drive)
    a=sub.add_parser('verify-delivery');a.set_defaults(func=verify_delivery)
    a=sub.add_parser('write-status');a.add_argument('--case',required=True);a.add_argument('--pilot',action='store_true')
    a.set_defaults(func=lambda args: write_status_json(args))
    for cmd,fn in [('gwh',gwh),('collect',collect),('prepare',prep),('relax',relax),('hess',hess),('refine',refine),('dedup',dedup)]:
        a=sub.add_parser(cmd);a.add_argument('--case',required=True,choices=sorted(cases()));a.add_argument('--pilot',action='store_true',help='keep small pilot in separate pilots/ folder')
        if cmd=='gwh':
            a.add_argument('--gwh',default=str(ROOT/'vendors/Glorious-Water-Hole/gloriouswaterhole.py'))
            a.add_argument('--nstruct',type=int,default=100)
            a.add_argument('--xtb',default=None)
        elif cmd=='relax':
            a.add_argument('--xtb',default=None);a.add_argument('--workers',type=int,default=2);a.add_argument('--limit',type=int)
        elif cmd=='hess':
            a.add_argument('--xtb',default=None);a.add_argument('--workers',type=int,default=1)
            a.add_argument('--limit',type=int);a.add_argument('--imag-thresh',type=float,default=5.0)
        elif cmd=='refine':
            a.add_argument('--xtb',default=None);a.add_argument('--limit',type=int)
            a.add_argument('--imag-thresh',type=float,default=5.0);a.add_argument('--step-A',type=float,default=0.15)
            a.add_argument('--max-iters',type=int,default=P.DEFAULT_MAX_REFINE_ITERS)
        elif cmd=='dedup':
            a.add_argument('--cocos',help='optional current source cocos.py; default bundled submitted snapshot')
            a.add_argument('--mdatom',help='path to installed MDAtomLocalEnvironment atom_local_environment_computer.py')
            a.add_argument('--seed',type=int,default=20261030);a.add_argument('--samples',type=int,default=50)
            a.add_argument('--imag-thresh',type=float,default=5.0);a.add_argument('--step-A',type=float,default=0.15)
            a.add_argument('--max-iters',type=int,default=P.DEFAULT_MAX_REFINE_ITERS)
            a.add_argument('--skip-robustness',action='store_true')
        a.set_defaults(func=fn)
    x=p.parse_args();x.func(x)

if __name__=='__main__':main()
