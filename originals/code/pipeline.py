#!/usr/bin/env python3
"""Reproducible Rb+(H2O)4 / ionic hydrate workflow.

This code DOES NOT approximate GFN2-xTB energies: genuine new energies require xtb.
It deliberately preserves all GWH constrained structures until after free relaxation.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment

HARTREE_TO_KCAL_MOL = 627.5094740631
ENERGY_RE = re.compile(r'\bE\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*Eh')
TOTAL_RE = re.compile(r'TOTAL ENERGY\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)', re.I)
OPT_RE = re.compile(r'GEOMETRY OPTIMIZATION CONVERGED', re.I)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_xyz(path):
    """Strict multi-frame XYZ reader with atom order and size retained."""
    p = Path(path)
    lines = p.read_text(encoding='utf-8').splitlines()
    i=0; frames=[]
    while i < len(lines):
        if not lines[i].strip():
            i+=1; continue
        try: n=int(lines[i].strip())
        except ValueError as ex: raise ValueError(f'{p}: bad atom count at line {i+1}') from ex
        if n<=0 or i+2+n>len(lines):raise ValueError(f'{p}: incomplete XYZ at {i+1}')
        comment=lines[i+1]; symbols=[]; xyz=[]
        for line in lines[i+2:i+2+n]:
            toks=line.split()
            if len(toks)<4:raise ValueError(f'{p}: invalid atom row: {line}')
            symbols.append(toks[0]);xyz.append([float(k) for k in toks[1:4]])
        frames.append({'symbols':symbols,'xyz':np.array(xyz,float),'comment':comment})
        i+=n+2
    if not frames: raise ValueError(f'No structures in {p}')
    return frames


def write_xyz(path, frames):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8') as f:
        for v in frames:
            f.write(f'{len(v["symbols"])}\n{v.get("comment", "")}\n')
            for symbol,r in zip(v['symbols'],v['xyz']):
                f.write(f'{symbol:2s} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}\n')


def energy_in_comment(comment):
    m=ENERGY_RE.search(comment)
    return float(m.group(1)) if m else None


def geometry_features(frame, central, shell_cut=None):
    symbols=np.array(frame['symbols']); points=frame['xyz']; center_ids=np.where(symbols==central)[0]
    if len(center_ids)!=1:raise ValueError(f'Expected one central {central}, got {len(center_ids)}')
    ci=center_ids[0]; center=points[ci]; o_ids=np.where(symbols=='O')[0]; h_ids=np.where(symbols=='H')[0]
    o_xyz=points[o_ids]; h_xyz=points[h_ids]
    mo=np.linalg.norm(o_xyz-center,axis=1)
    data={'n_water':len(o_ids),'n_hydrogen':len(h_ids),'distance_min_A':float(mo.min()),'distance_mean_A':float(mo.mean()),'distance_max_A':float(mo.max()),'distance_std_A':float(mo.std())}
    if shell_cut is not None:data['CN_at_cut']=int(np.sum(mo<shell_cut))
    if len(o_ids)*2 != len(h_ids):
        data['water_connectivity']='not_2H_per_O';data['orientation_cos_mean']=None
        return data
    # Stable globally minimum-distance assignment to two H per O; flag O-H dissociation.
    distance=np.linalg.norm(h_xyz[:,None,:]-np.repeat(o_xyz,2,axis=0)[None,:,:],axis=2)
    rows,cols=linear_sum_assignment(distance)
    if np.max(distance[rows,cols])>1.25:
        data['water_connectivity']='O_H_distance_above_1.25_A';data['orientation_cos_mean']=None
        return data
    pairs=[[] for _ in o_ids]
    for h_index,slot in zip(rows,cols):pairs[slot//2].append(h_xyz[h_index])
    cos=[]
    for oxygen,hpair in zip(o_xyz,pairs):
        mu=np.mean(hpair,axis=0)-oxygen
        radial=oxygen-center
        cos.append(float(np.dot(mu,radial)/(np.linalg.norm(mu)*np.linalg.norm(radial))))
    data['water_connectivity']='nominal_2H_per_O'
    data['orientation_cos_mean']=float(np.mean(cos))
    data['orientation_cos_min']=float(np.min(cos))
    data['orientation_cos_max']=float(np.max(cos))
    data['orientation_cos_values']=';'.join(f'{v:.8f}' for v in cos)
    return data


def write_csv(path,rows):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        w.writeheader();w.writerows(rows)


def metadata_rows(path, frames):
    if not path: return [None for _ in frames]
    with open(path, newline='', encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=len(frames):
        raise ValueError(f'Metadata rows {len(rows)} do not match XYZ frames {len(frames)}')
    for i,(r,frame) in enumerate(zip(rows,frames)):
        e=energy_in_comment(frame['comment'])
        if e is None or not np.isclose(float(r['E_final_Eh']),e,rtol=0,atol=1e-9):
            raise ValueError(f'energy mismatch at frame {i}: {r["E_final_Eh"]} vs {e}')
        if r.get('success','').lower() not in ('true','1','yes'):
            raise ValueError(f'source stage marked unsuccessful: frame {i}')
    return rows

def audit(args):
    frames=read_xyz(args.input)
    for fr in frames:
        if fr['symbols']!=frames[0]['symbols']:
            raise ValueError('Inconsistent atom order/formula across frames')
    metadata=metadata_rows(args.energies_csv,frames)
    energies=[energy_in_comment(v['comment']) for v in frames]
    if any(e is None for e in energies):
        raise ValueError('At least one input frame has no E = ... Eh energy in comment')
    min_e=min(energies)
    rows=[]
    for i,(fr,e) in enumerate(zip(frames,energies)):
        row={'frame_index':i,'id':f'frame_{i:03d}','energy_constrained_Eh':e,
             'delta_E_constrained_kcal_mol':(e-min_e)*HARTREE_TO_KCAL_MOL}
        if metadata[i]:
            row.update({'source_run_id':metadata[i]['id'],'source_batch':metadata[i]['batch'],
                        'source_structure_index':metadata[i]['structure'],'source_nfev':metadata[i]['nfev'],
                        'source_success':metadata[i]['success']})
        row.update(geometry_features(fr,args.central,args.shell_cut))
        rows.append(row)
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'constrained_audit.csv',rows)
    summary={'source':str(args.input),'sha256':sha256(args.input),'central':args.central,
             'frames':len(frames),'atoms_per_frame':len(frames[0]['symbols']),
             'formula_counts':{s:frames[0]['symbols'].count(s) for s in sorted(set(frames[0]['symbols']))},
             'lowest_constrained_energy_Eh':min_e,
             'highest_constrained_energy_Eh':max(energies),
             'energy_span_kcal_mol':(max(energies)-min_e)*HARTREE_TO_KCAL_MOL,
             'distance_min_A':min(x['distance_min_A'] for x in rows),
             'distance_max_A':max(x['distance_max_A'] for x in rows),
             'shell_cut_A':args.shell_cut,
             'status':'CONSTRAINED GWH DATA ONLY; no free xTB optimizations inferred'}
    (out/'constrained_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8.4,4.6))
    ax.hist([r['delta_E_constrained_kcal_mol'] for r in rows],bins=20,edgecolor='white')
    ax.set(xlabel='Constrained relative energy, kcal/mol',ylabel='Number of structures',title=f'{args.central}+ hydrated clusters: constrained GWH stage')
    ax.grid(axis='y',alpha=.18)
    fig.tight_layout();fig.savefig(out/'constrained_energy_distribution.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.4,4.6))
    x=[r['orientation_cos_mean'] for r in rows if r['orientation_cos_mean'] is not None]
    ax.hist(x,bins=18,edgecolor='white')
    ax.set(xlabel='Average cos(angle of O→H-midpoint vs ion→O)',ylabel='Number of structures',title='Water orientation in constrained geometries')
    ax.grid(axis='y',alpha=.18)
    fig.tight_layout();fig.savefig(out/'constrained_water_orientations.png',dpi=170);plt.close(fig)
    print(json.dumps(summary,indent=2,ensure_ascii=False))


def prepare(args):
    frames=read_xyz(args.input);metadata=metadata_rows(args.energies_csv,frames);out=Path(args.output)
    if out.exists() and any(out.iterdir()) and not args.allow_existing:
        raise RuntimeError('Output directory not empty; refusing to overwrite. Use --allow-existing to reuse existing prepared inputs.')
    out.mkdir(parents=True,exist_ok=True)
    entries=[]
    for i,fr in enumerate(frames):
        folder=out/f'frame_{i:03d}';folder.mkdir(exist_ok=True)
        inp=folder/'input.xyz'
        if inp.exists():
            if not args.allow_existing:raise FileExistsError(inp)
            if inp.read_text()!=render_xyz(fr):raise ValueError(f'Existing input differs: {inp}')
        else: inp.write_text(render_xyz(fr),encoding='utf-8')
        entries.append({'id':f'frame_{i:03d}','index':i,'constrained_energy_Eh':energy_in_comment(fr['comment']),
                        'source_run_id':metadata[i]['id'] if metadata[i] else None,
                        'input_sha256':sha256(inp)})
    manifest={'source_sha256':sha256(args.input),'source_energy_csv_sha256':sha256(args.energies_csv) if args.energies_csv else None,
              'source':os.path.relpath(Path(args.input).resolve(),out.resolve()),'charge':args.charge,
              'method':'GFN2-xTB','solvent':'ALPB(water)','optimization':'unconstrained tight',
              'central':args.central,'input_count':len(frames),'entries':entries}
    man=out/'manifest.json'
    if man.exists() and json.loads(man.read_text())!=manifest:raise RuntimeError('Existing manifest differs, refusing overwrite')
    man.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Prepared',len(frames),'folders in',out)


def render_xyz(fr):
    lines=[str(len(fr['symbols'])),fr.get('comment','')]
    lines += [f'{s:2s} {x:.10f} {y:.10f} {z:.10f}' for s,(x,y,z) in zip(fr['symbols'],fr['xyz'])]
    return '\n'.join(lines)+'\n'


def get_result(folder,expected_symbols):
    opt=folder/'xtbopt.xyz';log=folder/'xtb.out'
    if not opt.exists() or not log.exists():return None
    data=log.read_text(errors='replace')
    if not OPT_RE.search(data) or 'normal termination of xtb' not in data.lower():return None
    try:
        fr=read_xyz(opt)
        if len(fr)!=1 or fr[0]['symbols']!=expected_symbols:return None
    except Exception:return None
    m=TOTAL_RE.findall(data)
    # TOTAL ENERGY at end of successful run, numeric Hartree
    e=float(m[-1]) if m else None
    if e is None:
        # optional xtbopt comment may be formatted like 'energy: -123...'
        mm=re.search(r'energy\s*:\s*([-+\d.eE]+)',fr[0]['comment'],re.I)
        if mm:e=float(mm.group(1))
    if e is None or not math.isfinite(e):return None
    return fr[0],e


def relax(args):
    out=Path(args.output); manifest=json.loads((out/'manifest.json').read_text())
    if shutil.which(args.xtb) is None and not Path(args.xtb).exists():
        raise RuntimeError('xTB executable unavailable; install xTB first. No simulated energy will be generated.')
    frames=read_xyz(out/'frame_000'/'input.xyz')
    expected=frames[0]['symbols']
    ids=[e['id'] for e in manifest['entries']]
    if args.limit:ids=ids[:args.limit]
    env=os.environ.copy();env.update({'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
    def single(id):
        folder=out/id;existing=get_result(folder,expected)
        if existing is not None:return (id,'skipped_success',existing[1])
        cmd=[args.xtb,'input.xyz','--gfn','2','--chrg',str(manifest['charge']),'--alpb','water','--opt','tight']
        with (folder/'xtb.out').open('w') as log:
            proc=subprocess.run(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
        res=get_result(folder,expected)
        status='success' if proc.returncode==0 and res else f'FAILED_exit_{proc.returncode}'
        (folder/'execution.json').write_text(json.dumps({'cmd':cmd,'returncode':proc.returncode,'status':status},indent=2))
        return (id,status,res[1] if res else None)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(single,id) for id in ids]
        for future in as_completed(jobs):
            id,status,e=future.result(); print(id,status, e)
    print('Run complete; execute the dedup command after verifying logs.')


def dedup(args):
    out=Path(args.output); manifest=json.loads((out/'manifest.json').read_text())
    expected=read_xyz(out/'frame_000'/'input.xyz')[0]['symbols']
    good=[]; bad=[]
    for rec in manifest['entries']:
        folder=out/rec['id']; data=get_result(folder,expected)
        if not data:bad.append(rec['id']);continue
        frame,e=data
        good.append({'id':rec['id'],'source_run_id':rec.get('source_run_id'),'frame':frame,'energy_Eh':e,
                     'constrained_Eh':rec['constrained_energy_Eh']})
    if not good:raise RuntimeError('No converged xtb results; cannot run CoCoS or report minimum energies')
    post=out/'postprocess'; post.mkdir(exist_ok=True)
    (post/'missing_failed_ids.txt').write_text('\n'.join(bad)+'\n' if bad else '(none)\n')
    if bad and not args.allow_incomplete:
        raise RuntimeError(
            f'{len(bad)} of {len(manifest["entries"])} inputs are not verified as converged. '
            'CoCoS and energy-window selection blocked to prevent misleading final results. '
            'Inspect xtb_runs/postprocess/missing_failed_ids.txt, rerun failures; '
            'for diagnostics only, use --allow-incomplete.'
        )
    good.sort(key=lambda r:(r['energy_Eh'],r['id']))
    post=out/'postprocess';post.mkdir(exist_ok=True)
    ener0=good[0]['energy_Eh']; rows=[]
    for j,r in enumerate(good):
        rows.append({'sorted_index':j,'source_id':r['id'],'source_run_id':r.get('source_run_id'),'energy_free_Eh':r['energy_Eh'],
                     'delta_free_kcal_mol':(r['energy_Eh']-ener0)*HARTREE_TO_KCAL_MOL,
                     'energy_constrained_Eh':r['constrained_Eh'],
                     'geometry_change_A_unaligned':None})
    write_csv(post/'all_freely_optimized.csv',rows)
    # Scientific diagnostic: ranks are computed separately within this fixed stoichiometry.
    # This is *not* thermodynamic free energy; the variable energy_Eh is xTB total electronic energy.
    rank_constrained = {r['id']: i+1 for i,r in enumerate(sorted(good, key=lambda r:(r['constrained_Eh'],r['id'])))}
    min_constrained = min(r['constrained_Eh'] for r in good)
    rank_rows=[]
    for j,r in enumerate(good):
        start_rank=rank_constrained[r['id']]
        after_rank=j+1
        rank_rows.append({
            'source_id':r['id'], 'source_run_id':r.get('source_run_id'),
            'constrained_rank':start_rank, 'unconstrained_rank':after_rank,
            'rank_shift_after_minus_before':after_rank-start_rank,
            'constrained_Eh':r['constrained_Eh'], 'unconstrained_Eh':r['energy_Eh'],
            'delta_constrained_kcal_mol':(r['constrained_Eh']-min_constrained)*HARTREE_TO_KCAL_MOL,
            'delta_unconstrained_kcal_mol':(r['energy_Eh']-ener0)*HARTREE_TO_KCAL_MOL,
        })
    write_csv(post/'rank_shift_constrained_vs_unconstrained.csv',rank_rows)
    if len(good)>1:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6.2,5.8))
        ax.scatter([r['constrained_rank'] for r in rank_rows],
                   [r['unconstrained_rank'] for r in rank_rows],s=24,alpha=.75)
        ax.plot([1,len(good)],[1,len(good)],linestyle='--',linewidth=1,color='gray')
        ax.set(xlabel='Rank after constrained optimization',
               ylabel='Rank after unconstrained optimization',
               title='Energy-rank stability upon unconstrained relaxation')
        ax.grid(alpha=.15)
        fig.tight_layout();fig.savefig(post/'rank_shift_after_relaxation.png',dpi=180);plt.close(fig)
    (post/'COMPLETENESS.txt').write_text(
        f'Verified unconstrained optimizations: {len(good)}/{len(manifest["entries"])}\n'
        + ('INCOMPLETE: diagnostics only.\n' if bad else 'COMPLETE: all requested inputs checked.\n')
        + 'These energies are calculated total energies, NOT Gibbs free energies.\n'
    )
    tosave=[{'symbols':r['frame']['symbols'],'xyz':r['frame']['xyz'],
             'comment':f'id={r["id"]} E = {r["energy_Eh"]:.14f} Eh'} for r in good]
    write_xyz(post/'free_sorted_for_cocos.xyz',tosave)
    print(f'Converged = {len(good)}, not verified = {len(bad)}')
    if not args.cocos: return
    code=Path(args.cocos).resolve(); assert code.exists()
    wrapper=Path(__file__).with_name('seeded_cocos.py').resolve();assert wrapper.exists()
    cmd=[sys.executable,str(wrapper),str(args.seed),str(code),'-f','free_sorted_for_cocos.xyz','--DoNotUseGlobalOpt',
         '--NumSamples',str(args.samples),'--DoNotSaveDoubles']
    proc=subprocess.run(cmd,cwd=post,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=False)
    (post/'cocos_stdout.log').write_text(proc.stdout)
    if proc.returncode!=0:raise RuntimeError('CoCoS failed; see cocos_stdout.log')
    survivors=[int(x) for x in re.findall(r'added structure (\d+)\b',proc.stdout)]
    observed=list(post.glob('*_quick_cleanup.trj'))
    if len(observed)!=1:raise RuntimeError('Ambiguous/missing CoCoS output')
    kept=read_xyz(observed[0]);survivors=list(dict.fromkeys(survivors))
    if len(kept)!=len(survivors):
        raise RuntimeError(f'Cannot map CoCoS output: {len(kept)} frames but {len(survivors)} log indices')
    retained=[good[i] for i in survivors]
    write_xyz(post/'unique_lowest_energy_representatives.xyz',
              [{'symbols':r['frame']['symbols'],'xyz':r['frame']['xyz'],
                'comment':f'id={r["id"]} E = {r["energy_Eh"]:.14f} Eh'} for r in retained])
    write_csv(post/'unique_representatives.csv',[
        {'source_id':r['id'],'source_run_id':r.get('source_run_id'),'energy_free_Eh':r['energy_Eh'],
         'delta_free_kcal_mol':(r['energy_Eh']-ener0)*HARTREE_TO_KCAL_MOL,
         **geometry_features(r['frame'],manifest['central'],args.shell_cut)} for r in retained])
    (post/'cocos_parameters.json').write_text(json.dumps({'seed':args.seed,'samples':args.samples,
      'cocos_sha256':sha256(code),'cocos_cmd':cmd,'comment':
      'One representative per comparison class, sorted by free energy before CoCoS. Sensitivity to criterion and sampling must be tested.'},indent=2))
    print('CoCoS classes:',len(retained),'lowest-energy representatives saved; verify sensitivity before publication')
    if args.mdatom:
        mdscript=Path(args.mdatom).resolve()
        if not mdscript.exists():raise FileNotFoundError(mdscript)
        # -e (not -n) is important for the author's current implementation,
        # since its analysis functions use a global loop variable named idx.
        mdcmd=[sys.executable,str(mdscript),'-t','unique_lowest_energy_representatives.xyz',
               '-e',manifest['central'],'--RunWaterAnalysis','--RecalcWater']
        mdproc=subprocess.run(mdcmd,cwd=post,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (post/'MDAtomLocalEnvironment.log').write_text(mdproc.stdout)
        (post/'MDAtomLocalEnvironment_command.json').write_text(json.dumps(mdcmd,indent=2))
        if mdproc.returncode !=0:
            raise RuntimeError('MDAtomLocalEnvironment failed: see log')
        print('MDAtomLocalEnvironment completed; interpret its time column as frame index, not physical simulation time.')


def main():
    pa=argparse.ArgumentParser(description=__doc__);subs=pa.add_subparsers(dest='command',required=True)
    p=subs.add_parser('audit');p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--central',default='Rb');p.add_argument('--shell-cut',type=float);p.add_argument('--energies-csv');p.set_defaults(func=audit)
    p=subs.add_parser('prepare');p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--central',default='Rb');p.add_argument('--charge',type=int,default=1)
    p.add_argument('--allow-existing',action='store_true');p.add_argument('--energies-csv');p.set_defaults(func=prepare)
    p=subs.add_parser('relax');p.add_argument('--output',required=True);p.add_argument('--xtb',default='xtb')
    p.add_argument('--workers',type=int,default=2);p.add_argument('--limit',type=int);p.set_defaults(func=relax)
    p=subs.add_parser('dedup');p.add_argument('--output',required=True);p.add_argument('--cocos')
    p.add_argument('--seed',type=int,default=20261029);p.add_argument('--samples',type=int,default=50)
    p.add_argument('--shell-cut',type=float);p.add_argument('--mdatom')
    p.add_argument('--allow-incomplete', action='store_true', help='Diagnostics only: do not treat partial sample as final data')
    p.set_defaults(func=dedup)
    args=pa.parse_args();args.func(args)

if __name__=='__main__':main()
