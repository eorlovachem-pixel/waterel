#!/usr/bin/env python3
"""Reproducible Li+/Na+/K+/Rb+ and F-/Cl-/Br- hydrate workflow.

This code DOES NOT approximate GFN2-xTB energies: genuine new energies require xtb.
It deliberately preserves all GWH constrained structures until after free relaxation.
CoCoS receives only confirmed minima (all ν ≥ −5 cm⁻¹).
Cations use --chrg 1; anions use --chrg -1. CASE_STATUS still gates production launches.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

try:
    import qc_checks
    import freq_parser
    import mode_follow
    import shell_descriptors
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import qc_checks
    import freq_parser
    import mode_follow
    import shell_descriptors

HARTREE_TO_KCAL_MOL = 627.5094740631
SUPPORTED_IONS = qc_checks.SUPPORTED_IONS
PROTOCOL_CHARGE = qc_checks.PROTOCOL_CHARGE
NU_MIN_CM1 = -5.0
MODE_STEP_A = 0.15
DEFAULT_MAX_REFINE_ITERS = 3
DEFAULT_COCOS_SEED = 20261030
DEFAULT_COCOS_SAMPLES = 50
DEFAULT_ROBUST_SEEDS = (20261030, 101, 202, 303, 7)
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
    desc = shell_descriptors.cluster_descriptors(frame, central, shell_cut=shell_cut)
    keep = (
        'n_water', 'n_hydrogen', 'distance_min_A', 'distance_mean_A', 'distance_max_A',
        'distance_std_A', 'CN_at_cut', 'water_connectivity', 'orientation_cos_mean',
        'orientation_cos_min', 'orientation_cos_max', 'orientation_cos_values',
        'CN', 'R_mean_A', 'R_min_A', 'R_max_A', 'R_std_A',
        'radial_deformation_std_over_mean', 'P2_mean', 'n_OH_O_hbonds',
        'n_XH_O_contacts', 'proton_transfer_candidate', 'full_first_shell',
        'CN_cut_A', 'CN_cut_status', 'xyz_sha256',
    )
    return {k: desc[k] for k in keep if k in desc}


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


def require_xtb(xtb):
    try:
        return str(qc_checks.resolve_working_xtb(xtb))
    except RuntimeError as exc:
        raise RuntimeError(f'{exc} No simulated energy will be generated.') from exc


def xtb_thread_env():
    env=os.environ.copy()
    env.update({'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
    return env


def stamp_if_exists(path):
    p=Path(path)
    if p.exists():
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        return p.with_name(p.stem + '_' + stamp + p.suffix)
    return p


def refuse_overwrite(path, allow=False):
    p=Path(path)
    if p.exists() and not allow:
        raise RuntimeError(f'Refusing overwrite of existing file: {p}')
    return p


def prepare(args):
    qc_checks.assert_supported_ion(args.central)
    expected = qc_checks.protocol_charge(args.central)
    if int(args.charge) != expected:
        raise ValueError(f'protocol charge for {args.central} is {expected:+d}, got {args.charge}')
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
    try:
        frames=read_xyz(opt)
        if len(frames)!=1:return None
        e,_info=qc_checks.require_completed_opt(data, frames[0], expected_symbols)
    except Exception:
        return None
    return frames[0],e


def relax(args):
    out=Path(args.output); manifest=json.loads((out/'manifest.json').read_text())
    args.xtb=require_xtb(args.xtb)
    print('Using xTB', args.xtb)
    frames=read_xyz(out/'frame_000'/'input.xyz')
    expected=frames[0]['symbols']
    ids=[e['id'] for e in manifest['entries']]
    if args.limit:ids=ids[:args.limit]
    env=xtb_thread_env()
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


def hess(args):
    """Run xTB Hessian on already freely optimized geometries. Does not overwrite xtb.out."""
    out=Path(args.output); manifest=json.loads((out/'manifest.json').read_text())
    args.xtb=require_xtb(args.xtb)
    print('Using xTB', args.xtb)
    expected=read_xyz(out/'frame_000'/'input.xyz')[0]['symbols']
    ids=[e['id'] for e in manifest['entries']]
    if args.limit:ids=ids[:args.limit]
    env=xtb_thread_env()
    rows=[]
    for id in ids:
        folder=out/id
        opt_res=get_result(folder,expected)
        hdir=folder/'hess'
        hlog=hdir/'xtb_hess.out'
        if opt_res is None:
            rows.append({'id':id,'status':'OPT_NOT_VERIFIED','is_genuine_minimum_candidate':False})
            continue
        if hlog.is_file() and 'normal termination of xtb' in hlog.read_text(errors='replace').lower():
            try:
                parsed=freq_parser.require_completed_hessian(hlog, imag_thresh_cm1=args.imag_thresh)
                rows.append({'id':id,'status':'skipped_success','classification':parsed['classification'],
                             'n_frequencies_parsed':parsed['n_frequencies_parsed'],
                             'imaginary_cm1':';'.join(f'{x:.4f}' for x in parsed['imaginary_cm1']),
                             'is_genuine_minimum_candidate':parsed['is_genuine_minimum_candidate']})
                continue
            except Exception:
                pass
        hdir.mkdir(exist_ok=True)
        shutil.copy2(folder/'xtbopt.xyz', hdir/'xtbopt.xyz')
        cmd=[args.xtb,'xtbopt.xyz','--gfn','2','--chrg',str(manifest['charge']),'--alpb','water','--hess']
        with hlog.open('w') as log:
            proc=subprocess.run(cmd,cwd=hdir,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
        try:
            parsed=freq_parser.require_completed_hessian(hlog, imag_thresh_cm1=args.imag_thresh)
            status='success' if proc.returncode==0 else f'FAILED_exit_{proc.returncode}'
            rows.append({'id':id,'status':status,'classification':parsed['classification'],
                         'n_frequencies_parsed':parsed['n_frequencies_parsed'],
                         'imaginary_cm1':';'.join(f'{x:.4f}' for x in parsed['imaginary_cm1']),
                         'is_genuine_minimum_candidate':parsed['is_genuine_minimum_candidate']})
        except Exception as exc:
            rows.append({'id':id,'status':f'PARSE_FAIL_{proc.returncode}','error':str(exc),
                         'is_genuine_minimum_candidate':False})
        print(id, rows[-1]['status'], rows[-1].get('classification'))
    post=out/'postprocess'; post.mkdir(exist_ok=True)
    dest=post/'hessian_status.csv'
    if dest.exists():
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        dest=post/f'hessian_status_{stamp}.csv'
    write_csv(dest, rows)
    print('Wrote', dest)


def _branch_dir(folder, name, step_A):
    return folder/'mode_follow'/f'{name}_{float(step_A):.2f}A'


def _read_hess(hlog, imag_thresh):
    return freq_parser.require_completed_hessian(hlog, imag_thresh_cm1=imag_thresh)


def collect_confirmed_minima(out, imag_thresh=5.0, step_A=MODE_STEP_A,
                             max_iters=DEFAULT_MAX_REFINE_ITERS):
    """Collect geometries with all ν ≥ −imag_thresh cm⁻¹. Does not invent frequencies.

    Walks mode-follow ± branches up to max_iters. Saddles never go to CoCoS.
    """
    out=Path(out); manifest=json.loads((out/'manifest.json').read_text())
    expected=read_xyz(out/'frame_000'/'input.xyz')[0]['symbols']
    confirmed=[]; notes=[]
    for rec in manifest['entries']:
        _collect_from_folder(
            folder=out/rec['id'], rec=rec, node_id=rec['id'], origin='free_opt',
            is_root=True, depth=0, max_iters=int(max_iters),
            imag_thresh=imag_thresh, step_A=step_A, expected=expected,
            confirmed=confirmed, notes=notes,
        )
    confirmed.sort(key=lambda r:(r['energy_Eh'], r['id']))
    return confirmed, notes, expected, manifest


def _collect_from_folder(folder, rec, node_id, origin, is_root, depth, max_iters,
                         imag_thresh, step_A, expected, confirmed, notes):
    opt_dir = folder if is_root else folder/'opt'
    opt=get_result(opt_dir, expected)
    if opt is None:
        notes.append({'id':node_id,'status':'OPT_NOT_VERIFIED' if is_root else 'REFINE_INCOMPLETE',
                      'depth':depth})
        return
    frame,e=opt
    hlog=folder/'hess'/'xtb_hess.out'
    if not hlog.is_file():
        notes.append({'id':node_id,'status':'HESS_MISSING' if is_root else 'REFINE_INCOMPLETE',
                      'depth':depth})
        return
    try:
        parsed=_read_hess(hlog, imag_thresh)
    except Exception as exc:
        notes.append({'id':node_id,'status':'HESS_PARSE' if is_root else 'REFINE_HESS_PARSE',
                      'error':str(exc),'depth':depth})
        return
    nu_min=-abs(float(imag_thresh))
    if freq_parser.is_confirmed_minimum(parsed, nu_min_cm1=nu_min):
        confirmed.append({
            'id':node_id,'parent_id':rec['id'],'origin':origin,
            'frame':frame,'energy_Eh':e,
            'constrained_Eh':rec.get('constrained_energy_Eh'),
            'source_run_id':rec.get('source_run_id'),
            'n_imag':len(parsed['imaginary_cm1']),
            'lowest_freq_cm1':min(parsed['frequencies_cm1']),
            'xyz':str(opt_dir/'xtbopt.xyz'),'hess_log':str(hlog),
            'refine_depth':depth,
        })
        notes.append({'id':node_id,
                      'status':'CONFIRMED_FREE_OPT' if is_root else 'CONFIRMED_REFINE',
                      'depth':depth})
        return
    notes.append({'id':node_id,'status':'SADDLE','imaginary_cm1':parsed['imaginary_cm1'],
                  'depth':depth})
    if depth >= int(max_iters):
        notes.append({'id':node_id,'status':'MAX_ITERS_STILL_SADDLE','depth':depth})
        return
    for name in ('plus','minus'):
        _collect_from_folder(
            folder=_branch_dir(folder, name, step_A), rec=rec,
            node_id=f'{node_id}_{name}',
            origin=(f'mode_follow_{name}' if origin=='free_opt' else f'{origin}_{name}'),
            is_root=False, depth=depth+1, max_iters=max_iters,
            imag_thresh=imag_thresh, step_A=step_A, expected=expected,
            confirmed=confirmed, notes=notes,
        )


def _ensure_displacements(folder, imag_thresh, step_A, xyz_path=None, g98_path=None):
    mfdir=folder/'mode_follow'
    plus=mfdir/f'displaced_plus_{float(step_A):.2f}A.xyz'
    minus=mfdir/f'displaced_minus_{float(step_A):.2f}A.xyz'
    report=mfdir/'mode_follow.json'
    if plus.is_file() and minus.is_file() and report.is_file():
        return json.loads(report.read_text(encoding='utf-8')), 'skipped_existing'
    xyz_src=Path(xyz_path) if xyz_path else folder/'xtbopt.xyz'
    g98=Path(g98_path) if g98_path else folder/'hess'/'g98.out'
    if not g98.is_file():
        return None, 'NO_G98'
    if mfdir.exists() and any(mfdir.iterdir()):
        raise RuntimeError(f'Refusing to mix new mode-follow files into non-empty {mfdir}')
    mfdir.mkdir(parents=True, exist_ok=True)
    frame=mode_follow.read_xyz(xyz_src)
    modes=mode_follow.parse_g98_modes(g98)
    chosen=mode_follow.choose_imaginary_g98_mode(modes, imag_thresh)
    if chosen is None:
        return None, 'NO_IMAG_MODE_IN_G98'
    scaled, meta=mode_follow.scale_mode_max_atom(chosen['cartesian'], step_A)
    payload={
        'xyz':str(xyz_src.resolve()),'g98':str(g98.resolve()),
        'chosen_g98_frequency_cm1':chosen['frequency_cm1'],'scale':meta,'displaced':True,'files':{},
    }
    for sign, name in ((+1.0,'plus'),(-1.0,'minus')):
        tagged=mode_follow.displace(frame, scaled, sign,
            comment=f"displaced {name} along imag mode nu={chosen['frequency_cm1']:.4f} cm-1 max-atom {step_A} A")
        out_xyz=mfdir/f'displaced_{name}_{float(step_A):.2f}A.xyz'
        mode_follow.write_xyz(out_xyz, tagged)
        payload['files'][name]=str(out_xyz)
    report.write_text(json.dumps(payload, indent=2, default=str), encoding='utf-8')
    return payload, 'wrote'


def _run_opt_hess(folder, xyz_src, xtb, charge, expected, imag_thresh, env):
    opt=folder/'opt'; opt.mkdir(parents=True, exist_ok=True)
    inp=opt/'input.xyz'
    src_text=Path(xyz_src).read_text(encoding='utf-8')
    if inp.exists() and inp.read_text(encoding='utf-8')!=src_text:
        raise RuntimeError(f'Existing {inp} differs from {xyz_src}; refusing overwrite')
    if not inp.exists():
        inp.write_text(src_text, encoding='utf-8')
    existing=get_result(opt, expected)
    if existing is None:
        cmd=[xtb,'input.xyz','--gfn','2','--chrg',str(charge),'--alpb','water','--opt','tight']
        with (opt/'xtb.out').open('w') as log:
            proc=subprocess.run(cmd,cwd=opt,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
        (opt/'execution.json').write_text(json.dumps({'cmd':cmd,'returncode':proc.returncode},indent=2))
        existing=get_result(opt, expected)
        if existing is None:
            return {'status':f'OPT_FAILED_exit_{proc.returncode}'}
        opt_status='success'
    else:
        opt_status='skipped_success'
    hdir=folder/'hess'; hlog=hdir/'xtb_hess.out'
    if hlog.is_file() and 'normal termination of xtb' in hlog.read_text(errors='replace').lower():
        try:
            parsed=_read_hess(hlog, imag_thresh)
            return {'status':'skipped_success','opt':opt_status,'hess':parsed,
                    'is_confirmed':freq_parser.is_confirmed_minimum(parsed, nu_min_cm1=-abs(imag_thresh))}
        except Exception:
            pass
    hdir.mkdir(exist_ok=True)
    shutil.copy2(opt/'xtbopt.xyz', hdir/'xtbopt.xyz')
    cmd=[xtb,'xtbopt.xyz','--gfn','2','--chrg',str(charge),'--alpb','water','--hess']
    with hlog.open('w') as log:
        proc=subprocess.run(cmd,cwd=hdir,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
    try:
        parsed=_read_hess(hlog, imag_thresh)
        return {'status':'success' if proc.returncode==0 else f'FAILED_exit_{proc.returncode}',
                'opt':opt_status,'hess':parsed,
                'is_confirmed':freq_parser.is_confirmed_minimum(parsed, nu_min_cm1=-abs(imag_thresh))}
    except Exception as exc:
        return {'status':f'HESS_PARSE_{proc.returncode}','error':str(exc),'is_confirmed':False}


def _refine_from_folder(folder, node_id, is_root, depth, max_iters, xtb, charge,
                        expected, imag_thresh, step_A, env, rows):
    opt_dir = folder if is_root else folder/'opt'
    opt=get_result(opt_dir, expected)
    if opt is None:
        rows.append({'id':node_id,'status':'OPT_NOT_VERIFIED' if is_root else 'REFINE_INCOMPLETE',
                     'depth':depth})
        return
    hlog=folder/'hess'/'xtb_hess.out'
    if not hlog.is_file():
        rows.append({'id':node_id,'status':'HESS_MISSING' if is_root else 'REFINE_INCOMPLETE',
                     'depth':depth})
        return
    try:
        parsed=_read_hess(hlog, imag_thresh)
    except Exception as exc:
        rows.append({'id':node_id,'status':'HESS_PARSE','error':str(exc),'depth':depth})
        return
    if freq_parser.is_confirmed_minimum(parsed, nu_min_cm1=-abs(imag_thresh)):
        rows.append({'id':node_id,'status':'ALREADY_MINIMUM','n_imag':0,'depth':depth,
                     'is_confirmed':True})
        return
    if depth >= int(max_iters):
        rows.append({'id':node_id,'status':'MAX_ITERS_STILL_SADDLE','depth':depth,
                     'is_confirmed':False})
        print(node_id, 'MAX_ITERS_STILL_SADDLE', depth)
        return
    try:
        _rep, disp_status=_ensure_displacements(
            folder, imag_thresh, step_A,
            xyz_path=opt_dir/'xtbopt.xyz', g98_path=folder/'hess'/'g98.out',
        )
    except Exception as disp_exc:
        rows.append({'id':node_id,'status':'DISPLACE_REFUSED','error':str(disp_exc),'depth':depth})
        return
    if disp_status in ('NO_G98','NO_IMAG_MODE_IN_G98'):
        rows.append({'id':node_id,'status':disp_status,'depth':depth}); return
    for name in ('plus','minus'):
        xyz_in=folder/'mode_follow'/f'displaced_{name}_{float(step_A):.2f}A.xyz'
        bid=f'{node_id}_{name}'
        if not xyz_in.is_file():
            rows.append({'id':bid,'status':'DISPLACE_MISSING','depth':depth+1}); continue
        try:
            result=_run_opt_hess(_branch_dir(folder,name,step_A), xyz_in, xtb,
                                 charge, expected, imag_thresh, env)
        except Exception as run_exc:
            rows.append({'id':bid,'status':'REFINE_ERROR','error':str(run_exc),'depth':depth+1})
            continue
        rows.append({'id':bid,'status':result.get('status'),
                     'is_confirmed':result.get('is_confirmed', False),
                     'classification':(result.get('hess') or {}).get('classification'),
                     'depth':depth+1})
        print(bid, rows[-1]['status'], rows[-1].get('is_confirmed'))
        if result.get('is_confirmed'):
            continue
        status=str(result.get('status') or '')
        if status.startswith('OPT_FAILED') or status.startswith('HESS_PARSE'):
            continue
        _refine_from_folder(
            _branch_dir(folder,name,step_A), bid, is_root=False, depth=depth+1,
            max_iters=max_iters, xtb=xtb, charge=charge, expected=expected,
            imag_thresh=imag_thresh, step_A=step_A, env=env, rows=rows,
        )


def refine(args):
    """Displace saddles ±0.15 Å, reoptimize, re-check frequencies, repeat until
    confirmed or max_iters. Never overwrites xtb.out."""
    out=Path(args.output); manifest=json.loads((out/'manifest.json').read_text())
    args.xtb=require_xtb(args.xtb)
    print('Using xTB', args.xtb)
    expected=read_xyz(out/'frame_000'/'input.xyz')[0]['symbols']
    ids=[e['id'] for e in manifest['entries']]
    if args.limit: ids=ids[:args.limit]
    imag_thresh=float(args.imag_thresh); step_A=float(getattr(args,'step_A', MODE_STEP_A))
    max_iters=int(getattr(args,'max_iters', DEFAULT_MAX_REFINE_ITERS))
    env=xtb_thread_env(); rows=[]
    for id in ids:
        _refine_from_folder(
            out/id, id, is_root=True, depth=0, max_iters=max_iters,
            xtb=args.xtb, charge=manifest['charge'], expected=expected,
            imag_thresh=imag_thresh, step_A=step_A, env=env, rows=rows,
        )
    post=out/'postprocess'; post.mkdir(exist_ok=True)
    dest=stamp_if_exists(post/'refine_status.csv')
    write_csv(dest, rows)
    print('Wrote', dest)


def run_cocos_on_xyz(post, xyz_name, cocos, seed, samples):
    wrapper=Path(__file__).with_name('seeded_cocos.py').resolve()
    if not wrapper.exists():
        raise FileNotFoundError(wrapper)
    code=Path(cocos).resolve()
    if not code.exists():
        raise FileNotFoundError(code)
    cmd=[sys.executable,str(wrapper),str(seed),str(code),'-f',xyz_name,
         '--DoNotUseGlobalOpt','--NumSamples',str(samples),'--DoNotSaveDoubles']
    proc=subprocess.run(cmd,cwd=post,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=False)
    return cmd, proc


def parse_cocos_survivors(stdout_text):
    added=[int(x) for x in re.findall(r'added structure (\d+)\b', stdout_text)]
    removed=[{'structure':int(a),'copy_of':int(b)}
             for a,b in re.findall(r'Removed structure (\d+) as a copy of (\d+)', stdout_text)]
    return list(dict.fromkeys(added)), removed


def cocos_robustness(post, xyz_name, cocos, seeds, samples_list):
    post=Path(post)
    rows=[]
    for ns in samples_list:
        for seed in seeds:
            job=post/'cocos_sensitivity'/f'N{ns}_seed{seed}'
            if job.exists() and (job/'cocos_log.txt').is_file():
                text=(job/'cocos_log.txt').read_text(errors='replace')
                added, removed=parse_cocos_survivors(text)
                rows.append({'samples':ns,'seed':seed,'num_unique':len(added),
                             'survivor_sorted_indices':';'.join(map(str,added)),
                             'status':'skipped_existing'})
                continue
            job.mkdir(parents=True, exist_ok=True)
            shutil.copy2(post/xyz_name, job/xyz_name)
            cmd, proc=run_cocos_on_xyz(job, xyz_name, cocos, seed, ns)
            (job/'cocos_log.txt').write_text(proc.stdout)
            (job/'cocos_command.json').write_text(json.dumps({'cmd':cmd,'returncode':proc.returncode},indent=2))
            if proc.returncode!=0:
                rows.append({'samples':ns,'seed':seed,'num_unique':None,'status':f'FAILED_{proc.returncode}'})
                continue
            added, removed=parse_cocos_survivors(proc.stdout)
            rows.append({'samples':ns,'seed':seed,'num_unique':len(added),
                         'survivor_sorted_indices':';'.join(map(str,added)),
                         'removed':json.dumps(removed),'status':'success',
                         'cocos_sha256':sha256(cocos)})
    dest=stamp_if_exists(post/'cocos_sensitivity.csv')
    write_csv(dest, rows)
    nuniq={r['num_unique'] for r in rows if r.get('num_unique') is not None}
    added={r.get('survivor_sorted_indices') for r in rows if r.get('status') in ('success','skipped_existing')}
    stable=len(nuniq)==1 and len(added)==1
    (post/'cocos_robustness.json').write_text(json.dumps({
        'stable_grouping':stable,'n_unique_values':sorted(nuniq, key=lambda x: (x is None, x)),
        'added_index_sets':sorted(added, key=lambda x: x or ''),
        'note':'Stability on this seed set is not a proof of algorithmic completeness.',
    }, indent=2), encoding='utf-8')
    print('CoCoS robustness wrote', dest, 'stable_grouping=', stable)
    return rows, stable


def dedup(args):
    """Rank free opts, then CoCoS **only** on confirmed minima (ν ≥ −5 cm⁻¹)."""
    out=Path(args.output); imag_thresh=float(getattr(args,'imag_thresh',5.0))
    step_A=float(getattr(args,'step_A', MODE_STEP_A))
    max_iters=int(getattr(args,'max_iters', DEFAULT_MAX_REFINE_ITERS))
    confirmed, notes, expected, manifest=collect_confirmed_minima(
        out, imag_thresh, step_A, max_iters=max_iters)
    qc_checks.assert_supported_ion(manifest['central'])
    good=[]; bad=[]
    for rec in manifest['entries']:
        folder=out/rec['id']; data=get_result(folder,expected)
        if not data:bad.append(rec['id']);continue
        frame,e=data
        good.append({'id':rec['id'],'source_run_id':rec.get('source_run_id'),'frame':frame,'energy_Eh':e,
                     'constrained_Eh':rec['constrained_energy_Eh']})
    post=out/'postprocess'; post.mkdir(exist_ok=True)
    (post/'missing_failed_ids.txt').write_text('\n'.join(bad)+'\n' if bad else '(none)\n')
    (post/'confirmed_minima_notes.json').write_text(json.dumps(notes, indent=2, default=str), encoding='utf-8')
    if not good:raise RuntimeError('No converged xtb results; cannot run CoCoS or report minimum energies')
    if bad and not args.allow_incomplete:
        raise RuntimeError(
            f'{len(bad)} of {len(manifest["entries"])} inputs are not verified as converged. '
            'CoCoS and energy-window selection blocked to prevent misleading final results. '
            'Inspect xtb_runs/postprocess/missing_failed_ids.txt, rerun failures; '
            'for diagnostics only, use --allow-incomplete.'
        )
    missing_hess=[n for n in notes if n.get('status') in ('HESS_MISSING','HESS_PARSE')]
    if missing_hess and not args.allow_incomplete:
        raise RuntimeError(
            f'{len(missing_hess)} frames have no verified Hessian. Run hess (and refine for saddles) first.'
        )
    good.sort(key=lambda r:(r['energy_Eh'],r['id']))
    ener0=good[0]['energy_Eh']; rows=[]
    for j,r in enumerate(good):
        rows.append({'sorted_index':j,'source_id':r['id'],'source_run_id':r.get('source_run_id'),'energy_free_Eh':r['energy_Eh'],
                     'delta_free_kcal_mol':(r['energy_Eh']-ener0)*HARTREE_TO_KCAL_MOL,
                     'energy_constrained_Eh':r['constrained_Eh'],
                     'geometry_change_A_unaligned':None})
    dest_all=stamp_if_exists(post/'all_freely_optimized.csv') if (post/'all_freely_optimized.csv').exists() else post/'all_freely_optimized.csv'
    write_csv(dest_all, rows)
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
    dest_rank=stamp_if_exists(post/'rank_shift_constrained_vs_unconstrained.csv') if (post/'rank_shift_constrained_vs_unconstrained.csv').exists() else post/'rank_shift_constrained_vs_unconstrained.csv'
    write_csv(dest_rank, rank_rows)
    if len(good)>1 and not (post/'rank_shift_after_relaxation.png').exists():
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
    if not confirmed:
        raise RuntimeError(
            'No confirmed minima (all ν ≥ −5 cm⁻¹). Run hess and refine; saddles are excluded from CoCoS.'
        )
    emin=confirmed[0]['energy_Eh']
    conf_rows=[]
    tosave=[]
    for j,r in enumerate(confirmed):
        delta=(r['energy_Eh']-emin)*HARTREE_TO_KCAL_MOL
        conf_rows.append({
            'sorted_index':j,'id':r['id'],'parent_id':r['parent_id'],'origin':r['origin'],
            'energy_Eh':r['energy_Eh'],'delta_kcal_mol':delta,
            'n_significant_imaginary':r['n_imag'],'lowest_freq_cm1':r['lowest_freq_cm1'],
            'xyz':r['xyz'],'hess_log':r['hess_log'],
        })
        tosave.append({'symbols':r['frame']['symbols'],'xyz':r['frame']['xyz'],
                       'comment':(f'id={r["id"]} origin={r["origin"]} E = {r["energy_Eh"]:.14f} Eh '
                                  f'delta_E = {delta:.8f} kcal/mol n_imag = {r["n_imag"]}')})
    write_csv(post/'confirmed_minima.csv', conf_rows)
    write_xyz(post/'confirmed_minima_for_cocos.xyz', tosave)
    # Compatibility name used by older robustness helper: same confirmed-only content.
    write_xyz(post/'free_sorted_for_cocos.xyz', tosave)
    n_saddle=sum(1 for n in notes if n.get('status')=='SADDLE')
    (post/'COMPLETENESS.txt').write_text(
        f'Verified unconstrained optimizations: {len(good)}/{len(manifest["entries"])}\n'
        f'Confirmed minima (ν ≥ {NU_MIN_CM1} cm⁻¹), including mode-follow branches: {len(confirmed)}\n'
        f'Saddles excluded from CoCoS until refine confirms a minimum: {n_saddle}\n'
        + ('INCOMPLETE: diagnostics only.\n' if bad else 'COMPLETE: all requested inputs checked.\n')
        + 'These energies are calculated total energies, NOT Gibbs free energies.\n'
        + 'CoCoS input is confirmed minima only.\n'
    )
    print(f'Converged = {len(good)}, not verified = {len(bad)}, confirmed minima = {len(confirmed)}')
    input_sha=sha256(post/'confirmed_minima_for_cocos.xyz')
    (post/'confirmed_minima_sha256.txt').write_text(input_sha+'\n', encoding='utf-8')
    unique_csv=post/'unique_representatives.csv'
    unique_xyz=post/'unique_lowest_energy_representatives.xyz'
    skip_cocos=False
    if unique_csv.exists() or unique_xyz.exists():
        prev=(post/'cocos_input_sha256.txt').read_text().strip() if (post/'cocos_input_sha256.txt').exists() else None
        if prev==input_sha:
            print('CoCoS outputs exist for this confirmed-minima set; skipping CoCoS (resume).')
            skip_cocos=True
        else:
            raise RuntimeError(
                f'Refusing overwrite of {unique_csv} / {unique_xyz}. '
                'Archive them first; confirmed-minima SHA does not match the previous CoCoS input.'
            )
    if not args.cocos:
        return
    if not skip_cocos:
        cmd, proc=run_cocos_on_xyz(post, 'confirmed_minima_for_cocos.xyz', args.cocos, args.seed, args.samples)
        (post/'cocos_stdout.log').write_text(proc.stdout)
        (post/'cocos_input_sha256.txt').write_text(input_sha+'\n', encoding='utf-8')
        if proc.returncode!=0:raise RuntimeError('CoCoS failed; see cocos_stdout.log')
        survivors, removed=parse_cocos_survivors(proc.stdout)
        observed=list(post.glob('*_quick_cleanup.trj'))
        if len(observed)!=1:raise RuntimeError('Ambiguous/missing CoCoS output')
        kept=read_xyz(observed[0])
        if len(kept)!=len(survivors):
            raise RuntimeError(f'Cannot map CoCoS output: {len(kept)} frames but {len(survivors)} log indices')
        retained=[confirmed[i] for i in survivors]
        write_xyz(unique_xyz, [{'symbols':r['frame']['symbols'],'xyz':r['frame']['xyz'],
                  'comment':f'id={r["id"]} origin={r["origin"]} E = {r["energy_Eh"]:.14f} Eh'} for r in retained])
        write_csv(unique_csv,[
            {'source_id':r['id'],'parent_id':r['parent_id'],'origin':r['origin'],
             'source_run_id':r.get('source_run_id'),'energy_free_Eh':r['energy_Eh'],
             'delta_free_kcal_mol':(r['energy_Eh']-emin)*HARTREE_TO_KCAL_MOL,
             'n_significant_imaginary':r['n_imag'],
             **geometry_features(r['frame'],manifest['central'],args.shell_cut)} for r in retained])
        (post/'cocos_parameters.json').write_text(json.dumps({
            'seed':args.seed,'samples':args.samples,'cocos_sha256':sha256(args.cocos),'cocos_cmd':cmd,
            'removed':removed,
            'comment':'CoCoS on confirmed minima only (ν ≥ −5 cm⁻¹), energy-sorted. Saddles excluded.',
        },indent=2), encoding='utf-8')
        print('CoCoS classes:',len(retained),'lowest-energy representatives saved; verify sensitivity before publication')
    if not getattr(args,'skip_robustness', False):
        seeds=getattr(args,'robust_seeds', None) or DEFAULT_ROBUST_SEEDS
        samples_list=getattr(args,'robust_samples', None) or (args.samples, 100)
        cocos_robustness(post, 'confirmed_minima_for_cocos.xyz', args.cocos, seeds, samples_list)
    if args.mdatom:
        mdscript=Path(args.mdatom).resolve()
        if not mdscript.exists():raise FileNotFoundError(mdscript)
        mdlog=post/'MDAtomLocalEnvironment.log'
        if mdlog.exists():
            print('MDAtomLocalEnvironment.log exists; skipping (resume).')
            return
        mdcmd=[sys.executable,str(mdscript),'-t','unique_lowest_energy_representatives.xyz',
               '-e',manifest['central'],'--RunWaterAnalysis','--RecalcWater']
        mdproc=subprocess.run(mdcmd,cwd=post,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        mdlog.write_text(mdproc.stdout)
        (post/'MDAtomLocalEnvironment_command.json').write_text(json.dumps(mdcmd,indent=2))
        if mdproc.returncode !=0:
            raise RuntimeError('MDAtomLocalEnvironment failed: see log')
        print('MDAtomLocalEnvironment completed; interpret its time column as frame index, not physical simulation time.')


def main():
    pa=argparse.ArgumentParser(description=__doc__);subs=pa.add_subparsers(dest='command',required=True)
    p=subs.add_parser('audit');p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--central',default='Rb',choices=list(SUPPORTED_IONS));p.add_argument('--shell-cut',type=float);p.add_argument('--energies-csv');p.set_defaults(func=audit)
    p=subs.add_parser('prepare');p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--central',default='Rb',choices=list(SUPPORTED_IONS));p.add_argument('--charge',type=int,default=PROTOCOL_CHARGE)
    p.add_argument('--allow-existing',action='store_true');p.add_argument('--energies-csv');p.set_defaults(func=prepare)
    p=subs.add_parser('relax');p.add_argument('--output',required=True);p.add_argument('--xtb',default=None)
    p.add_argument('--workers',type=int,default=2);p.add_argument('--limit',type=int);p.set_defaults(func=relax)
    p=subs.add_parser('hess');p.add_argument('--output',required=True);p.add_argument('--xtb',default=None)
    p.add_argument('--workers',type=int,default=1);p.add_argument('--limit',type=int)
    p.add_argument('--imag-thresh',type=float,default=5.0);p.set_defaults(func=hess)
    p=subs.add_parser('refine');p.add_argument('--output',required=True);p.add_argument('--xtb',default=None)
    p.add_argument('--limit',type=int);p.add_argument('--imag-thresh',type=float,default=5.0)
    p.add_argument('--step-A',type=float,default=MODE_STEP_A)
    p.add_argument('--max-iters',type=int,default=DEFAULT_MAX_REFINE_ITERS);p.set_defaults(func=refine)
    p=subs.add_parser('dedup');p.add_argument('--output',required=True);p.add_argument('--cocos')
    p.add_argument('--seed',type=int,default=DEFAULT_COCOS_SEED);p.add_argument('--samples',type=int,default=DEFAULT_COCOS_SAMPLES)
    p.add_argument('--shell-cut',type=float);p.add_argument('--mdatom')
    p.add_argument('--imag-thresh',type=float,default=5.0);p.add_argument('--step-A',type=float,default=MODE_STEP_A)
    p.add_argument('--max-iters',type=int,default=DEFAULT_MAX_REFINE_ITERS)
    p.add_argument('--skip-robustness',action='store_true')
    p.add_argument('--allow-incomplete', action='store_true', help='Diagnostics only: do not treat partial sample as final data')
    p.set_defaults(func=dedup)
    args=pa.parse_args();args.func(args)

if __name__=='__main__':main()
