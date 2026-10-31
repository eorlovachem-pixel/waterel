#!/usr/bin/env python3
"""First-shell, H-bond, and connectivity descriptors for Ion(H2O)n.

Does not invent quantum-chemical energies. Geometry numbers come from XYZ
coordinates already on disk. CN cutoffs are provisional and marked
BLOCKING_UNCONFIRMED until Denis agrees the selection rule.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

try:
    import qc_checks
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import qc_checks

# Engineering defaults only. Not an agreed publication filter.
# K6 confirmed maxima are ~2.706 Å, so CN(K)=6 at 3.20 Å. Other ions have no
# calculated first-shell distribution in this kit.
CN_CUT_PROVISIONAL_A = {
    "Li": 2.80,
    "Na": 3.10,
    "K": 3.20,
    "Rb": 3.40,
    "F": 2.80,
    "Cl": 3.40,
    "Br": 3.50,
}
CN_CUT_STATUS = "BLOCKING_UNCONFIRMED"
OH_ASSIGN_MAX_A = 1.25
HB_H_ACC_MAX_A = 2.50
HB_OO_MAX_A = 3.50
HB_ANGLE_MIN_DEG = 130.0
ANION_XH_MAX_A = 2.50
ANION_OHX_ANGLE_MIN_DEG = 130.0
PROTON_TRANSFER_MARGIN_A = 0.15


def sha256_xyz_bytes(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_frame(frame) -> str:
    lines = [f'{len(frame["symbols"])}', frame.get("comment", "")]
    for symbol, r in zip(frame["symbols"], np.asarray(frame["xyz"], float)):
        lines.append(f"{symbol:2s} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    blob = ("\n".join(lines) + "\n").encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def unit(v):
    n = np.linalg.norm(v)
    if n < 1e-12:
        return v * 0.0
    return v / n


def angle_deg(a, b, c):
    """Angle ABC in degrees."""
    v1 = np.asarray(a, float) - np.asarray(b, float)
    v2 = np.asarray(c, float) - np.asarray(b, float)
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-12 or n2 < 1e-12:
        return float("nan")
    cos = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(np.degrees(np.arccos(cos)))


def assign_water_hydrogens(o_xyz, h_xyz):
    """Globally assign two H atoms to each O. Returns pairs or raises."""
    if len(h_xyz) != 2 * len(o_xyz):
        return None, "not_2H_per_O"
    if len(o_xyz) == 0:
        return [], "no_oxygen"
    slots = np.repeat(o_xyz, 2, axis=0)
    dist = np.linalg.norm(h_xyz[:, None, :] - slots[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(dist)
    if float(np.max(dist[rows, cols])) > OH_ASSIGN_MAX_A:
        return None, "O_H_distance_above_1.25_A"
    pairs = [[] for _ in range(len(o_xyz))]
    pair_idx = [[] for _ in range(len(o_xyz))]
    for h_index, slot in zip(rows, cols):
        pairs[slot // 2].append(h_xyz[h_index])
        pair_idx[slot // 2].append(int(h_index))
    return {"xyz": pairs, "h_indices": pair_idx}, "nominal_2H_per_O"


def cluster_descriptors(frame, ion, n_water=None, shell_cut=None, xyz_path=None):
    """Return a JSON-serialisable dict. Lists of per-water values kept as CSV strings too."""
    qc_checks.assert_supported_ion(ion)
    symbols = np.array(frame["symbols"])
    points = np.asarray(frame["xyz"], float)
    center_ids = np.where(symbols == ion)[0]
    if len(center_ids) != 1:
        raise ValueError(f"Expected one central {ion}, got {len(center_ids)}")
    ci = int(center_ids[0])
    center = points[ci]
    o_ids = np.where(symbols == "O")[0]
    h_ids = np.where(symbols == "H")[0]
    o_xyz = points[o_ids]
    h_xyz = points[h_ids]
    if n_water is None:
        n_water = int(len(o_ids))
    n_water = int(n_water)
    mo = np.linalg.norm(o_xyz - center, axis=1) if len(o_ids) else np.array([])
    cut = float(shell_cut) if shell_cut is not None else CN_CUT_PROVISIONAL_A.get(ion)
    first_shell = np.where(mo < cut)[0] if cut is not None and len(mo) else np.arange(len(o_ids))
    r_shell = mo[first_shell] if len(first_shell) else np.array([])
    cn = int(len(first_shell))
    r_mean = float(r_shell.mean()) if len(r_shell) else None
    r_min = float(r_shell.min()) if len(r_shell) else None
    r_max = float(r_shell.max()) if len(r_shell) else None
    r_std = float(r_shell.std(ddof=0)) if len(r_shell) else None
    radial_def = None
    radial_span = None
    if r_mean and r_mean > 0 and r_min is not None:
        radial_def = float(r_std / r_mean)
        radial_span = float((r_max - r_min) / r_mean)

    assigned, connectivity = assign_water_hydrogens(o_xyz, h_xyz)
    cos_vals = []
    p2_vals = []
    water_rows = []
    hbonds = []
    anion_xh = []
    proton_transfer = False
    oho_angles = []
    oo_distances = []

    if assigned is not None:
        for local_i, (oxygen, hpair, hloc) in enumerate(
            zip(o_xyz, assigned["xyz"], assigned["h_indices"])
        ):
            mu = np.mean(hpair, axis=0) - oxygen
            radial = oxygen - center
            c = float(np.dot(mu, radial) / (np.linalg.norm(mu) * np.linalg.norm(radial)))
            c = float(np.clip(c, -1.0, 1.0))
            p2 = 0.5 * (3.0 * c * c - 1.0)
            cos_vals.append(c)
            p2_vals.append(p2)
            in_shell = bool(local_i in set(int(x) for x in first_shell))
            water_rows.append({
                "o_index": int(o_ids[local_i]),
                "h_indices": [int(h_ids[j]) for j in hloc],
                "R_ion_O_A": float(mo[local_i]),
                "in_first_shell": in_shell,
                "cos_theta": c,
                "P2": float(p2),
            })
            if qc_checks.is_anion(ion):
                for hp, h_local in zip(hpair, hloc):
                    d_hx = float(np.linalg.norm(hp - center))
                    d_ho = float(np.linalg.norm(hp - oxygen))
                    ang = angle_deg(oxygen, hp, center)
                    if d_hx < ANION_XH_MAX_A and ang >= ANION_OHX_ANGLE_MIN_DEG:
                        anion_xh.append({
                            "h_index": int(h_ids[h_local]),
                            "o_index": int(o_ids[local_i]),
                            "H_X_A": d_hx,
                            "O_H_A": d_ho,
                            "angle_OHX_deg": ang,
                        })
                    if d_hx + PROTON_TRANSFER_MARGIN_A < d_ho:
                        proton_transfer = True

        for i in range(len(o_xyz)):
            for j in range(i + 1, len(o_xyz)):
                oo = float(np.linalg.norm(o_xyz[i] - o_xyz[j]))
                oo_distances.append(oo)
                oio = angle_deg(o_xyz[i], center, o_xyz[j])
                if i in first_shell and j in first_shell:
                    oho_angles.append(oio)
                # O-H···O hydrogen bonds
                for hp, h_local in zip(assigned["xyz"][i], assigned["h_indices"][i]):
                    hacc = float(np.linalg.norm(hp - o_xyz[j]))
                    ang = angle_deg(o_xyz[i], hp, o_xyz[j])
                    if hacc < HB_H_ACC_MAX_A and oo < HB_OO_MAX_A and ang >= HB_ANGLE_MIN_DEG:
                        hbonds.append({
                            "donor_O": int(o_ids[i]),
                            "acceptor_O": int(o_ids[j]),
                            "h_index": int(h_ids[h_local]),
                            "H_O_A": hacc,
                            "O_O_A": oo,
                            "angle_OHO_deg": ang,
                        })
                for hp, h_local in zip(assigned["xyz"][j], assigned["h_indices"][j]):
                    hacc = float(np.linalg.norm(hp - o_xyz[i]))
                    ang = angle_deg(o_xyz[j], hp, o_xyz[i])
                    if hacc < HB_H_ACC_MAX_A and oo < HB_OO_MAX_A and ang >= HB_ANGLE_MIN_DEG:
                        hbonds.append({
                            "donor_O": int(o_ids[j]),
                            "acceptor_O": int(o_ids[i]),
                            "h_index": int(h_ids[h_local]),
                            "H_O_A": hacc,
                            "O_O_A": oo,
                            "angle_OHO_deg": ang,
                        })

    full_shell = bool(cn == n_water)
    data = {
        "ion": ion,
        "charge": qc_checks.protocol_charge(ion),
        "n_water_nominal": n_water,
        "n_oxygen": int(len(o_ids)),
        "n_hydrogen": int(len(h_ids)),
        "central_index": ci,
        "xyz_sha256": sha256_xyz_bytes(xyz_path) if xyz_path else sha256_frame(frame),
        "xyz_path": str(xyz_path) if xyz_path else None,
        "CN_cut_A": cut,
        "CN_cut_status": CN_CUT_STATUS,
        "CN": cn,
        "first_shell_O_indices": ";".join(str(int(o_ids[k])) for k in first_shell),
        "ion_O_distances_A": ";".join(f"{v:.8f}" for v in mo),
        "R_mean_A": r_mean,
        "R_min_A": r_min,
        "R_max_A": r_max,
        "R_std_A": r_std,
        "radial_deformation_std_over_mean": radial_def,
        "radial_deformation_span_over_mean": radial_span,
        "full_first_shell": full_shell,
        "water_connectivity": connectivity,
        "orientation_cos_mean": float(np.mean(cos_vals)) if cos_vals else None,
        "orientation_cos_min": float(np.min(cos_vals)) if cos_vals else None,
        "orientation_cos_max": float(np.max(cos_vals)) if cos_vals else None,
        "orientation_cos_values": ";".join(f"{v:.8f}" for v in cos_vals),
        "P2_mean": float(np.mean(p2_vals)) if p2_vals else None,
        "P2_values": ";".join(f"{v:.8f}" for v in p2_vals),
        "O_ion_O_angle_count": len(oho_angles),
        "O_ion_O_angle_mean_deg": float(np.mean(oho_angles)) if oho_angles else None,
        "O_ion_O_angle_min_deg": float(np.min(oho_angles)) if oho_angles else None,
        "O_ion_O_angle_max_deg": float(np.max(oho_angles)) if oho_angles else None,
        "O_ion_O_angles_deg": ";".join(f"{v:.4f}" for v in oho_angles),
        "OO_distance_count": len(oo_distances),
        "OO_mean_A": float(np.mean(oo_distances)) if oo_distances else None,
        "OO_min_A": float(np.min(oo_distances)) if oo_distances else None,
        "n_OH_O_hbonds": len(hbonds),
        "n_XH_O_contacts": len(anion_xh),
        "proton_transfer_candidate": proton_transfer,
        "species_class": "anion" if qc_checks.is_anion(ion) else "cation",
        "CN_equals_nominal_n": full_shell,
        "descriptor_notes": (
            "cosθ = cos(angle between ion→O and O→H-midpoint). "
            "P2 = (3 cos²θ − 1)/2. "
            f"CN cutoff {cut} Å is {CN_CUT_STATUS}. "
            "Do not treat an unweighted mean over structures as a thermodynamic average."
        ),
        "hbonds_json": json.dumps(hbonds, ensure_ascii=False),
        "anion_xh_json": json.dumps(anion_xh, ensure_ascii=False),
        "waters_json": json.dumps(water_rows, ensure_ascii=False),
        # Keys expected by pipeline.geometry_features
        "distance_min_A": float(mo.min()) if len(mo) else None,
        "distance_mean_A": float(mo.mean()) if len(mo) else None,
        "distance_max_A": float(mo.std()) and float(mo.max()) if len(mo) else None,
        "distance_std_A": float(mo.std()) if len(mo) else None,
        "n_water": int(len(o_ids)),
        "CN_at_cut": cn if cut is not None else None,
    }
    if len(mo):
        data["distance_max_A"] = float(mo.max())
        data["distance_min_A"] = float(mo.min())
        data["distance_mean_A"] = float(mo.mean())
        data["distance_std_A"] = float(mo.std())
    return data


def csv_row(desc, extra=None):
    skip = {"hbonds_json", "anion_xh_json", "waters_json"}
    row = {k: v for k, v in desc.items() if k not in skip}
    if extra:
        row.update(extra)
    return row
