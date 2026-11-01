#!/usr/bin/env python3
"""Numerical publication plots from existing confirmed-minima tables.

These are Matplotlib figures, not Burette renders. Missing ion/n series are
omitted; no interpolated or invented points. Unweighted means over structures
are labeled as sample statistics, not thermodynamic averages.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PUB = ROOT / "figures" / "publication"
SRC = ROOT / "code" / "figures_matplotlib.py"


def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    return plt


def load_rows(csv_path):
    with open(csv_path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fnum(row, key):
    v = row.get(key)
    if v in (None, "", "None"):
        return None
    return float(v)


def plot_k06_energy(rows, dest: Path, plt):
    labels = [r["id"] for r in rows]
    de = [fnum(r, "delta_kcal_mol") for r in rows]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    x = np.arange(len(labels))
    ax.bar(x, de, color="#4c72b0", edgecolor="black", width=0.55)
    ax.set_xticks(x, labels, rotation=15, ha="right")
    ax.set_ylabel(r"$\Delta E$ / kcal mol$^{-1}$")
    ax.set_xlabel("Confirmed minimum (same CoCoS group 0)")
    ax.set_title(r"[K(H$_2$O)$_6$]$^+$ pilot: relative GFN2-xTB/ALPB(water) energies")
    ax.axhline(0.0, color="0.3", lw=0.8)
    note = (
        "Sample of 3 Hessian-confirmed geometries from a 2-start pilot. "
        "Not a thermodynamic average. All three belong to one CoCoS class."
    )
    ax.text(0.0, -0.22, note, transform=ax.transAxes, fontsize=8, va="top")
    fig.tight_layout()
    fig.savefig(dest.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(dest.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(dest.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_k06_distances(rows, dest: Path, plt):
    labels = [r["id"] for r in rows]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    x = np.arange(len(labels))
    mins = [fnum(r, "distance_min_A") for r in rows]
    means = [fnum(r, "distance_mean_A") for r in rows]
    maxs = [fnum(r, "distance_max_A") for r in rows]
    ax.errorbar(
        x, means,
        yerr=[np.array(means) - np.array(mins), np.array(maxs) - np.array(means)],
        fmt="o", color="#c44e52", capsize=4, label=r"$R_\mathrm{mean}$ (min–max)",
    )
    ax.set_xticks(x, labels, rotation=15, ha="right")
    ax.set_ylabel(r"K–O distance / Å")
    ax.set_title(r"[K(H$_2$O)$_6$]$^+$: ion–oxygen distances of confirmed minima")
    ax.legend(frameon=False)
    ax.text(
        0.0, -0.22,
        "Individual structures, not a Boltzmann average. CN cut 3.20 Å is BLOCKING_UNCONFIRMED.",
        transform=ax.transAxes, fontsize=8, va="top",
    )
    fig.tight_layout()
    fig.savefig(dest.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(dest.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(dest.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_k06_orientation(rows, dest: Path, plt):
    labels = [r.get("id") or r.get("label") for r in rows]
    cos = [fnum(r, "orientation_cos_mean") for r in rows]
    p2 = [fnum(r, "P2_mean") for r in rows]
    if any(v is None for v in cos):
        return False
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    x = np.arange(len(labels))
    w = 0.35
    ax.bar(x - w / 2, cos, width=w, color="#55a868", edgecolor="black", label=r"$\langle\cos\theta\rangle$")
    if not any(v is None for v in p2):
        ax.bar(x + w / 2, p2, width=w, color="#8172b3", edgecolor="black", label=r"$\langle P_2\rangle$")
    ax.set_xticks(x, labels, rotation=15, ha="right")
    ax.set_ylabel("orientation statistic (sample mean over waters)")
    ax.set_title(r"[K(H$_2$O)$_6$]$^+$: water orientation of confirmed minima")
    ax.legend(frameon=False)
    ax.set_ylim(-1.05, 1.05)
    ax.axhline(0.0, color="0.4", lw=0.6)
    ax.text(
        0.0, -0.22,
        r"$\cos\theta$: angle(ion→O, O→H-midpoint). $P_2=(3\cos^2\theta-1)/2$. Unweighted over waters of one geometry.",
        transform=ax.transAxes, fontsize=8, va="top",
    )
    fig.tight_layout()
    fig.savefig(dest.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(dest.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(dest.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return True


def plot_missing_series(dest: Path, plt):
    """Document that R(n), CN(n), <cosθ>(n) cannot be drawn from one n."""
    fig, ax = plt.subplots(figsize=(8.0, 3.6))
    ax.set_axis_off()
    ax.text(
        0.02, 0.95,
        "Series plots not generated\n\n"
        "Allowed calculated data in this kit: K+(H2O)6 pilot only "
        "(2 GWH starts, 3 Hessian-confirmed minima, 1 CoCoS group).\n"
        "R(n), <cosθ>(n) and CN(n) versus hydration number require more than one n.\n"
        "Li+, Na+, Rb+, F− production series are ON_HOLD / NOT_ASSIGNED / unclaimed.\n"
        "Cl−/Br− packings exist (n=4–14) but production GWH/xTB is NOT_ASSIGNED.\n"
        "No placeholder points were inserted.",
        va="top", ha="left", fontsize=11, family="DejaVu Sans",
        transform=ax.transAxes, wrap=True,
    )
    fig.tight_layout()
    fig.savefig(dest.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(dest.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(dest.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_source_copy():
    PUB.mkdir(parents=True, exist_ok=True)
    text = Path(__file__).read_text(encoding="utf-8")
    (PUB / "figures_matplotlib.py").write_text(text, encoding="utf-8")


def make_k06_publication_plots(table_csv: Path):
    plt = _setup()
    PUB.mkdir(parents=True, exist_ok=True)
    write_source_copy()
    rows = load_rows(table_csv)
    if not rows:
        raise RuntimeError(f"no rows in {table_csv}")
    plot_k06_energy(rows, PUB / "Fig_K06_DeltaE", plt)
    plot_k06_distances(rows, PUB / "Fig_K06_R_ion_O", plt)
    # orientation may need a richer CSV
    plot_k06_orientation(rows, PUB / "Fig_K06_orientation", plt)
    plot_missing_series(PUB / "Fig_Series_R_CN_cos_NOT_GENERATED", plt)
    manifest = {
        "software": "matplotlib",
        "not_burette": True,
        "dpi": 300,
        "source_csv": str(table_csv),
        "source_code": str(PUB / "figures_matplotlib.py"),
        "files": sorted(p.name for p in PUB.glob("Fig_*")),
        "definitions": {
            "delta_kcal_mol": "E - min(E) of this confirmed set, Hartree→kcal/mol 627.5094740631",
            "R": "Euclidean ion–O distances, Å",
            "cos_theta": "cos(angle between ion→O and O→H-midpoint)",
            "P2": "(3 cos^2 theta - 1)/2",
            "CN": "count of O with R < 3.20 Å (BLOCKING_UNCONFIRMED cutoff)",
        },
    }
    (PUB / "MATPLOTLIB_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return manifest
