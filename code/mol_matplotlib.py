#!/usr/bin/env python3
"""Ball-and-stick images with Matplotlib. NEVER label these as Burette renders."""
from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import pipeline as P
    import visualize as V
    import shell_descriptors as SD
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import pipeline as P
    import visualize as V
    import shell_descriptors as SD

COLORS = {
    "H": "#f0f0f0",
    "O": "#e74c3c",
    "K": "#8e44ad",
    "Li": "#9b59b6",
    "Na": "#3498db",
    "Rb": "#6c3483",
    "F": "#27ae60",
    "Cl": "#1e8449",
    "Br": "#922b21",
}
RADII = {"H": 0.22, "O": 0.38, "K": 0.62, "Li": 0.45, "Na": 0.52, "Rb": 0.70, "F": 0.45, "Cl": 0.58, "Br": 0.62}
COV = {"H": 0.31, "O": 0.66, "K": 2.03, "Li": 1.28, "Na": 1.66, "Rb": 2.20, "F": 0.64, "Cl": 1.02, "Br": 1.20}


def covalent_bonds(symbols, xyz, scale=1.20):
    """Only intramolecular O–H (checked 0.7–1.25 Å). Ion–O is coordination, not covalent."""
    bonds = []
    n = len(symbols)
    for i in range(n):
        for j in range(i + 1, n):
            pair = {symbols[i], symbols[j]}
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if pair == {"O", "H"} and 0.70 <= d <= 1.25:
                bonds.append((i, j, d))
    return bonds


def coordination_contacts(symbols, xyz, ion, cut):
    contacts = []
    try:
        ci = symbols.index(ion)
    except ValueError:
        return contacts
    for j, s in enumerate(symbols):
        if s != "O":
            continue
        d = float(np.linalg.norm(xyz[ci] - xyz[j]))
        if d < float(cut):
            contacts.append((ci, j, d))
    return contacts


def draw_frame(ax, frame, title, hbonds=None, xh=None, ion="K"):
    symbols = list(frame["symbols"])
    xyz = np.asarray(frame["xyz"], float)
    # depth sort
    order = np.argsort(xyz[:, 1])
    bonds = covalent_bonds(symbols, xyz)
    for i, j, d in bonds:
        a, b = xyz[i], xyz[j]
        ax.plot([a[0], b[0]], [a[2], b[2]], color="#333333", lw=1.8, zorder=1, solid_capstyle="round")
    cut = 3.20
    for i, j, d in coordination_contacts(symbols, xyz, ion, cut):
        a, b = xyz[i], xyz[j]
        ax.plot([a[0], b[0]], [a[2], b[2]], color="#7d3c98", lw=1.0, ls=(0, (3, 2)), zorder=1.5, alpha=0.85)
    if hbonds:
        for hb in hbonds:
            a = xyz[int(hb["h_index"])]
            b = xyz[int(hb["acceptor_O"])]
            ax.plot([a[0], b[0]], [a[2], b[2]], color="#3498db", lw=1.1, ls=(0, (2, 2)), zorder=2)
    if xh:
        ci = symbols.index(ion)
        for contact in xh:
            a = xyz[int(contact["h_index"])]
            b = xyz[ci]
            ax.plot([a[0], b[0]], [a[2], b[2]], color="#16a085", lw=1.1, ls=(0, (1, 2)), zorder=2)
    for i in order:
        s = symbols[i]
        ax.scatter(
            xyz[i, 0], xyz[i, 2],
            s=(RADII.get(s, 0.35) * 420) ** 1,
            c=COLORS.get(s, "#888888"),
            edgecolors="black",
            linewidths=0.45,
            zorder=3 + float(xyz[i, 1]),
        )
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(r"$x$ / Å")
    ax.set_ylabel(r"$z$ / Å")
    ax.tick_params(labelsize=8)


def panel_figure(entries, dest: Path, caption, filename):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import json
    n = len(entries)
    fig, axes = plt.subplots(1, n, figsize=(8.0 * max(n, 1) / 2.2, 5.4))
    if n == 1:
        axes = [axes]
    for ax, e in zip(axes, entries):
        raw = P.read_xyz(e["xyz_path"])[0]
        aligned = V.align_ion_and_oxygens(raw, e["ion"])
        desc = e.get("descriptors") or SD.cluster_descriptors(raw, e["ion"], n_water=e.get("n_water"), xyz_path=e["xyz_path"])
        hb = json.loads(desc.get("hbonds_json") or "[]")
        xh = json.loads(desc.get("anion_xh_json") or "[]")
        title = f"{e['id']}\nE = {e['energy_Eh']:.8f} Eh\nCoCoS {e.get('cocos_group')}"
        draw_frame(ax, aligned, title, hbonds=hb, xh=xh, ion=e["ion"])
    fig.suptitle(caption + "\nrenderer: matplotlib (not Burette.app)", fontsize=11)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    png = dest.with_suffix(".png")
    svg = dest.with_suffix(".svg")
    pdf = dest.with_suffix(".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [png, svg, pdf]


def k06_matplotlib_structures(entries):
    pub = Path(__file__).resolve().parents[1] / "figures" / "publication"
    files = []
    files += panel_figure(
        entries,
        pub / "Fig_K06_Pilot",
        r"[K(H$_2$O)$_6$]$^+$ pilot: 3 frequency-confirmed geometries, one CoCoS group",
        "Fig_K06_Pilot",
    )
    files += panel_figure(
        entries,
        pub / "Fig_Conformers",
        r"Lowest-energy structurally compared minima of [K(H$_2$O)$_6$]$^+$ (same CoCoS class)",
        "Fig_Conformers",
    )
    files += panel_figure(
        [min(entries, key=lambda e: e["energy_Eh"])],
        pub / "Fig_Coordination_Shells",
        r"First shell of the lowest-energy [K(H$_2$O)$_6$]$^+$ minimum (frame_001)",
        "Fig_Coordination_Shells",
    )
    files += panel_figure(
        [min(entries, key=lambda e: e["energy_Eh"])],
        pub / "Fig_Structures_Cations",
        r"Only cation hydrate calculated here: [K(H$_2$O)$_6$]$^+$. Li/Na/Rb series not run.",
        "Fig_Structures_Cations",
    )
    note = pub / "Fig_Structures_Anions_NOT_GENERATED.txt"
    note.write_text(
        "Fig_Structures_Anions was not generated: no F-/Cl-/Br- calculations exist in this kit.\n"
        "CASE_STATUS has no anion rows. No dummy structures were drawn.\n",
        encoding="utf-8",
    )
    return files
