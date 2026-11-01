# Method limitations

- ν ≥ −5 cm⁻¹ is a documented engineering threshold, not a validated spectroscopic cutoff.
- CN cuts (K 3.20 Å, Cl 3.40 Å, Br 3.50 Å) are `BLOCKING_UNCONFIRMED`.
- Electronic E is not a free energy. E(n) values of different n are not hydration energies.
- GWH `'zyz'` vs SciPy-intrinsic `'ZYZ'` is an open author question.
- Cl⁻/Br⁻ packing R = Shannon/Pauling VI sums (3.21 / 3.36 Å). Aqueous Cl–O first shells are often ~3.14–3.19 Å; free opt would relax R if it were run.
- Author K packings for n=6–14 were **not** replaced by Sloane codes. New Cl/Br directions are Sloane; RMSD vs K is in `inputs/ANION_PACKING_PROVENANCE.json`.
- Production 100-start chemistry is blocked until Denis assigns cases.
- Burette.app is not available on this Windows/WSL host; molecular images use xyzrender/Matplotlib and are labeled as such.
