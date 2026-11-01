# Computational Methods (draft)

Electronic energies are GFN2-xTB total energies in hartree, **not** Gibbs free energies. Relative energies inside one stoichiometry are \(\Delta E = (E-E_{\min})\times 627.5094740631\) kcal/mol.

## Programs (pinned)

- GWH `7411811374b6422a15f23d6d86a47f2145f8ac47`, `Rotation.from_euler('zyz')` (SciPy-extrinsic). Not rewritten to `'ZYZ'`.
- xTB v6.7.1 Intel-static (`edcfbbe`); `/usr/bin/xtb` is not used for `--opt`.
- CoCoS `47f6680da2eee90d9b18536fd94d2bc8f6d05aad`, seed 20261030, confirmed minima only (\(\nu \ge -5\) cm\(^{-1}\)).
- MDAtomLocalEnvironment `3c43de70abef1f64fb31a2fe2693e4e73b88083c` when a case is assigned.
- Packing radii for Cl⁻/Br⁻: mendeleev 1.3.0, X⁻ VI + O²⁻ VI, pm/100. Directions: Neil Sloane dim-3 spherical codes.

## Protocol

GWH constrained Euler Powell (rigid NIST water, fixed ion–O) → unconstrained GFN2-xTB/ALPB(water) of **all** valid GWH geometries (no CoCoS pre-filter) → Hessian → mode following \(\pm 0.15\) Å → tight reopt → final Hessian → CoCoS → MDAtom + Python descriptors.

Charge: K⁺ +1; Cl⁻/Br⁻ −1. Multiplicity 1.

## What has actually been computed

Only the K⁺(H₂O)₆ **pilot** (2 GWH starts). Production 100-start series for K₄–K₁₄, Cl₄–Cl₁₄, Br₄–Br₁₄ remain `ON_HOLD` / `NOT_ASSIGNED`.
