# Figure index

Generated: 2026-10-31T23:20:57.812080+00:00

## Renderer policy

- **Burette.app** (macOS Mol* workspace) is **not available** on this Windows/WSL host.
- Reason: Burette.app is a macOS 12+ Apple Silicon desktop application (Homebrew cask SergeiNikolenko/burette). This host is Windows/WSL. Finder Quick Look, the Codex plugin, and the iPhone target do not run in Cursor on Windows. The public browser demo is a preview surface, not a headless exporter; screenshots would include UI chrome, which is forbidden for publication figures.
- **xyzrender** is Burette's documented XYZ fallback renderer. Files under `figures/burette/` produced by the xyzrender CLI are labeled `renderer=xyzrender`, not Burette.app.
- **Matplotlib** ball-and-stick and numerical plots live under `figures/publication/` and are labeled as Matplotlib.
- No dummy Li/Na/Rb/F/Cl/Br structures were drawn.

## Atom colour legend (ball-and-stick)

| Element | Colour |
|---|---|
| K (central) | purple `#8e44ad` / xyzrender highlight `#5b2c6f` |
| O | red `#e74c3c` |
| H | white |
| covalent bonds | dark grey |
| O–H···O (if drawn) | dashed blue |
| anion X···H (if drawn) | dashed teal |

xyzrender exports found: 7 PNG, 7 SVG.
xyzrender on PATH at generation: `None`

## Geometry-backed figures

| Figure ID | Geometry ID | Source XYZ | Method | E / Eh | CoCoS | CN | Renderer | File |
|---|---|---|---|---|---|---|---|---|
| Fig_K06_Pilot | frame_001 + frame_000_minus + frame_000_plus | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib (structure); xyzrender if present under figures/burette/ | figures/publication/Fig_K06_Pilot.png |
| Fig_Conformers | frame_001 + frame_000_minus + frame_000_plus | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib | figures/publication/Fig_Conformers.png |
| Fig_Coordination_Shells | frame_001 | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib | figures/publication/Fig_Coordination_Shells.png |
| Fig_Structures_Cations | frame_001 | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib | figures/publication/Fig_Structures_Cations.png |
| Fig_K06_DeltaE | frame_001 | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib | figures/publication/Fig_K06_DeltaE.png |
| Fig_K06_R_ion_O | frame_001 | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib | figures/publication/Fig_K06_R_ion_O.png |
| Fig_K06_orientation | frame_001 | results/K_06/pilot_2/xtb_runs/.../xtbopt.xyz | GFN2-xTB/ALPB(water), Hessian-confirmed | -30.51713548722800 | 0 | 6 | matplotlib | figures/publication/Fig_K06_orientation.png |

## Not generated (no calculated structures)

- `Fig_Structures_Anions` — no F−/Cl−/Br− XYZ in this kit.
- `R(n)`, `CN(n)`, `<cosθ>(n)` series — only n = 6 exists.

## Per-geometry scenes

- `frame_001` SHA-256 `cbd2d5f4b8cfcc03bccf4e9f938707c7c1fe9d5725135f7c7f9de990e4f43450` E = -30.51713548722800 Eh CoCoS 0 CN 6 scene `figures/burette/scenes/frame_001/original.xyz`
- `frame_000_minus` SHA-256 `7497af13f32b6b578d9ee57b0b66d8c5c5e2caf13818615131e6971ca51a1688` E = -30.51710085691600 Eh CoCoS 0 CN 6 scene `figures/burette/scenes/frame_000_minus/original.xyz`
- `frame_000_plus` SHA-256 `e244b557a8b486765c47e91a7473df1eb45663e348a2fa2df4261927f1602ec4` E = -30.51706743837400 Eh CoCoS 0 CN 6 scene `figures/burette/scenes/frame_000_plus/original.xyz`
