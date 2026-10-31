Live pins from Ubuntu WSL on 2026-10-30. Machine-readable copy: `PINNED_VERSIONS.json`.

| Component | Revision |
| --- | --- |
| GWH | `7411811374b6422a15f23d6d86a47f2145f8ac47` — author commit “fixed Euler angles definition” (`from_euler('zyz')`; SciPy docs: lowercase=extrinsic, uppercase `ZYZ`=intrinsic; source not rewritten) |
| CoCoS | `47f6680da2eee90d9b18536fd94d2bc8f6d05aad` — file SHA256 `4c538fb8…` identical to uploaded snapshot |
| MDAtomLocalEnvironment | `3c43de70abef1f64fb31a2fe2693e4e73b88083c` |
| xTB source tag | v6.7.1 = `26b28010e805f7d1aeeef39813feb473e69cc4be` |
| xTB **runtime** | Intel-static `~/xtb-standalone/xtb-dist/bin/xtb`, reports **6.7.1 (edcfbbe)**, SHA256 `debf27a9e0fa4bfb5ca75aafe4b90d8211f08ec2f4a482f375a4987212eaa12a` |
| `/usr/bin/xtb` | 6.7.1 (unknown), libgfortran 16 — **ANCOPT crash, do not use for `--opt`** |
| Python | 3.14.4; numpy 2.3.5; scipy 1.16.3; matplotlib 3.10.7+dfsg1; pytest 9.2.0 |
| K_06 packing.vecfle | SHA256 `06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5` |

Existing `cases/Li_07` and `pilots/Li_07` were run with GWH `bfc474c` (old `zyx`). Do not mix or overwrite.

| Burette | https://github.com/SergeiNikolenko/Burette HEAD `ba8171ec5bc0315bb3ace11cf2374ee908b67728` — macOS 12+ Apple Silicon app / Quick Look / Codex plugin. **Not runnable on this Windows/WSL host.** Browser demo is not a headless exporter. Scene XYZ stored in `figures/burette/scenes/`. |
| xyzrender | 0.3.8 (`/home/antiz/.local/bin/xyzrender`). Burette’s documented XYZ fallback. Exports labeled `xyzrender`, not Burette.app. Canvas request 2400 px; observed PNG long edge 6207–7200 px. |
