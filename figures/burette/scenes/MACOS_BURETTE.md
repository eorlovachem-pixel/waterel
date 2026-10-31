# Open these scenes in Burette.app (macOS only)

This Windows/WSL host cannot run Burette.app (macOS 12+, Apple Silicon; Homebrew cask `SergeiNikolenko/burette`).

On a supported Mac:

1. `brew tap SergeiNikolenko/burette && brew install --cask burette`
2. Open Burette.app once.
3. `Cmd+O` each file `figures/burette/scenes/<id>/original.xyz` (coordinates identical to the frequency-confirmed xTB geometries).
4. Renderer: Mol* or **xyzrender** (Burette settings → Structure rendering). Camera JSON is `scene.json` in the same folder.
5. Do not export a screenshot that includes the Burette UI chrome. Use Burette/xyzrender export if the installed version supports it.

Pinned Burette HEAD at generation: `ba8171ec5bc0315bb3ace11cf2374ee908b67728`.
Headless figures in this kit were produced with **xyzrender 0.3.8** (`--config tube -S 2400 --hy --no-orient -B white --nci --hl M`) and are labeled xyzrender, not Burette.app.
