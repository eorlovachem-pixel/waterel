#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, sys

linux = Path("/home/antiz/hydration-research/Hydration_unassigned_2026-10-30")
win = list(Path("/mnt/c/Users/antiz/OneDrive").glob(
    "*/Cursor_K6_research_handoff/Cursor_research/Hydration_17_workflows_original/Hydration_unassigned_2026-10-30"
))[0]
for rel in [
    "PINNED_VERSIONS.json",
    "PINNED_VERSIONS.md",
    "figures/burette/scenes/MACOS_BURETTE.md",
    "figures/burette/scenes/XYZRENDER_VERSION.json",
    ".gitignore",
]:
    src = win / rel
    if src.is_file():
        dst = linux / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(src.read_bytes().replace(b"\r\n", b"\n"))
        print("synced", rel)

# copy linux-only version files back
for rel in [
    "SHA256SUMS.txt",
    "figures/burette/scenes/EXPORT_SIZES.json",
    "figures/burette/scenes/XYZRENDER_VERSION.json",
]:
    src = linux / rel
    if src.is_file():
        dst = win / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(src.read_bytes())
        print("copied back", rel)

os.chdir(linux)
env = os.environ.copy()
# Per-command identity only; does not write ~/.gitconfig
cmd_prefix = ["git", "-c", "user.name=Lily", "-c", "user.email=lily@users.noreply.github.com"]
subprocess.run(["git", "add", "-A"], cwd=linux)
st = subprocess.run(["git", "status", "--short"], cwd=linux, text=True, capture_output=True)
print(st.stdout[-4000:] if st.stdout else "no status")
msg = (
    "Add K6 pilot postprocess, anion-ready pipeline, figures, and docs\n\n"
    "Keep ON_HOLD. Publish existing Hessian-confirmed K+(H2O)6 pilot data, "
    "xyzrender/Matplotlib figures, and tests. No 100-start production run."
)
r = subprocess.run(cmd_prefix + ["commit", "-m", msg], cwd=linux)
print("COMMIT", r.returncode)
sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=linux, text=True, capture_output=True)
print("HEAD", (sha.stdout or "").strip())
r = subprocess.run(["git", "push", "-u", "origin", "HEAD:main"], cwd=linux)
print("PUSH", r.returncode)
sys.exit(0 if r.returncode == 0 else r.returncode)
