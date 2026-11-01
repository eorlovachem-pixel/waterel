#!/usr/bin/env python3
"""Rewrite SHA256SUMS.txt for published (non-gitignored) files."""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "SHA256SUMS.txt"


def main() -> None:
    p = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files"],
        capture_output=True, text=True, check=True,
    )
    paths = [ln.strip() for ln in p.stdout.splitlines() if ln.strip() and ln.strip() != "SHA256SUMS.txt"]
    lines = []
    for rel in sorted(paths):
        fp = ROOT / rel
        if not fp.is_file():
            continue
        digest = hashlib.sha256(fp.read_bytes()).hexdigest()
        lines.append(f"{digest}  {rel.replace(chr(92), '/')}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines)} hashes to {OUT}")


if __name__ == "__main__":
    main()
