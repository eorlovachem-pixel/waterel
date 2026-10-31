#!/usr/bin/env python3
"""Diagnose WSL/Linux programs, Python packages, vendor SHAs and network.

Writes environment_check_report.json next to this script. Does not run a
production GWH series and does not treat `xtb --version` as a chemistry test.
"""
from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
import qc_checks  # noqa: E402

VENDOR_REPOS = {
    "Glorious-Water-Hole": "https://gitlab.com/madschumacher/glorious-water-hole",
    "CoCoS": "https://gitlab.com/madschumacher/cocoso",
    "MDAtomLocalEnvironment": "https://gitlab.com/madschumacher/mdatomlocalenvironment",
    "xtb": "https://github.com/grimme-lab/xtb",
}
URLS = [
    "https://github.com",
    "https://gitlab.com",
    "https://pypi.org",
    "https://conda.anaconda.org",
    "https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz",
]
STANDALONE_XTB = Path.home() / "xtb-standalone" / "xtb-dist" / "bin" / "xtb"


def run(cmd, timeout=20):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
        return {
            "cmd": cmd,
            "returncode": p.returncode,
            "stdout": (p.stdout or "")[-4000:],
            "stderr": (p.stderr or "")[-2000:],
        }
    except FileNotFoundError as e:
        return {"cmd": cmd, "returncode": None, "error": f"FileNotFoundError: {e}"}
    except subprocess.TimeoutExpired:
        return {"cmd": cmd, "returncode": None, "error": "timeout"}


def git_info(path: Path):
    if not (path / ".git").exists():
        return {"path": str(path), "present": path.exists(), "git": False}
    return {
        "path": str(path),
        "present": True,
        "git": True,
        "HEAD": run(["git", "-C", str(path), "rev-parse", "HEAD"]).get("stdout", "").strip(),
        "describe": run(["git", "-C", str(path), "describe", "--tags", "--always"]).get("stdout", "").strip(),
        "log1": run(["git", "-C", str(path), "log", "-1", "--format=%H %ci %s"]).get("stdout", "").strip(),
        "remote": run(["git", "-C", str(path), "remote", "get-url", "origin"]).get("stdout", "").strip(),
        "status": run(["git", "-C", str(path), "status", "-sb"]).get("stdout", "").strip(),
    }


def http_head(url: str, timeout=15):
    try:
        req = Request(url, method="HEAD")
        with urlopen(req, timeout=timeout) as resp:
            return {"url": url, "ok": True, "status": getattr(resp, "status", None), "final": resp.geturl()}
    except Exception as e:
        return {"url": url, "ok": False, "error": f"{type(e).__name__}: {e}"}


def ping_host(host: str):
    try:
        socket.setdefaulttimeout(4)
        socket.getaddrinfo(host, 443)
        return {"host": host, "dns": True}
    except Exception as e:
        return {"host": host, "dns": False, "error": f"{type(e).__name__}: {e}"}


def pkg_version(name):
    try:
        mod = __import__(name)
        return getattr(mod, "__version__", "imported")
    except Exception as e:
        return f"MISSING {type(e).__name__}: {e}"


def which_xtb_details(path: Path):
    if not path.exists():
        return {"path": str(path), "exists": False}
    info = {
        "path": str(path.resolve()),
        "exists": True,
        "sha256": qc_checks.sha256_file(path),
        "size": path.stat().st_size,
        "version_cmd": run([str(path), "--version"]),
    }
    ldd = run(["ldd", str(path)])
    info["ldd_excerpt"] = (ldd.get("stdout") or ldd.get("error") or "")[:1500]
    info["statically_linked"] = "statically linked" in (run(["file", str(path)]).get("stdout") or "")
    return info


def main():
    report = {
        "timestamp_UTC": datetime.now(timezone.utc).isoformat(),
        "cwd": str(Path.cwd()),
        "root": str(ROOT),
        "platform": run(["uname", "-a"]),
        "os_release": Path("/etc/os-release").read_text() if Path("/etc/os-release").exists() else None,
        "user": os.environ.get("USER"),
        "home": str(Path.home()),
        "python": {
            "executable": sys.executable,
            "version": sys.version,
            "numpy": pkg_version("numpy"),
            "scipy": pkg_version("scipy"),
            "matplotlib": pkg_version("matplotlib"),
            "pytest": pkg_version("pytest"),
        },
        "compilers": {
            "gcc": run(["gcc", "--version"]),
            "gfortran": run(["gfortran", "--version"]),
            "cmake": run(["cmake", "--version"]),
            "meson": run(["meson", "--version"]),
            "ninja": run(["ninja", "--version"]),
            "make": run(["make", "--version"]),
        },
        "git": run(["git", "--version"]),
        "miniforge": {
            "dir_exists": (Path.home() / "miniforge3").is_dir(),
            "conda": run([str(Path.home() / "miniforge3" / "bin" / "conda"), "--version"])
            if (Path.home() / "miniforge3" / "bin" / "conda").exists()
            else {"present": False},
        },
        "xtb_on_path": shutil.which("xtb"),
        "system_xtb": which_xtb_details(Path("/usr/bin/xtb")) if Path("/usr/bin/xtb").exists() else {"exists": False},
        "standalone_xtb": which_xtb_details(STANDALONE_XTB),
        "libgfortran": run(["dpkg-query", "-W", "-f=${Package} ${Version}\\n", "libgfortran5"]),
        "vendors": {},
        "gwh_euler": None,
        "k06_vecfle": None,
        "network_dns": [ping_host(h) for h in ["github.com", "gitlab.com", "pypi.org", "conda.anaconda.org"]],
        "network_http": [http_head(u) for u in URLS],
        "warnings": [],
        "blockers": [],
    }

    for name in VENDOR_REPOS:
        report["vendors"][name] = git_info(ROOT / "vendors" / name)

    gwh = ROOT / "vendors" / "Glorious-Water-Hole" / "gloriouswaterhole.py"
    if gwh.is_file():
        report["gwh_euler"] = qc_checks.gwh_uses_intrinsic_zyz(gwh)
        if report["gwh_euler"].get("zyz_ZYZ_discrepancy"):
            report["blockers"].append(
                "GWH contains both from_euler('zyz') and from_euler('ZYZ'). "
                "SciPy lowercase is extrinsic, uppercase is intrinsic. Source was not rewritten."
            )
        elif not report["gwh_euler"].get("author_method_intact", report["gwh_euler"].get("intrinsic_zyz")):
            report["blockers"].append(
                "GWH does not contain the author from_euler('zyz') string "
                "(commit 7411811374b6422a15f23d6d86a47f2145f8ac47). "
                "Do not rewrite zyz to ZYZ without an author decision."
            )
        else:
            report["warnings"].append(
                "Author GWH uses from_euler('zyz'). SciPy docs treat lowercase as extrinsic "
                "and 'ZYZ' as intrinsic. Source was not rewritten; confirm convention with Denis."
            )
    else:
        report["blockers"].append("GWH source missing; run setup_vendors.sh")

    k06 = ROOT / "inputs" / "K_06" / "packing.vecfle"
    if k06.is_file():
        report["k06_vecfle"] = qc_checks.validate_vecfle_bytes(
            k06, "K", 6, 3.04, "06ff2c3bf4b43e0d0382a684539dd5b3edf8e618500c857b9d4aa045700a11e5"
        )
        report["k06_vecfle"]["expected_natoms"] = qc_checks.expected_natoms(6)
        report["k06_vecfle"]["expected_formula"] = qc_checks.expected_formula("K", 6)
        report["k06_vecfle"]["charge"] = 1
    else:
        report["blockers"].append("missing inputs/K_06/packing.vecfle")

    path_xtb = shutil.which("xtb")
    report["resolved_working_xtb"] = None
    try:
        resolved = qc_checks.resolve_working_xtb(os.environ.get("XTB_BIN"), kit=ROOT)
        report["resolved_working_xtb"] = str(resolved)
    except RuntimeError as exc:
        report["blockers"].append(str(exc))
    if path_xtb == "/usr/bin/xtb":
        if report["resolved_working_xtb"] and report["resolved_working_xtb"] != "/usr/bin/xtb":
            report["warnings"].append(
                "PATH still lists /usr/bin/xtb, but pipeline/workflow will use "
                f"{report['resolved_working_xtb']}. Source env_xtb.sh in interactive shells."
            )
        else:
            report["warnings"].append(
                "PATH resolves xtb to /usr/bin/xtb 6.7.1 (Ubuntu, libgfortran 16). "
                "ANCOPT is known to crash. Source env_xtb.sh so the static 6.7.1 (edcfbbe) binary is used."
            )
    if not report["standalone_xtb"].get("exists") and not report["resolved_working_xtb"]:
        report["warnings"].append(
            "Standalone Intel-static xtb not found at ~/xtb-standalone/xtb-dist/bin/xtb. "
            "Install from GitHub release v6.7.1 without replacing /usr/bin/xtb."
        )
    if not report["network_http"][0]["ok"] and not report["standalone_xtb"].get("exists"):
        report["blockers"].append("GitHub unreachable and no local xtb tarball/binary")

    preflight = ROOT / "preflight_runs" / "latest.json"
    if preflight.is_file():
        latest = json.loads(preflight.read_text(encoding="utf-8"))
        report["preflight_latest"] = {
            "file": str(preflight),
            "usable_for_free_optimization": latest.get("usable_for_free_optimization"),
            "timestamp_UTC": latest.get("timestamp_UTC"),
        }
        if not latest.get("usable_for_free_optimization"):
            report["blockers"].append("latest xTB preflight did not pass opt+hess")
    else:
        report["warnings"].append("No preflight_runs/latest.json yet. Run: python3 code/xtb_preflight.py")

    out = ROOT / "environment_check_report.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(
        {
            "report": str(out),
            "python": report["python"]["version"].split()[0],
            "path_xtb": path_xtb,
            "standalone_xtb": report["standalone_xtb"].get("path"),
            "resolved_working_xtb": report.get("resolved_working_xtb"),
            "gwh_zyz": None if report["gwh_euler"] is None else report["gwh_euler"].get("author_method_intact"),
            "scipy_interprets_author_zyz_as": None if report["gwh_euler"] is None else report["gwh_euler"].get("scipy_interprets_author_zyz_as"),
            "k06_ok": report["k06_vecfle"] is not None,
            "warnings": report["warnings"],
            "blockers": report["blockers"],
        },
        indent=2,
        ensure_ascii=False,
    ))
    if report["blockers"]:
        raise SystemExit("environment_check found blockers; see environment_check_report.json")


if __name__ == "__main__":
    main()
