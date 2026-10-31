#!/usr/bin/env python3
"""Copy selected kit files from the Windows OneDrive tree to the Linux ext4 kit."""
from pathlib import Path

win = (
    Path("/mnt/c/Users/antiz/OneDrive")
    / "Документы"
    / "Cursor_K6_research_handoff"
    / "Cursor_research"
    / "Hydration_17_workflows_original"
    / "Hydration_unassigned_2026-10-30"
)
kit = Path.home() / "hydration-research" / "Hydration_unassigned_2026-10-30"
FILES = [
    "code/pipeline.py",
    "code/qc_checks.py",
    "code/freq_parser.py",
    "code/mode_follow.py",
    "code/xtb_preflight.py",
    "code/shell_descriptors.py",
    "code/visualize.py",
    "code/figures_matplotlib.py",
    "code/export_canonical.py",
    "code/mol_matplotlib.py",
    "code/write_indexes.py",
    "workflow.py",
    "env_xtb.sh",
    "CASE_STATUS.csv",
    "CASE_STATUS.json",
    "RUN_ONE_APPROVED_CASE.sh",
    "environment_check.py",
    "check_cocos_robustness.py",
    "COMMAND_K06_SERIES_100.md",
    "METHODOLOGY_AUDIT.md",
    "ON_HOLD_AND_SCOPE.md",
    "MESSAGE_TO_DENIS.md",
    "ASSIGNMENT_AUDIT.md",
    "CHANGELOG.md",
    "tests/test_protocol_commands.py",
    "tests/test_euler_zyz.py",
    "tests/test_xyz_and_logs.py",
    "tests/test_pipeline_confirmed.py",
    "tests/test_inputs_k06.py",
    "tests/test_shell_and_anions.py",
    "tests/test_visualize.py",
    "README_WINDOWS_WSL.md",
    "PINNED_VERSIONS.md",
    "run_k06_pilot_2.py",
    "SHA256SUMS.txt",
    "_sync_pipeline_to_linux.py",
    "tools/sync_to_linux.py",
    "tools/apply_cleanup.py",
    ".gitignore",
]


def main():
    for rel in FILES:
        src = win / rel
        dst = kit / rel
        if not src.is_file():
            print("skip missing", rel)
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        data = src.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        dst.write_bytes(data)
        print("wrote", dst, len(data))


if __name__ == "__main__":
    main()
