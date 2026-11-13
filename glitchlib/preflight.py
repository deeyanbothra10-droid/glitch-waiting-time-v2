"""Preflight checks for the real-data scripts 07, 08 and 09. Each check returns a message on failure;
require_all stops the script with every failing reason listed."""
import json
import subprocess
import sys
from pathlib import Path

from glitchlib.guards import PeekingError, check_unblinding_allowed
from glitchlib.reproduce import reproducing

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def check_frozen(root):
    """prereg/FROZEN exists and final_dataset.csv matches prereg/data_freeze.json (glitchlib.guards)."""
    try:
        check_unblinding_allowed(root)
    except PeekingError as err:
        return str(err)
    return None


def check_git_clean(root):
    result = subprocess.run(["git", "status", "--porcelain"], cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        return f"git status failed: {result.stderr.strip()}"
    if result.stdout.strip():
        return "the git working tree is not clean (commit or remove the changes first)."
    return None


def check_tests_pass(root):
    result = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=root, capture_output=True, text=True)
    if result.returncode != 0:
        return "the test suite does not pass: " + result.stdout.strip().splitlines()[-1]
    return None


def check_validation_passed(root):
    path = Path(root) / "results" / "validation.json"
    if not path.exists():
        return "results/validation.json does not exist (run scripts/06_validate.py)."
    record = json.loads(path.read_text(encoding="utf-8"))
    failed = [c["name"] for c in record.get("criteria", []) if c.get("result") != "PASS"]
    if not record.get("criteria") or failed:
        return f"results/validation.json is not all PASS (failing: {failed or 'no criteria recorded'})."
    return None


def check_no_provisional(cfg):
    """D4p item 5 (07): no decision may still be tagged provisional."""
    tagged = cfg.provisional_items()
    if tagged:
        return f"decisions still tagged provisional: {', '.join(tagged)} (D4p: 07 refuses)."
    return None


def check_absent(root, output):
    if (Path(root) / output).exists():
        return f"{output} already exists: this analysis runs once."
    return None


def require_all(output, root=PROJECT_ROOT, run_tests=True, refuse_provisional_cfg=None):
    """Stop unless every precondition holds: frozen and checksum-matched data, a clean git tree, a passing
    test suite, an all-PASS validation, and no existing output file; with refuse_provisional_cfg (07, D4p),
    also no decision still tagged provisional in that config. In reproduction mode only (run_all.py
    --reproduce, GLITCH_REPRODUCE=1, in a fresh clone) the clean-tree and run-once checks are skipped; the
    freeze, checksum, tests, validation and provisional checks still apply."""
    repro = reproducing()
    checks = [check_frozen(root), None if repro else check_git_clean(root), check_tests_pass(root) if run_tests else None,
              check_validation_passed(root), None if repro else check_absent(root, output),
              check_no_provisional(refuse_provisional_cfg) if refuse_provisional_cfg is not None else None]
    failures = [c for c in checks if c]
    if failures:
        sys.exit("STOP: preflight failed:\n  - " + "\n  - ".join(failures))
