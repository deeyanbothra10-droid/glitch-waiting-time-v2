"""Run the whole pipeline in order, stopping at the first error.

Steps that need the freeze are skipped with a message until it exists: 03_intervals and 04_uncertainty
run after prereg/FROZEN exists (the data freeze follows the pre-registration freeze); 07-09 need the
full guard and their own preflight; 10_figures always runs and itself skips the real-data figures.
Writes results/environment.txt with the Python and package versions.

Reproduction mode (post-freeze tooling): python run_all.py --reproduce DIR clones this repository's HEAD
into DIR (which must not exist or be empty), runs every step there, 07-10 included, with GLITCH_REPRODUCE=1
so that the 07-09 preflight skips only its run-once and clean-tree checks and the test suite skips only the
data-freeze-commit assertion (the no-peeking guard and its checksum check stay in force), then compares
every file the run wrote with the committed version (glitchlib.reproduce) and writes
DIR/reproduction_report.json. Everything is written under DIR; this repository is only read.
"""
import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from glitchlib.guards import FROZEN_REL, unblinding_allowed  # noqa: E402
from glitchlib.reproduce import REPRODUCE_ENV, compare_tree  # noqa: E402

ALWAYS = ["00_provenance", "01_pulsar_params", "emulation_fixture", "02_clean", "method_counts", "collinearity",
          "05_power", "06_validate"]
AFTER_PREREG_FREEZE = ["03_intervals", "04_uncertainty"]
AFTER_GUARD = ["07_primary", "08_secondary", "09_robustness"]


def environment():
    import matplotlib
    import numpy
    import pandas
    import scipy
    import yaml
    lines = [f"python {platform.python_version()}", f"numpy {numpy.__version__}", f"pandas {pandas.__version__}",
             f"scipy {scipy.__version__}", f"matplotlib {matplotlib.__version__}", f"pyyaml {yaml.__version__}",
             f"platform {platform.platform()}"]
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "environment.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def step(name):
    print(f"=== {name}", flush=True)
    result = subprocess.run([sys.executable, str(ROOT / "scripts" / f"{name}.py")], cwd=ROOT)
    if result.returncode != 0:
        sys.exit(f"STOP: scripts/{name}.py failed (exit {result.returncode}).")


def pipeline():
    environment()
    for name in ALWAYS:
        step(name)
    if (ROOT / FROZEN_REL).exists():
        for name in AFTER_PREREG_FREEZE:
            step(name)
    else:
        print(f"Skipped {', '.join(AFTER_PREREG_FREEZE)}: {FROZEN_REL} does not exist (pre-registration not frozen).")
    if unblinding_allowed():
        for name in AFTER_GUARD:
            step(name)
    else:
        print(f"Skipped {', '.join(AFTER_GUARD)}: the guard does not pass (not frozen, or data checksum differs).")
    step("10_figures")


def reproduce(directory):
    directory = Path(directory).resolve()
    if directory.exists() and any(directory.iterdir()):
        sys.exit(f"STOP: {directory} exists and is not empty; reproduction needs a fresh folder.")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    subprocess.run(["git", "clone", "--quiet", str(ROOT), str(directory)], check=True)
    cloned_at = time.time()
    time.sleep(1.1)                       # every file the run writes is then strictly newer than the clone
    env = dict(os.environ, **{REPRODUCE_ENV: "1"})
    result = subprocess.run([sys.executable, str(directory / "run_all.py")], cwd=directory, env=env)
    report = compare_tree(directory, cloned_at)
    report.update({"commit": head, "pipeline_exit_code": result.returncode})
    (directory / "reproduction_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Reproduction of {head[:7]} in {directory}: pipeline exit {result.returncode}; "
          f"{len(report['reproduced'])} files reproduced, {len(report['differs'])} differ, "
          f"{len(report['not_regenerated'])} not regenerated, {len(report['new_files'])} new.")
    for rel, diffs in report["differs"].items():
        print(f"  DIFFERS {rel}: {len(diffs)} value(s); first: {diffs[0]}")
    for rel in report["not_regenerated"]:
        print(f"  not regenerated: {rel}")
    print(f"Wrote {directory / 'reproduction_report.json'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reproduce", type=Path, help="fresh folder: clone HEAD there, rerun everything, compare")
    args = parser.parse_args()
    if args.reproduce:
        reproduce(args.reproduce)
    else:
        pipeline()


if __name__ == "__main__":
    main()
