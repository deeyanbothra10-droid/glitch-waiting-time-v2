"""The one real run of the primary analysis -> results/primary.json.

Preflight (glitchlib.preflight): prereg/FROZEN exists and data/processed/final_dataset.csv matches
prereg/data_freeze.json; the git tree is clean; the tests pass; results/validation.json is all PASS;
results/primary.json does not exist yet; no decision is still tagged provisional (D4p). Then D3 models on the primary data (admissible intervals of
D13-included pulsars) and the likelihood-ratio statistics of T1 and T2 with their signed roots, once
(D4l: the decisions, against null references simulated at the fitted sigma, are made in 08). A fit
failure stops the run (D4k). Prints the numbers; no reading.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import models, preflight  # noqa: E402
from glitchlib.analysis import primary, primary_data  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_REL = "results/primary.json"


def main():
    cfg = load_config()
    if not cfg["D4p.preflight_07_refuses_if_any_provisional"]:
        sys.exit("STOP: only D4p's refusal of provisional decisions is implemented.")
    preflight.require_all(OUT_REL, refuse_provisional_cfg=cfg)
    checksum = check_unblinding_allowed()
    try:
        result = primary(primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL)), cfg)
    except models.FitError as err:
        sys.exit(f"STOP: a fit failed on the real data, so the run stops (D4k): {err}")
    record = {**result, "final_dataset_sha256": checksum, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (ROOT / OUT_REL).write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    for test in ("T1", "T2"):
        r = result[test]
        print(f"{test}: {r['parameter']} = {r['estimate']:.4f} (SE {r['se']:.4f}, {result['interval_level']:.0%} "
              f"Wald interval {r['interval'][0]:.4f} to {r['interval'][1]:.4f}, uncalibrated); LR {r['statistic']:.3f}, df {r['df']}, "
              f"signed root {r['signed_root']:.3f}; chi-square p = {r['p_value']:.4g} (descriptive)")
    print("The T1 (D4l) and T2 (D4n) decisions and the calibrated intervals (D4m) use references simulated at "
          "D4n's matched sigma (scripts/08_secondary.py).")
    print(f"Gate: {result['gate']}; sigma = {result['sigma']:.4f}. Wrote {OUT_REL}")


if __name__ == "__main__":
    main()
