"""EXPLORATORY, POST HOC, NOT PRE-REGISTERED (requested 4 Nov, after the frozen run 07-09)
-> results/exploratory/heterogeneity_separate_a.json. Uses real waiting times, so the guard must pass.

Heterogeneity of d with pulsar-specific size coefficients: on the pulsars with at least 6 primary intervals,
the LR test of M4 with a separate a and a separate d per pulsar against M4 with a separate a per pulsar and
one common d. Reports the statistic, df, p, the common d-hat with its SE (reduced model) and each pulsar's
d-hat with its SE (full model). Numbers only; decides nothing; changes no committed result.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import exploratory  # noqa: E402
from glitchlib.analysis import primary_data  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_DIR = ROOT / "results" / "exploratory"
MIN_INTERVALS_PER_PULSAR = 6          # as requested; exploratory


def main():
    cfg = load_config()
    checksum = check_unblinding_allowed()
    data = primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL))
    t = exploratory.per_pulsar_lr_test(data, cfg, MIN_INTERVALS_PER_PULSAR, ("a",), ("a", "d"))
    reduced, full = t["reduced"], t["full"]
    record = {"label": exploratory.LABEL, "min_intervals": MIN_INTERVALS_PER_PULSAR, "pulsars": t["pulsars"],
              "n_intervals": t["n_intervals"], "statistic": t["statistic"], "df": t["df"], "p_value": t["p_value"],
              "loglik_reduced": t["loglik_reduced"], "loglik_full": t["loglik_full"],
              "common_d": {"estimate": reduced["estimates"]["d"], "se": reduced["se"]["d"]},
              "per_pulsar_d": {p: {"estimate": full["estimates"][f"d[{p}]"], "se": full["se"][f"d[{p}]"]} for p in t["pulsars"]},
              "reduced_model": reduced, "full_model": full,
              "final_dataset_sha256": checksum, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "heterogeneity_separate_a.json").write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    print(exploratory.LABEL)
    print(f"Separate a and d per pulsar vs separate a and common d: {len(t['pulsars'])} pulsars, {t['n_intervals']} intervals; "
          f"LR = {t['statistic']:.3f}, df = {t['df']}, p = {t['p_value']:.4g}")
    print(f"Common d (separate a) = {record['common_d']['estimate']:.4f} (SE {record['common_d']['se']:.4f})")
    for p, v in record["per_pulsar_d"].items():
        print(f"  d[{p}] = {v['estimate']:.4f} (SE {v['se']:.4f})")
    print("Wrote results/exploratory/heterogeneity_separate_a.json")


if __name__ == "__main__":
    main()
