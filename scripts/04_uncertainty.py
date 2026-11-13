"""Analysis dataset and data freeze -> data/processed/final_dataset.csv and prereg/data_freeze.json.

From data/processed/intervals.csv (03_intervals.py): adds ln R (R = (dff/1000)/df1), sigma_lnR =
sqrt((dff_err/dff)^2 + (df1_err/df1)^2), ln dnu_frac, ln dnudot_frac, and tau_sd_days from
data/processed/pulsar_params.csv (01_pulsar_params.py); the Δν̇ method class is already there from the
source audit (02_clean.py). Then records the dataset's SHA-256 in prereg/data_freeze.json.
Prints counts per pulsar and per method class only.
"""
import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, write_freeze_record  # noqa: E402
from glitchlib.intervals import add_derived  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

INTERVALS_PATH = ROOT / "data" / "processed" / "intervals.csv"
PARAMS_PATH = ROOT / "data" / "processed" / "pulsar_params.csv"


def main():
    load_config()
    intervals = pd.read_csv(INTERVALS_PATH)
    params = pd.read_csv(PARAMS_PATH)[["name", "tau_sd_days"]].rename(columns={"name": "pulsar"})
    final = add_derived(intervals).merge(params, on="pulsar", how="left")
    out = ROOT / FINAL_DATASET_REL
    final.to_csv(out, index=False)
    record = write_freeze_record(date.today().isoformat(), git_commit(), ROOT)

    included = final[final["included_pulsar"]]
    print(f"Wrote {FINAL_DATASET_REL}: {len(final)} intervals, {len(included)} in included pulsars.")
    print(f"Data freeze recorded in prereg/data_freeze.json: SHA-256 {record['sha256']}")
    print(included.groupby("pulsar").size().rename("intervals").to_string())
    print(final.groupby("method_class").size().rename("intervals").to_string())


if __name__ == "__main__":
    main()
