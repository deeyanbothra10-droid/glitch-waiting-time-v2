"""Spin-down timescale per pulsar -> data/processed/pulsar_params.csv.

tau_sd = nu / (2 |nudot|), in days, from data/meta/pulsar_params.csv. A row with nudot >= 0 stops the
script. Rows whose nu or nudot are not filled in yet get NaN and are listed.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.cleaning import tau_sd_days  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import read_meta  # noqa: E402

OUT_PATH = ROOT / "data" / "processed" / "pulsar_params.csv"


def main():
    load_config()
    params = read_meta("pulsar_params")
    params["tau_sd_days"] = tau_sd_days(params["nu_hz"], params["nudot_hz_s"])
    params["tau_sd_years"] = params["tau_sd_days"] / 365.25
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    params.to_csv(OUT_PATH, index=False)
    missing = params.loc[params["tau_sd_days"].isna(), "name"].tolist()
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}: {len(params)} pulsars, {len(params) - len(missing)} with tau_sd.")
    if missing:
        print(f"nu or nudot not filled in yet for {len(missing)}: {', '.join(missing)}")


if __name__ == "__main__":
    main()
