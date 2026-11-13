"""Waiting-time intervals -> data/processed/intervals.csv (built, not analysed).

From data/processed/glitches_clean.csv (02_clean.py): within each pulsar, ordered by mjd, every usable
glitch starts an interval ending at the next glitch (usable or not); the last one is censored at the
end of monitoring (data/meta/obs_end.csv). D11: an interval spanning a listed gap is censored at the gap
start and flagged. D10: a pulsar with no end-of-monitoring date loses its last interval (logged to
logs/dropped_last_intervals.csv). D13: intervals touching a glitch with epoch error > 30 d or unknown
are flagged. Prints counts only; computes no statistic involving waiting times.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.config import load_config  # noqa: E402
from glitchlib.intervals import build_intervals, epoch_flags  # noqa: E402
from glitchlib.io import gap_starts, read_meta  # noqa: E402

IN_PATH = ROOT / "data" / "processed" / "glitches_clean.csv"
OUT_PATH = ROOT / "data" / "processed" / "intervals.csv"
DROPPED_PATH = ROOT / "logs" / "dropped_last_intervals.csv"


def main():
    cfg = load_config()
    glitches = pd.read_csv(IN_PATH)
    obs = read_meta("obs_end")
    obs_end = dict(zip(obs["name"], obs["obs_end_mjd"].astype(float)))
    gaps, undated = gap_starts(read_meta("monitoring_gaps"))
    if cfg["D11.interval_spanning_gap"] != "censor_at_gap_start":
        sys.exit("STOP: only D11 interval_spanning_gap = censor_at_gap_start is implemented here.")
    intervals, dropped = build_intervals(glitches, obs_end, gaps)
    if len(intervals):
        intervals["epoch_flag"] = epoch_flags(intervals, cfg["D13.flag_mjd_err_above_days"],
                                              cfg["D13.flag_mjd_err_unknown"])
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    intervals.to_csv(OUT_PATH, index=False)
    DROPPED_PATH.parent.mkdir(parents=True, exist_ok=True)
    dropped.to_csv(DROPPED_PATH, index=False)

    print(f"Wrote {OUT_PATH.relative_to(ROOT)}: {len(intervals)} intervals; "
          f"{len(dropped)} last intervals dropped for want of an end date ({DROPPED_PATH.relative_to(ROOT)}).")
    if len(undated):
        print(f"Monitoring gaps not applied (no start date yet): {', '.join(undated['name'])}")
    if len(intervals):
        counts = intervals.groupby("pulsar").agg(
            completed=("event", "sum"), censored=("event", lambda e: int((~e.astype(bool)).sum())),
            gap_flagged=("gap_flag", "sum"), epoch_flagged=("epoch_flag", "sum"),
            included=("included_pulsar", "first"))
        print(counts[counts["included"]].to_string())


if __name__ == "__main__":
    main()
