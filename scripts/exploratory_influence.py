"""EXPLORATORY, POST HOC, NOT PRE-REGISTERED (requested 1 Nov, after the frozen run 07-09)
-> results/exploratory/influence.json. Uses real waiting times, so the no-peeking guard must pass.

a. Leave-one-interval-out d-hat of M4 over the primary data (D13 count reapplied after each removal, as for
   every subset): the minimum, the maximum, and the 5 intervals whose removal moves d-hat most.
b. Heterogeneity of d: LR test, on the pulsars with at least 6 primary intervals, of M4 with a separate d per
   pulsar against M4 with a common d.
Numbers only; decides nothing; changes no committed result.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import exploratory, models  # noqa: E402
from glitchlib.analysis import primary_data  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_DIR = ROOT / "results" / "exploratory"
MOST_INFLUENTIAL = 5                  # as requested; exploratory
MIN_INTERVALS_PER_PULSAR = 6          # as requested; exploratory


def main():
    cfg = load_config()
    checksum = check_unblinding_allowed()
    data = primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL))
    full = float(models.fit(cfg["D3.T2_model"], data, cfg["D6.primary_errors_in_lnR"]).estimates[cfg["D2.T2.parameter"]])
    rows = exploratory.leave_one_out_d(data, cfg)
    loo = pd.DataFrame(rows)
    loo = loo.join(data[["pulsar", "mjd", "t", "event", "method_class"]].reset_index(drop=True), on="position")
    ok = loo[loo["d_without"].notna()].copy()
    ok["change"] = ok["d_without"] - full
    top = ok.reindex(ok["change"].abs().sort_values(ascending=False).index).head(MOST_INFLUENTIAL)
    part_a = {"label": exploratory.LABEL, "d_full": full, "fits": int(len(loo)), "failed": int(loo["d_without"].isna().sum()),
              "fits_where_D13_dropped_a_pulsar": int(ok["pulsars_dropped"].map(len).gt(0).sum()),
              "d_without_min": float(ok["d_without"].min()), "d_without_max": float(ok["d_without"].max()),
              "most_influential": [{"pulsar": r.pulsar, "glitch_mjd": float(r.mjd), "T_days": float(r.t),
                                    "completed": bool(r.event), "method_class": r.method_class,
                                    "d_without": float(r.d_without), "change": float(r.change),
                                    "pulsars_dropped": list(r.pulsars_dropped)} for r in top.itertuples()],
              "all": [{"pulsar": r.pulsar, "glitch_mjd": float(r.mjd), "d_without": None if pd.isna(r.d_without) else float(r.d_without)}
                      for r in loo.itertuples()]}
    try:
        part_b = {"label": exploratory.LABEL, "min_intervals": MIN_INTERVALS_PER_PULSAR,
                  **exploratory.separate_d_test(data, cfg, MIN_INTERVALS_PER_PULSAR)}
    except models.FitError as err:
        part_b = {"label": exploratory.LABEL, "status": f"fit failed: {err}"}
    record = {"label": exploratory.LABEL, "a_leave_one_out_d": part_a, "b_separate_d_per_pulsar": part_b,
              "final_dataset_sha256": checksum, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "influence.json").write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")

    print(exploratory.LABEL)
    print(f"a. d-hat (M4, all {len(data)} intervals) = {full:.4f}; leave-one-out over {part_a['fits']} fits "
          f"({part_a['failed']} failed; D13 dropped a pulsar in {part_a['fits_where_D13_dropped_a_pulsar']}): "
          f"min {part_a['d_without_min']:.4f}, max {part_a['d_without_max']:.4f}")
    for r in part_a["most_influential"]:
        print(f"   {r['pulsar']} glitch MJD {r['glitch_mjd']:.1f}: T = {r['T_days']:.1f} d "
              f"({'completed' if r['completed'] else 'censored'}), {r['method_class']}; d-hat without = {r['d_without']:.4f} "
              f"(change {r['change']:+.4f}){'; D13 dropped ' + ', '.join(r['pulsars_dropped']) if r['pulsars_dropped'] else ''}")
    if "statistic" in part_b:
        print(f"b. separate d per pulsar vs common d ({len(part_b['pulsars'])} pulsars with >= {MIN_INTERVALS_PER_PULSAR} "
              f"intervals, {part_b['n_intervals']} intervals): LR = {part_b['statistic']:.3f}, df = {part_b['df']}, "
              f"p = {part_b['p_value']:.4g}; common d = {part_b['common_d']:.4f}")
        for p, v in part_b["separate_d"].items():
            print(f"   d[{p}] = {v['estimate']:.4f} (SE {v['se']:.4f})")
    else:
        print(f"b. {part_b['status']}")
    print("Wrote results/exploratory/influence.json")


if __name__ == "__main__":
    main()
