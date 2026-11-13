"""EXPLORATORY, POST HOC, NOT PRE-REGISTERED (requested 1 Nov, after the frozen run 07-09)
-> results/exploratory/posthoc.json. Uses real waiting times, so the no-peeking guard must pass.

a. Among the primary data's short-window-measured intervals, the completed ones with T shorter than the
   source's post-glitch span w (D4f: 40 d for Espinoza+2011, 125 d for Wang+2000) are counted and
   excluded; M2 and M4 are refitted (D13 count applied to the subset, as for robustness cases) and b-hat
   and d-hat are reported with uncalibrated Wald and approximate calibrated intervals (primary factors).
b. d-hat (M4) per pulsar for every pulsar with at least 6 primary intervals, with intervals and method mix.
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
from glitchlib.analysis import primary_data, robustness_case  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_DIR = ROOT / "results" / "exploratory"
MIN_INTERVALS_PER_PULSAR = 6          # as requested for item b; exploratory, not a decision


def fmt(case):
    if case.get("status") != "run":
        return case.get("status")
    e = case["estimates_D4p"]
    return ", ".join(f"{k} = {e[k]['estimate']:.4f}, uncalibrated [{e[k]['uncalibrated'][0]:.4f}, {e[k]['uncalibrated'][1]:.4f}], "
                     f"approx. calibrated [{e[k]['calibrated_interval'][0]:.4f}, {e[k]['calibrated_interval'][1]:.4f}], "
                     f"signed root {e[k]['signed_root']:.3f}" for k in ("b", "d"))


def main():
    cfg = load_config()
    checksum = check_unblinding_allowed()
    data = primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL))
    prim = json.loads((ROOT / "results" / "primary.json").read_text(encoding="utf-8"))
    sec = json.loads((ROOT / "results" / "secondary.json").read_text(encoding="utf-8"))
    primary_estimates = {prim[t]["parameter"]: prim[t]["estimate"] for t in ("T1", "T2")}
    factors = {prim[t]["parameter"]: sec["intervals_b_d"][t]["factor"] for t in ("T1", "T2")}

    short = exploratory.short_window_shorter_than_span(data, cfg)
    sw = data["method_class"] == cfg["D4f.short_window.class"]
    spans = {k: v["w_days"] for k, v in cfg["D4f.short_window.sources"].items()}
    try:
        refit = robustness_case(data[~short], cfg, "short-window intervals with T < w excluded", None, factors,
                                primary_estimates)
    except models.FitError as err:
        refit = {"status": f"fit failed: {err}"}
    part_a = {"label": exploratory.LABEL, "spans_days": spans,
              "short_window_intervals": int(sw.sum()), "short_window_completed": int((sw & data["event"].astype(bool)).sum()),
              "excluded_T_below_span": int(short.sum()),
              "excluded_by_source": {k: int((short & (data["source_key"] == k)).sum()) for k in spans},
              "excluded_by_pulsar": {k: int(v) for k, v in data[short].groupby("pulsar").size().items()},
              "refit": refit}
    part_b = {"label": exploratory.LABEL, "min_intervals": MIN_INTERVALS_PER_PULSAR,
              "pulsars": exploratory.per_pulsar_cases(data, cfg, MIN_INTERVALS_PER_PULSAR, factors, primary_estimates)}
    record = {"label": exploratory.LABEL, "a_short_window_T_below_span": part_a, "b_d_hat_per_pulsar": part_b,
              "final_dataset_sha256": checksum, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "posthoc.json").write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")

    print(exploratory.LABEL)
    print(f"a. short-window intervals {part_a['short_window_intervals']} (completed {part_a['short_window_completed']}); "
          f"T < w excluded: {part_a['excluded_T_below_span']} {part_a['excluded_by_source']} {part_a['excluded_by_pulsar']}")
    if refit.get("status") == "run":
        print(f"   refit: n = {refit['n_intervals']} ({refit['n_completed']} completed, {refit['n_pulsars']} pulsars, "
              f"dropped {refit['pulsars_dropped']}); {fmt(refit)}")
    else:
        print(f"   refit: {refit['status']}")
    print(f"b. per pulsar (>= {MIN_INTERVALS_PER_PULSAR} intervals):")
    for row in part_b["pulsars"]:
        print(f"   {row['pulsar']}: n = {row['n_intervals']} ({row['n_completed']} completed), mix {row['method_mix']}; "
              f"{fmt(row['fit'])}")
    print("Wrote results/exploratory/posthoc.json")


if __name__ == "__main__":
    main()
