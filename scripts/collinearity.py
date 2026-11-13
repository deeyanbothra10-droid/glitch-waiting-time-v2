"""Collinearity of the two size terms (D15, D4b) -> results/collinearity.json. Allowed before freezing:
uses glitch sizes only, no waiting times.

D4b: Pearson correlation of ln(Δν/ν) and ln(Δν̇/ν̇) over usable glitches, for pulsars with at least
4 admissible usable glitches; the D15 rule uses the unweighted median |r| across those pulsars (T2 is
declared underpowered if it exceeds 0.9). Also reported: the pooled within-pulsar (demeaned)
correlation and the within-pulsar SD of ln R.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.config import load_config  # noqa: E402
from glitchlib.diagnostics import collinearity_summary  # noqa: E402
from glitchlib.intervals import add_derived  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

IN_PATH = ROOT / "data" / "processed" / "glitches_clean.csv"
OUT_PATH = ROOT / "results" / "collinearity.json"


def main():
    cfg = load_config()
    col = "D4b.diagnostics.collinearity"
    glitches = pd.read_csv(IN_PATH)
    usable = add_derived(glitches[glitches["usable"]])
    n_admissible = usable[usable["admissible"]].groupby("pulsar").size()
    minimum = cfg[f"{col}.min_admissible_usable_glitches_per_pulsar"]
    eligible = sorted(n_admissible[n_admissible >= minimum].index)
    data = usable[usable["pulsar"].isin(eligible)]
    if eligible:
        summary = collinearity_summary(data, cfg[f"{col}.correlation"], minimum)
        per = summary["per_pulsar"]
        median_abs_r, pooled = summary["median_abs_r"], summary["pooled_within_pulsar_r"]
    else:
        per, median_abs_r, pooled = pd.DataFrame(columns=["pulsar", "n", "r", "sd_ln_R"]), None, None
    threshold = cfg["D15.collinearity.underpowered_if_above"]
    record = {"correlation": cfg[f"{col}.correlation"], "over": cfg[f"{col}.over"],
              "min_admissible_usable_glitches_per_pulsar": minimum, "n_pulsars": len(eligible),
              "median_abs_r": median_abs_r, "pooled_within_pulsar_r": pooled,
              "D15_threshold": threshold,
              "median_abs_r_above_threshold": None if median_abs_r is None else bool(median_abs_r > threshold),
              "per_pulsar": per.to_dict("records"), "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}: {len(eligible)} pulsars with >= {minimum} admissible usable glitches.")
    if eligible:
        print(per.to_string(index=False))
        print(f"median |r| = {median_abs_r:.3f} (D15 threshold {threshold}); pooled within-pulsar r = {pooled:.3f}")
    else:
        print("No pulsar qualifies yet: admissibility needs the source audit's df1_method and exp_fitted (D9).")


if __name__ == "__main__":
    main()
