"""Post-freeze forward forecast v1 (NOT part of the pre-registration) -> results/forecast/forecast_v1.json and
results/forecast/README.md.

Snapshot: the catalogue downloaded 2026-10-30 (MJD 61343), data/raw/glitches.csv as cleaned by 02_clean.py.
For each D13-included pulsar, its last catalogued glitch g: if g has dff > 0 and df1 > 0, the wait to the next
glitch is forecast from the frozen fits in results/primary.json (M2: b-hat, that pulsar's c-hat and sigma-hat,
with ln R_g = ln((dff/1000)/df1); M0 as the baseline), conditioned on no glitch before MJD 61343
(glitchlib.forecast). Conditional median and the central 50% and 90% windows, as waits and calendar dates.
Otherwise the pulsar is listed as not forecastable, with the reason. The scoring rule is written into the
record and the README. Changes no committed result.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import forecast  # noqa: E402
from glitchlib.io import RAW_GLITCHES_PATH, git_commit, sha256_file  # noqa: E402

OUT_DIR = ROOT / "results" / "forecast"
SNAPSHOT_DATE, SNAPSHOT_MJD = "2026-10-30", 61343.0
SCORING_DATE, SCORING_MJD = "2028-11-01", 62076.0
PROBS = {"median": 0.5, "w50_lo": 0.25, "w50_hi": 0.75, "w90_lo": 0.05, "w90_hi": 0.95}

RULE = [
    "Label: post-freeze forward forecast v1; NOT part of the pre-registration (prereg-v2).",
    f"Snapshot: the Jodrell Bank glitch catalogue downloaded {SNAPSHOT_DATE} (MJD {SNAPSHOT_MJD:.0f}); for each of the "
    "21 included pulsars, its last catalogued glitch g. Pulsars whose g lacks dff > 0 or df1 > 0 are not forecastable.",
    "Forecast: ln T ~ Normal(mu, sigma), T = wait after g in days, from the frozen fits in results/primary.json: "
    "M2 mu = c_p + b ln R_g with R_g = (dff/1000)/df1; M0 (baseline) mu = c_p; each model's own sigma-hat. "
    f"Conditioned on T > t0 = {SNAPSHOT_MJD:.0f} - epoch(g). Median and central 50% (0.25-0.75) and 90% (0.05-0.95) "
    "windows of the conditional distribution.",
    f"Scoring, on {SCORING_DATE} (MJD {SCORING_MJD:.0f}): download the catalogue again. For each forecast pulsar the "
    "outcome is the first glitch after g: the first catalogued glitch of that pulsar whose epoch exceeds g's epoch by "
    "more than g's catalogued epoch uncertainty (0 if none is given).",
    f"Void: if that glitch's epoch is before MJD {SNAPSHOT_MJD:.0f}, the pulsar is void and not scored.",
    f"Score per pulsar: log predictive density of the observed wait (density of T in days) under the conditional "
    f"forecast, or, if no glitch by MJD {SCORING_MJD:.0f}, log of the conditional survival at the wait to that date; "
    "M2 minus M0. Total score = sum over scored pulsars; positive means M2 forecasts better.",
    "Coverage: for the 50% and 90% windows, the fraction of scored pulsars whose observed wait lies inside the "
    "window (ends included). A pulsar with no glitch by the scoring date counts as outside if the window ends before "
    "the scoring date, and as not determinable (excluded from that coverage) otherwise; the counts are reported.",
    "Functions: glitchlib.forecast (conditional_quantiles, log_predictive, coverage, outcome_for).",
]


def window(mu, sigma, t0, g_mjd):
    q = dict(zip(PROBS, forecast.conditional_quantiles(mu, sigma, t0, list(PROBS.values()))))
    return {"mu": mu, "sigma": sigma, "waits_days": q,
            "dates": {k: forecast.mjd_to_date(g_mjd + v) for k, v in q.items()}}


def main():
    prim = json.loads((ROOT / "results" / "primary.json").read_text(encoding="utf-8"))
    m2, m0 = prim["fits"]["M2"]["estimates"], prim["fits"]["M0"]["estimates"]
    glitches = pd.read_csv(ROOT / "data" / "processed" / "glitches_clean.csv")
    included = glitches[glitches["included_pulsar"]]
    last = included.sort_values("mjd", kind="mergesort").groupby("pulsar").tail(1).sort_values("pulsar")
    forecasts, not_forecastable = {}, {}
    for g in last.itertuples():
        missing = [name for name, v in (("dff", g.dff), ("df1", g.df1)) if not (pd.notna(v) and v > 0)]
        if missing:
            not_forecastable[g.pulsar] = {"g_mjd": float(g.mjd), "reason": f"last glitch lacks {' and '.join(missing)} > 0",
                                          "dff": None if pd.isna(g.dff) else float(g.dff), "df1": None if pd.isna(g.df1) else float(g.df1)}
            continue
        ln_r = float(np.log((g.dff / 1000.0) / g.df1))
        t0 = SNAPSHOT_MJD - float(g.mjd)
        forecasts[g.pulsar] = {
            "g": {"mjd": float(g.mjd), "date": forecast.mjd_to_date(g.mjd), "mjd_err": None if pd.isna(g.mjd_err) else float(g.mjd_err),
                  "dff": float(g.dff), "df1": float(g.df1), "ln_R": ln_r, "source_key": g.source_key, "method_class": g.method_class},
            "t0_days": t0,
            "M2": window(m2[f"c[{g.pulsar}]"] + m2["b"] * ln_r, m2["sigma"], t0, float(g.mjd)),
            "M0": window(m0[f"c[{g.pulsar}]"], m0["sigma"], t0, float(g.mjd))}
    record = {"label": "POST-FREEZE FORWARD FORECAST v1, NOT PART OF THE PRE-REGISTRATION", "rule": RULE,
              "snapshot": {"date": SNAPSHOT_DATE, "mjd": SNAPSHOT_MJD}, "scoring": {"date": SCORING_DATE, "mjd": SCORING_MJD},
              "frozen_fits": {"M2": {"b": m2["b"], "sigma": m2["sigma"]}, "M0": {"sigma": m0["sigma"]},
                              "source": "results/primary.json"},
              "inputs_sha256": {"data/raw/glitches.csv": sha256_file(RAW_GLITCHES_PATH),
                                "data/processed/glitches_clean.csv": sha256_file(ROOT / "data" / "processed" / "glitches_clean.csv"),
                                "results/primary.json": sha256_file(ROOT / "results" / "primary.json")},
              "forecasts": forecasts, "not_forecastable": not_forecastable, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "forecast_v1.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    lines = ["# Forward forecast v1 (post-freeze; NOT part of the pre-registration)", "",
             "Written by `scripts/forecast_v1.py`; the forecasts are in `forecast_v1.json`.", "", "## Rule (fixed before any outcome)", ""]
    lines += [f"{k}. {r}" for k, r in enumerate(RULE, 1)]
    lines += ["", "## Forecasts (conditional on no glitch before 2026-10-30)", "",
              "| Pulsar | Last glitch g | Model | Median | 50% window | 90% window |", "|---|---|---|---|---|---|"]
    for p, f in forecasts.items():
        for m in ("M2", "M0"):
            d = f[m]["dates"]
            lines.append(f"| {p} | {f['g']['date']} | {m} | {d['median']} | {d['w50_lo']} to {d['w50_hi']} | "
                         f"{d['w90_lo']} to {d['w90_hi']} |")
    lines += ["", "## Not forecastable", ""] + [f"- {p}: {v['reason']} (g at MJD {v['g_mjd']:.1f})" for p, v in not_forecastable.items()]
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for p, f in forecasts.items():
        for m in ("M2", "M0"):
            d = f[m]["dates"]
            print(f"{p:11s} g {f['g']['date']} {m}: median {d['median']}; 50% {d['w50_lo']} to {d['w50_hi']}; "
                  f"90% {d['w90_lo']} to {d['w90_hi']}")
    for p, v in not_forecastable.items():
        print(f"{p:11s} not forecastable: {v['reason']}")
    print("Wrote results/forecast/forecast_v1.json and results/forecast/README.md")


if __name__ == "__main__":
    main()
