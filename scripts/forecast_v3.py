"""Post-freeze forward forecast v3 (NOT part of the pre-registration) -> results/forecast/forecast_v3.json and
results/forecast/README_v3.md. Supersedes forecast v2 (kept unchanged; tag forecast-v2) before any outcome was
known; v1 is also kept unchanged.

As v2 (scripts/forecast_v2.py), except the scored set: a forecastable pulsar is scored only if
  (a) it is covered by the D4g Jodrell Bank rule,
  (b) it has no documented obs_end before MJD 61343 in data/meta/obs_end.csv (an obs_end equal to the last
      catalogued glitch epoch is the D4g placeholder, "no censoring information", not a documented end), and
  (c) MJD 61343 - epoch(g) is no longer than its longest catalogued wait between consecutive glitches
(glitchlib.forecast.scoring_eligibility). Scored pulsars are conditioned on no glitch before MJD 61343. Every
other forecastable pulsar is reported, not scored, with the rules it fails and their numbers, and conditioned as
in v2: on no glitch before its documented obs_end, or unconditioned if obs_end equals g. The scoring rule is
otherwise v2's, including its two added rules. Changes no committed result.
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
from glitchlib.cleaning import jodrell_bank_pulsars  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import RAW_GLITCHES_PATH, git_commit, read_meta, sha256_file  # noqa: E402

OUT_DIR = ROOT / "results" / "forecast"
SNAPSHOT_DATE, SNAPSHOT_MJD = "2026-10-30", 61343.0
SCORING_DATE, SCORING_MJD = "2028-11-01", 62076.0
PROBS = {"median": 0.5, "w50_lo": 0.25, "w50_hi": 0.75, "w90_lo": 0.05, "w90_hi": 0.95}
UNSCORED_LABEL = "reported, not scored"
NOTES = {"J0537-6910": "its glitch rate (one glitch every 100-200 d) makes a 580 d glitch-free stretch to the "
                       "2026-10-30 snapshot implausible, so the snapshot is not used for conditioning"}

RULE = [
    "Label: post-freeze forward forecast v3; NOT part of the pre-registration (prereg-v2). Supersedes forecast v2 "
    "(tag forecast-v2) and v1 (tag forecast-v1), both kept unchanged, before any outcome was known.",
    f"Snapshot: the Jodrell Bank glitch catalogue downloaded {SNAPSHOT_DATE} (MJD {SNAPSHOT_MJD:.0f}); for each of the "
    "21 included pulsars, its last catalogued glitch g. Pulsars whose g lacks dff > 0 or df1 > 0 are not forecastable.",
    "Forecast: ln T ~ Normal(mu, sigma), T = wait after g in days, from the frozen fits in results/primary.json: "
    "M2 mu = c_p + b ln R_g with R_g = (dff/1000)/df1; M0 (baseline) mu = c_p; each model's own sigma-hat. Median and "
    "central 50% (0.25-0.75) and 90% (0.05-0.95) windows of the conditional distribution.",
    "Scored set: a forecastable pulsar is scored only if (a) it is covered by the D4g Jodrell Bank rule; (b) it has no "
    f"documented obs_end before MJD {SNAPSHOT_MJD:.0f} in data/meta/obs_end.csv (an obs_end equal to the last "
    "catalogued glitch epoch is the D4g placeholder, not a documented end of monitoring); and (c) "
    f"{SNAPSHOT_MJD:.0f} - epoch(g) is no longer than its longest catalogued wait between consecutive glitches. "
    f"Scored pulsars are conditioned on T > {SNAPSHOT_MJD:.0f} - epoch(g).",
    "Reported, not scored: every other forecastable pulsar, with the rules it fails and their numbers; conditioned on "
    "T > documented obs_end - epoch(g), or unconditioned if obs_end equals epoch(g).",
    f"Scoring, on {SCORING_DATE} (MJD {SCORING_MJD:.0f}): download the catalogue again. For each scored pulsar the "
    "outcome is the first glitch after g: the first catalogued glitch of that pulsar whose epoch exceeds g's epoch by "
    "more than g's catalogued epoch uncertainty (0 if none is given).",
    f"Void: if that glitch's epoch is before MJD {SNAPSHOT_MJD:.0f}, the pulsar is void and not scored.",
    f"Score per scored pulsar: log predictive density of the observed wait (density of T in days) under the "
    f"conditional forecast, or, if no glitch by MJD {SCORING_MJD:.0f}, log of the conditional survival at the wait to "
    "that date; M2 minus M0. Total score = sum over scored pulsars; positive means M2 forecasts better.",
    "Coverage: for the 50% and 90% windows, the fraction of scored pulsars whose observed wait lies inside the "
    "window (ends included). A pulsar with no glitch by the scoring date counts as outside if the window ends before "
    "the scoring date, and as not determinable (excluded from that coverage) otherwise; the counts are reported.",
    "Functions: glitchlib.forecast (conditional_quantiles, log_predictive, coverage, outcome_for, forecast_window, "
    "scoring_eligibility).",
]


def main():
    cfg = load_config()
    prim = json.loads((ROOT / "results" / "primary.json").read_text(encoding="utf-8"))
    m2, m0 = prim["fits"]["M2"]["estimates"], prim["fits"]["M0"]["estimates"]
    glitches = pd.read_csv(ROOT / "data" / "processed" / "glitches_clean.csv")
    included = glitches[glitches["included_pulsar"]]
    jodrell = set(jodrell_bank_pulsars(included, cfg))
    obs = read_meta("obs_end").set_index("name")
    last = included.sort_values("mjd", kind="mergesort").groupby("pulsar").tail(1).sort_values("pulsar")
    scored, unscored, not_forecastable = {}, {}, {}
    for g in last.itertuples():
        missing = [name for name, v in (("dff", g.dff), ("df1", g.df1)) if not (pd.notna(v) and v > 0)]
        if missing:
            not_forecastable[g.pulsar] = {"g_mjd": float(g.mjd), "reason": f"last glitch lacks {' and '.join(missing)} > 0"}
            continue
        epochs = included.loc[included["pulsar"] == g.pulsar, "mjd"].to_numpy(float)
        obs_end = float(obs.loc[g.pulsar, "obs_end_mjd"])
        elig = forecast.scoring_eligibility(g.pulsar in jodrell, obs_end, float(epochs.max()), epochs, float(g.mjd),
                                            SNAPSHOT_MJD)
        if elig["scored"]:
            t0, conditioning, target = SNAPSHOT_MJD - float(g.mjd), f"no glitch before MJD {SNAPSHOT_MJD:.0f} (snapshot)", scored
        elif obs_end > float(g.mjd):
            t0, conditioning, target = obs_end - float(g.mjd), f"no glitch before documented obs_end MJD {obs_end:.1f}", unscored
        else:
            t0, conditioning, target = 0.0, "unconditioned (obs_end equals g)", unscored
        ln_r = float(np.log((g.dff / 1000.0) / g.df1))
        entry = {"g": {"mjd": float(g.mjd), "date": forecast.mjd_to_date(g.mjd),
                       "mjd_err": None if pd.isna(g.mjd_err) else float(g.mjd_err), "dff": float(g.dff), "df1": float(g.df1),
                       "ln_R": ln_r, "source_key": g.source_key, "method_class": g.method_class},
                 "jodrell_bank": g.pulsar in jodrell, "obs_end_mjd": obs_end, "obs_end_source": str(obs.loc[g.pulsar, "source"]),
                 "eligibility": elig, "conditioning": conditioning,
                 "M2": forecast.forecast_window(m2[f"c[{g.pulsar}]"] + m2["b"] * ln_r, m2["sigma"], t0, float(g.mjd), PROBS),
                 "M0": forecast.forecast_window(m0[f"c[{g.pulsar}]"], m0["sigma"], t0, float(g.mjd), PROBS)}
        if target is unscored:
            entry["label"] = f"{UNSCORED_LABEL}: fails " + "; ".join(f"({k}) {v}" for k, v in elig["fails"].items())
        if g.pulsar in NOTES:
            entry["note"] = NOTES[g.pulsar]
        target[g.pulsar] = entry
    record = {"label": "POST-FREEZE FORWARD FORECAST v3, NOT PART OF THE PRE-REGISTRATION", "supersedes": ["forecast-v2", "forecast-v1"],
              "rule": RULE, "snapshot": {"date": SNAPSHOT_DATE, "mjd": SNAPSHOT_MJD},
              "scoring": {"date": SCORING_DATE, "mjd": SCORING_MJD}, "n_scored": len(scored), "scored_pulsars": list(scored),
              "frozen_fits": {"M2": {"b": m2["b"], "sigma": m2["sigma"]}, "M0": {"sigma": m0["sigma"]}, "source": "results/primary.json"},
              "inputs_sha256": {"data/raw/glitches.csv": sha256_file(RAW_GLITCHES_PATH),
                                "data/processed/glitches_clean.csv": sha256_file(ROOT / "data" / "processed" / "glitches_clean.csv"),
                                "data/meta/obs_end.csv": sha256_file(ROOT / "data" / "meta" / "obs_end.csv"),
                                "results/primary.json": sha256_file(ROOT / "results" / "primary.json"),
                                "results/forecast/forecast_v1.json": sha256_file(OUT_DIR / "forecast_v1.json"),
                                "results/forecast/forecast_v2.json": sha256_file(OUT_DIR / "forecast_v2.json")},
              "scored": scored, "reported_not_scored": unscored, "not_forecastable": not_forecastable,
              "git_commit": git_commit(), "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (OUT_DIR / "forecast_v3.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    def table(entries, with_reason):
        head = "| Pulsar | Last glitch g | Conditioning | Model | Median | 50% window | 90% window |"
        rows = [head + (" Fails |" if with_reason else ""), "|---|---|---|---|---|---|---|" + ("---|" if with_reason else "")]
        for p, f in entries.items():
            for m in ("M2", "M0"):
                d = f[m]["dates"]
                reason = (" " + "; ".join(f"({k}) {v}" for k, v in f["eligibility"]["fails"].items()) + " |") if with_reason else ""
                rows.append(f"| {p} | {f['g']['date']} | {f['conditioning']} | {m} | {d['median']} | "
                            f"{d['w50_lo']} to {d['w50_hi']} | {d['w90_lo']} to {d['w90_hi']} |" + (reason if m == "M2" else (" |" if with_reason else "")))
        return rows
    lines = ["# Forward forecast v3 (post-freeze; NOT part of the pre-registration)", "",
             "Written by `scripts/forecast_v3.py`; the forecasts are in `forecast_v3.json`. Forecasts v1 and v2 "
             "(`forecast_v1.json`, `README.md`, `forecast_v2.json`, `README_v2.md`; tags `forecast-v1`, `forecast-v2`) "
             "are kept unchanged.", "",
             "## Why v2 was superseded", "",
             "v2 scored every Jodrell Bank pulsar conditioned on no glitch before the 2026-10-30 snapshot. For some of them "
             "that assumption is not supported: J1814-1744's last glitch is 5,717 d before the snapshot, 4.9x its longest "
             "catalogued wait (1,166 d); and J0205+6449, J0631+1036 and J1814-1744 have a documented obs_end in Dec 2018 "
             "(end of the Basu et al. 2022 Jodrell Bank sample). v3 adds rules (b) and (c) below. It was made before any "
             "outcome was known.", "",
             "## Scored set", "",
             f"{len(scored)} pulsars are scored: {', '.join(scored)}. With so few pulsars and windows this wide, the score "
             "is a weak test; it is reported whatever it shows.", "",
             "## Rule (fixed before any outcome)", ""]
    lines += [f"{k}. {r}" for k, r in enumerate(RULE, 1)]
    lines += ["", "## Scored forecasts", ""] + table(scored, False)
    lines += ["", "## Reported, not scored", ""] + table(unscored, True)
    lines += [""] + [f"Note, {p}: {f['note']}." for p, f in unscored.items() if "note" in f]
    lines += ["", "## Not forecastable", ""] + [f"- {p}: {v['reason']} (g at MJD {v['g_mjd']:.1f})" for p, v in not_forecastable.items()]
    (OUT_DIR / "README_v3.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("SCORED:", ", ".join(scored))
    for p, f in unscored.items():
        print(f"not scored {p}: " + "; ".join(f"({k}) {v}" for k, v in f["eligibility"]["fails"].items()) + f" [{f['conditioning']}]")
    print("Wrote results/forecast/forecast_v3.json and results/forecast/README_v3.md")


if __name__ == "__main__":
    main()
