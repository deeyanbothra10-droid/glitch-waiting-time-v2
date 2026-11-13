"""Secondary analyses and diagnostics on the real data -> results/secondary.json. Same preflight as 07
(with results/secondary.json absent), and results/primary.json must exist.

D4 secondary (implied beta from M1; within-pulsar Spearman, completed intervals only; D4a permutation tests
with D4b.simulation.real_run_n_permutations and the D4a seed; LOO-CV for M0-M4; D5 Poisson baseline), D15
placebo (fit, Spearman, permutation), D15/D4b missingness (conditional logit of "df1 missing" on ln of the
following wait over glitches of included pulsars whose following wait is completed; descriptive: df1 <= 0),
and the reference distributions (S1, S1m, S4-mix and the null runs S2, S3 and S3m, with the S3/S3m size
coefficient set to the real M3 estimate, 06's replicate counts): D4n item 1, each run at its own matched
generating sigma, the one at which its median fitted sigma-hat of the primary model equals the real
sigma-hat (bisection on D4p's range [0.01, 2]; reported, including any scenario run at the lower end and
any bisection that did not converge within D4p's 30 steps, which then uses the final bracket's midpoint). From them: the D4l T1 decision,
the D4n T2 decision (d-hat below the smaller 2.5th percentile of d-hat under S3 and S3m; the signed root
reported alongside, never deciding), the D4j outcome with the D4l bounds (powers: the S1m rejection rates)
with D4n's report-only short-window label, and, for missingness, the percentile of the real log odds ratio
in its S1m distribution (D4i), from estimable replicates only (D4k). A missingness check that is not estimable is
reported as such; any other fit failure stops the run (D4k).
Implied beta is reported per pulsar with its interval and D7 class, and decides nothing (D4i).
Prints numbers only.
"""
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import analysis, decision, diagnostics, models, preflight, validation  # noqa: E402
from glitchlib.analysis import primary_data, secondary, summarise_secondary  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.intervals import following_waits  # noqa: E402
from glitchlib.io import gap_starts, git_commit, read_meta  # noqa: E402

OUT_REL = "results/secondary.json"


def missingness(cfg):
    m = "D4b.diagnostics.missingness"
    glitches = pd.read_csv(ROOT / "data" / "processed" / "glitches_clean.csv")
    glitches = glitches[glitches["included_pulsar"]].copy()
    obs = read_meta("obs_end")
    gaps, _ = gap_starts(read_meta("monitoring_gaps"))
    waits = following_waits(glitches, dict(zip(obs["name"], obs["obs_end_mjd"].astype(float))), gaps)
    data = glitches.join(waits)
    data = data[data["following_t"].notna()].assign(ln_wait=lambda d: np.log(d["following_t"]),
                                                    df1_missing=lambda d: d["df1"].isna(),
                                                    df1_le_0=lambda d: d["df1"] <= 0)
    level = cfg["D4b.interval_estimates.level"]

    def check(frame, outcome):
        try:
            return diagnostics.missingness_check(frame, outcome, "ln_wait", "following_completed", level)
        except models.FitError as err:              # D4k: separation is reported, not a stop
            return {"status": cfg["D4k.missingness_not_estimable.reported_as"], "reason": str(err)}
    main = check(data, "df1_missing")
    if "p_value" in main:
        main["significant"] = bool(main["p_value"] < cfg[f"{m}.alpha"])
    return {"df1_missing": main, "df1_le_0_descriptive": check(data[data["df1"].notna()], "df1_le_0")}


def references(cfg, sigma_hat, s3_slope):
    """D4n item 1: every reference scenario at its own matched generating sigma (replaces D4k item 3), then
    run with 06's replicate count; emulated scenarios calibrated at their matched sigma (D4e)."""
    workers = os.cpu_count()
    templates, targets = validation.real_templates(cfg)
    index = validation.real_sigma_index(cfg)
    jobs = validation.reference_jobs(cfg) + validation.null_jobs(cfg, s3_slope)
    scenario_of = {job.name: job.name.removesuffix("-null") for job in jobs}
    if set(scenario_of.values()) != set(cfg["D4n.matched_sigma.scenarios"]):
        sys.exit(f"STOP: the reference jobs {sorted(scenario_of.values())} are not D4n's scenarios "
                 f"{cfg['D4n.matched_sigma.scenarios']}.")
    raw, matched, calibration = {}, {}, {}
    for job in jobs:
        try:
            match = validation.match_job_sigma(templates, targets, job, sigma_hat, cfg, workers)
        except validation.SigmaMatchError as err:
            sys.exit(f"STOP: D4n sigma matching failed for {job.name} (a search step had no analysed replicate; "
                     f"no decision covers this): {err}")
        matched[scenario_of[job.name]] = match
        if match["at_lower_end"]:
            print(f"D4p: {scenario_of[job.name]}'s fitted sigma-hat at sigma = {match['sigma']} already exceeds the real "
                  f"sigma-hat; run at sigma = {match['sigma']}", flush=True)
        if not match["converged"]:
            print(f"D4p: {scenario_of[job.name]}: bisection did not converge; run at the midpoint of the final bracket, "
                  f"sigma = {match['sigma']}", flush=True)
        noise = None
        if job.scenario.emulated:
            calibration[job.name] = validation.calibrate(templates, targets, cfg, match["sigma"], index, workers)
            noise = {c["pulsar"]: c["noise_sd"] for c in calibration[job.name]}
        raw[job.name] = validation.run_job(templates, job, cfg, match["sigma"], index, noise, workers)
    failures = {name: validation.failure_check(rec, cfg) for name, rec in raw.items()}
    nulls = {k: raw[k] for k in ("S2-null", "S3-null", "S3m-null")}
    factors = validation.interval_factors(nulls, cfg)
    raw, thr, refs, powers = validation.classified(raw, nulls, cfg)
    sigma_report = {scenario_of[name]: {"matched_generating_sigma": matched[scenario_of[name]]["sigma"],
                                        "run_at_lower_end": matched[scenario_of[name]]["at_lower_end"],
                                        "bisection_converged": matched[scenario_of[name]]["converged"],
                                        "bisection_steps": matched[scenario_of[name]]["steps"],
                                        "mean_fitted_sigma_hat": validation.sigma_hat_mean(rec),
                                        "median_fitted_sigma_hat": float(validation.ok(rec)["sigma_hat"].median())}
                    for name, rec in raw.items()}
    return (refs, thr, factors, calibration, powers, validation.missingness_reference(raw["S1m"]), failures,
            sigma_report)


def main():
    preflight.require_all(OUT_REL)
    primary_path = ROOT / "results" / "primary.json"
    if not primary_path.exists():
        sys.exit("STOP: results/primary.json does not exist; run scripts/07_primary.py first.")
    cfg = load_config()
    checksum = check_unblinding_allowed()
    prim = json.loads(primary_path.read_text(encoding="utf-8"))
    data = primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL))
    tau = dict(zip(data["pulsar"], data["tau_sd_days"]))
    try:
        sec = secondary(data, cfg, tau, cfg["D4b.simulation.real_run_n_permutations"], cfg["D4a.permutation.seed"], True)
    except models.FitError as err:
        sys.exit(f"STOP: a fit failed on the real data, so the run stops (D4k): {err}")
    refs, thresholds, factors, calibration, powers, s1m_missingness, reference_failures, sigma_report = references(
        cfg, prim["sigma"], prim["fits"]["M3"]["estimates"]["a"])                 # D4l/D4m: S3, S3m slope = real M3 estimate
    level = cfg["D4b.interval_estimates.level"]
    intervals = {test: {"calibrated": validation.calibrated_interval(prim[test]["estimate"], prim[test]["se"],
                                                                     factors[prim[test]["parameter"]], level),
                        "factor": factors[prim[test]["parameter"]],
                        cfg["D4m.calibrated_intervals.plain_wald_label"]: prim[test]["interval"]}
                 for test in ("T1", "T2")}
    decisions = {"T1": bool(prim["T1"]["signed_root"] > thresholds["T1"]),                      # D4l
                 "T2": bool(prim["T2"]["estimate"] < thresholds["T2_d_cut"]),                       # D4n
                 "T2_signed_root_D4m_report_only": bool(prim["T2"]["signed_root"]
                                                        > thresholds["T2_signed_root_D4m_report_only"])}
    try:
        d_no_sw = analysis.d_without_short_window(data, cfg)                                        # D4n item 3
    except models.FitError as err:
        sys.exit(f"STOP: a fit failed on the real data, so the run stops (D4k): {err}")
    outcome = decision.classify(decisions["T1"], decisions["T2"], prim["T1"]["estimate"], sec["placebo"]["b"], refs,
                                powers["T2"], cfg, d_no_sw)
    outcome["d_without_short_window"] = d_no_sw
    miss = missingness(cfg)
    miss["S1m_reference"] = {k: s1m_missingness[k] for k in ("replicates", "estimable", "not_estimable")}
    miss["S1m_reference"]["log_or_interval"] = validation.percentile_interval(s1m_missingness["values"],
                                                                               cfg["D4b.interval_estimates.level"])
    if "log_odds_ratio" in miss["df1_missing"]:
        miss["df1_missing"]["percentile_in_S1m"] = decision.percentile_of(miss["df1_missing"]["log_odds_ratio"],
                                                                          s1m_missingness["values"])
    glitches = pd.read_csv(ROOT / "data" / "processed" / "glitches_clean.csv")
    usable = glitches[glitches["usable"]]
    median_dnu = (usable["dff"] * 1e-9).groupby(usable["pulsar"]).median()
    beta = sec["implied_beta"]["table"]
    beta["D7_class"] = np.where(beta["pulsar"].map(median_dnu) >= cfg["D7.vela_like_min_dnu_frac"], "Vela-like", "Crab-like")
    record = {**summarise_secondary(sec), "missingness": miss, "references": decision.references_for_json(refs),
              "directional_thresholds": thresholds, "directional_decisions": decisions,
              "powers_at_matched_sigma": powers, "real_sigma_hat": prim["sigma"], "reference_sigma_D4n": sigma_report,
              "intervals_b_d": intervals,
              "reference_calibration": calibration, "reference_failures": reference_failures, "outcome": outcome,
              "final_dataset_sha256": checksum, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (ROOT / OUT_REL).write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    print(json.dumps({k: record[k] for k in ("directional_decisions", "placebo", "outcome", "powers_at_matched_sigma")},
                     indent=2, default=float))
    print(f"Wrote {OUT_REL}")


if __name__ == "__main__":
    main()
