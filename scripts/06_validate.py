"""Validation of the whole analysis on synthetic data with the real pulsar structure -> results/validation.json.
Synthetic data only: no real waiting times are read.

For every scenario (S1, S1 with measurement noise off, S2, S3, S1m, and the S4 variants S4w, S4l, S4g,
S4g0, S4sw for both short-window source types, and S4-mix) at each sigma of the validation grid (D14 plus
D4k), runs exactly the primary, secondary and diagnostic analyses of 07/08 (glitchlib.analysis) and
reports: mean slope, rejection rates of the D4l directional tests (thresholds: 97.5th percentiles of the
signed-root LR statistic under S2 for T1 and S3 for T2, 1,000 replicates per sigma with their own seeds),
the D4j outcome rates with the D4l bounds (references: S1, S1m and S4-mix at the same sigma; powers: the
S1m rejection rates), the D4l report-only diagnostics, the S1m missingness reference (estimable
replicates only, D4k), preference over size-only (AIC; LOO-CV on the first D4b replicates), forward and
placebo Spearman, permutation rejection rates (999 permutations on the first 200 replicates, D4b), and the
failed replicates with their reasons. Then checks the D14 criteria with the D4h false-positive bound and
the D4k failed-replicate limit, and prints PASS or FAIL for each.

Each finished job is saved to results/validation_partial/; a rerun resumes from there (refused if the
commit or config checksum differs). results/validation.json is written only when every job is done.
With --reuse DIR (repeatable), jobs are read read-only from earlier runs' partial directories, searched in
order, and only the missing jobs are run and saved in --store; records saved before D4l get the S1/S3/S1m
diagnostic reruns (same seeds, checked to reproduce the saved b-hat). The source of every job is recorded.
D4m adds S3m (S3 through the emulation) as an evaluation scenario and as a T2 null, the A2-rate check under
S3m, the interval calibration factors, and the S1m-versus-S3m fractions (report only).
D4n: T2 rejects when d-hat is below the smaller of the 2.5th percentiles of d-hat under S3-null and S3m-null
at the same sigma; the D4l/D4m signed-root T2 and the outcomes it would give are reported alongside (never
deciding). Every power and outcome rate is reported against the generating sigma and the mean fitted
sigma-hat of the primary model; A2 replicates get the "short-window-sensitive" label where the saved
records contain the refit without short-window glitches. With --saved-only, the run stops instead of
simulating anything that is not already saved.
"""
import argparse
import json
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import decision, validation  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_PATH = ROOT / "results" / "validation.json"
PARTIAL_DIR = ROOT / "results" / "validation_partial"


def summarise(rec, refs, cfg, t2_power, t2_power_d4m, sigma):
    good = validation.ok(rec)
    level = cfg["D4b.interval_estimates.level"]
    alpha = cfg["D4.per_test_alpha"]
    out = {"replicates": int(len(rec)), "analysed": int(len(good)), "discarded": int(rec["discarded"].sum()),
           "pulsars_dropped_total": int(rec["pulsars_dropped"].sum()), "failures": validation.failure_check(rec, cfg),
           "seconds_max": float(rec["seconds"].max()) if "seconds" in rec else None}
    if not len(good):
        return out
    results = validation.replicate_outcomes(good, refs, t2_power, cfg)
    old_t2 = [r["outcome"] for r in validation.replicate_outcomes(good, refs, t2_power_d4m, cfg, "T2_rejects_D4m")]
    a2 = [r for r in results if r["outcome"] == "A2"]
    out.update({
        "b_mean": float(good["b"].mean()), "b_median": float(good["b"].median()),
        "b_interval": validation.percentile_interval(good["b"], level),
        "d_mean": float(good["d"].mean()), "placebo_b_mean": float(good["placebo_b"].mean()),
        "T1_rejection_rate": float(good["T1_rejects"].mean()), "T2_rejection_rate": float(good["T2_rejects"].mean()),
        "T2_rejection_rate_signed_root_D4m_report_only": float(good["T2_rejects_D4m"].mean()),
        "T1_rejection_rate_chi2_old": float(good["T1_rejects_chi2"].mean()),
        "T2_rejection_rate_chi2_old": float(good["T2_rejects_chi2"].mean()),
        "outcomes": decision.rates([r["outcome"] for r in results]),
        "outcomes_with_D4m_T2_report_only": decision.rates(old_t2),
        "sigma": {"generating": sigma, "mean_fitted_sigma_hat": validation.sigma_hat_mean(rec)},     # D4n item 4
        "short_window_label": ({"A2_replicates": len(a2), "labelled": sum(r["short_window_label"] is not None for r in a2)}
                               if "d_without_short_window" in good else
                               "not available: these saved replicates predate D4n (needs new simulation)"),
        "attenuation_label_rate": float(np.mean([r["attenuation_label"] is not None for r in results])),
        "aic_prefers_M2_over_M3": float(good["aic_M2_below_M3"].mean()),
        "aic_prefers_M4_over_M3": float(good["aic_M4_below_M3"].mean()),
        "spearman_forward_mean": float(np.nanmean(good["spearman_forward"])),
        "spearman_placebo_mean": float(np.nanmean(good["spearman_placebo"])),
        "median_completed_intervals": float(good["n_completed"].median()), "median_pulsars": float(good["n_pulsars"].median())})
    old_t2 = good[good["T2_rejects_chi2"].astype(bool)]
    diagnostics = {"old_T2_rejections": int(len(old_t2)),
                   "old_T2_rejections_fraction_d_negative": float((old_t2["d"] < 0).mean()) if len(old_t2) else None,
                   "b_sd": float(good["b"].std(ddof=1)), "b_mean_se": float(good["b_se"].mean()),
                   "d_sd": float(good["d"].std(ddof=1))}
    if "d_se" in good:
        diagnostics["d_mean_se"] = float(good["d_se"].mean())
    slopes = {c.replace("meas_slope_", ""): float(good[c].mean()) for c in good.columns if c.startswith("meas_slope_")}
    if slopes:
        diagnostics["measured_on_true_ln_dnudot_slope_mean_by_class"] = slopes
    out["d4l_report_only"] = diagnostics
    if "perm_p" in good:
        perm = good.dropna(subset=["perm_p"])
        out.update({"permutation_replicates": int(len(perm)),
                    "permutation_rejection_rate": float((perm["perm_p"] < alpha).mean()),
                    "placebo_significantly_negative_rate": float(perm["placebo_significantly_negative"].astype(bool).mean())})
    if "missingness_log_or" in good:
        m = validation.missingness_reference(rec)
        out["missingness"] = {"estimable": m["estimable"], "not_estimable": m["not_estimable"],
                              "log_or_mean": float(np.mean(m["values"])) if m["estimable"] else None,
                              "log_or_interval": validation.percentile_interval(m["values"], level)}
    if "loo_M4" in good:
        loo = good.dropna(subset=["loo_M4"])
        out.update({"loo_replicates": int(len(loo)),
                    "loo_prefers_M4_over_M3": float((loo["loo_M4"] > loo["loo_M3"]).mean()),
                    "loo_prefers_M2_over_M3": float((loo["loo_M2"] > loo["loo_M3"]).mean())})
    return out


def criteria(per_sigma, cfg):
    """D14 criteria at each sigma with the D4h false-positive bound (rejection rate of the D4l directional
    tests <= alpha + 3 sqrt(alpha (1 - alpha) / n), n = analysed replicates), and the D4k failed-replicate
    limit for every job, null and diagnostic jobs included."""
    pc, fp = "D14.pass_criteria", "D4h.false_positive_bound"
    truth = cfg["D4b.simulation.S1.slope"]
    bound = lambda scen: validation.false_positive_bound(cfg[f"{fp}.alpha"], cfg[f"{fp}.mc_sd_multiplier"],
                                                         scen["analysed"])
    out = []
    for sigma, entry in per_sigma.items():
        s = entry["scenarios"]
        checks = [("S2 T1 false-positive rate", s["S2"]["T1_rejection_rate"], bound(s["S2"])),
                  ("S2 T2 false-positive rate", s["S2"]["T2_rejection_rate"], bound(s["S2"])),
                  ("S3 T2 false-positive rate", s["S3"]["T2_rejection_rate"], bound(s["S3"])),
                  ("S1 noise off: |median b - 1|", abs(s["S1-noise-off"]["b_median"] - truth),
                   cfg[f"{pc}.S1_bhat_median_max_abs_deviation_noise_free"])]
        checks.append(("S3m A2 rate (D4m)", s["S3m"]["outcomes"]["A2"],
                       validation.false_positive_bound(cfg[f"{fp}.alpha"], cfg[f"{fp}.mc_sd_multiplier"],
                                                       s["S3m"]["analysed"])))
        checks += [(f"{name}: failed-replicate fraction (D4k)", f["fraction"], f["limit"])
                   for name, f in entry["failures_all_jobs"].items()]
        for name, value, bound_value in checks:
            out.append({"name": f"{name} (sigma {sigma})", "value": value, "bound": bound_value,
                        "result": "PASS" if value <= bound_value else "FAIL"})
    return out


def run_or_load(store, templates, job, cfg, sigma, s_idx, noise, workers, label=None, reuse=(), sources=None,
                saved_only=False):
    """A job's records: from the first --reuse directory that has it (read-only), else from this run's
    store (resume), else run now and saved there. sources records where each job came from."""
    name = label or job.name
    saved, directory = validation.load_saved(reuse, sigma, name)
    if saved is not None:
        if sources is not None:
            sources[f"{sigma}/{name}"] = str(directory)
        return saved
    if sources is not None:
        sources[f"{sigma}/{name}"] = str(store.dir)
    saved = store.load_job(sigma, name)
    if saved is not None:
        print(f"sigma={sigma} {name}: resumed {len(saved)} replicates", flush=True)
        return saved
    if saved_only:
        raise SystemExit(f"STOP (--saved-only): sigma={sigma} {name} is not saved anywhere and would need new simulation.")
    rec = validation.run_job(templates, job, cfg, sigma, s_idx, noise, workers)
    store.save_job(sigma, name, rec)
    print(f"sigma={sigma} {name}: {len(rec)} replicates, {int(rec['failed'].sum())} failed; saved", flush=True)
    return rec


def merge_diagnostics(old, diag):
    """Reuse mode: the diagnostic rerun has the same seeds, so its datasets must be identical; check b-hat
    and take d-hat's SE and the measurement slopes from it."""
    a, b = validation.ok(old).set_index("replicate"), validation.ok(diag).set_index("replicate")
    common = a.index.intersection(b.index)
    if len(common) != len(a) or not np.allclose(a.loc[common, "b"], b.loc[common, "b"], rtol=0, atol=1e-9):
        raise SystemExit("STOP: the diagnostic rerun did not reproduce the saved replicates (b-hat differs).")
    extra = [c for c in diag.columns if c == "d_se" or c.startswith("meas_slope_")]
    return old.merge(diag[["replicate"] + extra], on="replicate", how="left")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse", type=Path, action="append", default=[],
                        help="an earlier run's partial directory to read jobs from (read-only); repeatable, "
                             "searched in order")
    parser.add_argument("--store", type=Path, default=PARTIAL_DIR,
                        help="where this run saves the jobs it runs (resumable)")
    parser.add_argument("--saved-only", action="store_true",
                        help="stop instead of simulating any job or calibration that is not already saved")
    args = parser.parse_args()
    cfg = load_config()
    workers = os.cpu_count()
    commit = git_commit()
    reuse = args.reuse
    store = validation.PartialStore(args.store, commit, cfg["meta.decisions_md_sha256"])
    reuse_meta = {str(d): json.loads((d / "meta.json").read_text(encoding="utf-8")) for d in reuse}
    for directory, meta in reuse_meta.items():
        print(f"Reusing jobs from {directory} (commit {meta['git_commit'][:7]}), read-only", flush=True)
    sources = {}
    templates, targets = validation.real_templates(cfg)
    grid = validation.validation_sigma_grid(cfg)
    per_sigma = {}
    for s_idx, sigma in enumerate(grid):
        calibration = validation.load_saved_calibration(reuse, sigma) or store.load_calibration(sigma)
        if calibration is None:
            if args.saved_only:
                raise SystemExit(f"STOP (--saved-only): the sigma={sigma} calibration is not saved; it would need new simulation.")
            calibration = validation.calibrate(templates, targets, cfg, sigma, s_idx, workers)
            store.save_calibration(sigma, calibration)
        noise = {c["pulsar"]: c["noise_sd"] for c in calibration}
        records, failures = {}, {}
        get = lambda job, label=None: run_or_load(store, templates, job, cfg, sigma, s_idx, noise, workers, label,
                                                  reuse, sources, args.saved_only)
        for job in validation.validation_jobs(cfg):
            rec = validation.from_saved_records(get(job), cfg)
            if "d_se" not in rec and job.name in ("S1", "S3", "S1m"):        # records saved before D4l
                diag = get(replace(job, secondary=False), f"{job.name}-diagnostic")
                failures[f"{job.name}-diagnostic"] = validation.failure_check(diag, cfg)
                rec = merge_diagnostics(rec, diag)
            records[job.name] = rec
        nulls = {job.name: get(job) for job in validation.null_jobs(cfg, cfg["D4b.simulation.S3.slope"])}
        factors = validation.interval_factors(nulls, cfg)
        s3m_null = validation.ok(nulls["S3m-null"])
        s1m_raw = validation.ok(records["S1m"])
        s1m_vs_s3m = {"S1m_b_above_S3m_97.5": float((s1m_raw["b"] > np.percentile(s3m_null["b"], 97.5)).mean()),
                      "S1m_d_below_S3m_2.5": float((s1m_raw["d"] < np.percentile(s3m_null["d"], 2.5)).mean())}
        records, thr, refs, powers = validation.classified(records, nulls, cfg)
        failures.update({name: validation.failure_check(rec, cfg) for name, rec in {**records, **nulls}.items()})
        t1_power, t2_power, t2_power_d4m = powers["T1"], powers["T2"], powers["T2_signed_root_D4m_report_only"]
        per_sigma[str(sigma)] = {"calibration": calibration,
                                 "calibration_outside_tolerance": [c["pulsar"] for c in calibration if not c["within_tolerance"]],
                                 "directional_thresholds": thr,
                                 "null_mean_fitted_sigma_hat": {n: validation.sigma_hat_mean(r) for n, r in nulls.items()},
                                 "interval_factors": factors, "d4m_report_only": s1m_vs_s3m,
                                 "references": decision.references_for_json(refs),
                                 "powers_from_S1m": {"T1": t1_power, "T2": t2_power,
                                                     "T2_signed_root_D4m_report_only": t2_power_d4m,
                                                     "generating_sigma": sigma,
                                                     "S1m_mean_fitted_sigma_hat": validation.sigma_hat_mean(records["S1m"])},
                                 "failures_all_jobs": failures,
                                 "scenarios": {name: summarise(rec, refs, cfg, t2_power, t2_power_d4m, sigma)
                                               for name, rec in records.items()}}
    checks = criteria(per_sigma, cfg)
    record = {"description": "Validation on synthetic data with the real pulsar structure (D14/D4h/D4k criteria, "
                             "D4j outcomes with the D4l T1, the D4n T2 and the D4l bounds)",
              "criteria": checks, "all_pass": all(c["result"] == "PASS" for c in checks), "sigma_grid": grid,
              "per_test_alpha": cfg["D4.per_test_alpha"], "seed": cfg["D4b.simulation.seed"],
              "replicate_time_limit_s": validation.REPLICATE_TIME_LIMIT_S,
              "pulsars": [tp.name for tp in templates], "config_sha256": cfg["meta.decisions_md_sha256"],
              "git_commit": commit, "reused_job_sources": reuse_meta, "job_source_by_sigma": sources,
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "per_sigma": per_sigma}
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    for c in checks:
        if c["result"] == "FAIL" or "failed-replicate" not in c["name"]:
            print(f"{c['result']}  {c['name']}: {c['value']:.4f} (bound {c['bound']:.4f})")
    print(f"Failed-replicate checks: {sum(c['result'] == 'PASS' for c in checks if 'failed-replicate' in c['name'])} PASS "
          f"of {sum('failed-replicate' in c['name'] for c in checks)}")
    print("D4n item 4: generating sigma -> mean fitted sigma-hat (M2); S1m powers; A2 rates (new T2 / D4m T2)")
    for sigma, entry in per_sigma.items():
        sc, pw = entry["scenarios"], entry["powers_from_S1m"]
        print(f"  sigma {sigma}: S1m sigma-hat {pw['S1m_mean_fitted_sigma_hat']:.3f}; T1 power {pw['T1']:.3f}, "
              f"T2 power {pw['T2']:.3f} (D4m T2 {pw['T2_signed_root_D4m_report_only']:.3f}); "
              + ", ".join(f"{n} A2 {v['outcomes']['A2']:.3f}/{v['outcomes_with_D4m_T2_report_only']['A2']:.3f}"
                          for n, v in sc.items() if "outcomes" in v))
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}; all criteria pass: {record['all_pass']}")


if __name__ == "__main__":
    main()
