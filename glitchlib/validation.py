"""Power and validation runs (scripts 05, 06): simulated replicates with the real pulsar structure, analysed
by exactly the functions of the real run (glitchlib.analysis).

Templates come from glitchlib.templates (the D13-included pulsars). Each replicate keeps the intervals
started by glitches marked admissible, applies the D4e inclusion rule, and is analysed. Seeds: replicate
r of job j at sigma index s uses SeedSequence((D4b seed, j, s, r)), so results do not depend on how the
work is split across processes.

Run control (not research settings): each replicate has a REPLICATE_TIME_LIMIT_S wall-clock limit and is
recorded as a failure (reason "timeout") if it exceeds it; failed fits and any other error in a replicate
are recorded with their reason, and never stop the run.
Progress is printed after every chunk of replicates.
"""
import json
import signal
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from glitchlib import analysis, diagnostics, models, simulate
from glitchlib.intervals import apply_replicate_inclusion, following_waits
from glitchlib.io import read_meta
from glitchlib.templates import build_templates

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPLICATE_TIME_LIMIT_S = 120


class ReplicateTimeout(Exception):
    pass


def _on_alarm(signum, frame):
    raise ReplicateTimeout


@dataclass
class Job:
    job_id: int
    name: str
    scenario: simulate.Scenario
    method_rule: str | None        # emulated jobs: "resampled" or one method; statistics jobs: None
    noise: bool                    # statistics jobs: measurement noise on (False only for the D14 check)
    short_window_days: float | None
    n_replicates: int
    secondary: bool                # run the full secondary analyses (06) or only what 05 needs


def real_templates(cfg, root=PROJECT_ROOT):
    """Templates for the D13-included pulsars from the processed glitch table (no waiting times)."""
    glitches = pd.read_csv(root / "data" / "processed" / "glitches_clean.csv")
    obs = read_meta("obs_end")
    params = pd.read_csv(root / "data" / "processed" / "pulsar_params.csv")
    pulsars = sorted(glitches.loc[glitches["included_pulsar"], "pulsar"].unique())
    return build_templates(glitches, dict(zip(obs["name"], obs["obs_end_mjd"].astype(float))),
                           read_meta("monitoring_gaps"), dict(zip(params["name"], params["tau_sd_days"])),
                           read_meta("source_audit"), cfg, pulsars, dict(zip(params["name"], params["nu_hz"])))


def _calibrate_one(args):
    template, scenario, sigma, settings, target, n_sequences, tolerance, entropy = args
    return simulate.calibrate_noise(template, scenario, sigma, "resampled", settings, target, n_sequences, tolerance,
                                    np.random.default_rng(np.random.SeedSequence(entropy)))


def calibrate(templates, targets, cfg, sigma, sigma_index, workers):
    """D4e: per pulsar, the ν̇ noise under S1m on the pulsar's own method mix, at this sigma."""
    scenario = simulate.scenario_from_config(cfg, cfg["D4e.calibration.scenario"])
    settings = simulate.emulation_from_config(cfg, cfg["D4c.mixtures.unknown_emulated_as"])
    jobs = [(tp, scenario, sigma, settings, targets[tp.name], cfg["D4d.answers.calibration_sequences_per_step"],
             cfg["D4d.answers.calibration_tolerance_fraction_of_target"], (cfg["D4b.simulation.seed"], 999, sigma_index, k))
            for k, tp in enumerate(templates)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(_calibrate_one, jobs))


def _replicate(args):
    """One replicate with the wall-clock limit: a timeout or failed fit is recorded, never raised."""
    rep = args[4]
    started = time.monotonic()
    previous = signal.signal(signal.SIGALRM, _on_alarm)
    signal.alarm(REPLICATE_TIME_LIMIT_S)
    try:
        record = _replicate_body(args)
    except ReplicateTimeout:
        record = {"replicate": rep, "pulsars_dropped": 0, "discarded": False, "failed": True,
                  "failure": f"timeout after {REPLICATE_TIME_LIMIT_S} s"}
    except Exception as err:            # any other error is recorded for this replicate, never fatal to the run
        record = {"replicate": rep, "pulsars_dropped": 0, "discarded": False, "failed": True,
                  "failure": f"error: {type(err).__name__}: {err}"}
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    record["seconds"] = round(time.monotonic() - started, 3)
    return record


def _replicate_body(args):
    templates, job, sigma, sigma_index, rep, cfg, noise_sd = args
    rng = np.random.default_rng(np.random.SeedSequence((cfg["D4b.simulation.seed"], job.job_id, sigma_index, rep)))
    if job.scenario.emulated:
        settings = replace(simulate.emulation_from_config(cfg, cfg["D4c.mixtures.unknown_emulated_as"]),
                           short_window_days=job.short_window_days)
        intervals, glitches = simulate.simulate_emulated_dataset(templates, job.scenario, sigma, job.method_rule,
                                                                 settings, noise_sd, rng)
    else:
        intervals, glitches = simulate.simulate_statistics_dataset(templates, job.scenario, sigma, job.noise, rng)
    analysed = intervals[intervals["admissible"].astype(bool)] if len(intervals) else intervals
    kept, dropped, discarded = apply_replicate_inclusion(
        analysed, cfg["D13.min_admissible_completed_intervals_per_pulsar"],
        cfg["D4e.replicate_inclusion.discard_replicate_if_fewer_pulsars_than"], [tp.name for tp in templates])
    record = {"replicate": rep, "pulsars_dropped": len(dropped), "discarded": bool(discarded), "failed": False,
              "failure": None}
    if discarded:
        return record
    try:
        prim = analysis.primary(kept, cfg)
        tau = {tp.name: tp.tau_sd_days for tp in templates}
        sim = "D4b.simulation"
        perm = cfg[f"{sim}.validation_permutation.n_permutations"] if (
            job.secondary and rep < cfg[f"{sim}.validation_permutation.first_n_replicates"]) else None
        loo = job.secondary and rep < cfg[f"{sim}.validation_loo_cv_first_n_replicates"]
        sec = analysis.secondary(kept, cfg, tau, perm, int(rng.integers(2**31)), loo)
    except models.FitError as err:
        record["failed"] = True
        record["failure"] = f"fit: {err}"
        return record
    record.update({
        "n_intervals": prim["n_intervals"], "n_completed": prim["n_completed"], "n_pulsars": prim["n_pulsars"],
        "gate_met": prim["gate"]["met"], "sigma_hat": prim["sigma"],
        "b": prim["T1"]["estimate"], "b_se": prim["T1"]["se"], "T1_p": prim["T1"]["p_value"],
        "T1_signed": prim["T1"]["signed_root"], "T1_rejects_chi2": prim["T1"]["rejects_chi2"],
        "d": prim["T2"]["estimate"], "d_se": prim["T2"]["se"], "T2_p": prim["T2"]["p_value"],
        "T2_signed": prim["T2"]["signed_root"], "T2_rejects_chi2": prim["T2"]["rejects_chi2"],
        "a_M4": prim["fits"]["M4"]["estimates"]["a"],
        "aic_M2_below_M3": prim["aic"]["M2"] < prim["aic"]["M3"], "aic_M4_below_M3": prim["aic"]["M4"] < prim["aic"]["M3"],
        "placebo_b": sec["placebo"]["b"], "implied_beta_all_positive": sec["implied_beta"]["all_positive"],
        "spearman_forward": sec["spearman_forward"]["pooled_r"], "spearman_placebo": sec["spearman_placebo"]["pooled_r"]})
    if "permutation" in sec:
        record.update({"perm_p": sec["permutation"]["p_value"], "placebo_perm_p": sec["placebo_permutation"]["p_value"],
                       "placebo_significantly_negative": sec["placebo"]["significantly_negative"]})
    if "loo" in sec:
        record.update({f"loo_{m}": v for m, v in sec["loo"].items()})
    try:                                     # D4n item 3, report only: a failed refit is recorded, not a failure
        record["d_without_short_window"] = analysis.d_without_short_window(kept, cfg)
    except models.FitError as err:
        record["d_without_short_window"] = None
        record["d_without_short_window_failure"] = str(err)
    if job.scenario.emulated:
        record["missingness_log_or"] = missingness_log_or(glitches, templates, set(kept["pulsar"]), cfg)
        record.update(measurement_slopes(glitches))
    return record


def measurement_slopes(glitches):
    """D4l item 4 (report only): per method class, the least-squares slope of measured on true
    ln(Δν̇/ν̇) over the usable glitches of one emulated catalogue."""
    g = glitches[glitches["usable"].astype(bool)]
    out = {}
    for cls, h in g.groupby("method_class"):
        x, y = np.log(h["dnudot_frac"].to_numpy(float)), np.log(h["df1"].to_numpy(float) * 1e-3)
        if len(h) >= 3 and np.ptp(x) > 0:
            out[f"meas_slope_{cls}"] = float(np.polyfit(x, y, 1)[0])
    return out


def missingness_log_or(glitches, templates, pulsars, cfg):
    """D15/D4b missingness on one simulated catalogue: conditional logit of "Δν̇ not measured" on ln of the
    completed following wait, over glitches of the analysed pulsars (None if it cannot be estimated)."""
    g = glitches[glitches["pulsar"].isin(pulsars)].reset_index(drop=True)
    waits = following_waits(g, {tp.name: tp.obs_end_mjd for tp in templates},
                            {tp.name: [start for start, _ in tp.gaps] for tp in templates if tp.gaps})
    data = g.join(waits)
    data = data[data["following_t"].notna()].assign(ln_wait=lambda x: np.log(x["following_t"]),
                                                    df1_missing=lambda x: x["df1"].isna())
    try:
        return diagnostics.missingness_check(data, "df1_missing", "ln_wait", "following_completed",
                                             cfg["D4b.interval_estimates.level"])["log_odds_ratio"]
    except (ValueError, models.FitError):
        return None


def _chunk(args_list):
    return [_replicate(args) for args in args_list]


def run_job(templates, job, cfg, sigma, sigma_index, noise_sd, workers, replicates=None):
    """All replicates of one job at one sigma (or only the listed replicate numbers, regenerated from their
    own seeds), in parallel chunks; prints progress after each chunk and returns a DataFrame of
    per-replicate records, ordered by replicate."""
    reps = range(job.n_replicates) if replicates is None else sorted(int(r) for r in replicates)
    args = [(templates, job, sigma, sigma_index, rep, cfg, noise_sd) for rep in reps]
    if not args:
        return pd.DataFrame(columns=["replicate", "pulsars_dropped", "discarded", "failed", "failure"])
    size = max(1, len(args) // (8 * workers))
    chunks = [args[i:i + size] for i in range(0, len(args), size)]
    started, done, failed, records = time.monotonic(), 0, 0, []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(_chunk, c) for c in chunks]):
            chunk = future.result()
            records.extend(chunk)
            done += len(chunk)
            failed += sum(r["failed"] for r in chunk)
            print(f"  [{job.name} sigma={sigma}] {done}/{len(args)} replicates, {failed} failed, "
                  f"{time.monotonic() - started:.0f} s", flush=True)
    return pd.DataFrame(records).sort_values("replicate").reset_index(drop=True)


def ok(records):
    """Replicates that were analysed (not discarded by D4e, no failed fit)."""
    return records[~records["discarded"] & ~records["failed"]] if len(records) else records


def rate(records, column):
    good = ok(records)
    return float(good[column].mean()) if len(good) else None


def false_positive_bound(alpha, multiplier, n):
    """D4h: a false-positive rate passes if it is at most alpha + multiplier * sqrt(alpha (1 - alpha) / n)."""
    return float(alpha + multiplier * np.sqrt(alpha * (1 - alpha) / n))


def percentile_interval(values, level):
    """D4b: percentile interval for a simulated quantity."""
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    lo, hi = 100 * (1 - level) / 2, 100 * (1 + level) / 2
    return [float(np.percentile(values, lo)), float(np.percentile(values, hi))] if values.size else [None, None]


# ---- D4l ---------------------------------------------------------------------------------------------
def null_jobs(cfg, s3_slope):
    """D4l/D4m: the null references of the directional tests: S2 for T1; S3 and S3m (measured) for T2, with
    size coefficient s3_slope (D4b's in 06, the real M3 estimate in 08). Job ids 41-43, so their seeds are
    disjoint from every evaluation job's."""
    n = cfg["D4l.directional_tests.null_replicates"]
    s2 = simulate.scenario_from_config(cfg, cfg["D4l.directional_tests.T1_null_scenario"])
    s3 = replace(simulate.scenario_from_config(cfg, cfg["D4l.directional_tests.T2_null_scenario"]), slope=s3_slope)
    s3m = replace(simulate.scenario_from_config(cfg, "S3m"), slope=s3_slope)
    return [Job(41, "S2-null", s2, None, True, None, n, False), Job(42, "S3-null", s3, None, True, None, n, False),
            Job(43, "S3m-null", s3m, "resampled", True, None, n, False)]


def from_saved_records(records, cfg):
    """Old saved replicates (before D4l): the signed roots recovered exactly from the stored chi-square
    p-values (LR = chi2.isf(p, 1)) and the signs of b-hat and d-hat; the old decisions renamed *_chi2."""
    from scipy import stats
    out = records.rename(columns={"T1_rejects": "T1_rejects_chi2", "T2_rejects": "T2_rejects_chi2"}).copy()
    for test, est in (("T1", "b"), ("T2", "d")):
        if f"{test}_signed" in out:
            continue
        direction = cfg[f"D4l.directional_tests.{test}_predicted_direction"]
        lr = stats.chi2.isf(out[f"{test}_p"].astype(float), 1)
        out[f"{test}_signed"] = [analysis.signed_root(s, e, direction) if np.isfinite(s) else np.nan
                                 for s, e in zip(lr, out[est].astype(float))]
    return out


def thresholds(null_records, cfg):
    """The decision thresholds from the null references (each at the sigma it was run at):
      T1 (D4l): the 97.5th percentile of the signed root under S2-null;
      T2 (D4n): the cut on d-hat, the smaller of the 2.5th percentiles of d-hat under S3-null and S3m-null.
    Also, report only, the superseded D4m signed-root threshold for T2 (the larger 97.5th percentile under
    S3 and S3m), so the old and new T2 can be shown side by side; it never changes an outcome."""
    from glitchlib import decision
    cut, cut_per_null = decision.t2_d_cut({name: ok(null_records[f"{name}-null"])["d"]
                                           for name in cfg["D4n.T2.null_scenarios"]}, cfg)
    old_per_null = {name: decision.null_threshold(ok(null_records[f"{name}-null"])["T2_signed"], cfg)
                    for name in cfg["D4m.T2_rejects_if_exceeds_null_percentile_under"]}
    return {"T1": decision.null_threshold(ok(null_records["S2-null"])["T1_signed"], cfg),
            "T2_d_cut": cut, "T2_d_cut_per_null": cut_per_null,
            "T2_signed_root_D4m_report_only": max(old_per_null.values()),
            "T2_signed_root_D4m_per_null_report_only": old_per_null}


def apply_directional(records, thr):
    """T1 (D4l) rejects when its signed root exceeds the S2 threshold; T2 (D4n) when d-hat is below the
    cut. T2_rejects_D4m: the superseded signed-root T2, report only."""
    out = records.copy()
    out["T1_rejects"] = out["T1_signed"] > thr["T1"]
    out["T2_rejects"] = out["d"] < thr["T2_d_cut"]
    out["T2_rejects_D4m"] = out["T2_signed"] > thr["T2_signed_root_D4m_report_only"]
    return out


def classified(records, nulls, cfg):
    """One sigma, as in 05/06/08: the thresholds from the null references, the T1/T2 decisions on every
    scenario's records, the D4j references (S1, S1m, S4-mix) and the powers (S1m rejection rates; the D4m
    signed-root T2 power report only). Returns (records with decisions, thresholds, refs, powers)."""
    from glitchlib import decision
    thr = thresholds(nulls, cfg)
    records = {name: apply_directional(rec, thr) for name, rec in records.items()}
    s1m, s4, s1 = (ok(records[n]) for n in ("S1m", "S4-mix", "S1"))
    refs = decision.references(s1m["b"], s1m["placebo_b"], s4["b"], s1["b"], cfg)
    powers = {"T1": float(s1m["T1_rejects"].mean()), "T2": float(s1m["T2_rejects"].mean()),
              "T2_signed_root_D4m_report_only": float(s1m["T2_rejects_D4m"].mean())}
    return records, thr, refs, powers


def replicate_outcomes(good, refs, t2_power, cfg, t2_column="T2_rejects"):
    """decision.classify for every analysed replicate (rows of ok(records), in order), with the D4n
    short-window refit where the records have it."""
    from glitchlib import decision
    sw = good["d_without_short_window"] if "d_without_short_window" in good else pd.Series(None, index=good.index)
    return [decision.classify(r.T1_rejects, getattr(r, t2_column), r.b, r.placebo_b, refs, t2_power, cfg,
                              None if pd.isna(d_sw) else float(d_sw))
            for r, d_sw in zip(good.itertuples(), sw)]


def sigma_hat_mean(records):
    """D4n item 4: mean fitted sigma-hat of the primary model (D3, M2) over the analysed replicates."""
    good = ok(records)
    return float(good["sigma_hat"].mean()) if len(good) else None


def interval_factors(null_records, cfg):
    """D4m: (SD of the estimate / mean model SE) from the matching null reference at the same sigma
    (S2 for b, S3m for d)."""
    ci = "D4m.calibrated_intervals"
    out = {}
    for est, null in (("b", cfg[f"{ci}.b_null"]), ("d", cfg[f"{ci}.d_null"])):
        good = ok(null_records[f"{null}-null"])
        out[est] = float(good[est].std(ddof=1) / good[f"{est}_se"].mean())
    return out


calibrated_interval = analysis.calibrated_interval     # D4m; shared with 09 (glitchlib.analysis)


# ---- D4k -----------------------------------------------------------------------------------------
def validation_sigma_grid(cfg):
    """D14 sigma grid plus the D4k additions (06 only)."""
    return sorted(set(cfg["D14.sigma_grid"]) | set(cfg["D4k.validation_sigma_grid_additional"]))


def failure_check(records, cfg):
    """D4k: failed replicates are excluded and listed with their reason; more than the allowed fraction
    of failures makes the job FAIL until the cause is explained in the findings log."""
    failed = records[records["failed"]] if len(records) else records
    fraction = float(len(failed) / len(records)) if len(records) else 0.0
    limit = cfg["D4k.failed_replicates.max_fraction_per_job"]
    return {"failed": int(len(failed)), "fraction": fraction, "limit": limit,
            "result": "PASS" if fraction <= limit else "FAIL",
            "reasons": failed[["replicate", "failure"]].to_dict("records") if len(failed) else []}


def validation_jobs(cfg):
    """The evaluation jobs of 06 (D14 scenarios, S1m, S3m and the S4 variants) with their job ids, which
    fix their seeds."""
    rep = "D4b.simulation.replicates"
    sc = lambda name: simulate.scenario_from_config(cfg, name)
    s4 = sc("S4")
    out = [Job(1, "S1", sc("S1"), None, True, None, cfg[f"{rep}.S1"], True),
           Job(2, "S2", sc("S2"), None, True, None, cfg[f"{rep}.S2"], True),
           Job(5, "S3", sc("S3"), None, True, None, cfg[f"{rep}.S3"], True),
           Job(6, "S1-noise-off", sc("S1"), None, False, None, cfg[f"{rep}.S1"], True),
           Job(7, "S3m", sc("S3m"), "resampled", True, None, cfg[f"{rep}.S3"], True),     # D4m (reading: S3's count)
           Job(3, "S1m", sc("S1m"), "resampled", True, None, cfg["D4h.confirmed.S1m_replicates_per_sigma"], True),
           Job(4, "S4-mix", s4, "resampled", True, None, cfg[f"{rep}.S4_per_variant"], True)]
    for k, (method, name) in enumerate([("window", "S4w"), ("local-relaxation", "S4l"), ("global", "S4g"),
                                        ("global-noexp", "S4g0")]):
        out.append(Job(10 + k, name, s4, method, True, None, cfg[f"{rep}.S4_per_variant"], True))
    for k, source in enumerate(cfg["D4f.short_window.sources"].values()):
        out.append(Job(20 + k, f"S4sw-{source['source_type']}", s4, "short-window", True,
                       float(source["w_days"]), cfg[f"{rep}.S4_per_variant"], True))
    return out


def load_saved(directories, sigma, name):
    """A saved job from the first of the partial-store directories that has it: (records, directory), or
    (None, None)."""
    for directory in directories:
        path = Path(directory) / f"sigma_{sigma}__{name}.pkl"
        if path.exists():
            return pd.read_pickle(path), Path(directory)
    return None, None


def load_saved_calibration(directories, sigma):
    """A saved D4e calibration for this sigma from the first directory that has it, or None."""
    for directory in directories:
        path = Path(directory) / f"sigma_{sigma}__calibration.json"
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
    return None


def detection_counts(glitches):
    """D4o item 5 (report only): per exponential group, over the emulated glitches of one catalogue: all
    glitches, those with an exponential, those whose exponential was detected (fitted, D4d), and
    detections where no exponential was present."""
    out = {}
    for group, g in glitches.groupby("exp_group"):
        present, fitted = g["exp_present"].astype(bool), g["exp_fitted"].astype(bool)
        out[group] = {"glitches": int(len(g)), "with_exponential": int(present.sum()),
                      "detected": int((fitted & present).sum()), "detected_without_exponential": int((fitted & ~present).sum())}
    return out


def _emulation_counts(args):
    templates, job, sigma, sigma_index, rep, cfg, noise_sd, settings, label_spec = args
    rng = np.random.default_rng(np.random.SeedSequence((cfg["D4b.simulation.seed"], job.job_id, sigma_index, rep)))
    _, glitches = simulate.simulate_emulated_dataset(templates, job.scenario, sigma, job.method_rule, settings, noise_sd, rng)
    if label_spec is not None:              # count by another spec's groups (e.g. the placeholder by D4o group)
        glitches = glitches.copy()
        for pulsar, g in glitches.groupby("pulsar"):
            glitches.loc[g.index, "exp_group"] = simulate.exponential_ranges(g, pulsar, label_spec)[0]
    return detection_counts(glitches)


def emulated_detection_counts(templates, job, cfg, sigma, sigma_index, noise_sd, settings, workers, label_spec=None):
    """D4o item 5: detection counts summed over every replicate of an emulated job, from the same seeds as
    the job's replicates (the emulation only; nothing is fitted). settings: the emulation settings to use
    (the D4o ones, or the D14 placeholder for the report-only comparison); label_spec: count by this
    spec's groups instead of the drawing spec's."""
    args = [(templates, job, sigma, sigma_index, rep, cfg, noise_sd, settings, label_spec) for rep in range(job.n_replicates)]
    total = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for counts in pool.map(_emulation_counts, args, chunksize=max(1, len(args) // (8 * workers))):
            for group, c in counts.items():
                acc = total.setdefault(group, dict.fromkeys(c, 0))
                for k, v in c.items():
                    acc[k] += v
    return total


def floor_sigma_index(cfg):
    """Seed index of the report-only sigma-hat floor runs (1 Nov request): three past the validation grid,
    so their seeds differ from the grid, the real-sigma runs and the D4n search steps."""
    return real_sigma_index(cfg) + 2


def reference_jobs(cfg):
    """D4k/D4l: the S1, S1m and S4-mix reference runs for the real classification, with 06's job ids,
    replicate counts and seeding (job ids 1, 3 and 4)."""
    return [Job(1, "S1", simulate.scenario_from_config(cfg, "S1"), None, True, None,
                cfg["D4b.simulation.replicates.S1"], False),
            Job(3, "S1m", simulate.scenario_from_config(cfg, "S1m"), "resampled", True, None,
                cfg["D4h.confirmed.S1m_replicates_per_sigma"], False),
            Job(4, "S4-mix", simulate.scenario_from_config(cfg, "S4"), "resampled", True, None,
                cfg["D4b.simulation.replicates.S4_per_variant"], False)]


def real_sigma_index(cfg):
    """Seed index for runs at the real fitted sigma: one past the validation grid, so its seeds differ
    from every grid sigma's (seeding otherwise as in 06)."""
    return len(validation_sigma_grid(cfg))


def missingness_reference(records):
    """D4k: S1m missingness reference from the estimable replicates only, with the counts."""
    good = ok(records)
    values = good["missingness_log_or"].astype(float) if "missingness_log_or" in good else pd.Series(dtype=float)
    estimable = values.dropna()
    return {"replicates": int(len(good)), "estimable": int(len(estimable)),
            "not_estimable": int(len(good) - len(estimable)), "values": estimable.to_numpy()}


# ---- D4n: sigma matched on fitted scatter (08) ------------------------------------------------------
class SigmaMatchError(Exception):
    pass


def match_sigma(median_sigma_hat_at, target, cfg):
    """D4n item 1 with D4p item 1: the generating sigma on D4p's search range at which
    median_sigma_hat_at(sigma) (the median fitted sigma-hat of the primary model over D4n's replicates per
    step) is within D4n's tolerance of target, by bisection. If the fitted sigma-hat at the lower end
    already exceeds the target, returns the lower end (at_lower_end = True). If no step is within the
    tolerance after D4p's maximum number of bisection steps, returns the midpoint of the final bracket
    (converged = False). Returns {sigma, at_lower_end, converged, steps: [[sigma, median sigma-hat], ...]}."""
    m, p = "D4n.matched_sigma", "D4p.matched_sigma_search"
    if cfg[f"{m}.search"] != "bisection":
        raise ValueError(f"D4n.matched_sigma.search = {cfg[f'{m}.search']} is not implemented.")
    lo, hi = (float(x) for x in cfg[f"{p}.search_range"])
    tol = cfg[f"{m}.tolerance_in_sigma_hat"]
    steps = []

    def evaluate(sigma):
        value = float(median_sigma_hat_at(sigma))
        steps.append([sigma, value])
        return value

    at_lo = evaluate(lo)
    if at_lo > target + tol:
        return {"sigma": lo, "at_lower_end": True, "converged": True, "steps": steps}
    if abs(at_lo - target) <= tol:
        return {"sigma": lo, "at_lower_end": False, "converged": True, "steps": steps}
    for _ in range(cfg[f"{p}.max_bisection_steps"]):
        mid = (lo + hi) / 2
        value = evaluate(mid)
        if abs(value - target) <= tol:
            return {"sigma": mid, "at_lower_end": False, "converged": True, "steps": steps}
        lo, hi = (mid, hi) if value < target else (lo, mid)
    return {"sigma": (lo + hi) / 2, "at_lower_end": False, "converged": False, "steps": steps}


def matched_sigma_index(cfg):
    """Seed index of the D4n search steps: two past the validation grid, so their seeds differ from every
    grid sigma's and from the full runs at the matched sigma (real_sigma_index). The same replicate seeds
    are used at every step."""
    return real_sigma_index(cfg) + 1


def match_job_sigma(templates, targets, job, target, cfg, workers):
    """D4n item 1 for one reference job: bisection on the generating sigma, each step running D4n's
    replicates per step (emulated jobs: noise calibrated at that sigma, D4e)."""
    index = matched_sigma_index(cfg)
    search_job = replace(job, n_replicates=cfg["D4n.matched_sigma.replicates_per_step"], secondary=False)

    def median_at(sigma):
        noise = None
        if job.scenario.emulated:
            noise = {c["pulsar"]: c["noise_sd"] for c in calibrate(templates, targets, cfg, sigma, index, workers)}
        good = ok(run_job(templates, search_job, cfg, sigma, index, noise, workers))
        if not len(good):
            raise SigmaMatchError(f"{job.name}: no analysed replicate at sigma {sigma}")
        return float(good["sigma_hat"].median())

    return match_sigma(median_at, target, cfg)


# ---- partial saves and resume (06) -----------------------------------------------------------------
class PartialStore:
    """Finished jobs saved under results/validation_partial/, so an interrupted 06 resumes. Resuming is
    refused if the saved commit or config checksum differs from the current ones."""

    def __init__(self, directory, commit, config_sha256):
        self.dir = Path(directory)
        self.meta = {"git_commit": commit, "config_sha256": config_sha256}
        meta_path = self.dir / "meta.json"
        if meta_path.exists():
            saved = json.loads(meta_path.read_text(encoding="utf-8"))
            if saved != self.meta:
                raise SystemExit(f"STOP: {self.dir} holds partial results from {saved}, but the current run is "
                                 f"{self.meta}. Resume refused; remove {self.dir} to start from scratch.")
        elif self.dir.exists() and any(self.dir.iterdir()):
            raise SystemExit(f"STOP: {self.dir} holds files but no meta.json; remove it to start from scratch.")
        if commit is None:
            raise SystemExit("STOP: no git commit found; partial results could not be tied to the code.")
        self.dir.mkdir(parents=True, exist_ok=True)
        meta_path.write_text(json.dumps(self.meta, indent=2) + "\n", encoding="utf-8")

    def _path(self, sigma, name, kind):
        return self.dir / f"sigma_{sigma}__{name}.{kind}"

    def load_job(self, sigma, name):
        path = self._path(sigma, name, "pkl")
        return pd.read_pickle(path) if path.exists() else None

    def save_job(self, sigma, name, records):
        tmp = self._path(sigma, name, "pkl.tmp")
        records.to_pickle(tmp)
        tmp.replace(self._path(sigma, name, "pkl"))           # written whole or not at all

    def load_calibration(self, sigma):
        path = self._path(sigma, "calibration", "json")
        return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None

    def save_calibration(self, sigma, calibration):
        self._path(sigma, "calibration", "json").write_text(json.dumps(calibration, indent=2, default=float) + "\n",
                                                            encoding="utf-8")
