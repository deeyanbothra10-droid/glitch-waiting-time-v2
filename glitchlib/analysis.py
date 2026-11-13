"""The analyses of the real run (scripts 07, 08, 09), imported unchanged by the validation (05, 06).

Every function takes analysis-ready intervals (from final_dataset.csv, or a simulated replicate built by
the same interval code) and the validated config; nothing here decides what data are analysed.
"""
import numpy as np
import pandas as pd
from scipy import stats

from glitchlib import diagnostics, models

PRIMARY_MODELS = ("M0", "M2", "M3", "M4")


def _with_ln_t(intervals):
    return intervals.assign(ln_t=np.log(intervals["t"].astype(float)))


def gate(intervals, cfg):
    """D12 gate: admissible completed intervals and pulsars, against the D12 minimums."""
    completed = intervals[intervals["event"].astype(bool)]
    n_completed, n_pulsars = int(len(completed)), int(completed["pulsar"].nunique())
    return {"admissible_completed": n_completed, "pulsars": n_pulsars,
            "met": bool(n_completed >= cfg["D12.gate.min_admissible_completed_intervals"]
                        and n_pulsars >= cfg["D12.gate.min_pulsars"])}


DIRECTION = {"b_positive": 1.0, "d_negative": -1.0}


def signed_root(statistic, estimate, direction):
    """D4l: sqrt of the LR statistic, positive when the estimate lies in the predicted direction."""
    return float(np.sign(estimate) * DIRECTION[direction] * np.sqrt(max(statistic, 0.0)))


def primary(intervals, cfg):
    """D3 models and the likelihood-ratio statistics of T1 (M2 vs M0) and T2 (M4 vs M3), with D4b Wald
    intervals. D4l: the decision uses the signed root of each LR statistic against its simulated null
    reference (glitchlib.decision), so here each test reports its signed root; the chi-square decision at
    the D4 per-test alpha is kept as rejects_chi2 for reference only."""
    me = cfg["D6.primary_errors_in_lnR"]
    level = cfg["D4b.interval_estimates.level"]
    alpha = cfg["D4.per_test_alpha"]
    fits = {m: models.fit(m, intervals, me) for m in PRIMARY_MODELS}
    out = {"n_intervals": int(len(intervals)), "n_completed": int(intervals["event"].astype(bool).sum()),
           "n_pulsars": int(intervals["pulsar"].nunique()), "gate": gate(intervals, cfg),
           "fits": {m: f.summary() for m, f in fits.items()}, "aic": {m: float(f.aic) for m, f in fits.items()},
           "per_test_alpha": alpha, "interval_level": level}
    for test, param in (("T1", cfg["D2.T1.parameter"]), ("T2", cfg["D2.T2.parameter"])):
        comparison = cfg[f"D4.{test}_comparison"]
        full = fits[comparison["full"]]
        lr = models.lr_test(full, fits[comparison["reduced"]])
        estimate = float(full.estimates[param])
        out[test] = {**lr, "parameter": param, "estimate": estimate,
                     "se": float(full.se[param]), "interval": list(full.intervals(level)[param]),
                     "signed_root": signed_root(lr["statistic"], estimate,
                                                cfg[f"D4l.directional_tests.{test}_predicted_direction"]),
                     "rejects_chi2": bool(lr["p_value"] < alpha)}
    out["sigma"] = float(fits[cfg["D3.T1_model"]].estimates["sigma"])
    return out


def d_without_short_window(intervals, cfg):
    """D4n item 3 (report only): d-hat of the T2 model refitted without the intervals started by
    short-window glitches (method class of D4f), with the D13 inclusion count applied to that subset as
    for every robustness subset (09). None if no pulsar is left. A FitError is not caught here (D4k: on
    real data it stops the run)."""
    from glitchlib.intervals import apply_inclusion
    subset = intervals[intervals["method_class"] != cfg["D4f.short_window.class"]]
    kept, _ = apply_inclusion(subset, cfg["D13.min_admissible_completed_intervals_per_pulsar"])
    if kept.empty:
        return None
    fit = models.fit(cfg["D3.T2_model"], kept, cfg["D6.primary_errors_in_lnR"])
    return float(fit.estimates[cfg["D2.T2.parameter"]])


def implied_beta_check(intervals, cfg, tau_sd_days):
    """M1 fit and implied beta per pulsar with its Wald interval (D4i: reported, decides nothing), and,
    for reference only, whether beta > 0 for every pulsar with at least 5 admissible intervals (the old
    D12 A condition)."""
    fit = models.fit("M1", intervals, cfg["D6.primary_errors_in_lnR"])
    table = diagnostics.implied_beta(fit, tau_sd_days)
    z = stats.norm.ppf(0.5 + cfg["D4i.implied_beta.interval_level"] / 2)
    table["beta_low"], table["beta_high"] = table["beta"] - z * table["se_beta"], table["beta"] + z * table["se_beta"]
    counts = intervals.groupby("pulsar").size()
    a = "D12.outcomes.A_support"
    table["n_intervals"] = table["pulsar"].map(counts)
    checked = table[table["n_intervals"] >= cfg[f"{a}.implied_beta_min_admissible_intervals"]]
    return {"table": table, "all_positive": bool((checked["beta"] > cfg[f"{a}.implied_beta_min"]).all()),
            "n_checked": int(len(checked))}


def secondary(intervals, cfg, tau_sd_days, n_permutations, seed, loo):
    """D4 secondary analyses and D15 placebo: implied beta, Spearman (completed intervals only), the
    placebo fit, the D4a permutation tests (primary and placebo; skipped if n_permutations is None),
    LOO-CV for M0-M4 (if loo), and the D5 per-pulsar Poisson baseline."""
    me = cfg["D6.primary_errors_in_lnR"]
    level = cfg["D4b.interval_estimates.level"]
    sp = "D4a.spearman"
    data = _with_ln_t(intervals)
    placebo = diagnostics.placebo_data(data)
    placebo = placebo.assign(ln_t=np.log(placebo["t"].astype(float)))
    placebo_fit = models.fit("M2", placebo, me)
    out = {"implied_beta": implied_beta_check(intervals, cfg, tau_sd_days),
           "spearman_forward": diagnostics.spearman_within(data, "ln_R", "ln_t", cfg[f"{sp}.per_pulsar_min_completed_admissible"],
                                                           cfg[f"{sp}.label"], level),
           "spearman_placebo": diagnostics.spearman_within(placebo, "ln_R", "ln_t",
                                                           cfg[f"{sp}.per_pulsar_min_completed_admissible"],
                                                           cfg[f"{sp}.label"], level),
           "placebo": {"b": float(placebo_fit.estimates["b"]), "se": float(placebo_fit.se["b"]),
                       "interval": list(placebo_fit.intervals(level)["b"]), "n": int(len(placebo))},
           "poisson": models.fit_poisson_per_pulsar(intervals)}
    if n_permutations is not None:
        intervals_rule = cfg["D4a.permutation.intervals"]
        out["permutation"] = diagnostics.permutation_test(intervals, intervals_rule, n_permutations, seed, me)
        out["placebo_permutation"] = diagnostics.permutation_test(placebo, intervals_rule, n_permutations, seed, me)
        alpha = cfg["D4b.diagnostics.placebo.alpha"]
        p, b = out["placebo_permutation"]["p_value"], out["placebo_permutation"]["b_obs"]
        out["placebo"]["significantly_negative"] = bool(p < alpha and b < 0)      # D4b; reported (D4d)
        out["placebo"]["significantly_positive"] = bool(p < alpha and b > 0)
    if loo:
        out["loo"] = {m: models.loo_cv(m, intervals, me)["mean_log_score"] for m in ("M0", "M1", "M2", "M3", "M4")}
    return out


def kaplan_meier(t, event):
    """Kaplan-Meier survival estimate: (event times, survival just after each)."""
    t, event = np.asarray(t, float), np.asarray(event, bool)
    times = np.unique(t[event])
    survival, s = [], 1.0
    for time in times:
        at_risk = np.sum(t >= time)
        s *= 1.0 - np.sum((t == time) & event) / at_risk
        survival.append(s)
    return times, np.array(survival)


def summarise_secondary(sec):
    """Flatten the secondary results into JSON-ready values."""
    out = {"implied_beta": sec["implied_beta"]["table"].to_dict("records"),
           "implied_beta_all_positive": sec["implied_beta"]["all_positive"],
           "placebo": sec["placebo"], "poisson": sec["poisson"].to_dict("records")}
    for key in ("spearman_forward", "spearman_placebo"):
        s = dict(sec[key])
        s["per_pulsar"] = s["per_pulsar"].to_dict("records")
        out[key] = s
    for key in ("permutation", "placebo_permutation"):
        if key in sec:
            out[key] = {k: v for k, v in sec[key].items() if k != "b_perm"}
    if "loo" in sec:
        out["loo"] = sec["loo"]
    return out




def primary_data(final):
    """The primary analysis data: admissible intervals of D13-included pulsars (final_dataset.csv rows)."""
    return final[final["admissible"].astype(bool) & final["included_pulsar"].astype(bool)].reset_index(drop=True)


def calibrated_interval(estimate, se, factor, level):
    """D4m: estimate +- z * SE * factor (the plain Wald interval when factor = 1)."""
    half = stats.norm.ppf(0.5 + level / 2) * se * factor
    return [float(estimate - half), float(estimate + half)]


def robustness_estimates(result, cfg, factors, primary_estimates):
    """D4p item 3: for b-hat and d-hat of one robustness case, the estimate, its calibrated interval
    (model SE times the D4m factor; D4p clarification: the primary's factors from 08, labelled
    "approximate (primary factors)"), the plain Wald interval (uncalibrated), the signed root, and whether
    it has the primary's sign. No outcome classification. factors and primary_estimates: {"b": ., "d": .}."""
    level = cfg["D4b.interval_estimates.level"]
    out = {}
    for test in ("T1", "T2"):
        r = result[test]
        param = r["parameter"]
        out[param] = {"estimate": r["estimate"],
                      "calibrated_interval": calibrated_interval(r["estimate"], r["se"], factors[param], level),
                      "calibrated_interval_label": cfg["D4p.robustness_09.interval_label"],
                      cfg["D4m.calibrated_intervals.plain_wald_label"]: r["interval"], "signed_root": r["signed_root"],
                      "same_sign_as_primary": bool(np.sign(r["estimate"]) == np.sign(primary_estimates[param]))}
    return out


def robustness_case(data, cfg, label, exploratory_below, factors, primary_estimates):
    """One robustness case (09): the D13 inclusion count applied to the data being analysed (confirmed
    reading), then the primary analysis, with D4p's report for b-hat and d-hat (robustness_estimates).
    No data left: reported as such. A fit failure on real data stops the run (D4k): FitError is not
    caught here."""
    from glitchlib.intervals import apply_inclusion
    kept, dropped = apply_inclusion(data, cfg["D13.min_admissible_completed_intervals_per_pulsar"])
    if kept.empty:
        return {"case": label, "status": "no data: no pulsar has enough completed intervals in this subset",
                "pulsars_dropped": dropped}
    result = primary(kept, cfg)
    out = {"case": label, "status": "run", "n_intervals": result["n_intervals"], "n_completed": result["n_completed"],
           "n_pulsars": result["n_pulsars"], "pulsars_dropped": dropped, "T1": result["T1"], "T2": result["T2"],
           "sigma": result["sigma"], "estimates_D4p": robustness_estimates(result, cfg, factors, primary_estimates)}
    if exploratory_below is not None:
        out["exploratory"] = bool(result["n_intervals"] < exploratory_below)
    return out
