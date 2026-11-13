"""EXPLORATORY, POST HOC, NOT PRE-REGISTERED analyses (requested 1 Nov, after the frozen run 07-09).

Nothing here is part of the pre-registration or changes any committed result; scripts using it write only
to results/exploratory/. Fits reuse glitchlib.analysis.robustness_case (D13 count applied to the subset,
the primary models, b-hat and d-hat with calibrated intervals labelled "approximate (primary factors)").
"""
from dataclasses import replace

import numpy as np
from scipy import stats

from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from glitchlib import models, simulate
from glitchlib.analysis import robustness_case
from glitchlib.intervals import apply_inclusion

LABEL = "EXPLORATORY, POST HOC, NOT PRE-REGISTERED"


def short_window_shorter_than_span(intervals, cfg):
    """Completed intervals started by a short-window glitch whose T is shorter than its source's
    post-glitch span w (D4f.short_window.sources: w_days). Censored intervals are never flagged (T unknown)."""
    spans = {key: float(src["w_days"]) for key, src in cfg["D4f.short_window.sources"].items()}
    w = intervals["source_key"].map(spans)
    return ((intervals["method_class"] == cfg["D4f.short_window.class"]) & intervals["event"].astype(bool)
            & w.notna() & (intervals["t"] < w))


def per_pulsar_cases(data, cfg, min_intervals, factors, primary_estimates):
    """For every pulsar with at least min_intervals intervals in data: its method mix and its own fit
    (robustness_case on that pulsar alone). A fit failure is recorded for that pulsar, not raised."""
    counts = data.groupby("pulsar").size()
    out = []
    for pulsar in counts[counts >= min_intervals].index:
        d = data[data["pulsar"] == pulsar]
        row = {"pulsar": pulsar, "n_intervals": int(len(d)), "n_completed": int(d["event"].astype(bool).sum()),
               "method_mix": {k: int(v) for k, v in d["method_class"].value_counts().sort_index().items()}}
        try:
            row["fit"] = robustness_case(d, cfg, f"{pulsar} alone", None, factors, primary_estimates)
        except models.FitError as err:
            row["fit"] = {"status": f"fit failed: {err}"}
        out.append(row)
    return out


def _t2_fit(data, cfg):
    """d-hat of the T2 model (M4) on data after the D13 count (as for every subset): (estimate, pulsars dropped)."""
    kept, dropped = apply_inclusion(data, cfg["D13.min_admissible_completed_intervals_per_pulsar"])
    fit = models.fit(cfg["D3.T2_model"], kept, cfg["D6.primary_errors_in_lnR"])
    return float(fit.estimates[cfg["D2.T2.parameter"]]), dropped


def leave_one_out_d(data, cfg):
    """d-hat of M4 with each interval left out in turn (D13 count reapplied to what remains). Returns one row
    per left-out interval: its index in data, d-hat without it, the pulsars the D13 count dropped, or the
    fit failure."""
    rows = []
    for i in range(len(data)):
        try:
            d, dropped = _t2_fit(data.drop(index=data.index[i]), cfg)
            rows.append({"position": i, "d_without": d, "pulsars_dropped": dropped})
        except models.FitError as err:
            rows.append({"position": i, "d_without": None, "pulsars_dropped": [], "failure": str(err)})
    return rows


def _fit_per_pulsar(design, params):
    """Maximise the T2 model's likelihood with each coefficient in params replaced by one per pulsar.
    Returns (estimates, standard errors, loglik), both dicts keyed by parameter name."""
    if not params:
        theta, loglik, cov = models.maximise(design)
        names = design.names
    else:
        keep, extra, extra_names = [], [], []
        for j, name in enumerate(design.names[:-1]):
            if name in params:
                block = np.zeros((len(design.y), len(design.pulsars)))
                block[np.arange(len(design.y)), design.pulsar_index] = design.X[:, j]
                extra.append(block)
                extra_names += [f"{name}[{p}]" for p in design.pulsars]
            else:
                keep.append(j)
        X = np.hstack([design.X[:, keep]] + extra)
        names = [design.names[j] for j in keep] + extra_names + ["log_sigma"]
        theta, loglik, cov = models.maximise(replace(design, model=f"M4_separate_{'_'.join(params)}", X=X, names=names))
    se = np.sqrt(np.diag(cov))
    return dict(zip(names, map(float, theta))), dict(zip(names, map(float, se))), float(loglik)


def per_pulsar_lr_test(data, cfg, min_intervals, reduced_separate, full_separate):
    """LR test, on the pulsars with at least min_intervals intervals, of the T2 model (M4) with the
    coefficients in full_separate given one value per pulsar against the same model with only those in
    reduced_separate per pulsar (same intervals; df = extra parameters)."""
    counts = data.groupby("pulsar").size()
    sub = data[data["pulsar"].isin(counts[counts >= min_intervals].index)]
    design = models.build_design(cfg["D3.T2_model"], sub)
    est_r, se_r, ll_r = _fit_per_pulsar(design, tuple(reduced_separate))
    est_f, se_f, ll_f = _fit_per_pulsar(design, tuple(full_separate))
    statistic = float(2 * (ll_f - ll_r))
    df = len(est_f) - len(est_r)
    return {"pulsars": list(design.pulsars), "n_intervals": int(len(sub)), "reduced_separate": list(reduced_separate),
            "full_separate": list(full_separate), "loglik_reduced": ll_r, "loglik_full": ll_f,
            "statistic": statistic, "df": int(df), "p_value": float(stats.chi2.sf(statistic, df)),
            "reduced": {"estimates": est_r, "se": se_r}, "full": {"estimates": est_f, "se": se_f}}


def separate_d_test(data, cfg, min_intervals):
    """LR test, on the pulsars with at least min_intervals intervals, of M4 with a separate d per pulsar
    against M4 with a common d (same intervals; df = number of pulsars - 1)."""
    t = per_pulsar_lr_test(data, cfg, min_intervals, (), ("d",))
    return {"pulsars": t["pulsars"], "n_intervals": t["n_intervals"], "common_d": t["reduced"]["estimates"]["d"],
            "loglik_common": t["loglik_reduced"], "loglik_separate": t["loglik_full"], "statistic": t["statistic"],
            "df": t["df"], "p_value": t["p_value"],
            "separate_d": {p: {"estimate": t["full"]["estimates"][f"d[{p}]"], "se": t["full"]["se"][f"d[{p}]"]}
                           for p in t["pulsars"]}}


# ---- single-pulsar simulations ---------------------------------------------------------------------
def _single_pulsar_replicate(args):
    """One emulated replicate of one pulsar: admissible intervals, the D13 count, then M4 alone.
    Returns d-hat and M4's sigma-hat, or the reason it was discarded or failed."""
    template, scenario, sigma, method_rule, settings, noise_sd, entropy, cfg = args
    rng = np.random.default_rng(np.random.SeedSequence(entropy))
    intervals, _ = simulate.simulate_emulated_dataset([template], scenario, sigma, method_rule, settings,
                                                      {template.name: noise_sd}, rng)
    analysed = intervals[intervals["admissible"].astype(bool)] if len(intervals) else intervals
    kept, _ = apply_inclusion(analysed, cfg["D13.min_admissible_completed_intervals_per_pulsar"]) if len(analysed) else (analysed, [])
    if kept.empty:
        return {"rep": entropy[-1], "status": "discarded"}
    try:
        fit = models.fit(cfg["D3.T2_model"], kept, cfg["D6.primary_errors_in_lnR"])
    except models.FitError as err:
        return {"rep": entropy[-1], "status": f"failed: {err}"}
    return {"rep": entropy[-1], "status": "ok", "d": float(fit.estimates[cfg["D2.T2.parameter"]]),
            "sigma_hat": float(fit.estimates["sigma"]), "n_intervals": int(len(kept))}


def run_single_pulsar(template, scenario, method_rule, settings, sigma, noise_sd, n, key, cfg, workers):
    """n single-pulsar replicates; replicate r uses SeedSequence(key + (r,)). Returns a DataFrame by replicate."""
    args = [(template, scenario, sigma, method_rule, settings, noise_sd, tuple(key) + (r,), cfg) for r in range(n)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(_single_pulsar_replicate, args, chunksize=max(1, n // (8 * workers))))
    return pd.DataFrame(rows).sort_values("rep").reset_index(drop=True)
