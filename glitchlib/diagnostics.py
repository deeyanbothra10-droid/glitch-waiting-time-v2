"""Secondary analyses and diagnostics (D4, D4a, D15) built on glitchlib.models.

Every choice is an argument with no default; scripts pass the values from config/decisions.yaml.
These functions accept any data; whether real waiting times may be used is decided by the caller
through glitchlib.guards before it calls them.
"""
import numpy as np
import pandas as pd
from scipy import stats

from glitchlib import models

# Relative tolerance when comparing |b_perm| with |b_obs|, so a permutation equal to the observed
# order is not undercounted because the optimiser stops a hair away from the same maximum.
_PERMUTATION_TIE_RTOL = 1e-8


def spearman_within(data, x, y, min_completed, label, level):
    """Per-pulsar Spearman correlation of columns x and y over completed intervals only (D4a, D4b).

    Pulsars with fewer than min_completed completed intervals are skipped. The pooled value is the
    Fisher-z weighted mean with weights n_p - 3, back-transformed, with interval
    tanh(z +- z_q / sqrt(sum(n_p - 3))). Pulsars with |r| = 1 (infinite Fisher z) or an undefined r
    are flagged in the per-pulsar table and left out of the pooled value. Descriptive: no p-values.
    """
    completed = data[data["event"]]
    rows = []
    for pulsar, g in completed.groupby("pulsar"):
        if len(g) < min_completed:
            continue
        rows.append({"pulsar": pulsar, "n": len(g), "r": float(stats.spearmanr(g[x], g[y]).statistic)})
    per = pd.DataFrame(rows, columns=["pulsar", "n", "r"])
    # scipy can return 0.9999999999999999 for a perfect rank match, so compare with a rounding tolerance.
    per["abs_r_one"] = np.isclose(per["r"].abs(), 1.0, rtol=0, atol=1e-12)
    per["undefined"] = ~np.isfinite(per["r"])
    pooled = per[~per["abs_r_one"] & ~per["undefined"]]
    weights = pooled["n"] - 3
    if len(pooled):
        z_bar = float(np.sum(weights * np.arctanh(pooled["r"])) / np.sum(weights))
        half = stats.norm.ppf(0.5 + level / 2) / np.sqrt(np.sum(weights))
        pooled_r, interval = float(np.tanh(z_bar)), [float(np.tanh(z_bar - half)), float(np.tanh(z_bar + half))]
    else:
        pooled_r, interval = np.nan, [np.nan, np.nan]
    return {"label": label, "x": x, "y": y, "per_pulsar": per, "pooled_r": pooled_r,
            "interval": interval, "level": level, "n_pulsars_pooled": int(len(pooled)),
            "n_intervals_pooled": int(pooled["n"].sum()),
            "excluded_abs_r_one": sorted(per.loc[per["abs_r_one"], "pulsar"]),
            "excluded_undefined": sorted(per.loc[per["undefined"], "pulsar"])}


def _permute_within(values, groups, rng):
    out = values.copy()
    for idx in groups:
        out[idx] = values[rng.permutation(idx)]
    return out


def permutation_test(data, intervals, n_permutations, seed, measurement_error):
    """Stratified permutation test of the M2 slope (D4a).

    ln R is shuffled within each pulsar across all its intervals; t, event and pulsar stay fixed.
    Statistic: b-hat from the full censored M2 likelihood. Two-sided
    p = (1 + #{|b_perm| >= |b_obs|}) / (1 + N).
    """
    if intervals != "all":
        raise ValueError(f"D4a: permutation over intervals={intervals!r} is not implemented; only 'all'.")
    data = data.reset_index(drop=True)
    observed = models.fit("M2", data, measurement_error)
    b_obs = observed.estimates["b"]
    design = models.build_design("M2", data)
    ln_r = design.X[:, 0].copy()          # column 0 of the M2 design is ln R
    groups = [np.flatnonzero(design.pulsar_index == p) for p in range(len(design.pulsars))]
    rng = np.random.default_rng(seed)
    b_perm = np.empty(n_permutations)
    for k in range(n_permutations):
        design.X[:, 0] = _permute_within(ln_r, groups, rng)
        theta, _, _ = models.maximise(design, start=observed.theta)
        b_perm[k] = theta[0]
    exceed = np.abs(b_perm) >= np.abs(b_obs) * (1 - _PERMUTATION_TIE_RTOL)
    return {"b_obs": float(b_obs), "p_value": float((1 + exceed.sum()) / (1 + n_permutations)),
            "n_permutations": int(n_permutations), "seed": int(seed), "b_perm": b_perm}


def placebo_data(data):
    """Rows with the previous wait as outcome: t <- prev_t, event <- prev_event (D15).

    Intervals with no previous wait (a pulsar's first glitch) are dropped, since they have no outcome.
    """
    for col in ("prev_t", "prev_event"):
        if col not in data.columns:
            raise ValueError(f"placebo_data needs column {col!r}.")
    out = data[data["prev_t"].notna()].copy()
    out["t"] = out["prev_t"].astype(float)
    out["event"] = out["prev_event"].astype(bool)
    return out.reset_index(drop=True)


def placebo_fit(data, measurement_error):
    """M2 fitted with the previous wait T_(k-1) as outcome (D15)."""
    return models.fit("M2", placebo_data(data), measurement_error)


def _stratum_terms(beta, x, y):
    """Conditional log-likelihood of one stratum with its first and second derivatives in beta.

    Given that m of the stratum's n glitches have y = 1, the probability that it is exactly these m is
        exp(beta * sum_{y=1} x) / e_m,   e_m = sum over all m-subsets S of exp(beta * sum_{j in S} x_j),
    the elementary symmetric polynomial of degree m in a_j = exp(beta x_j). e_m and its beta-derivatives
    are built one glitch at a time with the recursion e_k <- e_k + a_j e_(k-1). x is centred first;
    this multiplies numerator and e_m by the same factor, so the likelihood is unchanged.
    """
    x = x - x.mean()
    m = int(y.sum())
    e = np.zeros(m + 1)
    d1 = np.zeros(m + 1)
    d2 = np.zeros(m + 1)
    e[0] = 1.0
    for xj in x:
        aj = np.exp(beta * xj)
        e_new = e.copy()
        e_new[1:] += aj * e[:-1]
        d1_new = d1.copy()
        d1_new[1:] += aj * (xj * e[:-1] + d1[:-1])
        d2_new = d2.copy()
        d2_new[1:] += aj * (xj ** 2 * e[:-1] + 2 * xj * d1[:-1] + d2[:-1])
        e, d1, d2 = e_new, d1_new, d2_new
    mean_s = d1[m] / e[m]                      # expected sum of x over the chosen subset
    loglik = beta * x[y].sum() - np.log(e[m])
    score = x[y].sum() - mean_s
    information = d2[m] / e[m] - mean_s ** 2   # variance of that sum
    return loglik, score, information


def conditional_logit(strata, y, x):
    """Conditional (stratum-matched) logistic regression with one slope; stratum intercepts are
    conditioned out. Returns beta-hat, its SE (inverse observed information) and the log-likelihood."""
    y = np.asarray(y, bool)
    x = np.asarray(x, float)
    groups = [np.flatnonzero(strata == s) for s in pd.unique(strata)]
    groups = [g for g in groups if 0 < y[g].sum() < len(g)]     # only strata where y varies inform beta
    if not groups:
        raise ValueError("conditional_logit: no stratum has both outcomes.")
    # Separation: if in every informative stratum the y = 1 rows have the largest possible sum of x (or in
    # every stratum the smallest), the likelihood keeps rising as beta -> +inf (or -inf): no finite estimate.
    def observed_is(extreme):
        return all(np.isclose(x[g][y[g]].sum(), extreme(np.sort(x[g]), int(y[g].sum())), rtol=0, atol=1e-12)
                   for g in groups)
    if observed_is(lambda xs, m: xs[-m:].sum()) or observed_is(lambda xs, m: xs[:m].sum()):
        raise models.FitError("conditional_logit: complete separation, the estimate is infinite.")

    def total(beta):
        parts = np.array([_stratum_terms(beta, x[g], y[g]) for g in groups])
        return parts.sum(axis=0)

    beta = 0.0
    for _ in range(100):
        loglik, score, info = total(beta)
        step = score / info
        # Halve the Newton step until the log-likelihood does not fall.
        while total(beta + step)[0] < loglik - 1e-12 and abs(step) > 1e-12:
            step /= 2
        beta += step
        if score * score / info < 1e-14:      # Newton decrement
            break
    else:
        raise models.FitError("conditional_logit did not converge (the estimate may be infinite).")
    loglik, score, info = total(beta)
    if not (np.isfinite(beta) and np.isfinite(info) and info > 0):
        raise models.FitError("conditional_logit: no finite estimate or standard error.")
    return float(beta), float(1 / np.sqrt(info)), float(loglik)


def missingness_check(data, outcome_col, ln_wait_col, following_completed_col, level):
    """Conditional (pulsar-stratified) logistic regression of an indicator on ln(following wait) (D15, D4b).

    Rows whose following wait is censored are excluded and counted. Pulsars where the indicator never
    varies carry no information about the slope and are listed. Returns the odds ratio per unit
    ln wait with its Wald interval at the given level, the Wald z and the two-sided p.
    """
    completed = data[following_completed_col].astype(bool)
    used = data[completed]
    no_variation = sorted(p for p, g in used.groupby("pulsar") if g[outcome_col].nunique() < 2)
    beta, se, loglik = conditional_logit(used["pulsar"].to_numpy(), used[outcome_col].to_numpy(bool),
                                         used[ln_wait_col].to_numpy(float))
    z = beta / se
    half = stats.norm.ppf(0.5 + level / 2) * se
    return {"outcome": outcome_col, "log_odds_ratio": beta, "se_log_odds_ratio": se,
            "odds_ratio": float(np.exp(beta)),
            "interval": [float(np.exp(beta - half)), float(np.exp(beta + half))], "level": level,
            "z": float(z), "p_value": float(2 * stats.norm.sf(abs(z))), "loglik": loglik,
            "n": int(len(used)), "n_excluded_censored_following_wait": int((~completed).sum()),
            "n_pulsars": int(used["pulsar"].nunique()), "pulsars_without_variation": no_variation}


def implied_beta(m1_fit, tau_sd_days):
    """Implied beta per pulsar from the M1 (slope-1) intercepts: beta_p = 2e-3 tau_sd / exp(c_p) - 0.5.

    The SE uses the delta method: d beta / d c_p = -2e-3 tau_sd / exp(c_p).
    tau_sd_days maps pulsar -> spin-down timescale in days.
    """
    if m1_fit.model != "M1":
        raise ValueError(f"implied_beta needs the M1 fit (ln T = ln R + c_p), got {m1_fit.model}.")
    missing = [p for p in m1_fit.pulsars if p not in tau_sd_days]
    if missing:
        raise ValueError(f"implied_beta: no tau_sd for {missing}.")
    est, se = m1_fit.estimates, m1_fit.se
    rows = []
    for p in m1_fit.pulsars:
        c, scale = est[f"c[{p}]"], 2e-3 * tau_sd_days[p]
        rows.append({"pulsar": p, "c_p": c, "se_c_p": se[f"c[{p}]"], "tau_sd_days": tau_sd_days[p],
                     "beta": scale / np.exp(c) - 0.5, "se_beta": scale / np.exp(c) * se[f"c[{p}]"]})
    return pd.DataFrame(rows)


def collinearity_summary(data, method, min_n):
    """Within-pulsar correlation of ln dnu_frac and ln dnudot_frac, and the within-pulsar SD of ln R (D15, D4b).

    For pulsars with at least min_n rows: the per-pulsar correlation, its unweighted median |r|
    (the D15 rule), and the pooled within-pulsar correlation (both variables demeaned within pulsar,
    then one Pearson correlation over all rows).
    """
    correlate = {"pearson": stats.pearsonr, "spearman": stats.spearmanr}
    if method not in correlate:
        raise ValueError(f"collinearity_summary: method must be one of {sorted(correlate)}.")
    rows, kept = [], []
    for pulsar, g in data.groupby("pulsar"):
        if len(g) < min_n:
            continue
        kept.append(g)
        r = float(correlate[method](g["ln_dnu_frac"], g["ln_dnudot_frac"]).statistic)
        rows.append({"pulsar": pulsar, "n": len(g), "r": r, "sd_ln_R": float(g["ln_R"].std(ddof=1))})
    per = pd.DataFrame(rows, columns=["pulsar", "n", "r", "sd_ln_R"])
    pooled_r = np.nan
    if kept:
        both = pd.concat(kept)
        demeaned = both[["ln_dnu_frac", "ln_dnudot_frac"]] - both.groupby("pulsar")[
            ["ln_dnu_frac", "ln_dnudot_frac"]].transform("mean")
        pooled_r = float(stats.pearsonr(demeaned["ln_dnu_frac"], demeaned["ln_dnudot_frac"]).statistic)
    return {"method": method, "per_pulsar": per, "median_abs_r": float(per["r"].abs().median()),
            "pooled_within_pulsar_r": pooled_r, "n_pulsars": int(len(per))}
