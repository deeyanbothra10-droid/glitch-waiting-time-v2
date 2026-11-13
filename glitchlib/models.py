"""Lognormal accelerated-failure-time (AFT) models M0-M4 with right censoring (D3).

Every model has one intercept c_p per pulsar and one shared sigma:
    ln T = c_p + (covariate terms) + sigma * eps,   eps ~ N(0, 1)
    M0: c_p                     M1: ln R + c_p (ln R enters with slope fixed at 1)
    M2: b ln R + c_p            M3: a ln(dnu_frac) + c_p
    M4: a ln(dnu_frac) + d ln(dnudot_frac) + c_p
A completed interval contributes log f(T), the lognormal density of T in days; a censored one
contributes log S(c) = log P(T > c). The fit maximises the sum with step-halving Newton steps using the
exact gradient and Hessian (capped at MAX_NEWTON_STEPS); standard errors come from the inverse of the
negative Hessian.

Input data: a DataFrame with columns pulsar, t (days: T if completed, c if censored), event
(True = completed) and whichever of ln_R, ln_dnu_frac, ln_dnudot_frac the model uses.
"""
import hashlib
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import log_ndtr

MODELS = {
    "M0": {"covariates": {}, "offset": None},
    "M1": {"covariates": {}, "offset": "ln_R"},
    "M2": {"covariates": {"b": "ln_R"}, "offset": None},
    "M3": {"covariates": {"a": "ln_dnu_frac"}, "offset": None},
    "M4": {"covariates": {"a": "ln_dnu_frac", "d": "ln_dnudot_frac"}, "offset": None},
}


class FitError(RuntimeError):
    """The maximum-likelihood fit did not converge or the data cannot identify the model."""


@dataclass
class Design:
    model: str
    X: np.ndarray            # covariate columns, then one indicator column per pulsar
    offset: np.ndarray
    y: np.ndarray            # ln t
    event: np.ndarray        # True = completed
    names: list              # parameter names, ending with log_sigma
    pulsars: list
    pulsar_index: np.ndarray


@dataclass
class FitResult:
    model: str
    measurement_error: str
    names: list
    theta: np.ndarray        # estimates in the order of names (last is log_sigma)
    cov: np.ndarray
    loglik: float
    n: int
    n_completed: int
    pulsars: list
    data_fingerprint: str

    @property
    def n_params(self):
        return len(self.theta)

    @property
    def aic(self):
        return 2 * self.n_params - 2 * self.loglik

    @property
    def estimates(self):
        out = dict(zip(self.names, self.theta))
        out["sigma"] = float(np.exp(self.theta[-1]))
        return out

    @property
    def se(self):
        se = np.sqrt(np.diag(self.cov))
        out = dict(zip(self.names, se))
        out["sigma"] = float(np.exp(self.theta[-1]) * se[-1])   # delta method
        return out

    def intervals(self, level):
        """Wald intervals estimate +- z_q * SE from the inverse Hessian (D4b), for every fitted
        parameter; sigma's interval is exp() of the log_sigma interval, the scale it was fitted on."""
        z_q = stats.norm.ppf(0.5 + level / 2)
        se = np.sqrt(np.diag(self.cov))
        out = {n: (float(t - z_q * s), float(t + z_q * s)) for n, t, s in zip(self.names, self.theta, se)}
        out["sigma"] = tuple(float(np.exp(v)) for v in out["log_sigma"])
        return out

    def intercepts(self):
        return {p: self.estimates[f"c[{p}]"] for p in self.pulsars}

    def summary(self):
        return {
            "model": self.model, "measurement_error": self.measurement_error,
            "estimates": {k: float(v) for k, v in self.estimates.items()},
            "se": {k: float(v) for k, v in self.se.items()},
            "loglik": float(self.loglik), "aic": float(self.aic), "n_params": self.n_params,
            "n": self.n, "n_completed": self.n_completed,
        }


def _fingerprint(data):
    cols = data[["pulsar", "t", "event"]].astype(str).to_csv(index=False).encode()
    return hashlib.sha256(cols).hexdigest()


def build_design(model, data):
    if model not in MODELS:
        raise ValueError(f"Unknown model {model!r}; expected one of {sorted(MODELS)}.")
    spec = MODELS[model]
    needed = ["pulsar", "t", "event"] + list(spec["covariates"].values())
    needed += [spec["offset"]] if spec["offset"] else []
    missing = [c for c in needed if c not in data.columns]
    if missing:
        raise ValueError(f"{model} needs columns {missing}.")
    t = data["t"].to_numpy(float)
    if not np.all(np.isfinite(t)) or np.any(t <= 0):
        raise ValueError("Every t must be a finite positive number of days.")
    numeric = [c for c in needed if c not in ("pulsar", "t", "event")]
    values = data[numeric].to_numpy(float) if numeric else np.empty((len(data), 0))
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{model}: non-finite values in {numeric}.")
    pulsars = sorted(data["pulsar"].unique())
    pulsar_index = pd.Categorical(data["pulsar"], categories=pulsars).codes
    indicators = np.zeros((len(data), len(pulsars)))
    indicators[np.arange(len(data)), pulsar_index] = 1.0
    covariates = data[list(spec["covariates"].values())].to_numpy(float)
    offset = data[spec["offset"]].to_numpy(float) if spec["offset"] else np.zeros(len(data))
    names = list(spec["covariates"]) + [f"c[{p}]" for p in pulsars] + ["log_sigma"]
    return Design(model, np.hstack([covariates, indicators]), offset, np.log(t),
                  data["event"].to_numpy(bool), names, pulsars, pulsar_index)


def _pointwise(theta, X, offset, y, event):
    """Per-interval log-likelihood: log f(T) if completed, log S(c) if censored."""
    log_sigma = theta[-1]
    z = (y - (X @ theta[:-1] + offset)) / np.exp(log_sigma)
    return np.where(event, stats.norm.logpdf(z) - log_sigma - y, log_ndtr(-z))


def _loglik_grad_hess(theta, X, offset, y, event):
    """Log-likelihood with its exact gradient and Hessian in (coefficients, log_sigma)."""
    log_sigma = theta[-1]
    sigma = np.exp(log_sigma)
    z = (y - (X @ theta[:-1] + offset)) / sigma
    loglik = np.where(event, stats.norm.logpdf(z) - log_sigma - y, log_ndtr(-z)).sum()
    # h = phi(z) / Phi(-z), the inverse Mills ratio, computed in logs so it never overflows.
    h = np.exp(stats.norm.logpdf(z) - log_ndtr(-z))
    k = h * (h - z)                                   # dh/dz
    g_mu = np.where(event, z / sigma, h / sigma)       # d loglik_i / d mu_i
    g_s = np.where(event, z ** 2 - 1, h * z)           # d loglik_i / d log_sigma
    h_mumu = np.where(event, -1 / sigma ** 2, -k / sigma ** 2)
    h_mus = np.where(event, -2 * z / sigma, -(z * k + h) / sigma)
    h_ss = np.where(event, -2 * z ** 2, -z * (z * k + h))
    grad = np.append(X.T @ g_mu, g_s.sum())
    p = X.shape[1]
    hess = np.empty((p + 1, p + 1))
    hess[:p, :p] = X.T @ (h_mumu[:, None] * X)
    hess[:p, p] = hess[p, :p] = X.T @ h_mus
    hess[p, p] = h_ss.sum()
    return loglik, grad, hess


def _start(design):
    """Least squares on ln t, ignoring censoring: a starting point only."""
    coef, *_ = np.linalg.lstsq(design.X, design.y - design.offset, rcond=None)
    resid = design.y - design.offset - design.X @ coef
    return np.append(coef, np.log(max(resid.std(), 1e-3)))


MAX_NEWTON_STEPS = 200      # iteration cap: a fit that has not converged by then is a failure (FitError)


class _NonFinite(Exception):
    pass


def _evaluate(theta, args):
    with np.errstate(all="ignore"):         # overflow is caught below, not printed
        loglik, grad, hess = _loglik_grad_hess(theta, *args)
    if not (np.isfinite(loglik) and np.all(np.isfinite(grad)) and np.all(np.isfinite(hess))):
        raise _NonFinite
    return loglik, grad, hess


def _ascent_step(grad, hess):
    """Newton step g' (-H)^-1 when -H is positive definite; otherwise the step for -H + mu I with mu just
    large enough to make it positive definite, so the step always points uphill."""
    neg = -hess
    try:
        np.linalg.cholesky(neg)
        return np.linalg.solve(neg, grad), True
    except np.linalg.LinAlgError:
        lowest = float(np.linalg.eigvalsh(neg)[0])
        mu = -lowest * 1.5 + 1e-8 * max(1.0, float(np.max(np.abs(np.diag(neg)))))
        return np.linalg.solve(neg + mu * np.eye(len(grad)), grad), False


def _newton(theta, args, max_steps):
    """Newton's method with step halving. Converged when -H is positive definite and the Newton
    decrement g' (-H)^-1 g, the predicted gain in log-likelihood, is below 1e-12 (a rule that does not
    depend on the scale of the covariates). Returns (theta, loglik, cov), or None if non-finite values
    appear or the cap is reached."""
    try:
        loglik, grad, hess = _evaluate(theta, args)
        for _ in range(max_steps):
            step, newton = _ascent_step(grad, hess)
            if newton and float(grad @ step) < 1e-12:
                return theta, loglik, np.linalg.inv(-hess)
            scale = 1.0
            while True:
                trial = theta + scale * step
                try:
                    new = _evaluate(trial, args)
                    if new[0] >= loglik - 1e-12:
                        break
                except _NonFinite:                 # an overshooting trial step is simply rejected
                    pass
                scale /= 2
                if scale < 1e-10:
                    return None
            theta, (loglik, grad, hess) = trial, new
    except (_NonFinite, np.linalg.LinAlgError):
        return None
    return None


def maximise(design, start=None):
    """Return (theta, loglik, cov) at the maximum likelihood; raise FitError if it is not reached.

    Starts from start (warm start: permutations, leave-one-out) or from least squares on ln t. A warm
    start that fails is retried from least squares. Non-finite values or no convergence within
    MAX_NEWTON_STEPS raise FitError; nothing loops without limit."""
    if not design.event.any():
        raise FitError(f"{design.model}: no completed intervals.")
    for p in range(len(design.pulsars)):
        if not design.event[design.pulsar_index == p].any():
            raise FitError(f"{design.model}: pulsar {design.pulsars[p]} has no completed interval, "
                           "so its intercept has no finite maximum-likelihood estimate.")
    args = (design.X, design.offset, design.y, design.event)
    starts = ([np.asarray(start, float)] if start is not None else []) + [_start(design)]
    for theta0 in starts:
        result = _newton(theta0, args, MAX_NEWTON_STEPS)
        if result is not None:
            return result
    raise FitError(f"{design.model}: no convergence within {MAX_NEWTON_STEPS} Newton steps, or non-finite "
                   "log-likelihood (e.g. extreme covariate values).")


def _fit_ignore_errors(model, data, start):
    design = build_design(model, data)
    theta, loglik, cov = maximise(design, start)
    return FitResult(model, "ignore", design.names, theta, cov, float(loglik), len(data),
                     int(design.event.sum()), design.pulsars, _fingerprint(data))


# Measurement-error handling, keyed by the D6.primary_errors_in_lnR values the schema allows.
MEASUREMENT_ERROR_HANDLERS = {"ignore": _fit_ignore_errors}


def fit(model, data, measurement_error, start=None):
    """Fit one model. measurement_error comes from config D6.primary_errors_in_lnR (no default)."""
    if measurement_error not in MEASUREMENT_ERROR_HANDLERS:
        raise ValueError(f"D6: measurement-error handling {measurement_error!r} is not implemented; "
                         f"implemented options are {sorted(MEASUREMENT_ERROR_HANDLERS)}.")
    return MEASUREMENT_ERROR_HANDLERS[measurement_error](model, data, start)


def log_score(fit_result, data):
    """log f(T) or log S(c) of each row of data under a fitted model (pulsars must be in the fit)."""
    unknown = set(data["pulsar"]) - set(fit_result.pulsars)
    if unknown:
        raise ValueError(f"Cannot score pulsars absent from the fit (no intercept): {sorted(unknown)}.")
    spec = MODELS[fit_result.model]
    est = fit_result.estimates
    mu = np.array([est[f"c[{p}]"] for p in data["pulsar"]], dtype=float)
    for name, column in spec["covariates"].items():
        mu = mu + est[name] * data[column].to_numpy(float)
    if spec["offset"]:
        mu = mu + data[spec["offset"]].to_numpy(float)
    X = np.zeros((len(data), 0))
    theta = np.array([est["log_sigma"]])
    return _pointwise(theta, X, mu, np.log(data["t"].to_numpy(float)), data["event"].to_numpy(bool))


def lr_test(full, reduced):
    """Likelihood-ratio test of a reduced model nested in a full one, fitted to the same data."""
    if full.data_fingerprint != reduced.data_fingerprint:
        raise ValueError("lr_test: the two fits were made on different data.")
    if full.measurement_error != reduced.measurement_error:
        raise ValueError("lr_test: the two fits used different measurement-error handling.")
    f_spec, r_spec = MODELS[full.model], MODELS[reduced.model]
    nested = (f_spec["offset"] is None and r_spec["offset"] is None
              and set(r_spec["covariates"].items()) < set(f_spec["covariates"].items()))
    if not nested:
        raise ValueError(f"lr_test: {reduced.model} is not nested in {full.model}.")
    statistic = max(2.0 * (full.loglik - reduced.loglik), 0.0)
    df = full.n_params - reduced.n_params
    return {"full": full.model, "reduced": reduced.model, "statistic": statistic, "df": df,
            "p_value": float(stats.chi2.sf(statistic, df))}


def loo_cv(model, data, measurement_error):
    """Leave one interval out, refit, score it (log f(T) or log S(c)); mean over all intervals (D4a)."""
    data = data.reset_index(drop=True)
    counts = data["pulsar"].value_counts()
    if (counts < 2).any():
        raise ValueError(f"loo_cv: pulsars with a single interval cannot be scored when it is left out: "
                         f"{sorted(counts[counts < 2].index)}.")
    full = fit(model, data, measurement_error)
    scores = np.empty(len(data))
    for i in range(len(data)):
        train = data.drop(index=i)
        held_out = data.iloc[[i]]
        scores[i] = log_score(fit(model, train, measurement_error, start=full.theta), held_out)[0]
    return {"model": model, "mean_log_score": float(scores.mean()), "n": len(data), "scores": scores}


def fit_poisson_per_pulsar(data):
    """D5 baseline: a homogeneous Poisson process (exponential waits) per pulsar, reported descriptively.

    The rate MLE is completed intervals / total time at risk; completed intervals contribute
    log f(T) = log(rate) - rate*T and censored ones log S(c) = -rate*c.
    """
    rows = []
    for pulsar, g in data.groupby("pulsar"):
        d = int(g["event"].sum())
        exposure = float(g["t"].sum())
        rate = d / exposure
        loglik = d * np.log(rate) - rate * exposure if d else 0.0
        rows.append({"pulsar": pulsar, "n": len(g), "n_completed": d, "exposure_days": exposure,
                     "rate_per_day": rate, "se_log_rate": 1 / np.sqrt(d) if d else np.inf,
                     "loglik": loglik})
    return pd.DataFrame(rows)
