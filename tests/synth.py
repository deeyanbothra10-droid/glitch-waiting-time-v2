"""Synthetic interval data with known parameters, for tests only (no real waiting times).

ln T = c_p + covariate terms + sigma * eps. Censoring is non-informative: each interval's censoring
time is fixed given its covariates at ln C = mu + sigma * z_(1 - censor_frac), so on average
censor_frac of intervals are censored.
"""
import numpy as np
import pandas as pd
from scipy import stats


def simulate_aft(model, coefs, intercepts, sigma, n_per, censor_frac, rng):
    pulsars = list(intercepts)
    n = n_per * len(pulsars)
    pulsar = np.repeat(pulsars, n_per)
    ln_dnu = rng.normal(np.log(1e-6), 1.0, n)
    ln_dnudot = rng.normal(np.log(1e-2), 0.7, n)
    ln_r = ln_dnu - ln_dnudot + np.log(1e3)        # R = (dff/1000)/df1 in catalogue units
    mu = np.array([intercepts[p] for p in pulsar], dtype=float)
    if model == "M1":
        mu = mu + ln_r
    mu = mu + coefs.get("b", 0.0) * ln_r if model == "M2" else mu
    if model in ("M3", "M4"):
        mu = mu + coefs["a"] * ln_dnu
    if model == "M4":
        mu = mu + coefs["d"] * ln_dnudot
    ln_t = mu + sigma * rng.normal(size=n)
    ln_c = mu + sigma * stats.norm.ppf(1 - censor_frac) if censor_frac > 0 else np.full(n, np.inf)
    return pd.DataFrame({
        "pulsar": pulsar, "t": np.exp(np.minimum(ln_t, ln_c)), "event": ln_t <= ln_c,
        "ln_R": ln_r, "ln_dnu_frac": ln_dnu, "ln_dnudot_frac": ln_dnudot,
    })
