"""Post-freeze forward forecast (NOT part of the pre-registration; results/forecast/ only).

The frozen lognormal AFT fits of results/primary.json give, for the wait T (days) after a pulsar's last
catalogued glitch g, ln T ~ Normal(mu, sigma) with mu = c_p + b ln R_g (M2) or mu = c_p (M0, baseline). Given
no glitch before the snapshot (T > t0), the forecast is the conditional distribution
    S(t | T > t0) = S(t) / S(t0),   S(t) = Phi(-(ln t - mu) / sigma).
Scoring (fixed with forecast v1): per pulsar, the log predictive density of the observed wait under the
conditional distribution (or its log survival at the scoring date if no glitch by then), M2 minus M0.
"""
import numpy as np
from scipy import stats

MJD_EPOCH = np.datetime64("1858-11-17")


def mjd_to_date(mjd):
    return str(MJD_EPOCH + np.timedelta64(int(np.floor(mjd)), "D"))


def conditional_quantiles(mu, sigma, t0, probs):
    """Waits t (days) with P(T <= t | T > t0) = q for each q in probs, for ln T ~ Normal(mu, sigma)."""
    log_s0 = stats.norm.logsf((np.log(t0) - mu) / sigma) if t0 > 0 else 0.0
    out = []
    for q in probs:
        z = stats.norm.isf(np.exp(np.log1p(-q) + log_s0))       # S(t) = (1 - q) S(t0)
        out.append(float(np.exp(mu + sigma * z)))
    return out


def log_predictive(mu, sigma, t0, t, observed):
    """Log predictive score of an outcome under the conditional forecast: log f(t) - log S(t0) (density of T
    in days) for an observed wait t, or log S(t) - log S(t0) if no glitch occurred by wait t."""
    log_s0 = stats.norm.logsf((np.log(t0) - mu) / sigma) if t0 > 0 else 0.0
    z = (np.log(t) - mu) / sigma
    if observed:
        return float(stats.norm.logpdf(z) - np.log(sigma) - np.log(t) - log_s0)
    return float(stats.norm.logsf(z) - log_s0)


def coverage(window, t, observed):
    """Whether an outcome falls in a forecast window [lo, hi] of waits: True or False when determinable;
    None for an unobserved wait t (no glitch by the scoring date) unless the window ends before t (a miss)."""
    lo, hi = window
    if observed:
        return bool(lo <= t <= hi)
    return False if hi < t else None


def outcome_for(next_epoch, g_mjd, snapshot_mjd, scoring_mjd):
    """Classify one pulsar's outcome from the scoring catalogue: 'void' if the next glitch's epoch is before the
    snapshot, 'observed' with its wait if it is at or before the scoring date, else 'censored' at the
    scoring date (next_epoch None: no later glitch catalogued)."""
    if next_epoch is not None and next_epoch < snapshot_mjd:
        return {"status": "void"}
    if next_epoch is not None and next_epoch <= scoring_mjd:
        return {"status": "observed", "t": float(next_epoch - g_mjd)}
    return {"status": "censored", "t": float(scoring_mjd - g_mjd)}


def forecast_window(mu, sigma, t0, g_mjd, probs):
    """The conditional quantiles for probs ({label: q}) as waits after g (days) and as calendar dates."""
    waits = dict(zip(probs, conditional_quantiles(mu, sigma, t0, list(probs.values()))))
    return {"mu": float(mu), "sigma": float(sigma), "t0_days": float(t0), "waits_days": waits,
            "dates": {k: mjd_to_date(g_mjd + v) for k, v in waits.items()}}


def scoring_eligibility(jodrell_bank, obs_end_mjd, last_glitch_mjd, glitch_epochs, g_mjd, snapshot_mjd):
    """Forecast v3: whether a pulsar is scored, with the rules it fails and their numbers.
    (a) covered by the D4g Jodrell Bank rule; (b) no documented obs_end before the snapshot (an obs_end equal
    to the last catalogued glitch epoch is the D4g placeholder, not a documented end of monitoring);
    (c) snapshot - epoch(g) no longer than the longest catalogued wait between consecutive glitches."""
    epochs = np.sort(np.asarray(glitch_epochs, float))
    longest = float(np.max(np.diff(epochs))) if len(epochs) > 1 else None
    gap = float(snapshot_mjd - g_mjd)
    documented = obs_end_mjd != last_glitch_mjd
    fails = {}
    if not jodrell_bank:
        fails["a"] = "not covered by the D4g Jodrell Bank rule"
    if documented and obs_end_mjd < snapshot_mjd:
        fails["b"] = f"documented obs_end MJD {obs_end_mjd:.1f} is before MJD {snapshot_mjd:.0f}"
    if longest is None or gap > longest:
        fails["c"] = (f"snapshot - g = {gap:.1f} d exceeds its longest catalogued wait {longest:.1f} d "
                      f"({gap / longest:.2f}x)" if longest else "fewer than two catalogued glitches")
    return {"scored": not fails, "fails": fails, "snapshot_minus_g_days": gap, "longest_wait_days": longest,
            "ratio": gap / longest if longest else None, "documented_obs_end": bool(documented)}
