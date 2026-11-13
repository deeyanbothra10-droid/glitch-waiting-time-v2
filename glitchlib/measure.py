"""Emulate how papers measure a glitch's Δν̇/ν̇, so the simulations see the same biases as the catalogue.

The synthetic signal is the fractional spin-down change s(t) = ν̇(t)/ν̇_0 - 1 (positive = faster
spin-down). Each glitch k at t_k adds, for t >= t_k,
    J_k * max(0, 1 - (t - t_k)/t_rec_k)          a step that relaxes linearly to zero over t_rec_k
  + A_k * exp(-(t - t_k)/tau_k)                  an exponential, A_k = Q_k Δν_k / (tau_k |ν̇|)
Each window between glitches is sampled from its start at its own cadence (D4d: the cadence of the
source of the glitch that starts it), with Gaussian white noise. Δν̇/ν̇ is then measured per glitch
by one of four methods:
  window          one constant level per inter-glitch window; step = level after - level before.
                  No recovery or exponential terms.
  local-relaxation per window: level + linear recovery slope (+ that window's own exponential if
                  detected). Total = this window's fit at t_k - previous window's fit at t_k.
  global          one fit over the whole series: constant + one global slope (ν̈) + a permanent step
                  per glitch (+ its exponential if detected). Total = step + exponential amplitude.
  global-noexp    as global without exponential terms. Total = step.
  short-window    (D4f) a constant level over w days before and w days after the glitch, each span cut
                  at the neighbouring glitch; step = level after - level before. No recovery or
                  exponential terms. w is set per glitch (by its source).
Exponential detection (D4d): in each window after a glitch, the exponential is included only if it
lowers the AIC of the level + slope fit of that window by at least min_aic_drop; the same per-window
decision is used by the local and global fits. Measurability (D4c): a glitch needs min_after samples
after it (window and short-window: also min_before before it); otherwise its Δν̇ is missing.
Limits: works in the ν̇ domain, not by phase-coherent fitting of times of arrival; no red timing
noise; exponential timescales are known to the fit (only amplitudes are fitted), so fits are linear;
glitch epochs are known exactly.
"""
import numpy as np

METHODS = ("window", "local-relaxation", "global", "global-noexp", "short-window")


def exponential_amplitude(Q, dnu_frac, tau_d_days, tau_sd_days):
    """A = Q Δν / (tau_d |ν̇|) as a fraction of ν̇, using ν/|ν̇| = 2 tau_sd."""
    return Q * dnu_frac * 2.0 * tau_sd_days / tau_d_days


def true_signal(t, glitch_times, J, t_rec, A, tau_d):
    s = np.zeros_like(t, dtype=float)
    for tk, jk, rk, ak, dk in zip(glitch_times, J, t_rec, A, tau_d):
        dt = t - tk
        after = dt >= 0
        s[after] += jk * np.clip(1.0 - dt[after] / rk, 0.0, None) + ak * np.exp(-dt[after] / dk)
    return s


def sample_times(t_start, t_end, glitch_times, cadence):
    """Sample epochs. cadence is one number, or one per window: window 0 runs from t_start to the first
    glitch, window j >= 1 from glitch j-1 to glitch j (or t_end). Each window is sampled from its start."""
    edges = np.concatenate([[t_start], np.asarray(glitch_times, float), [t_end]])
    cadence = np.broadcast_to(np.asarray(cadence, float), (len(edges) - 1,))
    return np.concatenate([np.arange(edges[j], edges[j + 1], cadence[j]) for j in range(len(edges) - 1)])


def sample_series(t_start, t_end, glitch_times, cadence, J, t_rec, A, tau_d, noise_sd, rng):
    """Samples of s(t) with white noise of standard deviation noise_sd: (t, noisy, clean)."""
    t = sample_times(t_start, t_end, glitch_times, cadence)
    clean = true_signal(t, glitch_times, J, t_rec, A, tau_d)
    return t, clean + noise_sd * rng.normal(size=t.size), clean


def _segments(t, glitch_times):
    """Sample indices of each window: before the first glitch, then after each glitch up to the next."""
    edges = np.concatenate([[-np.inf], np.asarray(glitch_times, float), [np.inf]])
    return [np.flatnonzero((t >= edges[j]) & (t < edges[j + 1])) for j in range(len(edges) - 1)]


def _lstsq(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    s2 = resid @ resid / (len(y) - X.shape[1])
    return beta, s2 * np.linalg.inv(X.T @ X), resid @ resid


def measure_short_window(t, y, glitch_times, w, min_after, min_before):
    """D4f: level over [t_k - w_k, t_k) minus level over [t_k, t_k + w_k), each cut at the neighbouring glitch."""
    total = np.full(len(glitch_times), np.nan)
    se = np.full(len(glitch_times), np.nan)
    edges = np.concatenate([[-np.inf], glitch_times, [np.inf]])
    for k, tk in enumerate(glitch_times):
        pre = y[(t >= max(tk - w[k], edges[k])) & (t < tk)]
        post = y[(t >= tk) & (t < min(tk + w[k], edges[k + 2]))]
        if len(post) < max(min_after, 2) or len(pre) < max(min_before, 2):
            continue
        total[k] = post.mean() - pre.mean()
        se[k] = np.sqrt(post.var(ddof=1) / len(post) + pre.var(ddof=1) / len(pre))
    return total, se


def _local_basis(tt, origin, tau):
    """Level + linear slope from origin, plus an exponential from origin when tau is given."""
    cols = [np.ones_like(tt), tt - origin]
    if tau is not None:
        cols.append(np.exp(-(tt - origin) / tau))
    return np.column_stack(cols)


def aic_drop(tt, yy, origin, tau):
    """AIC(level + slope) - AIC(level + slope + exponential) for one window, Gaussian errors with
    unknown variance: n ln(RSS_without / RSS_with) - 2."""
    rss0 = _lstsq(_local_basis(tt, origin, None), yy)[2]
    rss1 = _lstsq(_local_basis(tt, origin, tau), yy)[2]
    spread = np.sum((yy - yy.mean()) ** 2)
    if rss0 <= 1e-20 * spread:      # level + slope already fits to rounding error: nothing to explain
        return -np.inf
    if rss1 <= 0:
        return np.inf
    return len(yy) * np.log(rss0 / rss1) - 2.0


def detect_exponentials(t, y, glitch_times, tau_d, min_aic_drop):
    """Per glitch: is its exponential detected in the window after it (D4d)? Needs 4 samples (3 + 1)."""
    detected = np.zeros(len(glitch_times), bool)
    for k, idx in enumerate(_segments(t, glitch_times)[1:]):
        if len(idx) >= 4:
            detected[k] = aic_drop(t[idx], y[idx], glitch_times[k], tau_d[k]) >= min_aic_drop
    return detected


def measure_window(t, y, glitch_times, min_after, min_before):
    segs = _segments(t, glitch_times)
    counts = np.array([len(idx) for idx in segs])
    level = np.array([y[idx].mean() if len(idx) >= 2 else np.nan for idx in segs])
    var = np.array([y[idx].var(ddof=1) / len(idx) if len(idx) >= 2 else np.nan for idx in segs])
    total, se = level[1:] - level[:-1], np.sqrt(var[1:] + var[:-1])
    ok = (counts[1:] >= min_after) & (counts[:-1] >= min_before)
    return np.where(ok, total, np.nan), np.where(ok, se, np.nan)


def measure_local(t, y, glitch_times, tau_d, fit_exp, min_after):
    segs = _segments(t, glitch_times)
    fits = []
    for j, idx in enumerate(segs):
        # Window 0 lies before the first glitch: level + slope only. Window j >= 1 starts at glitch j-1
        # and includes that glitch's exponential when it was detected.
        if j == 0:
            origin, tau = (t[idx][0] if len(idx) else 0.0), None
        else:
            origin, tau = glitch_times[j - 1], (tau_d[j - 1] if fit_exp[j - 1] else None)
        n_par = 2 if tau is None else 3
        if len(idx) < n_par + 1:          # a fit needs at least one degree of freedom for its error
            fits.append(None)
            continue
        beta, cov, _ = _lstsq(_local_basis(t[idx], origin, tau), y[idx])
        fits.append((origin, tau, beta, cov))
    counts = np.array([len(idx) for idx in segs])
    total = np.full(len(glitch_times), np.nan)
    se = np.full(len(glitch_times), np.nan)
    for k, tk in enumerate(glitch_times):
        before, after = fits[k], fits[k + 1]
        if counts[k + 1] < min_after or before is None or after is None:
            continue
        x_pre = _local_basis(np.array([tk]), before[0], before[1])[0]    # previous window's fit at t_k
        x_post = _local_basis(np.array([tk]), after[0], after[1])[0]     # this window's fit at t_k
        total[k] = x_post @ after[2] - x_pre @ before[2]
        se[k] = np.sqrt(x_post @ after[3] @ x_post + x_pre @ before[3] @ x_pre)
    return total, se


def measure_global(t, y, glitch_times, tau_d, fit_exp, min_after):
    segs = _segments(t, glitch_times)
    counts = np.array([len(idx) for idx in segs])
    # A glitch needs min_after samples after it; otherwise its step cannot be told apart from the next
    # one, so its columns are left out and it is reported missing.
    measurable = counts[1:] >= min_after
    mid = 0.5 * (t[0] + t[-1])
    columns = [np.ones_like(t), (t - mid) / (t[-1] - t[0])]
    step_col, exp_col = {}, {}
    for k, tk in enumerate(glitch_times):
        if not measurable[k]:
            continue
        after = t >= tk
        step_col[k] = len(columns)
        columns.append(after.astype(float))
        if fit_exp[k]:
            exp_col[k] = len(columns)
            columns.append(after * np.exp(-np.clip(t - tk, 0.0, None) / tau_d[k]))
    beta, cov, _ = _lstsq(np.column_stack(columns), y)
    total = np.full(len(glitch_times), np.nan)
    se = np.full(len(glitch_times), np.nan)
    for k in step_col:
        c = np.zeros(len(beta))
        c[step_col[k]] = 1.0
        if k in exp_col:
            c[exp_col[k]] = 1.0
        total[k] = c @ beta
        se[k] = np.sqrt(c @ cov @ c)
    return total, se


def measure(method, t, y, glitch_times, tau_d, min_after, min_before_window, min_aic_drop, short_window_days):
    """Measured Δν̇/ν̇ (fraction), its standard error, and whether each glitch's exponential was fitted.

    NaN where the method cannot measure a glitch. Window, global-noexp and short-window never fit
    exponentials. short_window_days: the span w per glitch, needed only by short-window (else None).
    """
    glitch_times = np.asarray(glitch_times, float)
    tau_d = np.asarray(tau_d, float)
    none = np.zeros(len(glitch_times), bool)
    if method == "window":
        return (*measure_window(t, y, glitch_times, min_after, min_before_window), none)
    if method == "global-noexp":
        return (*measure_global(t, y, glitch_times, tau_d, none, min_after), none)
    if method == "short-window":
        w = np.broadcast_to(np.asarray(short_window_days, float), (len(glitch_times),))
        if not np.all(np.isfinite(w)) or np.any(w <= 0):
            raise ValueError("short-window needs a positive span w for every glitch.")
        return (*measure_short_window(t, y, glitch_times, w, min_after, min_before_window), none)
    if method in ("local-relaxation", "global"):
        fit_exp = detect_exponentials(t, y, glitch_times, tau_d, min_aic_drop)
        fitter = measure_local if method == "local-relaxation" else measure_global
        return (*fitter(t, y, glitch_times, tau_d, fit_exp, min_after), fit_exp)
    raise ValueError(f"Unknown measurement method {method!r}; expected one of {METHODS}.")
