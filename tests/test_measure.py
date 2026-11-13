"""The Δν̇ measurement emulator on synthetic spin-down series with known answers (D4d correctness tests)."""
import numpy as np
import pytest

from glitchlib import measure

GLITCHES = np.array([1000.0, 2300.0, 3100.0, 4500.0])
J = np.array([5e-3, 2e-3, 8e-3, 4e-3])
A = np.array([3e-3, 1e-3, 2e-3, 5e-3])
TAU = np.full(4, 20.0)
NEVER = np.full(4, 1e12)          # recovery time so long the steps never relax
AIC = 10.0                        # detection threshold used in these tests


def _series(t_rec, amp, noise_sd=0.0, seed=0, cadence=5.0):
    rng = np.random.default_rng(seed)
    t, y, _ = measure.sample_series(0.0, 6000.0, GLITCHES, cadence, J, t_rec, amp, TAU, noise_sd, rng)
    return t, y


def _measure(method, t, y, glitches=GLITCHES, tau=TAU, threshold=AIC):
    return measure.measure(method, t, y, glitches, tau, 4, 4, threshold, 40.0)


# ---- D4d (a): window, no exponentials, T <= t_rec ------------------------------------------------
@pytest.mark.parametrize("ratio", [0.2, 0.5, 1.0])
def test_a_window_linear_recovery(ratio):
    # One glitch at 1000 after a flat baseline; the window after it has length T = 1000 and the
    # step relaxes over t_rec = T / ratio. Expected J (1 - T / (2 t_rec)). With samples every c days
    # from the glitch, the sample mean is J (1 - (T - c) / (2 t_rec)); c = 0.01 makes that negligible.
    T, j = 1000.0, 4e-3
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 1000.0 + T, [1000.0], 0.01, [j], [T / ratio], [0.0], [10.0], 0.0, rng)
    value, _, fitted = measure.measure("window", t, y, [1000.0], [10.0], 4, 4, AIC, None)
    assert value[0] == pytest.approx(j * (1 - ratio / 2), rel=1e-4)
    assert not fitted[0]


# ---- D4d (b): window, one exponential, constant step ---------------------------------------------
@pytest.mark.parametrize("T,tau", [(500.0, 20.0), (100.0, 50.0), (60.0, 100.0)])
def test_b_window_exponential_with_constant_step(T, tau):
    # Δν̇_p + (Δν_d / T)(1 - e^(-T/tau_d)); as a fraction of ν̇ the second term is A tau_d / T (1 - e^(-T/tau_d)).
    j, a = 4e-3, 3e-3
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 1000.0 + T, [1000.0], 0.01, [j], [1e12], [a], [tau], 0.0, rng)
    value, _, _ = measure.measure("window", t, y, [1000.0], [tau], 4, 4, AIC, None)
    assert value[0] == pytest.approx(j + a * tau / T * (1 - np.exp(-T / tau)), rel=1e-4)


# ---- D4d (c): global and local exact when t_rec = T and every exponential is fitted --------------
def test_c_fitted_methods_exact_when_recovery_ends_at_next_glitch():
    # Vortex-creep toy: each step recovers at a common rate nu_ddot, so t_rec = T = J / nu_ddot for every
    # glitch, and the window before the first glitch is the end of an earlier (unmeasured) recovery at the
    # same rate. Then one global slope is exact, and so is a per-window line.
    nu_ddot = 5e-6
    j = np.array([5e-3, 2e-3, 8e-3, 4e-3])
    T = j / nu_ddot                                     # 1000, 400, 1600, 800 d
    glitches = 1000.0 + np.concatenate([[0.0], np.cumsum(T[:-1])])
    end = glitches[-1] + T[-1]
    earlier = 1000.0 - 3e-3 / nu_ddot                  # unmeasured glitch whose recovery ends at 1000
    all_t = np.concatenate([[earlier], glitches])
    amp = np.array([0.0, 3e-3, 1e-3, 2e-3, 5e-3])
    tau = np.full(5, 20.0)
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(700.0, end, all_t[1:], 2.0, np.concatenate([[3e-3], j])[1:], T, amp[1:],
                                    tau[1:], 0.0, rng)
    y = y + measure.true_signal(t, [earlier], [3e-3], [3e-3 / nu_ddot], [0.0], [20.0])
    for method in ("local-relaxation", "global"):
        value, _, fitted = measure.measure(method, t, y, glitches, tau[1:], 4, 4, AIC, None)
        assert fitted.all(), method
        assert value == pytest.approx(j + amp[1:], rel=1e-6), method


# ---- D4d (d): an exponential below the detection rule is never fitted -----------------------------
def test_d_undetected_exponential_is_never_fitted():
    # Two glitches; the first exponential is far below the noise, the second far above it.
    glitches = np.array([1000.0, 2000.0])
    tau = np.array([30.0, 30.0])
    rng = np.random.default_rng(3)
    t, y, _ = measure.sample_series(0.0, 3000.0, glitches, 1.0, [5e-3, 5e-3], [1e12, 1e12], [1e-7, 5e-3], tau,
                                    1e-4, rng)
    for method in ("local-relaxation", "global"):
        _, _, fitted = measure.measure(method, t, y, glitches, tau, 4, 4, AIC, None)
        assert list(fitted) == [False, True], method
    # Left out means exactly the fit without that term: the first glitch's local value is the same as
    # when no exponential is ever fitted.
    with_rule, _, _ = measure.measure("local-relaxation", t, y, glitches, tau, 4, 4, AIC, None)
    never, _, _ = measure.measure("local-relaxation", t, y, glitches, tau, 4, 4, np.inf, None)
    assert with_rule[0] == never[0]


def test_aic_drop_known_answer():
    # n ln(RSS_without / RSS_with) - 2, checked against fits done here.
    rng = np.random.default_rng(4)
    tt = np.arange(0.0, 200.0, 2.0)
    yy = 1e-3 - 1e-6 * tt + 2e-4 * np.exp(-tt / 25.0) + 5e-5 * rng.normal(size=tt.size)
    x0 = np.column_stack([np.ones_like(tt), tt])
    x1 = np.column_stack([x0, np.exp(-tt / 25.0)])
    rss = [np.sum((yy - x @ np.linalg.lstsq(x, yy, rcond=None)[0]) ** 2) for x in (x0, x1)]
    assert measure.aic_drop(tt, yy, 0.0, 25.0) == pytest.approx(tt.size * np.log(rss[0] / rss[1]) - 2)


# ---- other properties ------------------------------------------------------------------------
@pytest.mark.parametrize("method", measure.METHODS)
def test_permanent_steps_are_recovered_exactly_by_every_method(method):
    t, y = _series(NEVER, np.zeros(4))
    value, _, fitted = _measure(method, t, y)
    assert value == pytest.approx(J, rel=1e-8)
    assert not fitted.any()


@pytest.mark.parametrize("method", ["local-relaxation", "global"])
def test_methods_with_exponentials_report_the_total(method):
    t, y = _series(NEVER, A)
    value, _, fitted = _measure(method, t, y)
    assert fitted.all()
    assert value == pytest.approx(J + A, rel=1e-8)


@pytest.mark.parametrize("method", ["window", "global-noexp"])
def test_methods_without_exponentials_miss_most_of_them(method):
    # Only the exponential's average over a window leaks into the step (and the previous glitch's
    # tail leaks into the "before" level), so the result stays far from the total J + A.
    t, y = _series(NEVER, A)
    value, _, _ = _measure(method, t, y)
    assert np.all(np.abs(value - J) < 0.5 * A)


def test_exponential_amplitude_known_answer():
    # A = Q dnu_frac 2 tau_sd / tau_d = 0.02 * 1e-6 * 2e6 / 10 = 4e-3
    assert measure.exponential_amplitude(0.02, 1e-6, 10.0, 1e6) == pytest.approx(4e-3)


def _scatter_and_errors(amp):
    values, errors, fitted = [], [], []
    for seed in range(300):
        t, y = _series(NEVER, amp, noise_sd=1e-3, seed=seed)
        v, se, f = _measure("local-relaxation", t, y)
        values.append(v)
        errors.append(se)
        fitted.append(f)
    return np.std(values, axis=0) / np.mean(errors, axis=0), np.array(fitted)


def test_reported_errors_match_scatter_when_the_model_is_fixed():
    # Exponentials 10x the noise: detected in every realisation, so each glitch is always fitted by the
    # same model, and the least-squares error bars match the scatter.
    ratio, fitted = _scatter_and_errors(10 * A)
    assert fitted.all()
    assert np.all(np.abs(ratio - 1) < 0.15)


def test_marginal_detection_adds_scatter_beyond_the_error_bars():
    # Exponentials comparable to the noise: detected in some realisations only. The value then jumps
    # between two models, so its scatter exceeds the reported (fixed-model) error bar. Recorded in
    # logs/emulation_findings.md; this test pins the behaviour down.
    ratio, fitted = _scatter_and_errors(A)
    marginal = (fitted.mean(axis=0) > 0.05) & (fitted.mean(axis=0) < 0.95)
    assert marginal.any()
    assert np.all(ratio[marginal] > 1.5)


def test_each_window_has_its_own_cadence():
    t = measure.sample_times(0.0, 100.0, [30.0, 60.0], [10.0, 5.0, 2.0])
    assert list(t[t < 30]) == [0.0, 10.0, 20.0]
    assert list(t[(t >= 30) & (t < 60)]) == [30.0, 35.0, 40.0, 45.0, 50.0, 55.0]
    assert np.allclose(np.diff(t[t >= 60]), 2.0) and t[t >= 60][0] == 60.0


def test_glitch_without_enough_samples_is_unmeasured():
    glitches = np.array([1000.0, 1003.0, 3000.0])       # second glitch 3 d after the first, cadence 5 d
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 5000.0, glitches, 5.0, [5e-3, 2e-3, 4e-3], np.full(3, 1e12),
                                    np.zeros(3), np.full(3, 20.0), 0.0, rng)
    for method in measure.METHODS:
        value, _, _ = measure.measure(method, t, y, glitches, np.full(3, 20.0), 4, 4, AIC, 40.0)
        assert np.isnan(value[0]) and np.isfinite(value[2]), method


def test_window_needs_samples_before_but_other_methods_do_not():
    # Second glitch 30 d after the first at a 10-d cadence: 3 samples between them. Window cannot
    # measure the second glitch (3 before < 4) nor the first (3 after < 4); local and global can
    # measure the second (enough after) but not the first.
    glitches = np.array([1000.0, 1030.0])
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 3000.0, glitches, 10.0, [5e-3, 3e-3], np.full(2, 1e12),
                                    np.zeros(2), [20.0, 20.0], 0.0, rng)
    for method, expected in [("window", [False, False]), ("local-relaxation", [False, True]),
                             ("global", [False, True]), ("global-noexp", [False, True])]:
        value, _, _ = measure.measure(method, t, y, glitches, [20.0, 20.0], 4, 4, AIC, None)
        assert list(np.isfinite(value)) == expected, method


def test_unknown_method_is_refused():
    t, y = _series(NEVER, A)
    with pytest.raises(ValueError, match="Unknown measurement method"):
        _measure("unknown", t, y)


# ---- short-window (D4f) -----------------------------------------------------------------------
def test_short_window_known_answer_with_linear_recovery_and_exponential():
    # One glitch at 1000, flat before. Over [1000, 1000 + w) the level averages
    # J (1 - w / (2 t_rec)) + A tau / w (1 - e^(-w/tau)); the level before is 0.
    j, a, tau, t_rec, w = 4e-3, 3e-3, 20.0, 2000.0, 40.0
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 3000.0, [1000.0], 0.001, [j], [t_rec], [a], [tau], 0.0, rng)
    value, _, fitted = measure.measure("short-window", t, y, [1000.0], [tau], 4, 4, AIC, [w])
    expected = j * (1 - w / (2 * t_rec)) + a * tau / w * (1 - np.exp(-w / tau))
    assert value[0] == pytest.approx(expected, rel=1e-4)
    assert not fitted[0]


def test_short_window_spans_are_cut_at_neighbouring_glitches():
    # Glitches at 1000, 1020 and 1100 with permanent steps; w = 40 at a 1-d cadence. The second glitch's
    # "before" span is cut at 1000 (20 samples: the first step only) and its "after" span is a full 40 d,
    # cut at nothing; the first glitch's "after" span is cut at 1020. Permanent steps: every value = J.
    glitches = np.array([1000.0, 1020.0, 1100.0])
    j = np.array([5e-3, 2e-3, 3e-3])
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 2000.0, glitches, 1.0, j, np.full(3, 1e12), np.zeros(3), np.full(3, 20.0),
                                    0.0, rng)
    value, _, _ = measure.measure("short-window", t, y, glitches, np.full(3, 20.0), 4, 4, AIC, np.full(3, 40.0))
    assert value == pytest.approx(j, rel=1e-10)
    # A recovering first step shows the cut: the second glitch's "before" level is the first step's
    # average over [1000, 1020) only, not over [980, 1020).
    t, y, _ = measure.sample_series(0.0, 2000.0, glitches, 1.0, j, [100.0, 1e12, 1e12], np.zeros(3),
                                    np.full(3, 20.0), 0.0, rng)
    value, _, _ = measure.measure("short-window", t, y, glitches, np.full(3, 20.0), 4, 4, AIC, np.full(3, 40.0))
    before_2 = np.mean(5e-3 * (1 - np.arange(0, 20) / 100.0))
    after_2 = np.mean(5e-3 * np.clip(1 - np.arange(20, 60) / 100.0, 0, None)) + 2e-3
    assert value[1] == pytest.approx(after_2 - before_2, rel=1e-10)


def test_short_window_needs_samples_on_both_sides_and_a_span():
    glitches = np.array([1000.0, 1030.0])
    rng = np.random.default_rng(0)
    t, y, _ = measure.sample_series(0.0, 3000.0, glitches, 10.0, [5e-3, 3e-3], np.full(2, 1e12), np.zeros(2),
                                    [20.0, 20.0], 0.0, rng)
    value, _, _ = measure.measure("short-window", t, y, glitches, [20.0, 20.0], 4, 4, AIC, [40.0, 40.0])
    assert list(np.isfinite(value)) == [False, False]       # 3 samples between the glitches
    with pytest.raises(ValueError, match="span"):
        measure.measure("short-window", t, y, glitches, [20.0, 20.0], 4, 4, AIC, None)
