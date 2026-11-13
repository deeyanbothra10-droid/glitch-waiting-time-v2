"""Models M0-M4 on synthetic data with known parameters."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from synth import simulate_aft

from glitchlib import models

INTERCEPTS = {"P0": 7.0, "P1": 8.0, "P2": 9.0}
TRUTH = {
    "M0": {},
    "M1": {},
    "M2": {"b": 0.8},
    "M3": {"a": 0.5},
    "M4": {"a": 0.6, "d": -0.9},
}
SIGMA = 0.5


def _data(model, seed, n_per=40, censor_frac=0.10, coefs=None):
    rng = np.random.default_rng(seed)
    return simulate_aft(model, TRUTH[model] if coefs is None else coefs, INTERCEPTS, SIGMA, n_per,
                        censor_frac, rng)


@pytest.mark.parametrize("model", list(models.MODELS))
def test_recovers_parameters_within_3_se(model):
    data = _data(model, seed=20261031)
    assert 0.03 <= (~data["event"]).mean() <= 0.20          # roughly 10% censored
    result = models.fit(model, data, "ignore")
    truth = dict(TRUTH[model], sigma=SIGMA, **{f"c[{p}]": c for p, c in INTERCEPTS.items()})
    for name, value in truth.items():
        assert abs(result.estimates[name] - value) < 3 * result.se[name], (model, name)


def test_m2_on_m1_data_gives_slope_one():
    result = models.fit("M2", _data("M1", seed=11), "ignore")
    assert abs(result.estimates["b"] - 1.0) < 3 * result.se["b"]


def test_m2_on_m0_data_gives_slope_zero():
    result = models.fit("M2", _data("M0", seed=12), "ignore")
    assert abs(result.estimates["b"]) < 3 * result.se["b"]


def _lifelines_frame(data, covariates):
    frame = pd.get_dummies(data["pulsar"], drop_first=True, dtype=float)
    for col in covariates:
        frame[col] = data[col]
    frame["t"] = data["t"]
    frame["event"] = data["event"].astype(int)
    return frame


@pytest.mark.parametrize("model,covariates", [("M0", []), ("M2", ["ln_R"])])
def test_matches_lifelines(model, covariates):
    lifelines = pytest.importorskip("lifelines")
    data = _data(model, seed=13)
    ours = models.fit(model, data, "ignore")
    theirs = lifelines.LogNormalAFTFitter().fit(_lifelines_frame(data, covariates), "t", "event")
    assert ours.loglik == pytest.approx(theirs.log_likelihood_, abs=1e-4)
    assert ours.estimates["sigma"] == pytest.approx(np.exp(theirs.params_["sigma_"]["Intercept"]), rel=1e-4)
    assert ours.estimates["c[P0]"] == pytest.approx(theirs.params_["mu_"]["Intercept"], abs=1e-4)
    if covariates:
        assert ours.estimates["b"] == pytest.approx(theirs.params_["mu_"]["ln_R"], abs=1e-4)
        assert ours.se["b"] == pytest.approx(theirs.standard_errors_["mu_"]["ln_R"], rel=1e-3)


def test_pointwise_likelihood_known_answer():
    # One completed (T = 100 d) and one censored (c = 300 d) interval, mu = ln 150, sigma = 0.7.
    mu, sigma = np.log(150.0), 0.7
    expected = [stats.lognorm.logpdf(100.0, s=sigma, scale=150.0),
                stats.lognorm.logsf(300.0, s=sigma, scale=150.0)]
    got = models._pointwise(np.array([np.log(sigma)]), np.zeros((2, 0)), np.full(2, mu),
                            np.log([100.0, 300.0]), np.array([True, False]))
    assert got == pytest.approx(expected, rel=1e-12)


def test_gradient_and_hessian_match_finite_differences():
    design = models.build_design("M4", _data("M4", seed=14))
    theta = models.maximise(design)[0] + 0.05
    args = (design.X, design.offset, design.y, design.event)
    _, grad, hess = models._loglik_grad_hess(theta, *args)
    h = 1e-6
    num_grad = [(models._loglik_grad_hess(theta + h * e, *args)[0]
                 - models._loglik_grad_hess(theta - h * e, *args)[0]) / (2 * h) for e in np.eye(len(theta))]
    num_hess = [(models._loglik_grad_hess(theta + h * e, *args)[1]
                 - models._loglik_grad_hess(theta - h * e, *args)[1]) / (2 * h) for e in np.eye(len(theta))]
    assert np.allclose(grad, num_grad, atol=1e-5)
    assert np.allclose(hess, np.array(num_hess), atol=1e-4)


def test_wald_intervals_known_answer():
    result = models.fit("M2", _data("M2", seed=25), "ignore")
    lo, hi = result.intervals(0.95)["b"]
    assert lo == pytest.approx(result.estimates["b"] - 1.959964 * result.se["b"], rel=1e-6)
    assert hi == pytest.approx(result.estimates["b"] + 1.959964 * result.se["b"], rel=1e-6)
    s_lo, s_hi = result.intervals(0.95)["sigma"]
    assert s_lo < result.estimates["sigma"] < s_hi


def test_aic_counts_every_parameter():
    result = models.fit("M4", _data("M4", seed=15), "ignore")
    assert result.n_params == 2 + 3 + 1          # a, d, three intercepts, log_sigma
    assert result.aic == pytest.approx(2 * 6 - 2 * result.loglik)


def test_lr_test_size_under_null():
    # M0 data: the LR test of M2 vs M0 should reject at 5% about 5% of the time.
    p_values = []
    for seed in range(300):
        data = _data("M0", seed=1000 + seed)
        p_values.append(models.lr_test(models.fit("M2", data, "ignore"),
                                       models.fit("M0", data, "ignore"))["p_value"])
    p_values = np.array(p_values)
    assert stats.binomtest(int((p_values < 0.05).sum()), len(p_values), 0.05).pvalue > 0.001
    assert stats.kstest(p_values, "uniform").pvalue > 0.001


def test_lr_test_detects_strong_effect():
    data = _data("M2", seed=16)
    test = models.lr_test(models.fit("M2", data, "ignore"), models.fit("M0", data, "ignore"))
    assert test["df"] == 1 and test["p_value"] < 1e-6


def test_lr_test_refuses_mismatched_fits():
    a, b = _data("M0", seed=17), _data("M0", seed=18)
    with pytest.raises(ValueError, match="different data"):
        models.lr_test(models.fit("M2", a, "ignore"), models.fit("M0", b, "ignore"))
    with pytest.raises(ValueError, match="not nested"):
        models.lr_test(models.fit("M3", a, "ignore"), models.fit("M2", a, "ignore"))


def test_loo_cv_matches_closed_form_without_censoring():
    # M0 without censoring: the ML intercept is the pulsar mean of ln T and sigma^2 the mean squared
    # residual, so each leave-one-out score can be computed by hand.
    data = _data("M0", seed=19, n_per=8, censor_frac=0.0)
    ln_t = np.log(data["t"].to_numpy())
    expected = []
    for i in range(len(data)):
        keep = np.arange(len(data)) != i
        train, pul = ln_t[keep], data["pulsar"].to_numpy()[keep]
        means = {p: train[pul == p].mean() for p in INTERCEPTS}
        sigma = np.sqrt(np.mean((train - np.array([means[p] for p in pul])) ** 2))
        mu = means[data["pulsar"].iloc[i]]
        expected.append(stats.norm.logpdf(ln_t[i], mu, sigma) - ln_t[i])
    result = models.loo_cv("M0", data, "ignore")
    assert result["scores"] == pytest.approx(expected, rel=1e-6)
    assert result["mean_log_score"] == pytest.approx(np.mean(expected), rel=1e-6)


def test_loo_cv_prefers_true_model():
    data = _data("M1", seed=20, n_per=15)
    assert (models.loo_cv("M1", data, "ignore")["mean_log_score"]
            > models.loo_cv("M0", data, "ignore")["mean_log_score"])


def test_loo_cv_refuses_single_interval_pulsar():
    data = _data("M0", seed=21, n_per=5)
    data = pd.concat([data, data.iloc[[0]].assign(pulsar="lonely")], ignore_index=True)
    with pytest.raises(ValueError, match="single interval"):
        models.loo_cv("M0", data, "ignore")


def test_unknown_measurement_error_names_d6():
    with pytest.raises(ValueError, match="D6"):
        models.fit("M0", _data("M0", seed=22), "deming")


def test_pulsar_with_only_censored_intervals_is_refused():
    data = _data("M0", seed=23)
    data.loc[data["pulsar"] == "P2", "event"] = False
    with pytest.raises(models.FitError, match="P2"):
        models.fit("M0", data, "ignore")


def test_fit_does_not_modify_data():
    data = _data("M4", seed=24)
    before = data.copy()
    models.fit("M4", data, "ignore")
    pd.testing.assert_frame_equal(data, before)


def test_poisson_per_pulsar_known_answer():
    data = pd.DataFrame({"pulsar": ["A", "A", "A"], "t": [100.0, 300.0, 200.0],
                         "event": [True, True, False]})
    row = models.fit_poisson_per_pulsar(data).iloc[0]
    assert row["rate_per_day"] == pytest.approx(2 / 600)
    assert row["loglik"] == pytest.approx(2 * np.log(2 / 600) - 2)
    assert row["se_log_rate"] == pytest.approx(1 / np.sqrt(2))


def _time_limit(seconds):
    import signal

    def handler(signum, frame):
        raise TimeoutError(f"fit did not finish within {seconds} s")
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(seconds)


def test_fit_that_hung_in_scipy_trust_exact_now_converges():
    # The M3 fit of replicate 370 of S4w at sigma 0.2 (06_validate, 31 Oct) looped for hours inside
    # scipy's trust-exact (NaN in its subproblem). The capped Newton method must finish quickly and reach
    # the same maximum as an independent fit (lifelines).
    import pathlib
    import signal
    lifelines = pytest.importorskip("lifelines")
    data = pd.read_csv(pathlib.Path(__file__).parent / "data" / "hang_s4w_sigma0.2_rep370_M3.csv")
    _time_limit(10)
    try:
        ours = models.fit("M3", data, "ignore")
    finally:
        signal.alarm(0)
    frame = _lifelines_frame(data, ["ln_dnu_frac"])
    theirs = lifelines.LogNormalAFTFitter().fit(frame, "t", "event")
    assert ours.loglik == pytest.approx(theirs.log_likelihood_, abs=1e-3)


def test_non_finite_likelihood_is_a_fit_failure_not_a_loop():
    import signal
    data = _data("M3", seed=60)
    data["ln_dnu_frac"] = data["ln_dnu_frac"] * 1e154          # exp overflows: non-finite from the start
    _time_limit(10)
    try:
        with pytest.raises(models.FitError, match="non-finite"):
            models.fit("M3", data, "ignore")
    finally:
        signal.alarm(0)
