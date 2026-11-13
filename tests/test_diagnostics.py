"""Diagnostics on synthetic data with known answers."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from synth import simulate_aft

from glitchlib import diagnostics, models

INTERCEPTS = {"P0": 7.0, "P1": 8.0, "P2": 9.0}


def test_spearman_within_known_answer():
    # A: d = [-1, 1, -1, 1, 0] -> r = 1 - 6*4/(5*24) = 0.8. B: sum d^2 = 68 -> r = 1 - 6*68/(6*35).
    # A's censored row and C (only 3 completed) must be ignored.
    data = pd.DataFrame({
        "pulsar": ["A"] * 6 + ["B"] * 6 + ["C"] * 4,
        "event": [True] * 5 + [False] + [True] * 6 + [True, True, True, False],
        "x": [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4],
        "y": [2, 1, 4, 3, 5, -9, 6, 5, 4, 3, 1, 2, 1, 2, 3, 4],
    })
    result = diagnostics.spearman_within(data, "x", "y", min_completed=4, label="completed intervals only",
                                         level=0.95)
    per = result["per_pulsar"].set_index("pulsar")["r"]
    r_b = 1 - 6 * 68 / (6 * 35)
    assert per["A"] == pytest.approx(0.8) and per["B"] == pytest.approx(r_b)
    assert "C" not in per.index
    z_bar = (2 * np.arctanh(0.8) + 3 * np.arctanh(r_b)) / 5
    assert result["pooled_r"] == pytest.approx(np.tanh(z_bar))
    half = 1.959963984540054 / np.sqrt(5)                     # Fisher z +- 1.96 / sqrt(sum(n_p - 3))
    assert result["interval"] == pytest.approx([np.tanh(z_bar - half), np.tanh(z_bar + half)])
    assert result["label"] == "completed intervals only" and result["n_intervals_pooled"] == 11


def test_spearman_abs_r_one_is_flagged_and_excluded_from_pooled():
    data = pd.DataFrame({"pulsar": ["A"] * 5 + ["B"] * 5, "event": True,
                         "x": [1, 2, 3, 4, 5] * 2, "y": [2, 1, 4, 3, 5, 1, 2, 3, 4, 5]})
    result = diagnostics.spearman_within(data, "x", "y", min_completed=4, label="completed intervals only",
                                         level=0.95)
    assert result["excluded_abs_r_one"] == ["B"]
    assert bool(result["per_pulsar"].set_index("pulsar").loc["B", "abs_r_one"])
    assert result["pooled_r"] == pytest.approx(0.8) and result["n_pulsars_pooled"] == 1


def test_permute_within_keeps_each_pulsars_values():
    rng = np.random.default_rng(1)
    values = np.arange(10.0)
    groups = [np.arange(0, 4), np.arange(4, 10)]
    out = diagnostics._permute_within(values, groups, rng)
    assert sorted(out[:4]) == [0, 1, 2, 3] and sorted(out[4:]) == [4, 5, 6, 7, 8, 9]
    assert not np.array_equal(out, values)


def test_permutation_keeps_t_event_and_pulsar_fixed():
    data = simulate_aft("M0", {}, INTERCEPTS, 0.5, 10, 0.1, np.random.default_rng(2))
    before = data.copy()
    diagnostics.permutation_test(data, "all", 9, 1, "ignore")
    pd.testing.assert_frame_equal(data, before)


def test_permutation_strong_effect_gives_smallest_p():
    data = simulate_aft("M2", {"b": 1.0}, INTERCEPTS, 0.3, 30, 0.1, np.random.default_rng(3))
    result = diagnostics.permutation_test(data, "all", 199, 20261031, "ignore")
    assert result["p_value"] == pytest.approx(1 / 200)
    assert result["b_obs"] == pytest.approx(models.fit("M2", data, "ignore").estimates["b"])


def test_permutation_is_reproducible_with_seed():
    data = simulate_aft("M0", {}, INTERCEPTS, 0.5, 10, 0.1, np.random.default_rng(4))
    a = diagnostics.permutation_test(data, "all", 19, 7, "ignore")
    b = diagnostics.permutation_test(data, "all", 19, 7, "ignore")
    assert np.array_equal(a["b_perm"], b["b_perm"]) and a["p_value"] == b["p_value"]


def test_permutation_refuses_other_interval_handling():
    data = simulate_aft("M0", {}, INTERCEPTS, 0.5, 10, 0.1, np.random.default_rng(5))
    with pytest.raises(ValueError, match="D4a"):
        diagnostics.permutation_test(data, "completed_only", 9, 1, "ignore")


def test_permutation_p_values_uniform_under_null():
    # R irrelevant (M0 data, 10% censored). With N = 49, P(p <= k/50) = k/50 exactly under the null.
    n_perm, p_values = 49, []
    for seed in range(150):
        data = simulate_aft("M0", {}, INTERCEPTS, 0.5, 12, 0.1, np.random.default_rng(5000 + seed))
        p_values.append(diagnostics.permutation_test(data, "all", n_perm, seed, "ignore")["p_value"])
    p_values = np.array(p_values)
    for level in (0.1, 0.5):
        assert stats.binomtest(int((p_values <= level).sum()), len(p_values), level).pvalue > 0.001
    assert stats.kstest(p_values, "uniform").pvalue > 0.001


def test_placebo_data_uses_previous_wait():
    data = pd.DataFrame({"pulsar": ["A", "A", "A"], "t": [10.0, 20.0, 30.0], "event": [True, True, False],
                         "prev_t": [np.nan, 10.0, 20.0], "prev_event": [np.nan, True, True]})
    out = diagnostics.placebo_data(data)
    assert list(out["t"]) == [10.0, 20.0] and list(out["event"]) == [True, True]


def test_placebo_fit_recovers_negative_backward_slope():
    rng = np.random.default_rng(6)
    data = simulate_aft("M0", {}, INTERCEPTS, 0.5, 40, 0.1, rng)
    c = data["pulsar"].map(INTERCEPTS).to_numpy()
    data["prev_t"] = np.exp(c - 0.7 * data["ln_R"] + 0.5 * rng.normal(size=len(data)))
    data["prev_event"] = True
    result = diagnostics.placebo_fit(data, "ignore")
    assert abs(result.estimates["b"] + 0.7) < 3 * result.se["b"]


def test_conditional_logit_single_stratum_hand_worked():
    # x = 0, 1, 2 with the middle glitch missing: L = e^b / (1 + e^b + e^2b), maximised at b = 0.
    # Information there = variance of x under equal weights on {0, 1, 2} = 2/3, so SE = sqrt(3/2).
    beta, se, loglik = diagnostics.conditional_logit(np.array(["A"] * 3), [False, True, False], [0, 1, 2])
    assert beta == pytest.approx(0.0, abs=1e-10)
    assert se == pytest.approx(np.sqrt(1.5))
    assert loglik == pytest.approx(-np.log(3))


def test_conditional_logit_two_strata_hand_worked():
    # Add stratum B: x = 0, 1 with the second missing, L_B = e^b / (1 + e^b). With u = e^b the score
    # (1 - u^2)/(1 + u + u^2) + 1/(1 + u) = 0 reduces to u^3 - 2u - 2 = 0 (Cardano's formula).
    q = np.sqrt(1 - 8 / 27)
    u = np.cbrt(1 + q) + np.cbrt(1 - q)
    info = ((u + 4 * u ** 2) / (1 + u + u ** 2) - ((u + 2 * u ** 2) / (1 + u + u ** 2)) ** 2
            + u / (1 + u) ** 2)
    strata = np.array(["A", "A", "A", "B", "B"])
    beta, se, _ = diagnostics.conditional_logit(strata, [False, True, False, False, True], [0, 1, 2, 0, 1])
    assert beta == pytest.approx(np.log(u), rel=1e-9)
    assert se == pytest.approx(1 / np.sqrt(info), rel=1e-9)


def test_stratum_recursion_matches_brute_force_enumeration():
    from itertools import combinations
    x = np.array([0.3, -1.2, 2.0, 0.7, -0.4, 1.1])
    y = np.array([True, False, True, False, False, True])
    beta = 0.37
    xc = x - x.mean()
    sums = np.array([xc[list(s)].sum() for s in combinations(range(6), 3)])
    weights = np.exp(beta * sums)
    mean = np.sum(weights * sums) / weights.sum()
    expected = (beta * xc[y].sum() - np.log(weights.sum()), xc[y].sum() - mean,
                np.sum(weights * sums ** 2) / weights.sum() - mean ** 2)
    assert diagnostics._stratum_terms(beta, x, y) == pytest.approx(expected, rel=1e-12)


def test_missingness_recovers_known_odds_ratio():
    rng = np.random.default_rng(7)
    rows = []
    for p, alpha in enumerate([-4.0, -5.0, -3.5, -4.5, -6.0]):
        ln_w = rng.normal(6.0, 1.0, 150)
        missing = rng.random(150) < 1 / (1 + np.exp(-(alpha + 0.8 * ln_w)))
        rows.append(pd.DataFrame({"pulsar": f"P{p}", "ln_wait": ln_w, "missing": missing, "done": True}))
    rows.append(pd.DataFrame({"pulsar": "never", "ln_wait": rng.normal(6, 1, 20), "missing": False,
                              "done": True}))
    rows.append(pd.DataFrame({"pulsar": "P0", "ln_wait": [9.0, 9.5], "missing": [True, False],
                              "done": False}))        # censored following waits: excluded and counted
    result = diagnostics.missingness_check(pd.concat(rows), "missing", "ln_wait", "done", level=0.95)
    assert abs(result["log_odds_ratio"] - 0.8) < 3 * result["se_log_odds_ratio"]
    assert result["pulsars_without_variation"] == ["never"]
    assert result["n_excluded_censored_following_wait"] == 2
    assert result["odds_ratio"] == pytest.approx(np.exp(result["log_odds_ratio"]))
    lo, hi = result["interval"]
    assert lo == pytest.approx(np.exp(result["log_odds_ratio"] - 1.959964 * result["se_log_odds_ratio"]))
    assert lo < result["odds_ratio"] < hi


def _m1_fit(c_values, se_c):
    names = [f"c[{p}]" for p in c_values] + ["log_sigma"]
    theta = np.array(list(c_values.values()) + [0.0])
    cov = np.diag([se_c ** 2] * len(c_values) + [0.01])
    return models.FitResult("M1", "ignore", names, theta, cov, 0.0, 10, 10, list(c_values), "x")


def test_implied_beta_known_answer():
    # beta = 2e-3 tau / exp(c) - 0.5. Choose c so beta = 1 for tau = 1e6 d and beta = 0 for tau = 5e5 d.
    tau = {"A": 1e6, "B": 5e5}
    c = {"A": np.log(2e-3 * 1e6 / 1.5), "B": np.log(2e-3 * 5e5 / 0.5)}
    result = diagnostics.implied_beta(_m1_fit(c, 0.1), tau).set_index("pulsar")
    assert result.loc["A", "beta"] == pytest.approx(1.0) and result.loc["B", "beta"] == pytest.approx(0.0)
    assert result.loc["A", "se_beta"] == pytest.approx(1.5 * 0.1)


def test_implied_beta_refuses_other_models():
    fit = _m1_fit({"A": 1.0}, 0.1)
    fit.model = "M2"
    with pytest.raises(ValueError, match="M1"):
        diagnostics.implied_beta(fit, {"A": 1e6})


def test_collinearity_known_answer():
    rng = np.random.default_rng(8)
    x = rng.normal(size=500)
    data = pd.concat([
        pd.DataFrame({"pulsar": "A", "ln_dnu_frac": [1.0, 2, 3, 4], "ln_dnudot_frac": [3.0, 5, 7, 9],
                      "ln_R": [1.0, 2, 3, 2]}),
        pd.DataFrame({"pulsar": "B", "ln_dnu_frac": x, "ln_dnudot_frac": rng.normal(size=500),
                      "ln_R": rng.normal(size=500)}),
        pd.DataFrame({"pulsar": "C", "ln_dnu_frac": [1.0, 2], "ln_dnudot_frac": [1.0, 2], "ln_R": [0.0, 1]}),
    ])
    result = diagnostics.collinearity_summary(data, method="pearson", min_n=3)
    per = result["per_pulsar"].set_index("pulsar")
    assert per.loc["A", "r"] == pytest.approx(1.0) and abs(per.loc["B", "r"]) < 0.15
    assert "C" not in per.index
    assert per.loc["A", "sd_ln_R"] == pytest.approx(np.std([1, 2, 3, 2], ddof=1))
    assert result["median_abs_r"] == pytest.approx((1.0 + abs(per.loc["B", "r"])) / 2)
    with pytest.raises(ValueError):
        diagnostics.collinearity_summary(data, method="kendall", min_n=3)


def test_pooled_within_pulsar_correlation_removes_pulsar_offsets():
    # Within each pulsar y = x exactly, but pulsar offsets run the other way, so the raw pooled
    # correlation is negative while the demeaned one is exactly 1.
    data = pd.DataFrame({"pulsar": ["A"] * 4 + ["B"] * 4,
                         "ln_dnu_frac": [0.0, 1, 2, 3, 10, 11, 12, 13],
                         "ln_dnudot_frac": [10.0, 11, 12, 13, 0, 1, 2, 3], "ln_R": 0.0})
    result = diagnostics.collinearity_summary(data, method="pearson", min_n=4)
    assert stats.pearsonr(data["ln_dnu_frac"], data["ln_dnudot_frac"]).statistic < 0
    assert result["pooled_within_pulsar_r"] == pytest.approx(1.0)


def test_collinearity_refuses_unknown_method():
    data = pd.DataFrame({"pulsar": "A", "ln_dnu_frac": [1.0, 2, 3], "ln_dnudot_frac": [1.0, 3, 2],
                         "ln_R": 0.0})
    with pytest.raises(ValueError):
        diagnostics.collinearity_summary(data, method="kendall", min_n=3)


@pytest.mark.parametrize("missing", [[True, True, False, False, False, False], [False, False, False, False, True, True]])
def test_conditional_logit_refuses_complete_separation(missing):
    # The missing glitches have the two shortest (or longest) waits: the estimate is infinite.
    strata = np.array(["A"] * 6 + ["B"] * 4)
    x = np.r_[np.arange(1.0, 7.0), np.arange(1.0, 5.0)]
    y = np.r_[missing, [True, False, False, False] if missing[0] else [False, False, False, True]]
    with pytest.raises(models.FitError, match="separation"):
        diagnostics.conditional_logit(strata, y, x)
    # One stratum that is not separated makes the estimate finite again.
    y_mixed = np.r_[missing, [False, True, False, False]]
    beta, se, _ = diagnostics.conditional_logit(strata, y_mixed, x)
    assert np.isfinite(beta) and 0 < se < np.inf


def test_missingness_check_reports_separation_as_a_fit_failure():
    data = pd.DataFrame({"pulsar": ["A"] * 6, "ln_wait": np.arange(1.0, 7.0),
                         "missing": [True, True, False, False, False, False], "done": True})
    with pytest.raises(models.FitError):
        diagnostics.missingness_check(data, "missing", "ln_wait", "done", 0.95)
