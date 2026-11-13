"""The real-run analyses on synthetic data with known answers."""
import numpy as np
import pandas as pd
import pytest
from synth import simulate_aft

from glitchlib import analysis
from glitchlib.config import load_config

INTERCEPTS = {f"P{i}": 6.5 + 0.2 * i for i in range(6)}


def _data(model, coefs, seed, n_per=30):
    data = simulate_aft(model, coefs, INTERCEPTS, 0.4, n_per, 0.05, np.random.default_rng(seed))
    data["prev_t"] = np.r_[np.nan, data["t"].to_numpy()[:-1]]
    data.loc[data["pulsar"] != data["pulsar"].shift(), "prev_t"] = np.nan
    data["prev_event"] = data["prev_t"].notna()
    return data


def test_primary_rejects_under_the_model_and_reports_the_right_numbers():
    cfg = load_config()
    out = analysis.primary(_data("M4", {"a": 1.0, "d": -1.0}, 1), cfg)
    assert out["T1"]["rejects_chi2"] and out["T2"]["rejects_chi2"]
    assert out["T1"]["signed_root"] == pytest.approx(np.sqrt(out["T1"]["statistic"]))     # b > 0: positive
    assert out["T2"]["signed_root"] == pytest.approx(np.sqrt(out["T2"]["statistic"]))     # d < 0: positive
    assert out["T1"]["parameter"] == "b" and out["T2"]["parameter"] == "d"
    assert abs(out["T2"]["estimate"] + 1) < 3 * out["T2"]["se"]
    lo, hi = out["T2"]["interval"]
    assert lo < out["T2"]["estimate"] < hi and out["T1"]["df"] == 1 and out["per_test_alpha"] == 0.025
    assert out["gate"]["met"] and out["n_pulsars"] == 6


def test_primary_mostly_accepts_under_the_null():
    cfg = load_config()
    rejects = [analysis.primary(_data("M0", {}, seed), cfg)["T1"]["rejects_chi2"] for seed in range(100, 160)]
    assert np.mean(rejects) < 0.12


def test_gate_known_answer():
    cfg = load_config()
    data = pd.DataFrame({"pulsar": ["A"] * 20 + ["B"] * 20, "event": [True] * 39 + [False]})
    assert analysis.gate(data, cfg) == {"admissible_completed": 39, "pulsars": 2, "met": False}


def test_secondary_contains_every_analysis():
    cfg = load_config()
    data = _data("M1", {}, 2)
    tau = {p: 1e6 for p in INTERCEPTS}
    out = analysis.secondary(data, cfg, tau, 19, 5, True)
    assert set(out) >= {"implied_beta", "spearman_forward", "spearman_placebo", "placebo", "poisson",
                        "permutation", "placebo_permutation", "loo"}
    assert out["permutation"]["p_value"] == pytest.approx(1 / 20)
    assert out["loo"]["M1"] > out["loo"]["M0"]
    assert "significantly_negative" in out["placebo"]
    flat = analysis.summarise_secondary(out)
    assert "b_perm" not in flat["permutation"]


def test_implied_beta_check_uses_pulsars_with_enough_intervals():
    cfg = load_config()
    data = _data("M1", {}, 3)
    check = analysis.implied_beta_check(data, cfg, {p: 1e9 for p in INTERCEPTS})   # huge tau: beta > 0
    assert check["all_positive"] and check["n_checked"] == 6
    check = analysis.implied_beta_check(data, cfg, {p: 1.0 for p in INTERCEPTS})   # tiny tau: beta < 0
    assert not check["all_positive"]


def test_kaplan_meier_known_answer():
    # Times 1, 2, 2, 3+ (censored), 4: S(1) = 4/5, S(2) = 4/5 * 2/4 = 0.4, S(4) = 0.4 * 0/1 = 0.
    times, surv = analysis.kaplan_meier([1, 2, 2, 3, 4], [True, True, True, False, True])
    assert list(times) == [1, 2, 4] and surv == pytest.approx([0.8, 0.4, 0.0])


def test_d4k_5_robustness_case_reports_no_data_and_lets_fit_failures_through():
    from glitchlib import models
    cfg = load_config()
    data = _data("M2", {"b": 1.0}, 4)
    f, prim = {"b": 1.0, "d": 1.0}, {"b": 1.0, "d": -1.0}
    out = analysis.robustness_case(data[data["pulsar"] == "nobody"], cfg, "empty", None, f, prim)
    assert out["status"].startswith("no data")
    few = data.groupby("pulsar").head(3)                                        # 3 completed per pulsar < 4
    assert analysis.robustness_case(few, cfg, "too few", None, f, prim)["status"].startswith("no data")
    ok = analysis.robustness_case(data, cfg, "all", 40, f, prim)
    assert ok["status"] == "run" and ok["exploratory"] is False
    bad = data.assign(ln_R=data["ln_R"] * 1e154)                                # non-finite: a fit failure
    with pytest.raises(models.FitError):
        analysis.robustness_case(bad, cfg, "broken", None, f, prim)


def test_d4p_robustness_estimates_known_answer():
    cfg = load_config()
    result = {"T1": {"parameter": "b", "estimate": 0.5, "se": 0.1, "interval": [0.3, 0.7], "signed_root": 4.0},
              "T2": {"parameter": "d", "estimate": 0.2, "se": 0.05, "interval": [0.1, 0.3], "signed_root": -3.0}}
    out = analysis.robustness_estimates(result, cfg, {"b": 2.0, "d": 1.0}, {"b": 0.8, "d": -0.4})
    assert out["b"]["calibrated_interval"] == pytest.approx([0.5 - 1.959964 * 0.2, 0.5 + 1.959964 * 0.2])
    assert out["d"]["calibrated_interval"] == pytest.approx([0.2 - 1.959964 * 0.05, 0.2 + 1.959964 * 0.05])
    assert out["b"]["uncalibrated"] == [0.3, 0.7] and out["d"]["signed_root"] == -3.0
    assert out["b"]["calibrated_interval_label"] == "approximate (primary factors)"            # D4p clarification
    assert out["b"]["same_sign_as_primary"] and not out["d"]["same_sign_as_primary"]      # d: + here, - in primary


def test_signed_root_follows_the_predicted_direction():
    assert analysis.signed_root(4.0, +0.3, "b_positive") == 2.0
    assert analysis.signed_root(4.0, -0.3, "b_positive") == -2.0
    assert analysis.signed_root(4.0, -0.3, "d_negative") == 2.0
    assert analysis.signed_root(4.0, +0.3, "d_negative") == -2.0


def test_d4n_d_without_short_window_drops_those_intervals_and_refits_m4():
    from glitchlib import models
    cfg = load_config()
    kept = _data("M4", {"a": 1.0, "d": -1.0}, 5).assign(method_class="local-relaxation")
    sw = _data("M4", {"a": 1.0, "d": +1.0}, 6).assign(method_class="short-window")       # opposite d
    both = pd.concat([kept, sw], ignore_index=True)
    assert models.fit("M4", both, "ignore").estimates["d"] > -0.5                          # contaminated
    d = analysis.d_without_short_window(both, cfg)
    assert d == pytest.approx(models.fit("M4", kept, "ignore").estimates["d"], rel=1e-9)
    assert abs(d + 1) < 0.2
    assert analysis.d_without_short_window(sw, cfg) is None                                 # nothing left


def test_d4n_d_without_short_window_applies_the_d13_count_to_the_subset():
    from glitchlib import models
    cfg = load_config()
    data = _data("M4", {"a": 1.0, "d": -1.0}, 7).assign(method_class="global")
    p0 = data["pulsar"] == "P0"
    data.loc[p0, "method_class"] = "short-window"
    data.loc[p0 & (data.groupby("pulsar").cumcount() < 3), "method_class"] = "global"    # P0: 3 left, < 4
    expected = models.fit("M4", data[data["pulsar"] != "P0"], "ignore").estimates["d"]
    assert analysis.d_without_short_window(data, cfg) == pytest.approx(expected, rel=1e-9)
