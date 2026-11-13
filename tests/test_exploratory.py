"""Exploratory (post hoc) helpers on synthetic data with known answers."""
import numpy as np
import pandas as pd
import pytest
from synth import simulate_aft

from glitchlib import exploratory
from glitchlib.config import load_config


def test_short_window_shorter_than_span_known_answer():
    cfg = load_config()
    data = pd.DataFrame({
        "source_key": ["elsk11", "elsk11", "wmp+00", "wmp+00", "aeka18", "elsk11"],
        "method_class": ["short-window", "short-window", "short-window", "short-window", "local-relaxation",
                         "short-window"],
        "event": [True, True, True, True, True, False],
        "t": [39.0, 40.0, 124.0, 300.0, 10.0, 5.0]})
    # elsk11 w = 40 d, wmp+00 w = 125 d (D4f); not short-window or censored: never flagged
    assert exploratory.short_window_shorter_than_span(data, cfg).tolist() == [True, False, True, False, False, False]


def test_per_pulsar_cases_selects_pulsars_and_reports_method_mix():
    cfg = load_config()
    rng = np.random.default_rng(5)
    data = simulate_aft("M4", {"a": 1.0, "d": -1.0}, {"P0": 7.0, "P1": 7.5, "P2": 8.0}, 0.4, 30, 0.05, rng)
    data = data[~((data["pulsar"] == "P2") & (data.groupby("pulsar").cumcount() >= 5))].copy()   # P2: 5 intervals
    data["method_class"] = np.where(np.arange(len(data)) % 3 == 0, "global", "short-window")
    data["prev_t"], data["prev_event"] = np.nan, False
    out = exploratory.per_pulsar_cases(data, cfg, 6, {"b": 1.0, "d": 1.0}, {"b": 1.0, "d": -1.0})
    assert [r["pulsar"] for r in out] == ["P0", "P1"]                                  # P2 has 5 < 6
    p0 = data[data["pulsar"] == "P0"]
    assert out[0]["method_mix"] == p0["method_class"].value_counts().sort_index().to_dict()
    assert out[0]["fit"]["status"] == "run"
    assert abs(out[0]["fit"]["estimates_D4p"]["d"]["estimate"] + 1) < 0.5            # true d = -1


def test_per_pulsar_cases_records_a_fit_failure():
    cfg = load_config()
    rng = np.random.default_rng(6)
    data = simulate_aft("M4", {"a": 1.0, "d": -1.0}, {"P0": 7.0}, 0.4, 10, 0.0, rng)
    data["ln_dnudot_frac"] = data["ln_dnudot_frac"] * 1e154                           # non-finite: a fit failure
    data["method_class"] = "global"
    out = exploratory.per_pulsar_cases(data, cfg, 6, {"b": 1.0, "d": 1.0}, {"b": 1.0, "d": -1.0})
    assert out[0]["fit"]["status"].startswith("fit failed")


def _two_slope_data(d_values, seed, n_per=25):
    rng = np.random.default_rng(seed)
    parts = []
    for k, d in enumerate(d_values):
        part = simulate_aft("M4", {"a": 1.0, "d": d}, {f"P{k}": 7.0 + 0.3 * k}, 0.3, n_per, 0.0, rng)
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def test_leave_one_out_d_matches_a_direct_refit():
    from glitchlib import models
    cfg = load_config()
    data = _two_slope_data([-1.0, -1.0], 7, n_per=12)
    rows = exploratory.leave_one_out_d(data, cfg)
    assert len(rows) == len(data)
    direct = models.fit("M4", data.drop(index=data.index[5]), "ignore").estimates["d"]
    assert rows[5]["d_without"] == pytest.approx(direct, rel=1e-8) and rows[5]["pulsars_dropped"] == []


def test_separate_d_test_known_answers():
    cfg = load_config()
    same = exploratory.separate_d_test(_two_slope_data([-1.0, -1.0, -1.0], 8), cfg, 6)
    differ = exploratory.separate_d_test(_two_slope_data([-2.0, 0.0, 2.0], 9), cfg, 6)
    assert same["df"] == 2 and differ["df"] == 2 and same["statistic"] >= -1e-8
    assert differ["p_value"] < 1e-6 and differ["statistic"] > same["statistic"]
    assert differ["separate_d"]["P0"]["estimate"] == pytest.approx(-2.0, abs=0.5)
    assert differ["separate_d"]["P2"]["estimate"] == pytest.approx(2.0, abs=0.5)
    one = exploratory.separate_d_test(_two_slope_data([-1.0], 10), cfg, 6)          # one pulsar: identical models
    assert one["df"] == 0 and one["statistic"] == pytest.approx(0.0, abs=1e-6)


def test_single_pulsar_replicates_are_reproducible_and_fit_m4_alone():
    from test_simulate import S1M, SETTINGS, _templates
    cfg = load_config()
    template = _templates(80)[0]
    a = exploratory.run_single_pulsar(template, S1M, "resampled", SETTINGS, 0.4, 1e-3, 4, (1, 2, 3), cfg, 2)
    b = exploratory.run_single_pulsar(template, S1M, "resampled", SETTINGS, 0.4, 1e-3, 4, (1, 2, 3), cfg, 2)
    pd.testing.assert_frame_equal(a, b)
    assert list(a["rep"]) == [0, 1, 2, 3]
    assert all(s in ("ok", "discarded") or s.startswith("failed") for s in a["status"])
    ok = a[a["status"] == "ok"]
    assert len(ok) and (ok["sigma_hat"] > 0).all()
    # replicate 2 directly: the same seed gives the same M4 fit
    one = exploratory._single_pulsar_replicate((template, S1M, 0.4, "resampled", SETTINGS, 1e-3, (1, 2, 3, 2), cfg))
    assert one == {k: v for k, v in a.iloc[2].to_dict().items() if not (isinstance(v, float) and v != v)}


def test_single_pulsar_replicate_discarded_when_too_few_completed():
    import dataclasses
    from test_simulate import S1M, SETTINGS, _templates
    cfg = load_config()
    template = dataclasses.replace(_templates(81)[0], admissible_fraction=0.01)
    out = exploratory._single_pulsar_replicate((template, S1M, 0.4, "resampled", SETTINGS, 1e-3, (9, 0), cfg))
    assert out["status"] == "discarded"


def _per_pulsar_data(a_values, d_values, seed, n_per=30):
    rng = np.random.default_rng(seed)
    return pd.concat([simulate_aft("M4", {"a": a, "d": d}, {f"P{k}": 7.0 + 0.3 * k}, 0.3, n_per, 0.0, rng)
                      for k, (a, d) in enumerate(zip(a_values, d_values))], ignore_index=True)


def test_separate_d_given_separate_a_known_answers():
    cfg = load_config()
    common_d = _per_pulsar_data([0.2, 1.0, 1.8], [-1.0, -1.0, -1.0], 11)     # a differs, d common
    t = exploratory.per_pulsar_lr_test(common_d, cfg, 6, ("a",), ("a", "d"))
    assert t["df"] == 2 and t["p_value"] > 0.001
    assert t["reduced"]["estimates"]["d"] == pytest.approx(-1.0, abs=0.3)
    assert t["reduced"]["estimates"]["a[P0]"] == pytest.approx(0.2, abs=0.4)
    assert t["reduced"]["estimates"]["a[P2]"] == pytest.approx(1.8, abs=0.4)
    differ = exploratory.per_pulsar_lr_test(_per_pulsar_data([1.0, 1.0, 1.0], [-2.0, 0.0, 2.0], 12), cfg, 6,
                                            ("a",), ("a", "d"))
    assert differ["p_value"] < 1e-6 and differ["full"]["estimates"]["d[P2]"] == pytest.approx(2.0, abs=0.5)
    # the per-pulsar fits use the same likelihood: with every coefficient separate on one pulsar, nothing changes
    one = exploratory.per_pulsar_lr_test(_per_pulsar_data([1.0], [-1.0], 13), cfg, 6, (), ("a", "d"))
    assert one["df"] == 0 and one["statistic"] == pytest.approx(0.0, abs=1e-6)
