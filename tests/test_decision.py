"""D4j order with the D4l tests and bounds: one test per branch, known answers."""
import numpy as np
import pytest

from glitchlib import decision
from glitchlib.config import load_config


def _refs():
    s1m = np.linspace(0.0, 2.0, 1001)          # 2.5%: 0.05; 97.5%: 1.95
    placebo = np.linspace(-1.0, 1.0, 1001)     # 2.5%: -0.95
    s4 = np.linspace(-0.2, 0.2, 1001)          # 97.5%: 0.19
    s1 = np.linspace(0.0, 1.0, 1001)           # 97.5%: 0.975
    return decision.references(s1m, placebo, s4, s1, load_config())


def _out(t1, t2, b, placebo=0.0, t2_power=0.34):
    return decision.classify(t1, t2, b, placebo, _refs(), t2_power, load_config())


def test_reference_percentiles():
    refs = _refs()
    assert refs["s1m_b_lower"] == pytest.approx(0.05) and refs["s1m_b_refutation_bound"] == pytest.approx(0.05)
    assert refs["s1m_b_label_bound"] == pytest.approx(1.95) and refs["s1_b_upper"] == pytest.approx(0.975)
    assert refs["s4mix_b_support_bound"] == pytest.approx(0.19) and refs["s1m_placebo_bound"] == pytest.approx(-0.95)
    assert "s1m_b" not in decision.references_for_json(refs)


def test_null_threshold_is_the_97_5th_percentile():
    assert decision.null_threshold(np.linspace(-3, 3, 1001), load_config()) == pytest.approx(2.85)


def test_placebo_artefact_first():
    assert _out(True, True, 0.8, placebo=-0.99)["outcome"] == "D_placebo"


def test_c_only_by_low_slope_no_power_clause():
    assert _out(True, True, 0.02)["outcome"] == "C"
    assert _out(False, False, 0.5)["outcome"] == "inconclusive"                  # no power clause any more


def test_a2_b_a1():
    assert _out(True, True, 0.8)["outcome"] == "A2"
    assert _out(True, False, 0.8, t2_power=0.85)["outcome"] == "B"
    out = _out(True, False, 0.8)
    assert out["outcome"] == "A1" and out["T2_power"] == 0.34                  # D4p: the power itself, no fixed label


def test_a_upper_bound_is_s1s_97_5th_percentile():
    assert _out(True, True, 0.97)["outcome"] == "A2"
    assert _out(True, True, 0.98)["outcome"] == "inconclusive"                  # above S1's 97.5th percentile


def test_signal_within_bias_range():
    assert _out(True, True, 0.15)["outcome"] == "D_bias_range"
    assert _out(True, True, 0.19)["outcome"] == "D_bias_range"


def test_report_only_percentile_and_attenuation_label():
    out = _out(True, True, 0.8)
    assert out["b_percentile_in_S1m"] == pytest.approx(100 * 401 / 1001) and out["attenuation_label"] is None
    refs = decision.references(np.linspace(0.0, 0.5, 1001), np.linspace(-1, 1, 1001), np.linspace(-0.2, 0.2, 1001),
                               np.linspace(0.0, 1.0, 1001), load_config())
    labelled = decision.classify(True, True, 0.8, 0.0, refs, 0.34, load_config())
    assert labelled["outcome"] == "A2" and labelled["attenuation_label"] == "less attenuated than emulated"


def test_rates_are_exclusive():
    r = decision.rates(["A1", "A2", "C", "D_placebo", "D_bias_range", "inconclusive", "B", "A2"])
    assert sum(r[k] for k in decision.OUTCOMES) == pytest.approx(1.0) and r["A_total"] == pytest.approx(0.5)


def test_percentile_of():
    assert decision.percentile_of(0.5, np.linspace(0, 1, 101)) == pytest.approx(100 * 51 / 101)   # 51 of 101 <= 0.5


def test_d4m_a1_label_and_a2_only_evidence():
    out = _out(True, False, 0.8)
    assert out["outcome"] == "A1" and not out["creep_specific_evidence"]
    assert out["outcome_label"] == ("association consistent with the scaling, not specific to vortex creep; "
                                    "Δν̇ term not established (T2 power at the matched σ reported)")       # D4p item 4
    a2 = _out(True, True, 0.8)
    assert a2["creep_specific_evidence"] and a2["outcome_label"] is None


# ---- D4n -----------------------------------------------------------------------------------------
def test_d4n_t2_cut_is_the_smaller_2_5th_percentile():
    cut, per_null = decision.t2_d_cut({"S3": np.linspace(-2, 2, 1001), "S3m": np.linspace(-1, 1, 1001)}, load_config())
    assert per_null == pytest.approx({"S3": -1.9, "S3m": -0.95}) and cut == pytest.approx(-1.9)
    cut, _ = decision.t2_d_cut({"S3": np.linspace(-0.5, 0.5, 1001), "S3m": np.linspace(-1, 1, 1001)}, load_config())
    assert cut == pytest.approx(-0.95)                                       # now S3m's is the smaller


def test_d4n_short_window_label_only_on_a2_and_never_changes_the_outcome():
    cfg = load_config()
    a2 = decision.classify(True, True, 0.8, 0.0, _refs(), 0.34, cfg, d_without_short_window=0.0)
    assert a2["outcome"] == "A2" and a2["short_window_label"] == "short-window-sensitive"      # >= 0 counts
    assert decision.classify(True, True, 0.8, 0.0, _refs(), 0.34, cfg, -0.01)["short_window_label"] is None
    assert decision.classify(True, True, 0.8, 0.0, _refs(), 0.34, cfg)["short_window_label"] is None   # not computed
    a1 = decision.classify(True, False, 0.8, 0.0, _refs(), 0.34, cfg, 0.5)
    assert a1["outcome"] == "A1" and a1["short_window_label"] is None
