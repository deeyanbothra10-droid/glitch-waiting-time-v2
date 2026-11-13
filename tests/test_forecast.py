"""Forward-forecast maths (post-freeze, not pre-registered) on known answers."""
import numpy as np
import pytest
from scipy import stats

from glitchlib import forecast


def test_unconditional_quantiles_are_lognormal_quantiles():
    mu, sigma = 6.0, 0.7
    got = forecast.conditional_quantiles(mu, sigma, 0.0, [0.05, 0.5, 0.95])
    assert got == pytest.approx([float(stats.lognorm.ppf(q, s=sigma, scale=np.exp(mu))) for q in (0.05, 0.5, 0.95)])


def test_conditional_quantiles_known_answer():
    mu, sigma, t0 = 6.0, 0.7, 500.0
    s0 = stats.lognorm.sf(t0, s=sigma, scale=np.exp(mu))
    med = forecast.conditional_quantiles(mu, sigma, t0, [0.5])[0]
    assert stats.lognorm.sf(med, s=sigma, scale=np.exp(mu)) == pytest.approx(0.5 * s0)
    lo, hi = forecast.conditional_quantiles(mu, sigma, t0, [0.05, 0.95])
    assert t0 < lo < med < hi
    far = forecast.conditional_quantiles(mu, sigma, 1e6, [0.5])[0]                # deep tail stays finite
    assert np.isfinite(far) and far > 1e6


def test_log_predictive_known_answers():
    mu, sigma, t0 = 6.0, 0.7, 200.0
    dist = stats.lognorm(s=sigma, scale=np.exp(mu))
    assert forecast.log_predictive(mu, sigma, t0, 600.0, True) == pytest.approx(dist.logpdf(600.0) - dist.logsf(t0))
    assert forecast.log_predictive(mu, sigma, t0, 900.0, False) == pytest.approx(dist.logsf(900.0) - dist.logsf(t0))


def test_coverage_and_outcomes():
    assert forecast.coverage((100, 300), 200, True) is True and forecast.coverage((100, 300), 400, True) is False
    assert forecast.coverage((100, 300), 400, False) is False                     # window over, no glitch: a miss
    assert forecast.coverage((100, 500), 400, False) is None                       # not determinable
    assert forecast.outcome_for(61300.0, 61000.0, 61343.0, 62076.0) == {"status": "void"}
    assert forecast.outcome_for(61500.0, 61000.0, 61343.0, 62076.0) == {"status": "observed", "t": 500.0}
    assert forecast.outcome_for(None, 61000.0, 61343.0, 62076.0) == {"status": "censored", "t": 1076.0}
    assert forecast.mjd_to_date(61343.0) == "2026-10-30" and forecast.mjd_to_date(62076.0) == "2028-11-01"


def test_forecast_window_dates_follow_the_quantiles():
    w = forecast.forecast_window(6.0, 0.7, 0.0, 61000.0, {"median": 0.5, "lo": 0.05})
    assert w["waits_days"]["median"] == pytest.approx(float(np.exp(6.0)))
    assert w["dates"]["median"] == forecast.mjd_to_date(61000.0 + np.exp(6.0)) and w["t0_days"] == 0.0


def test_scoring_eligibility_rules():
    epochs = [60000.0, 60400.0, 61000.0]                        # longest wait 600 d
    ok = forecast.scoring_eligibility(True, 61000.0, 61000.0, epochs, 61000.0, 61343.0)   # placeholder obs_end
    assert ok["scored"] and ok["longest_wait_days"] == 600.0 and ok["snapshot_minus_g_days"] == 343.0
    b = forecast.scoring_eligibility(True, 61100.0, 61000.0, epochs, 61000.0, 61343.0)    # documented, before
    assert not b["scored"] and set(b["fails"]) == {"b"}
    c = forecast.scoring_eligibility(True, 60400.0, 60400.0, epochs[:2], 60400.0, 61343.0)  # 943 d > 400 d
    assert set(c["fails"]) == {"c"} and c["ratio"] == pytest.approx(943 / 400)
    a = forecast.scoring_eligibility(False, 61000.0, 61000.0, epochs, 61000.0, 61343.0)
    assert set(a["fails"]) == {"a"}
