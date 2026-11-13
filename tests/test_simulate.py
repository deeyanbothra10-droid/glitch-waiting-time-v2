"""Simulation scenarios on a default-free fixture: every parameter is defined here, not in config.

S4 bias directions and sizes are not tested (D4d): they are measured by scripts/emulation_fixture.py.
"""
import dataclasses

import numpy as np
import pandas as pd
import pytest

from glitchlib import simulate
from glitchlib.config import load_config

# ---- fixture (test-only values; the emulation settings mirror D4c/D4d) ---------------------------
N_PULSARS, N_POOL, N_REAL_GLITCHES = 10, 15, 25
FIRST, OBS_END = 50000.0, 75000.0                       # W_p = 25000 / 25 = 1000 d
_ONE_RANGE = simulate.ExponentialGroup("test", (0.005, 0.05), (1.0, 100.0), None)
SETTINGS = simulate.EmulationSettings(pre_glitch_days=365.0, min_samples_after=4, min_samples_before_window=4,
                                      min_aic_drop=10.0,
                                      exponentials=simulate.ExponentialSpec({}, None, _ONE_RANGE, _ONE_RANGE),
                                      unknown_as="window", short_window_days=None)
S1 = simulate.Scenario("S1", "ln_R", 1.0, False, None, None)
S2 = simulate.Scenario("S2", None, 0.0, False, None, None)
S1M = simulate.Scenario("S1m", "ln_R", 1.0, True, "equals_T", None)
S4 = simulate.Scenario("S4", None, 0.0, True, "lognormal_W", 0.3)
SIGMA = 0.4
METHOD_CLASSES = ["window", "local-relaxation", "global", "global-noexp", "unknown"]


def _pool(rng, n):
    dnu = np.exp(rng.normal(np.log(1e-6), 0.8, n))
    dnudot = np.exp(rng.normal(np.log(5e-3), 0.6, n))
    return pd.DataFrame({"dnu_frac": dnu, "dnu_frac_err": rng.uniform(0.005, 0.05, n) * dnu,
                         "dnudot_frac": dnudot, "dnudot_frac_err": rng.uniform(0.02, 0.2, n) * dnudot,
                         "method_class": rng.choice(METHOD_CLASSES, n), "cadence_days": rng.choice([7.0, 14.0, 30.0], n)})


def _templates(seed, n_pool=N_POOL):
    rng = np.random.default_rng(seed)
    return [simulate.PulsarTemplate(f"F{p}", FIRST, OBS_END, N_REAL_GLITCHES, _pool(rng, n_pool), 1e6, (), False, 1.0)
            for p in range(N_PULSARS)]


# ---- basic properties ----------------------------------------------------------------------------
def test_typical_wait_known_answer():
    assert simulate.typical_wait(50000.0, 75000.0, 25) == 1000.0


def test_sequence_stays_in_the_real_window():
    seq = simulate.simulate_sequence(_templates(1)[0], S4, SIGMA, np.random.default_rng(2))
    assert seq["mjd"].iloc[0] == FIRST
    assert (seq["mjd"] < OBS_END).all()
    assert seq["mjd"].iloc[-1] + seq["T_true"].iloc[-1] >= OBS_END   # the last wait runs past obs_end


def test_mean_wait_and_glitch_count_match_the_real_ones():
    # D4c: c_p = ln W_p - sigma^2/2, so the mean wait is W_p and counts match on average.
    template = _templates(3)[0]
    rng = np.random.default_rng(4)
    for sigma in (0.2, 0.6):
        seqs = [simulate.simulate_sequence(template, S2, sigma, rng) for _ in range(400)]
        waits = np.concatenate([s["T_true"].to_numpy() for s in seqs])
        assert abs(waits.mean() / 1000.0 - 1) < 0.02, sigma
        assert abs(np.mean([len(s) for s in seqs]) / N_REAL_GLITCHES - 1) < 0.03, sigma


def test_s1_waits_follow_the_model():
    template = _templates(5)[0]
    rng = np.random.default_rng(6)
    seq = pd.concat([simulate.simulate_sequence(template, S1, SIGMA, rng) for _ in range(200)])
    ln_r_pool = np.log(template.pool["dnu_frac"] / template.pool["dnudot_frac"] * 1e3)
    ln_r = np.log(seq["dnu_frac"] / seq["dnudot_frac"] * 1e3)
    centred_pool = ln_r_pool - ln_r_pool.mean()
    c_p = np.log(1000.0) - SIGMA ** 2 / 2 - np.log(np.mean(np.exp(centred_pool)))     # D4h intercept
    resid = np.log(seq["T_true"]) - c_p - (ln_r - ln_r_pool.mean())
    assert abs(resid.mean()) < 0.03 and abs(resid.std() - SIGMA) < 0.03
    assert seq["t_rec"].isna().all()                     # statistics layer: no recovery time


def test_statistics_layer_noise_uses_the_resampled_glitchs_relative_errors():
    template = _templates(7)[0]
    rng = np.random.default_rng(8)
    t = pd.concat([simulate.measure_statistics_layer(simulate.simulate_sequence(template, S2, SIGMA, rng), rng, True)
                   for _ in range(300)])
    for true, err, measured in (("dnudot_frac", "dnudot_frac_err", t["df1"] * 1e-3),
                                ("dnu_frac", "dnu_frac_err", t["dff"] * 1e-9)):
        z = np.log(measured / t[true]) / (t[err] / t[true])
        assert abs(z.mean()) < 0.03 and abs(z.std() - 1) < 0.03, true


def test_emulated_layer_adds_dnu_noise_too():
    template = _templates(9)[0]
    rng = np.random.default_rng(10)
    tables = [simulate.measure_emulated_layer(simulate.simulate_sequence(template, S4, SIGMA, rng), template,
                                              "global", SETTINGS, 1e-4, rng) for _ in range(40)]
    t = pd.concat(tables)
    z = np.log(t["dff"] * 1e-9 / t["dnu_frac"]) / (t["dnu_frac_err"] / t["dnu_frac"])
    assert abs(z.mean()) < 0.1 and abs(z.std() - 1) < 0.1


def test_s1m_recovery_ends_at_next_glitch_and_s4_is_independent():
    template = _templates(11)[0]
    rng = np.random.default_rng(12)
    s1m = pd.concat([simulate.simulate_sequence(template, S1M, SIGMA, rng) for _ in range(50)])
    assert np.allclose(s1m["t_rec"], s1m["T_true"])
    s4 = pd.concat([simulate.simulate_sequence(template, S4, SIGMA, rng) for _ in range(200)])
    assert abs(np.corrcoef(np.log(s4["t_rec"]), np.log(s4["T_true"]))[0, 1]) < 0.05
    assert abs(np.median(np.log(s4["t_rec"])) - np.log(1000.0)) < 0.03
    assert abs(np.log(s4["t_rec"]).std() - 0.3) < 0.02


def test_methods_come_from_the_resampled_glitch():
    template = _templates(13)[0]
    seq = simulate.simulate_sequence(template, S1M, SIGMA, np.random.default_rng(14))
    for unknown_as in ("window", "global"):
        classes, emulated = simulate.emulation_methods(seq, "resampled", unknown_as)
        assert list(classes) == list(template.pool["method_class"].iloc[seq["pool_index"]])
        assert all(e == unknown_as for c, e in zip(classes, emulated) if c == "unknown")
        assert all(e == c for c, e in zip(classes, emulated) if c != "unknown")


def test_window_cadence_follows_the_glitch_that_starts_it():
    seq = pd.DataFrame({"cadence_days": [7.0, 30.0, 14.0]})
    # Window 0 (before the first glitch) uses the first glitch's cadence; window j uses glitch j-1's.
    assert list(simulate.window_cadences(seq)) == [7.0, 7.0, 30.0, 14.0]


def test_missing_errors_are_filled_and_flagged():
    rng = np.random.default_rng(15)
    pools = {"A": _pool(rng, 6), "B": _pool(rng, 4)}
    pools["A"].loc[0, "dnudot_frac_err"] = np.nan                 # A keeps 5 known values: its own median
    pools["B"].loc[[0, 1], "dnu_frac_err"] = np.nan               # B keeps 2 known values: global median
    known_a = (pools["A"]["dnudot_frac_err"] / pools["A"]["dnudot_frac"]).dropna()
    rel_dnu = pd.concat([(p["dnu_frac_err"] / p["dnu_frac"]) for p in pools.values()]).dropna()
    filled = simulate.fill_missing_errors(pools, min_values=3)
    a, b = filled["A"], filled["B"]
    assert a.loc[0, "dnudot_frac_err"] / a.loc[0, "dnudot_frac"] == pytest.approx(np.median(known_a))
    assert list(a["dnudot_frac_err_imputed"]) == [True] + [False] * 5
    for i in (0, 1):
        assert b.loc[i, "dnu_frac_err"] / b.loc[i, "dnu_frac"] == pytest.approx(np.median(rel_dnu))
    assert list(b["dnu_frac_err_imputed"]) == [True, True, False, False]
    assert np.isnan(pools["A"].loc[0, "dnudot_frac_err"])         # the input is not modified


def test_pool_with_missing_errors_is_refused():
    template = _templates(16)[0]
    template.pool.loc[0, "dnudot_frac_err"] = np.nan
    with pytest.raises(ValueError, match="fill_missing_errors"):
        simulate.simulate_sequence(template, S2, SIGMA, np.random.default_rng(17))


def test_target_relative_error_ignores_imputed_values_and_falls_back_to_global_median():
    templates = _templates(18)
    templates[0].pool["dnudot_frac_err_imputed"] = [True] * 13 + [False] * 2     # only 2 real values
    real = [(tp.pool["dnudot_frac_err"] / tp.pool["dnudot_frac"]).to_numpy() for tp in templates]
    real[0] = real[0][13:]
    targets = simulate.target_relative_errors(templates, min_values=3)
    assert targets["F0"] == pytest.approx(np.median(np.concatenate(real)))
    assert targets["F1"] == pytest.approx(np.median(real[1]))


@pytest.mark.parametrize("sigma", [0.2, 0.6])
def test_noise_calibration_under_s1m_meets_its_tolerance_for_local_relaxation(sigma):
    # D4e: calibrated under S1m (t_rec = T). Local-relaxation has no model misfit there, so the
    # bisection reaches the target.
    template = _templates(19)[0]
    out = simulate.calibrate_noise(template, S1M, sigma, "local-relaxation", SETTINGS, 0.10, 4, 0.05,
                                   np.random.default_rng(20))
    assert out["within_tolerance"] and not out["at_floor"] and out["noise_sd"] > 0
    assert abs(out["achieved"] - 0.10) <= 0.05 * 0.10


def test_methods_with_misfit_under_s1m_can_sit_at_the_floor():
    # Window, global (one nu_ddot for windows whose recovery slopes J_k/T_k differ) and global-noexp
    # misfit even when t_rec = T. With them in the mix, zero noise can already exceed the target; the
    # calibration then reports at_floor with noise 0 (recorded in logs/emulation_findings.md).
    template = _templates(19)[0]
    out = simulate.calibrate_noise(template, S1M, 0.4, "window", SETTINGS, 0.10, 4, 0.05, np.random.default_rng(21))
    assert out["at_floor"] and out["noise_sd"] == 0.0 and not out["within_tolerance"]


def test_scenarios_and_settings_from_config():
    cfg = load_config()
    s = {name: simulate.scenario_from_config(cfg, name) for name in ("S1", "S2", "S3", "S1m", "S4")}
    assert (s["S1"].covariate, s["S1"].slope, s["S1"].emulated) == ("ln_R", 1, False)
    assert (s["S2"].covariate, s["S2"].emulated) == (None, False)
    assert (s["S3"].covariate, s["S3"].slope) == ("ln_dnu_frac", 1)
    assert (s["S1m"].emulated, s["S1m"].t_rec) == (True, "equals_T")
    assert (s["S4"].t_rec, s["S4"].t_rec_log_sd) == ("lognormal_W", 0.3)
    settings = simulate.emulation_from_config(cfg, "window")
    assert (settings.pre_glitch_days, settings.min_samples_after, settings.min_aic_drop) == (365, 4, 10)
    with pytest.raises(ValueError, match="unknown_as"):
        simulate.emulation_from_config(cfg, "local-relaxation")


def test_one_replicate_of_each_layer_runs_end_to_end():
    templates = _templates(21)[:3]
    rng = np.random.default_rng(22)
    intervals, glitches = simulate.simulate_statistics_dataset(templates, S1, SIGMA, True, rng)
    assert {"ln_R", "t", "event", "prev_t"} <= set(intervals.columns) and intervals["usable"].all()
    noise = {tp.name: 1e-4 for tp in templates}
    intervals, glitches = simulate.simulate_emulated_dataset(templates, S1M, SIGMA, "resampled", SETTINGS, noise, rng)
    assert {"method_class", "emulated_method", "exp_fitted"} <= set(glitches.columns)
    assert len(intervals) == int(glitches["usable"].sum())


def test_short_window_variant_and_mix_use_their_spans():
    import dataclasses
    template = _templates(23)[0]
    rng = np.random.default_rng(24)
    seq = simulate.simulate_sequence(template, S4, SIGMA, rng)
    fixed = dataclasses.replace(SETTINGS, short_window_days=40.0)
    table = simulate.measure_emulated_layer(seq, template, "short-window", fixed, 1e-4, rng)
    assert (table["short_window_days"] == 40.0).all() and table["usable"].any()
    pool = template.pool.assign(method_class="short-window", short_window_days=125.0)
    mixed = simulate.PulsarTemplate(template.name, FIRST, OBS_END, N_REAL_GLITCHES, pool, 1e6, (), False, 1.0)
    seq = simulate.simulate_sequence(mixed, S4, SIGMA, rng)
    table = simulate.measure_emulated_layer(seq, mixed, "resampled", SETTINGS, 1e-4, rng)
    assert (table["short_window_days"] == 125.0).all() and (table["emulated_method"] == "short-window").all()


def test_short_window_admission_rule():
    passing = [{"w": 40, "sigma": 0.2, "mean": 0.005, "se_mean": 0.005},     # 0.005 + 1.96*0.005 = 0.0148
               {"w": 150, "sigma": 0.6, "mean": 0.0, "se_mean": 0.01}]      # 0.0196
    assert simulate.short_window_admission(passing, 0.02, 1.96)["admitted"]
    failing = passing + [{"w": 125, "sigma": 0.4, "mean": 0.01, "se_mean": 0.006}]   # 0.02176
    verdict = simulate.short_window_admission(failing, 0.02, 1.96)
    assert not verdict["admitted"] and [c["passes"] for c in verdict["checks"]] == [True, True, False]


def test_admission_fails_when_a_span_gives_no_estimate():
    entries = [{"w": 40, "sigma": 0.2, "mean": None, "se_mean": None},
               {"w": 150, "sigma": 0.2, "mean": 0.0, "se_mean": 0.001}]
    verdict = simulate.short_window_admission(entries, 0.02, 1.96)
    assert not verdict["admitted"] and verdict["checks"][0]["upper_95"] is None


def test_glitches_inside_a_gap_go_unobserved_and_spanning_intervals_are_censored():
    import dataclasses
    template = dataclasses.replace(_templates(30)[0], gaps=((55000.0, 60000.0),))
    rng = np.random.default_rng(31)
    for _ in range(20):
        seq = simulate.simulate_sequence(template, S2, SIGMA, rng)
        assert not ((seq["mjd"] > 55000.0) & (seq["mjd"] < 60000.0)).any()
    intervals, _ = simulate.simulate_statistics_dataset([template], S2, SIGMA, True, rng)
    spanning = intervals[(intervals["mjd"] < 55000.0) & (intervals["mjd"] + intervals["t"] >= 55000.0)]
    assert len(spanning) == 1 and spanning["gap_flag"].all() and not spanning["event"].any()
    assert spanning["t"].iloc[0] == pytest.approx(55000.0 - spanning["mjd"].iloc[0])


def test_last_censored_interval_dropped_when_obs_end_is_the_last_real_glitch():
    import dataclasses
    rng = np.random.default_rng(32)
    keep = _templates(33)[0]
    drop = dataclasses.replace(keep, name="D", drop_censored_last=True)
    intervals, _ = simulate.simulate_statistics_dataset([keep, drop], S2, SIGMA, True, rng)
    assert (~intervals.loc[intervals["pulsar"] == keep.name, "event"].astype(bool)).sum() == 1
    assert intervals.loc[intervals["pulsar"] == "D", "event"].astype(bool).all()


def test_admissible_fraction_marks_glitches():
    import dataclasses
    template = dataclasses.replace(_templates(34)[0], admissible_fraction=0.3)
    rng = np.random.default_rng(35)
    seq = pd.concat([simulate.simulate_sequence(template, S2, SIGMA, rng) for _ in range(200)])
    assert abs(seq["admissible"].mean() - 0.3) < 0.02
    full = simulate.simulate_sequence(_templates(34)[0], S2, SIGMA, rng)
    assert full["admissible"].all()


def test_measurement_noise_off_gives_true_values():
    template = _templates(36)[0]
    rng = np.random.default_rng(37)
    table = simulate.measure_statistics_layer(simulate.simulate_sequence(template, S1, SIGMA, rng), rng, False)
    assert np.allclose(table["df1"], table["dnudot_frac"] * 1e3) and np.allclose(table["dff"], table["dnu_frac"] * 1e9)
    assert (table["df1_err"] == 0).all()


@pytest.mark.parametrize("scenario", [S1, simulate.Scenario("S3", "ln_dnu_frac", 1.0, False, None, None)])
def test_d4h_intercept_makes_the_mean_wait_W_with_a_covariate(scenario):
    # D4h: c_p = ln W - sigma^2/2 - ln m_p, m_p = pool mean of exp(slope * centred covariate).
    template = _templates(50)[0]
    rng = np.random.default_rng(51)
    waits = np.concatenate([simulate.simulate_sequence(template, scenario, 0.4, rng)["T_true"].to_numpy()
                            for _ in range(600)])
    assert abs(waits.mean() / 1000.0 - 1) < 0.03


# ---- D4o: literature-based exponential recoveries ------------------------------------------------------
def test_d4o_spec_from_config_known_answer():
    spec = simulate.exponential_spec_from_config(load_config())
    j, crab = spec.named["J0537-6910"], spec.named["B0531+21"]
    assert (j.Q, j.tau_d_days, j.min_dnu_microhz) == ((0.002, 0.04), (4.0, 37.0), 20)
    assert (crab.Q, crab.tau_d_days, crab.min_dnu_microhz) == ((0.58, 1.0), (2.0, 21.1), None)
    assert spec.others_split == 1e-7
    assert (spec.others_at_least.Q, spec.others_at_least.tau_d_days) == ((0.002, 0.06), (5.8, 488.0))
    assert (spec.others_below.Q, spec.others_below.tau_d_days) == ((0.10, 0.15), (54.0, 2340.0))


def test_d4o_exponential_ranges_by_group():
    spec = simulate.exponential_spec_from_config(load_config())
    seq = pd.DataFrame({"dnu_frac": [5e-8, 2e-7, 1e-7], "dnu_microhz": [10.0, 25.0, 20.0]})
    groups, q, tau, present = simulate.exponential_ranges(seq, "J0537-6910", spec)
    assert list(present) == [False, True, True]                             # Δν >= 20 µHz only
    assert (q[1].tolist(), tau[1].tolist()) == ([0.002, 0.04], [4.0, 37.0])
    groups, q, tau, present = simulate.exponential_ranges(seq, "B0531+21", spec)
    assert present.all() and set(groups) == {"crab"} and q[0].tolist() == [0.58, 1.0]
    groups, q, tau, present = simulate.exponential_ranges(seq, "J1841-0524", spec)
    assert list(groups) == ["others_below_split", "others_at_least_split", "others_at_least_split"]   # split at 1e-7
    assert tau[0].tolist() == [54.0, 2340.0] and q[2].tolist() == [0.002, 0.06]
    with pytest.raises(ValueError, match="dnu_microhz"):
        simulate.exponential_ranges(seq.drop(columns="dnu_microhz"), "J0537-6910", spec)


def test_d4o_draws_lie_in_their_ranges_and_absent_exponentials_have_q_zero():
    cfg = load_config()
    spec = simulate.exponential_spec_from_config(cfg)
    template = dataclasses.replace(_templates(70)[0], name="J0537-6910")
    template.pool["dnu_microhz"] = np.where(np.arange(len(template.pool)) % 2 == 0, 10.0, 30.0)
    seq = simulate.simulate_sequence(template, S1M, 0.4, np.random.default_rng(1))
    settings = dataclasses.replace(SETTINGS, exponentials=spec)
    draws = simulate._emulation_draws(seq, template, "resampled", settings, np.random.default_rng(2))[0]
    big = draws["dnu_microhz"] >= 20
    assert (draws.loc[~big, "Q"] == 0).all() and not draws.loc[~big, "exp_present"].any()
    assert draws.loc[big, "Q"].between(0.002, 0.04).all() and draws["tau_d"].between(4.0, 37.0).all()


def test_placeholder_spec_reproduces_the_old_d14_draws():
    # The superseded draw: tau_d = exp(U(ln 1, ln 100)) for every glitch, then Q = exp(U(ln 0.005, ln 0.05)).
    cfg = load_config()
    template = _templates(71)[0]
    seq = simulate.simulate_sequence(template, S1M, 0.4, np.random.default_rng(3))
    settings = dataclasses.replace(SETTINGS, exponentials=simulate.placeholder_exponential_spec(cfg))
    draws = simulate._emulation_draws(seq, template, "resampled", settings, np.random.default_rng(4))[0]
    rng = np.random.default_rng(4)
    n = len(seq)
    tau = np.exp(rng.uniform(np.log(1.0), np.log(100.0), n))
    q = np.exp(rng.uniform(np.log(0.005), np.log(0.05), n))
    assert np.allclose(draws["tau_d"], tau, rtol=1e-12) and np.allclose(draws["Q"], q, rtol=1e-12)
