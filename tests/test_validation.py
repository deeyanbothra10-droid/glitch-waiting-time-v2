"""The validation engine: replicate records, seeds and D4e inclusion (synthetic templates)."""
import dataclasses

import numpy as np
import pytest
import pandas as pd

from glitchlib import validation
from glitchlib.config import load_config
from test_simulate import S2, _templates


def _job(n=3):
    return validation.Job(7, "S2", S2, None, True, None, n, True)


def test_replicate_is_reproducible_and_complete():
    cfg = load_config()
    templates = _templates(40)[:4]
    args = (templates, _job(), 0.4, 1, 2, cfg, None)
    a, b = validation._replicate(args), validation._replicate(args)
    a.pop("seconds"), b.pop("seconds")                       # wall-clock time differs between runs
    assert a == b
    assert {"b", "T1_p", "T2_p", "T1_signed", "T2_signed", "d_se", "placebo_b", "perm_p", "loo_M4", "n_completed",
            "sigma_hat", "d_without_short_window"} <= set(a)


def test_replicate_with_too_few_pulsars_is_discarded():
    cfg = load_config()
    template = dataclasses.replace(_templates(41)[0], admissible_fraction=0.01)
    record = validation._replicate(([template], _job(), 0.4, 1, 0, cfg, None))
    assert record["discarded"] and record["pulsars_dropped"] == 1


def test_run_job_and_rates():
    cfg = load_config()
    records = validation.run_job(_templates(42)[:4], dataclasses.replace(_job(4), secondary=False), cfg, 0.4, 1,
                                 None, 2)
    assert len(records) == 4 and list(records["replicate"]) == [0, 1, 2, 3]
    assert 0 <= validation.rate(records, "T1_rejects_chi2") <= 1
    assert validation.percentile_interval([1, 2, 3, 4, 5], 0.5) == [2.0, 4.0]
    assert pd.isna(records.get("perm_p", pd.Series([np.nan]))).all()


def test_false_positive_bound_known_answer():
    # D4h: 0.025 + 3 sqrt(0.025 * 0.975 / 1000) = 0.0398
    assert validation.false_positive_bound(0.025, 3, 1000) == pytest.approx(0.039811, abs=1e-6)


def test_emulated_replicate_records_missingness():
    from test_simulate import S1M, SETTINGS
    cfg = load_config()
    job = validation.Job(8, "S1m", S1M, "resampled", True, None, 1, False)
    templates = _templates(43)[:4]
    record = validation._replicate((templates, job, 0.4, 1, 0, cfg, {tp.name: 1e-3 for tp in templates}))
    assert "missingness_log_or" in record and record["failure"] is None
    assert any(k.startswith("meas_slope_") for k in record)                     # D4l item 4


def test_replicate_over_the_time_limit_is_a_failure(monkeypatch):
    monkeypatch.setattr(validation, "REPLICATE_TIME_LIMIT_S", 1)
    monkeypatch.setattr(validation, "_replicate_body", lambda args: __import__("time").sleep(5))
    record = validation._replicate((None, None, 0.4, 1, 7, None, None))
    assert record["failed"] and record["failure"].startswith("timeout") and record["replicate"] == 7
    assert record["seconds"] < 3


def test_unexpected_error_in_a_replicate_is_recorded_not_raised(monkeypatch):
    def boom(args):
        raise ZeroDivisionError("float division by zero")
    monkeypatch.setattr(validation, "_replicate_body", boom)
    record = validation._replicate((None, None, 0.4, 1, 3, None, None))
    assert record["failed"] and record["failure"] == "error: ZeroDivisionError: float division by zero"


# ---- D4k -----------------------------------------------------------------------------------------
def test_d4k_1_failure_limit():
    cfg = load_config()
    rec = pd.DataFrame({"replicate": range(200), "failed": [True] * 2 + [False] * 198,
                        "failure": ["fit: x"] * 2 + [None] * 198})
    check = validation.failure_check(rec, cfg)
    assert check["result"] == "PASS" and check["failed"] == 2 and check["fraction"] == 0.01
    assert check["reasons"] == [{"replicate": 0, "failure": "fit: x"}, {"replicate": 1, "failure": "fit: x"}]
    rec.loc[2, ["failed", "failure"]] = [True, "timeout after 120 s"]
    assert validation.failure_check(rec, cfg)["result"] == "FAIL"               # 1.5% > 1%


def test_d4k_1_failed_replicates_are_excluded_from_rates():
    rec = pd.DataFrame({"discarded": [False, False, False], "failed": [False, True, False],
                        "T1_rejects": [True, True, False]})
    assert validation.rate(rec, "T1_rejects") == 0.5


def test_d4k_2_validation_grid():
    assert validation.validation_sigma_grid(load_config()) == [0.2, 0.4, 0.6, 0.8, 1.2]


def test_d4k_3_reference_jobs_and_seeds_as_in_06():
    cfg = load_config()
    s1, s1m, s4 = validation.reference_jobs(cfg)
    assert (s1.job_id, s1.name, s1.n_replicates) == (1, "S1", 1000)
    assert (s1m.job_id, s1m.name, s1m.n_replicates, s1m.method_rule) == (3, "S1m", 1000, "resampled")
    assert (s4.job_id, s4.name, s4.n_replicates, s4.scenario.t_rec) == (4, "S4-mix", 500, "lognormal_W")
    assert validation.real_sigma_index(cfg) == 5                                # past the 5 grid indices


def test_d4k_4_missingness_reference_counts_estimable_only():
    rec = pd.DataFrame({"discarded": [False] * 4, "failed": [False, False, False, True],
                        "missingness_log_or": [0.5, None, 1.5, 9.0]})
    ref = validation.missingness_reference(rec)
    assert (ref["replicates"], ref["estimable"], ref["not_estimable"]) == (3, 2, 1)
    assert list(ref["values"]) == [0.5, 1.5]


def test_partial_store_saves_resumes_and_refuses_a_different_run(tmp_path):
    store = validation.PartialStore(tmp_path / "partial", "abc", "cfg1")
    rec = pd.DataFrame({"replicate": [0, 1], "b": [0.1, 0.2]})
    store.save_job(0.4, "S2", rec)
    store.save_calibration(0.4, [{"pulsar": "A", "noise_sd": 1e-3}])
    again = validation.PartialStore(tmp_path / "partial", "abc", "cfg1")
    pd.testing.assert_frame_equal(again.load_job(0.4, "S2"), rec)
    assert again.load_job(0.4, "S3") is None and again.load_calibration(0.4)[0]["pulsar"] == "A"
    with pytest.raises(SystemExit, match="Resume refused"):
        validation.PartialStore(tmp_path / "partial", "def", "cfg1")           # different commit
    with pytest.raises(SystemExit, match="Resume refused"):
        validation.PartialStore(tmp_path / "partial", "abc", "cfg2")           # different config


# ---- D4l -----------------------------------------------------------------------------------------
def test_null_jobs_have_their_own_seeds_and_slope():
    cfg = load_config()
    s2, s3, s3m = validation.null_jobs(cfg, 0.73)
    assert (s2.job_id, s2.name, s2.n_replicates, s2.scenario.covariate) == (41, "S2-null", 1000, None)
    assert (s3.job_id, s3.name, s3.scenario.covariate, s3.scenario.slope) == (42, "S3-null", "ln_dnu_frac", 0.73)
    assert (s3m.job_id, s3m.name, s3m.method_rule, s3m.scenario.slope, s3m.scenario.emulated) == (
        43, "S3m-null", "resampled", 0.73, True)
    evaluation_ids = {1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 20, 21}
    assert not {s2.job_id, s3.job_id, s3m.job_id} & evaluation_ids


def test_saved_p_values_give_back_the_signed_roots():
    from glitchlib import analysis
    cfg = load_config()
    stats_ = [0.5, 3.0, 25.0, 400.0]
    rec = pd.DataFrame({"T1_p": [__import__("scipy").stats.chi2.sf(s, 1) for s in stats_], "b": [0.2, -0.1, 0.9, 1.0],
                        "T2_p": [__import__("scipy").stats.chi2.sf(s, 1) for s in stats_], "d": [0.1, -0.1, -0.8, 0.2],
                        "T1_rejects": [False, False, True, True], "T2_rejects": [False, False, True, True]})
    out = validation.from_saved_records(rec, cfg)
    assert out["T1_signed"].tolist() == pytest.approx([analysis.signed_root(s, b, "b_positive") for s, b in zip(stats_, rec["b"])])
    assert out["T2_signed"].tolist() == pytest.approx([analysis.signed_root(s, d, "d_negative") for s, d in zip(stats_, rec["d"])])
    assert "T1_rejects_chi2" in out and "T1_rejects" not in out


def test_directional_rejection_and_thresholds():
    cfg = load_config()
    nulls = {"S2-null": pd.DataFrame({"discarded": False, "failed": False, "T1_signed": np.linspace(-3, 3, 1001)}),
             "S3-null": pd.DataFrame({"discarded": False, "failed": False, "T2_signed": np.linspace(-2, 2, 1001),
                                      "d": np.linspace(-0.2, 0.2, 1001)}),
             "S3m-null": pd.DataFrame({"discarded": False, "failed": False, "T2_signed": np.linspace(-1, 1, 1001),
                                       "d": np.linspace(-0.1, 0.1, 1001)})}
    thr = validation.thresholds(nulls, cfg)
    assert thr["T1"] == pytest.approx(2.85)
    assert thr["T2_d_cut"] == pytest.approx(-0.19)                          # D4n: the smaller 2.5th percentile
    assert thr["T2_d_cut_per_null"] == pytest.approx({"S3": -0.19, "S3m": -0.095})
    assert thr["T2_signed_root_D4m_report_only"] == pytest.approx(1.9)      # D4m: the larger 97.5th, report only
    nulls["S3m-null"]["d"] = np.linspace(-0.4, 0.4, 1001)
    assert validation.thresholds(nulls, cfg)["T2_d_cut"] == pytest.approx(-0.38)     # now S3m's is the smaller
    rec = validation.apply_directional(pd.DataFrame({"T1_signed": [2.9, -5.0], "T2_signed": [1.8, 2.0],
                                                     "d": [-0.30, -0.10]}), thr)
    assert rec["T1_rejects"].tolist() == [True, False]
    assert rec["T2_rejects"].tolist() == [True, False]                       # decided by d-hat only
    assert rec["T2_rejects_D4m"].tolist() == [False, True]                   # the old signed-root test, report only


# ---- D4m -----------------------------------------------------------------------------------------
def test_interval_factors_and_calibrated_intervals():
    cfg = load_config()
    nulls = {"S2-null": pd.DataFrame({"discarded": False, "failed": False, "b": [-2.0, 0.0, 2.0], "b_se": [1.0, 1.0, 1.0]}),
             "S3m-null": pd.DataFrame({"discarded": False, "failed": False, "d": [-1.0, 0.0, 1.0], "d_se": [0.5, 0.5, 0.5]})}
    factors = validation.interval_factors(nulls, cfg)
    assert factors == pytest.approx({"b": 2.0, "d": 2.0})                   # SD 2 / SE 1; SD 1 / SE 0.5
    assert validation.calibrated_interval(1.0, 0.1, 2.0, 0.95) == pytest.approx([1 - 1.959964 * 0.2, 1 + 1.959964 * 0.2])
    assert validation.calibrated_interval(1.0, 0.1, 1.0, 0.95) == pytest.approx([1 - 1.959964 * 0.1, 1 + 1.959964 * 0.1])


def test_s3m_is_s3_through_s4_mix_emulation():
    from glitchlib import simulate
    cfg = load_config()
    s3m, s3, s4 = (simulate.scenario_from_config(cfg, n) for n in ("S3m", "S3", "S4"))
    assert (s3m.covariate, s3m.slope) == (s3.covariate, s3.slope) and s3m.emulated
    assert (s3m.t_rec, s3m.t_rec_log_sd) == (s4.t_rec, s4.t_rec_log_sd)


# ---- D4n -----------------------------------------------------------------------------------------
def test_sigma_hat_mean_uses_analysed_replicates_only():
    rec = pd.DataFrame({"discarded": [False, False, True], "failed": [False, False, False], "sigma_hat": [0.8, 1.0, 9.0]})
    assert validation.sigma_hat_mean(rec) == pytest.approx(0.9)


def test_match_sigma_bisection_known_answer():
    cfg = load_config()
    fitted = lambda s: np.sqrt(s ** 2 + 0.5 ** 2)                          # fitted scatter = sqrt(sigma^2 + 0.25)
    out = validation.match_sigma(fitted, 0.8, cfg)
    assert out["converged"] and not out["at_lower_end"] and abs(fitted(out["sigma"]) - 0.8) <= 0.01
    assert abs(out["sigma"] - np.sqrt(0.8 ** 2 - 0.25)) < 0.02
    assert out["steps"][0][0] == 0.01 and out["steps"][1][0] == pytest.approx(1.005)    # D4p: lower end 0.01


def test_match_sigma_lower_end_and_no_convergence():
    cfg = load_config()
    fitted = lambda s: np.sqrt(s ** 2 + 0.5 ** 2)
    out = validation.match_sigma(fitted, 0.4, cfg)                            # sigma = 0.01 already gives ~0.5 > 0.4
    assert out["sigma"] == 0.01 and out["at_lower_end"] and len(out["steps"]) == 1
    assert validation.match_sigma(fitted, 0.505, cfg)["sigma"] == 0.01         # within 0.01 at the lower end
    far = validation.match_sigma(fitted, 5.0, cfg)                             # beyond sigma = 2: never within 0.01
    assert not far["converged"] and len(far["steps"]) == 1 + 30                # D4p: 30 bisection steps
    lo, hi = 0.01, 2.0
    for s_, _ in far["steps"][1:]:
        lo = s_                                                                # every step moves the lower end up
    assert far["sigma"] == pytest.approx((lo + hi) / 2)                        # midpoint of the final bracket


def test_match_job_sigma_runs_200_replicates_per_step_with_fixed_search_seeds(monkeypatch):
    cfg = load_config()
    calls = []

    def fake_run_job(templates, job, cfg_, sigma, sigma_index, noise, workers):
        calls.append((job.n_replicates, job.secondary, sigma_index, noise))
        n = job.n_replicates
        return pd.DataFrame({"discarded": False, "failed": False, "sigma_hat": np.full(n, np.sqrt(sigma ** 2 + 0.36))})
    monkeypatch.setattr(validation, "run_job", fake_run_job)
    out = validation.match_job_sigma([], {}, _job(1000), 1.0, cfg, 1)
    assert abs(np.sqrt(out["sigma"] ** 2 + 0.36) - 1.0) <= 0.01
    assert {c[0] for c in calls} == {200} and not any(c[1] for c in calls)
    assert {c[2] for c in calls} == {validation.real_sigma_index(cfg) + 1} and all(c[3] is None for c in calls)


# ---- shared by 05, 06, 08 and the report-only scripts ------------------------------------------------
def test_run_job_listed_replicates_regenerate_the_same_records():
    cfg = load_config()
    job = dataclasses.replace(_job(4), secondary=False)
    templates = _templates(44)[:4]
    full = validation.run_job(templates, job, cfg, 0.4, 1, None, 2).drop(columns="seconds")
    some = validation.run_job(templates, job, cfg, 0.4, 1, None, 2, replicates=[3, 1]).drop(columns="seconds")
    assert list(some["replicate"]) == [1, 3]
    pd.testing.assert_frame_equal(some.reset_index(drop=True), full.iloc[[1, 3]].reset_index(drop=True))
    assert validation.run_job(templates, job, cfg, 0.4, 1, None, 2, replicates=[]).empty


def test_validation_jobs_keep_06s_ids_and_seed_indices_are_disjoint():
    cfg = load_config()
    ids = {job.name: job.job_id for job in validation.validation_jobs(cfg)}
    assert ids == {"S1": 1, "S2": 2, "S3": 5, "S1-noise-off": 6, "S3m": 7, "S1m": 3, "S4-mix": 4, "S4w": 10, "S4l": 11,
                   "S4g": 12, "S4g0": 13, "S4sw-espinoza": 20, "S4sw-wang": 21}
    grid = range(len(validation.validation_sigma_grid(cfg)))
    special = [validation.real_sigma_index(cfg), validation.matched_sigma_index(cfg), validation.floor_sigma_index(cfg)]
    assert len(set(special)) == 3 and not set(special) & set(grid)


def test_load_saved_searches_directories_in_order(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir(), b.mkdir()
    pd.DataFrame({"x": [1]}).to_pickle(b / "sigma_0.2__S1.pkl")
    pd.DataFrame({"x": [2]}).to_pickle(a / "sigma_0.4__S1.pkl")
    pd.DataFrame({"x": [3]}).to_pickle(b / "sigma_0.4__S1.pkl")
    assert validation.load_saved([a, b], 0.2, "S1")[0]["x"].tolist() == [1]
    rec, where = validation.load_saved([a, b], 0.4, "S1")
    assert rec["x"].tolist() == [2] and where == a
    assert validation.load_saved([a, b], 0.6, "S1") == (None, None)
    (b / "sigma_0.2__calibration.json").write_text('[{"pulsar": "P", "noise_sd": 0.1}]')
    assert validation.load_saved_calibration([a, b], 0.2) == [{"pulsar": "P", "noise_sd": 0.1}]


def test_classified_and_replicate_outcomes_known_answer():
    cfg = load_config()
    base = lambda **cols: pd.DataFrame({"discarded": False, "failed": False, **cols})
    lin = lambda lo, hi: np.linspace(lo, hi, 1001)
    nulls = {"S2-null": base(T1_signed=lin(-3, 3)), "S3-null": base(T2_signed=lin(-2, 2), d=lin(-0.2, 0.2)),
             "S3m-null": base(T2_signed=lin(-2, 2), d=lin(-0.1, 0.1))}            # T1 2.85; d cut -0.19; D4m 1.9
    ref = lambda b: base(T1_signed=np.full(1001, 3.0), T2_signed=np.zeros(1001), d=np.zeros(1001), b=b,
                         placebo_b=lin(-1, 1))
    test = base(replicate=[0, 1, 2], T1_signed=[3.0, 3.0, 1.0], T2_signed=[2.0, 0.0, 0.0], d=[0.0, -0.5, -0.5],
                b=[0.8, 0.8, 0.8], placebo_b=[0.0, 0.0, 0.0], d_without_short_window=[np.nan, 0.1, np.nan])
    records = {"S1m": ref(lin(0, 2)), "S4-mix": ref(lin(-0.2, 0.2)), "S1": ref(lin(0, 1)), "test": test}
    out, thr, refs, powers = validation.classified(records, nulls, cfg)
    assert thr["T2_d_cut"] == pytest.approx(-0.19) and powers["T2"] == 0.0 and powers["T2_signed_root_D4m_report_only"] == 0.0
    res = validation.replicate_outcomes(validation.ok(out["test"]), refs, powers["T2"], cfg)
    assert [r["outcome"] for r in res] == ["A1", "A2", "inconclusive"]        # T2 by d-hat; replicate 2: T1 fails
    assert res[1]["short_window_label"] == "short-window-sensitive" and res[0]["short_window_label"] is None
    old = validation.replicate_outcomes(validation.ok(out["test"]), refs, 0.0, cfg, "T2_rejects_D4m")
    assert [r["outcome"] for r in old] == ["A2", "A1", "inconclusive"]       # the signed-root T2, report only


def test_d4o_detection_counts_known_answer():
    glitches = pd.DataFrame({"exp_group": ["a", "a", "a", "b"], "exp_present": [True, True, False, True],
                             "exp_fitted": [True, False, True, False]})
    assert validation.detection_counts(glitches) == {
        "a": {"glitches": 3, "with_exponential": 2, "detected": 1, "detected_without_exponential": 1},
        "b": {"glitches": 1, "with_exponential": 1, "detected": 0, "detected_without_exponential": 0}}


def test_d4o_emulated_detection_counts_match_the_replicates_emulation():
    from test_simulate import S1M, SETTINGS
    from glitchlib import simulate
    cfg = load_config()
    job = validation.Job(8, "S1m", S1M, "resampled", True, None, 3, False)
    templates = _templates(45)[:3]
    noise = {tp.name: 1e-3 for tp in templates}
    total = validation.emulated_detection_counts(templates, job, cfg, 0.4, 1, noise, SETTINGS, 2)
    expected = {}
    for rep in range(3):
        rng = np.random.default_rng(np.random.SeedSequence((cfg["D4b.simulation.seed"], 8, 1, rep)))
        _, g = simulate.simulate_emulated_dataset(templates, S1M, 0.4, "resampled", SETTINGS, noise, rng)
        for k, v in validation.detection_counts(g)["test"].items():
            expected[k] = expected.get(k, 0) + v
    assert total == {"test": expected} and expected["glitches"] > 0
