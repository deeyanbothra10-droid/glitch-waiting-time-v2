"""Interval construction on synthetic glitch tables."""
import numpy as np
import pandas as pd
import pytest

from glitchlib.intervals import add_derived, build_intervals


def _example():
    # The Stage 4 worked example: glitches at 51000 (2100, 6.0), 52100 (30, missing),
    # 52900 (1800, 4.0), 53950 (2500, 7.5); obs_end 55000.
    return pd.DataFrame({"pulsar": "P", "mjd": [52900.0, 51000.0, 53950.0, 52100.0],
                         "dff": [1800.0, 2100.0, 2500.0, 30.0], "df1": [4.0, 6.0, 7.5, np.nan],
                         "usable": [True, True, True, False]})


def test_worked_example():
    intervals, dropped = build_intervals(_example(), {"P": 55000.0}, {})
    assert list(intervals["mjd"]) == [51000.0, 52900.0, 53950.0]      # 52100 starts none
    assert list(intervals["t"]) == [1100.0, 1050.0, 1050.0]
    assert list(intervals["event"]) == [True, True, False]
    assert np.isnan(intervals["prev_t"].iloc[0])
    assert list(intervals["prev_t"].iloc[1:]) == [800.0, 1050.0]      # previous glitch, usable or not
    assert dropped.empty


def test_orders_by_mjd_not_by_row_or_number():
    shuffled = _example().sample(frac=1.0, random_state=3)
    a, _ = build_intervals(shuffled, {"P": 55000.0}, {})
    b, _ = build_intervals(_example(), {"P": 55000.0}, {})
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))


def test_missing_obs_end_drops_last_interval_and_logs_it():
    intervals, dropped = build_intervals(_example(), {}, {})
    assert list(intervals["mjd"]) == [51000.0, 52900.0]
    assert list(dropped["mjd"]) == [53950.0]


def test_obs_end_before_last_glitch_is_refused():
    with pytest.raises(ValueError, match="end of monitoring"):
        build_intervals(_example(), {"P": 53000.0}, {})


def test_obs_end_equal_to_last_glitch_drops_and_logs_the_last_interval():
    # D4g: obs_end = last glitch epoch means no censoring information. The last interval is dropped
    # and logged, never entered with censoring time 0; the other intervals are unchanged.
    intervals, dropped = build_intervals(_example(), {"P": 53950.0}, {})
    assert list(intervals["mjd"]) == [51000.0, 52900.0] and (intervals["t"] > 0).all()
    assert list(intervals["event"]) == [True, True]
    assert dropped.to_dict("records") == [{"pulsar": "P", "mjd": 53950.0,
                                           "reason": "obs_end equals the last glitch epoch: no censoring information (D4g)"}]


def test_interval_status_known_answer():
    from glitchlib.intervals import interval_status
    data = _example()
    by_mjd = lambda status: dict(zip(data["mjd"], status))
    assert by_mjd(interval_status(data, {"P": 55000.0}, {})) == {
        51000.0: "completed", 52100.0: "completed", 52900.0: "completed", 53950.0: "censored_end"}
    assert by_mjd(interval_status(data, {"P": 53950.0}, {"P": [52000.0]})) == {
        51000.0: "censored_gap", 52100.0: "completed", 52900.0: "completed", 53950.0: "dropped_end_at_last_glitch"}
    assert by_mjd(interval_status(data, {}, {}))[53950.0] == "dropped_no_end"


def test_add_derived_known_answer():
    data = pd.DataFrame({"dff": [2100.0], "df1": [6.0], "dff_err": [21.0], "df1_err": [0.3]})
    out = add_derived(data).iloc[0]
    assert out["ln_R"] == pytest.approx(np.log(2.1 / 6.0))           # R = (dff/1000)/df1
    assert out["ln_dnu_frac"] == pytest.approx(np.log(2.1e-6))
    assert out["ln_dnudot_frac"] == pytest.approx(np.log(6.0e-3))
    assert out["sigma_lnR"] == pytest.approx(np.sqrt(0.01 ** 2 + 0.05 ** 2))


def _toy(counts):
    # counts: pulsar -> (completed, censored)
    rows = []
    for pulsar, (done, cens) in counts.items():
        rows += [{"pulsar": pulsar, "event": True}] * done + [{"pulsar": pulsar, "event": False}] * cens
    return pd.DataFrame(rows)


def test_inclusion_counts_completed_intervals_only():
    from glitchlib.intervals import apply_inclusion
    kept, dropped = apply_inclusion(_toy({"A": (4, 1), "B": (3, 2), "C": (5, 0)}), 4)
    assert sorted(kept["pulsar"].unique()) == ["A", "C"] and dropped == ["B"]
    assert len(kept) == 10


def test_replicate_kept_with_two_pulsars_discarded_with_one():
    from glitchlib.intervals import apply_replicate_inclusion
    kept, dropped, discarded = apply_replicate_inclusion(_toy({"A": (4, 0), "B": (4, 0), "C": (1, 1)}), 4, 2,
                                                         ["A", "B", "C", "D"])
    assert not discarded and dropped == ["C", "D"] and sorted(kept["pulsar"].unique()) == ["A", "B"]   # D: no interval
    kept, dropped, discarded = apply_replicate_inclusion(_toy({"A": (4, 0), "B": (3, 0)}), 4, 2, ["A", "B"])
    assert discarded and kept is None and dropped == ["B"]


def test_interval_spanning_a_gap_is_censored_at_the_gap_start():
    # D11: a gap starting at 53000 lies inside the interval 52900 -> 53950, which is censored at 53000.
    intervals, _ = build_intervals(_example(), {"P": 55000.0}, {"P": [53000.0]})
    row = intervals.set_index("mjd").loc[52900.0]
    assert row["t"] == 100.0 and not row["event"] and row["gap_flag"]
    other = intervals.set_index("mjd").loc[[51000.0, 53950.0]]
    assert list(other["t"]) == [1100.0, 1050.0] and not other["gap_flag"].any()


def test_gap_in_last_interval_censors_before_obs_end():
    intervals, _ = build_intervals(_example(), {"P": 55000.0}, {"P": [54500.0]})
    last = intervals.set_index("mjd").loc[53950.0]
    assert last["t"] == 550.0 and not last["event"] and last["gap_flag"]


def test_end_epoch_error_comes_from_the_glitch_that_ends_the_interval():
    data = _example().assign(mjd_err=[2.0, 1.0, 3.0, 45.0])      # rows: 52900, 51000, 53950, 52100
    intervals, _ = build_intervals(data, {"P": 55000.0}, {})
    by_start = intervals.set_index("mjd")
    assert by_start.loc[51000.0, "mjd_err"] == 1.0 and by_start.loc[51000.0, "mjd_err_end"] == 45.0
    assert by_start.loc[52900.0, "mjd_err_end"] == 3.0
    assert np.isnan(by_start.loc[53950.0, "mjd_err_end"])            # censored: no ending glitch


def test_completed_flags_use_order_only():
    from glitchlib.intervals import completed_flags
    data = _example()
    flags = completed_flags(data, {"P": [53000.0]})
    by_mjd = dict(zip(data["mjd"], flags))
    assert by_mjd == {51000.0: True, 52100.0: True, 52900.0: False, 53950.0: False}


def test_epoch_flags_known_answer():
    from glitchlib.intervals import epoch_flags
    data = pd.DataFrame({"mjd_err": [1.0, 40.0, 2.0, np.nan, 3.0],
                         "mjd_err_end": [2.0, 1.0, 31.0, 1.0, np.nan],
                         "event": [True, True, True, True, False]})
    assert list(epoch_flags(data, 30, True)) == [False, True, True, True, False]
    assert list(epoch_flags(data, 30, False)) == [False, True, True, False, False]


def test_no_usable_glitch_gives_an_empty_table_with_columns():
    data = _example().assign(usable=False, dff_err=1.0, df1_err=0.1)
    intervals, _ = build_intervals(data, {"P": 55000.0}, {})
    assert intervals.empty and {"t", "event", "prev_t", "dff", "df1"} <= set(intervals.columns)
    assert add_derived(intervals).empty


def test_following_waits_cover_unusable_glitches():
    from glitchlib.intervals import following_waits
    out = following_waits(_example(), {"P": 55000.0}, {})
    by_mjd = out.assign(mjd=_example()["mjd"]).set_index("mjd")
    assert by_mjd.loc[52100.0, "following_t"] == 800.0 and by_mjd.loc[52100.0, "following_completed"]
    assert by_mjd.loc[53950.0, "following_t"] == 1050.0 and not by_mjd.loc[53950.0, "following_completed"]
    dropped = following_waits(_example(), {"P": 53950.0}, {})
    assert np.isnan(dropped.loc[_example()["mjd"] == 53950.0, "following_t"]).all()


def test_epoch_uncertainty_exclusion_known_answer():
    from glitchlib.intervals import epoch_uncertainty_exclusion
    data = pd.DataFrame({"mjd_err": [3.0, 6.0, np.nan, 1.0, 8.0, 2.0],
                         "mjd_err_end": [4.0, 8.0, 1.0, np.nan, np.nan, 1.0],
                         "event": [True, True, True, True, False, True],
                         "t": [60.0, 99.0, 500.0, 500.0, 70.0, 10.0]})
    # 5 <= 6 kept; 10 > 9.9 excluded; start missing; end missing (completed); censored 8 > 7; sqrt(5) > 1.
    assert list(epoch_uncertainty_exclusion(data, 0.1)) == [False, True, True, True, True, True]
