"""Glitches -> waiting-time intervals, shared by the real-data scripts and the simulations.

Within each pulsar, glitches are ordered by mjd (D13). Each usable glitch starts an interval that ends
at the NEXT glitch of the same pulsar, usable or not: T = mjd_next - mjd, completed. The last glitch's
interval is right-censored at the pulsar's end of monitoring: c = obs_end - mjd. If obs_end is unknown
(D10) or equals the last glitch epoch (D4g: no censoring information), that last interval is dropped
and logged, never entered with c = 0. The previous wait (for the D15 placebo) is mjd - mjd_previous,
always completed, and missing for a pulsar's first glitch.
Monitoring gaps (D11): an interval that spans a listed gap (the gap starts after the interval's
start and before its end) is censored at the gap start and flagged.
"""
import numpy as np
import pandas as pd

DROP_REASONS = {"dropped_no_end": "no end-of-monitoring date (D10)",
                "dropped_end_at_last_glitch": "obs_end equals the last glitch epoch: no censoring information (D4g)"}


def _first_gap_inside(gaps, start, end):
    """Earliest gap start strictly inside (start, end), or None."""
    inside = [g for g in gaps if start < g < end]
    return min(inside) if inside else None


def _interval_end(pulsar, mjd, i, end, gaps):
    """Where the interval started by glitch i stops, and how: (stop epoch or None, status). Status is
    completed, censored_end, censored_gap, dropped_no_end or dropped_end_at_last_glitch. Uses
    comparisons only."""
    if i + 1 < len(mjd):
        stop, status = mjd[i + 1], "completed"
    elif not np.isfinite(end):
        return None, "dropped_no_end"
    elif end == mjd[i]:
        return None, "dropped_end_at_last_glitch"
    elif end < mjd[i]:
        raise ValueError(f"{pulsar}: end of monitoring {end} is before its last glitch {mjd[i]}.")
    else:
        stop, status = end, "censored_end"
    gap = _first_gap_inside(gaps, mjd[i], stop)
    return (gap, "censored_gap") if gap is not None else (stop, status)


def interval_status(glitches, obs_end, gaps):
    """Per glitch (same index): the status of the interval it would start (see _interval_end).
    obs_end: dict pulsar -> end-of-monitoring mjd; gaps: dict pulsar -> list of gap-start mjds."""
    status = pd.Series("", index=glitches.index, dtype=object)
    for pulsar, g in glitches.groupby("pulsar", sort=True):
        g = g.sort_values("mjd", kind="mergesort")
        mjd = g["mjd"].to_numpy(float)
        end = obs_end.get(pulsar, np.nan)
        for i, idx in enumerate(g.index):
            status[idx] = _interval_end(pulsar, mjd, i, end, gaps.get(pulsar, []))[1]
    return status


def completed_flags(glitches, gaps):
    """Per glitch (same index): does it have a next glitch of the same pulsar with no listed gap
    starting in between? Uses ordering and comparisons only, no time differences."""
    return interval_status(glitches, {}, gaps) == "completed"


def build_intervals(glitches, obs_end, gaps):
    """glitches: DataFrame with pulsar, mjd, usable and any per-glitch columns (carried along; an
    mjd_err column, if present, also gives mjd_err_end, the error of the glitch that ends the interval).
    obs_end: dict pulsar -> end-of-monitoring mjd (NaN or absent = unknown).
    gaps: dict pulsar -> list of gap-start mjds (D11).

    Returns (intervals, dropped): intervals started by usable glitches, and a table of last intervals
    dropped (D10: no end-of-monitoring date; D4g: obs_end equals the last glitch epoch).
    """
    rows, dropped = [], []
    for pulsar, g in glitches.groupby("pulsar", sort=True):
        g = g.sort_values("mjd", kind="mergesort").reset_index(drop=True)
        mjd = g["mjd"].to_numpy(float)
        mjd_err = g["mjd_err"].to_numpy(float) if "mjd_err" in g.columns else np.full(len(g), np.nan)
        end = obs_end.get(pulsar, np.nan)
        for i in range(len(g)):
            if not bool(g.loc[i, "usable"]):
                continue
            stop, status = _interval_end(pulsar, mjd, i, end, gaps.get(pulsar, []))
            if stop is None:
                dropped.append({"pulsar": pulsar, "mjd": mjd[i], "reason": DROP_REASONS[status]})
                continue
            row = g.loc[i].to_dict()
            row["prev_t"] = mjd[i] - mjd[i - 1] if i > 0 else np.nan
            row["prev_event"] = True if i > 0 else np.nan
            row["event"] = status == "completed"
            row["gap_flag"] = status == "censored_gap"
            row["mjd_err_end"] = mjd_err[i + 1] if status == "completed" else np.nan
            row["t"] = stop - mjd[i]
            rows.append(row)
    added = ["prev_t", "prev_event", "event", "mjd_err_end", "gap_flag", "t"]
    intervals = pd.DataFrame(rows, columns=list(dict.fromkeys(list(glitches.columns) + added)))
    if len(intervals) and (intervals["t"] <= 0).any():
        bad = intervals.loc[intervals["t"] <= 0, ["pulsar", "mjd"]].to_dict("records")
        raise ValueError(f"Intervals of zero or negative length (glitches on the same mjd): {bad}.")
    return intervals, pd.DataFrame(dropped, columns=["pulsar", "mjd", "reason"])


def following_waits(glitches, obs_end, gaps):
    """D15 missingness: for EVERY glitch (usable or not), the wait to the next glitch and whether it is
    completed, by the same rule as build_intervals; NaN where the last interval is dropped (D10, D4g).
    Uses real waiting times: call only after the guard (glitchlib.guards)."""
    t = pd.Series(np.nan, index=glitches.index)
    completed = pd.Series(False, index=glitches.index)
    for pulsar, g in glitches.groupby("pulsar", sort=True):
        g = g.sort_values("mjd", kind="mergesort")
        mjd = g["mjd"].to_numpy(float)
        end = obs_end.get(pulsar, np.nan)
        for i, idx in enumerate(g.index):
            stop, status = _interval_end(pulsar, mjd, i, end, gaps.get(pulsar, []))
            if stop is not None:
                t[idx], completed[idx] = stop - mjd[i], status == "completed"
    return pd.DataFrame({"following_t": t, "following_completed": completed})


def epoch_flags(intervals, max_err_days, flag_unknown):
    """D13: flag an interval that touches a glitch whose epoch error exceeds max_err_days (or is unknown,
    if flag_unknown): the glitch that starts it, or the glitch that ends it when it is completed."""
    def bad(err):
        err = err.astype(float)
        return (err > max_err_days) | (err.isna() if flag_unknown else False)
    return bad(intervals["mjd_err"]) | (intervals["event"].astype(bool) & bad(intervals["mjd_err_end"]))


def epoch_uncertainty_exclusion(intervals, max_fraction_of_T):
    """D8 item 9 with D4h: exclude an interval if sqrt(mjd_err_start^2 + mjd_err_end^2) exceeds
    max_fraction_of_T * T, or if an epoch error is missing. A censored interval has no ending glitch,
    so only its start error counts (its missing end error is not an X in the catalogue)."""
    start = intervals["mjd_err"].astype(float)
    end = intervals["mjd_err_end"].astype(float)
    completed = intervals["event"].astype(bool)
    missing = start.isna() | (completed & end.isna())
    combined = np.sqrt(start ** 2 + np.where(completed, end, 0.0) ** 2)
    return missing | (combined > max_fraction_of_T * intervals["t"].astype(float))


def apply_inclusion(intervals, min_completed):
    """D13 inclusion: keep pulsars with at least min_completed completed intervals in the data being
    analysed (for the primary analysis, the caller passes admissible intervals only).
    Returns (kept intervals, sorted list of dropped pulsars)."""
    completed = intervals[intervals["event"].astype(bool)].groupby("pulsar").size()
    keep = set(completed[completed >= min_completed].index)
    dropped = sorted(set(intervals["pulsar"]) - keep)
    return intervals[intervals["pulsar"].isin(keep)], dropped


def apply_replicate_inclusion(intervals, min_completed, min_pulsars, pulsars):
    """D4e: the D13 rule applied to one simulated replicate. Failing pulsars are dropped from the fit;
    the replicate is discarded (kept = None) only if fewer than min_pulsars pulsars remain.
    pulsars: every simulated pulsar, so one with no analysed interval at all also counts as dropped."""
    kept, _ = apply_inclusion(intervals, min_completed)
    remaining = set(kept["pulsar"])
    dropped = sorted(set(pulsars) - remaining)
    discarded = len(remaining) < min_pulsars
    return (None if discarded else kept), dropped, discarded


def add_derived(intervals):
    """Add ln R, ln dnu_frac, ln dnudot_frac and sigma_lnR from catalogue-unit columns
    dff (1e-9), df1 (1e-3), dff_err, df1_err. R = (dff/1000)/df1."""
    out = intervals.copy()
    out["ln_R"] = np.log((out["dff"] / 1000.0) / out["df1"])
    out["ln_dnu_frac"] = np.log(out["dff"] * 1e-9)
    out["ln_dnudot_frac"] = np.log(out["df1"] * 1e-3)
    out["sigma_lnR"] = np.sqrt((out["dff_err"] / out["dff"]) ** 2 + (out["df1_err"] / out["df1"]) ** 2)
    return out
