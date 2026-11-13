"""Synthetic glitch sequences for the validation and power scenarios (D14, D4b, D4c, D4d).

Each pulsar is simulated inside its real observing window [first glitch, obs_end]. Allowed real inputs
(real glitch structure only, never real waiting times; the D4b exception): glitch count, first-glitch date, end of monitoring, the
Δν/ν and Δν̇/ν̇ values, errors, measurement methods and source cadences of its usable glitches, and
its spin-down timescale. Real waiting times are never used.

    W_p = (obs_end - first glitch) / number of glitches                       typical wait (D4b)
    ln T = c_p + slope * (x - mean_p x) + sigma * eps,  c_p = ln W_p - sigma^2/2 - ln m_p   (D4h)
        m_p = mean over the pool of exp(slope * (x - mean_p x)); 1 without a covariate. The mean wait is W_p.
        S1, S1m: x = ln R (slope 1)   S2, S4: no covariate   S3: x = ln dnu_frac (slope 1)
Glitch properties are resampled with replacement from the pulsar's real usable glitches. The first
glitch sits on the real first-glitch date; glitches follow until the next would fall after obs_end,
whose interval is then censored there.

Measured Δν/ν = true * exp(e), SD(e) = real dff_err/dff of the resampled glitch, in every layer (D4d).
Two layers for Δν̇ (D4c):
  statistics (S1, S2, S3): measured Δν̇ = true * exp(e), SD(e) = real df1_err/df1. No emulation.
  emulated (S1m with t_rec = T; S4 with t_rec ~ lognormal(median W_p, log-sd 0.3), independent of T):
      a synthetic ν̇ series is measured by glitchlib.measure with each glitch's method; each window
      is sampled at the cadence of the source of the glitch that starts it (D4d). Each glitch's
      exponential recovery (Q, τ_d) is drawn log-uniform within its D4o group (J0537−6910, the Crab,
      all others split at Δν/ν = 1e-7); J0537−6910 glitches below D4o's Δν get none (Q = 0).
Either way the catalogue goes through the same interval code as the real data.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from glitchlib import measure
from glitchlib.cleaning import fill_relative_errors
from glitchlib.intervals import add_derived, build_intervals

POOL_COLUMNS = ["dnu_frac", "dnu_frac_err", "dnudot_frac", "dnudot_frac_err", "method_class", "cadence_days"]


@dataclass
class PulsarTemplate:
    name: str
    first_mjd: float
    obs_end_mjd: float
    n_glitches: int          # all real glitches of the pulsar, usable or not
    pool: pd.DataFrame       # real usable glitches, columns POOL_COLUMNS (+ error-imputed flags)
    tau_sd_days: float
    gaps: tuple              # monitoring gaps ((start, end), ...): glitches inside one go unobserved (D11)
    drop_censored_last: bool  # D4g: obs_end = last real glitch epoch, so the censored last interval is dropped
    admissible_fraction: float  # share of the real glitches that are admissible; each simulated glitch is
                                # marked admissible with this probability (1.0: every glitch)

    @property
    def typical_wait(self):
        return typical_wait(self.first_mjd, self.obs_end_mjd, self.n_glitches)


@dataclass
class Scenario:
    name: str
    covariate: str | None    # None, "ln_R" or "ln_dnu_frac"
    slope: float
    emulated: bool
    t_rec: str | None        # emulated only: "equals_T" or "lognormal_W"
    t_rec_log_sd: float | None


@dataclass(frozen=True)
class ExponentialGroup:
    name: str
    Q: tuple                         # log-uniform range
    tau_d_days: tuple                # log-uniform range
    min_dnu_microhz: float | None    # exponential only for glitches with Δν at least this (None: always)


@dataclass(frozen=True)
class ExponentialSpec:
    """Which group, and so which Q and τ_d ranges, each simulated glitch's exponential recovery uses."""
    named: dict                      # pulsar -> ExponentialGroup
    others_split: float | None       # Δν/ν split for every other pulsar (None: one range)
    others_at_least: ExponentialGroup
    others_below: ExponentialGroup


@dataclass
class EmulationSettings:
    pre_glitch_days: float
    min_samples_after: int
    min_samples_before_window: int
    min_aic_drop: float      # exponential detection threshold (D4d)
    exponentials: ExponentialSpec   # D4o groups (the D14 placeholder only for report-only comparisons)
    unknown_as: str          # how a glitch of class "unknown" is emulated: "window" or "global"
    short_window_days: float | None   # w for the single-method short-window variant (D4f); None otherwise,
                                      # when short-window glitches take w from their pool row


def typical_wait(first_mjd, obs_end_mjd, n_glitches):
    return (obs_end_mjd - first_mjd) / n_glitches


def scenario_from_config(cfg, name):
    """Scenario settings from D4b and D4c."""
    sim = "D4b.simulation"
    if name in ("S1", "S2", "S3"):
        if name not in cfg["D4c.layers.statistics_only"]:
            raise ValueError(f"D4c: {name} is not listed as a statistics-only scenario.")
        covariate = cfg[f"{sim}.{name}.covariate"]
        covariate = None if covariate == "none" else covariate
        slope = cfg[f"{sim}.{name}.slope"] if covariate else 0.0
        return Scenario(name, covariate, slope, False, None, None)
    if name == "S1m":
        return Scenario("S1m", cfg[f"{sim}.S1.covariate"], cfg[f"{sim}.S1.slope"], True,
                        cfg["D4c.layers.S1m_t_rec"], None)
    if name == "S3m":
        # D4m: S3's size-only waits through the measurement emulation, as for S4-mix (reading: S4's t_rec).
        if cfg["D4m.S3m.t_rec"] != "as_S4_mix" or cfg["D4m.S3m.base"] != "S3":
            raise ValueError("D4m: only S3m = S3 with S4-mix's emulation is implemented.")
        s3, s4 = scenario_from_config(cfg, "S3"), scenario_from_config(cfg, "S4")
        return Scenario("S3m", s3.covariate, s3.slope, True, s4.t_rec, s4.t_rec_log_sd)
    if name == "S4":
        if cfg[f"{sim}.S4_t_rec.median"] != "W_p" or cfg["D4c.layers.S4_t_rec"] != "independent_of_T":
            raise ValueError("D4b/D4c: S4 t_rec must be lognormal with median W_p, independent of T.")
        return Scenario("S4", None, 0.0, True, "lognormal_W", cfg[f"{sim}.S4_t_rec.log_sd"])
    raise ValueError(f"Unknown scenario {name!r}.")


def emulation_from_config(cfg, unknown_as):
    """Emulation settings from D4c, D4d and D14. unknown_as is D4c.mixtures.unknown_emulated_as for the
    all-data mixture or D4c.mixtures.unknown_sensitivity_emulated_as for its sensitivity run."""
    em = "D4c.emulation"
    allowed = (cfg["D4c.mixtures.unknown_emulated_as"], cfg["D4c.mixtures.unknown_sensitivity_emulated_as"])
    if unknown_as not in allowed:
        raise ValueError(f"D4c: unknown_as must be one of {allowed}.")
    if cfg["D4o.detection"] != "D4d_aic_only" or cfg["D4o.draw"] != "log_uniform_within_group":
        raise ValueError("D4o: only log-uniform draws within groups with D4d's AIC detection are implemented.")
    return EmulationSettings(cfg[f"{em}.pre_glitch_monitoring_days"], cfg[f"{em}.min_samples_after"],
                             cfg[f"{em}.min_samples_before_window"], cfg["D4d.exponential_detection.min_aic_drop"],
                             exponential_spec_from_config(cfg), unknown_as, None)


def exponential_spec_from_config(cfg):
    """D4o: J0537−6910 (exponential only if Δν >= its threshold), the Crab, and every other included pulsar
    with ranges split at Δν/ν = D4o's split."""
    g = "D4o.groups"
    named = {}
    thresholds = {"J0537-6910": cfg[f"{g}.J0537-6910.exponential_only_if_dnu_microhz_at_least"], "crab": None}
    for key, threshold in thresholds.items():
        group = ExponentialGroup(key, tuple(cfg[f"{g}.{key}.Q"]), tuple(cfg[f"{g}.{key}.tau_d_days"]), threshold)
        named.update({pulsar: group for pulsar in cfg[f"{g}.{key}.pulsars"]})
    if cfg[f"{g}.others.pulsars"] != "all_other_included":
        raise ValueError("D4o: the third group must be every other included pulsar.")
    side = lambda name, key: ExponentialGroup(name, tuple(cfg[f"{g}.others.{key}.Q"]),
                                              tuple(cfg[f"{g}.others.{key}.tau_d_days"]), None)
    return ExponentialSpec(named, cfg[f"{g}.others.dnu_frac_split"], side("others_at_least_split", "at_least_split"),
                           side("others_below_split", "below_split"))


def placeholder_exponential_spec(cfg):
    """The superseded D14 placeholder (one range for every glitch), for report-only comparisons with the
    placeholder record; never used by the analysis."""
    q, tau = "D14.exponentials.Q", "D14.exponentials.tau_d_days"
    group = ExponentialGroup("D14_placeholder", (cfg[f"{q}.low"], cfg[f"{q}.high"]),
                             (cfg[f"{tau}.low"], cfg[f"{tau}.high"]), None)
    return ExponentialSpec({}, None, group, group)


def exponential_ranges(seq, pulsar, spec):
    """Per glitch of one pulsar's sequence: group name, Q range, τ_d range, and whether it has an
    exponential at all (D4o: J0537−6910 needs Δν >= the threshold, from the pool's dnu_microhz)."""
    n = len(seq)
    if pulsar in spec.named:
        group = spec.named[pulsar]
        groups = [group] * n
        if group.min_dnu_microhz is not None:
            if "dnu_microhz" not in seq.columns:
                raise ValueError(f"{pulsar}: D4o's Δν threshold needs the pool's dnu_microhz (from dff and ν).")
            present = seq["dnu_microhz"].to_numpy(float) >= group.min_dnu_microhz
        else:
            present = np.ones(n, bool)
    else:
        if spec.others_split is None:
            groups = [spec.others_at_least] * n
        else:
            at_least = seq["dnu_frac"].to_numpy(float) >= spec.others_split
            groups = [spec.others_at_least if a else spec.others_below for a in at_least]
        present = np.ones(n, bool)
    q = np.array([gr.Q for gr in groups], float).reshape(n, 2)
    tau = np.array([gr.tau_d_days for gr in groups], float).reshape(n, 2)
    return np.array([gr.name for gr in groups], dtype=object), q, tau, present


def fill_missing_errors(pools, min_values):
    """D4d: fill missing errors of the real glitch pools (glitchlib.cleaning.fill_relative_errors, the
    same rule as for the real data). pools: dict pulsar -> DataFrame with dnu_frac, dnu_frac_err,
    dnudot_frac, dnudot_frac_err. Returns a new dict with flags dnu_frac_err_imputed, dnudot_frac_err_imputed."""
    table = pd.concat([pool.assign(pulsar=name) for name, pool in pools.items()], keys=list(pools))
    filled = fill_relative_errors(table, [("dnu_frac", "dnu_frac_err"), ("dnudot_frac", "dnudot_frac_err")],
                                  np.ones(len(table), bool), min_values)
    return {name: filled.loc[name].drop(columns="pulsar") for name in pools}


def _covariate(pool, name):
    if name == "ln_R":
        return np.log(pool["dnu_frac"]) - np.log(pool["dnudot_frac"]) + np.log(1e3)
    if name == "ln_dnu_frac":
        return np.log(pool["dnu_frac"])
    raise ValueError(f"Unknown covariate {name!r}.")


def _check_pool(template):
    missing = [c for c in POOL_COLUMNS if c not in template.pool.columns]
    if missing:
        raise ValueError(f"{template.name}: pool lacks columns {missing}.")
    values = template.pool[["dnu_frac", "dnu_frac_err", "dnudot_frac", "dnudot_frac_err", "cadence_days"]]
    values = values.to_numpy(float)
    if not np.all(np.isfinite(values)) or np.any(values[:, [0, 2, 4]] <= 0):
        raise ValueError(f"{template.name}: pool has missing or non-positive values, errors or cadences "
                         "(fill missing errors with fill_missing_errors first).")


def simulate_sequence(template, scenario, sigma, rng):
    """True glitch sequence: epochs, resampled properties, waits and (emulated layer) recovery times."""
    _check_pool(template)
    pool = template.pool.reset_index(drop=True)
    c_p = np.log(template.typical_wait) - sigma ** 2 / 2
    centred = None
    if scenario.covariate:
        x = _covariate(pool, scenario.covariate).to_numpy()
        centred = x - x.mean()
        c_p -= np.log(np.mean(np.exp(scenario.slope * centred)))     # D4h: ln m_p
    rows, t = [], template.first_mjd
    while t < template.obs_end_mjd:
        i = rng.integers(len(pool))
        ln_t = c_p + sigma * rng.normal()
        if centred is not None:
            ln_t += scenario.slope * centred[i]
        T = np.exp(ln_t)
        row = {"pulsar": template.name, "mjd": t, "pool_index": i, "T_true": T, "t_rec": np.nan}
        if scenario.emulated:
            if scenario.t_rec == "equals_T":
                row["t_rec"] = T
            elif scenario.t_rec == "lognormal_W":
                row["t_rec"] = np.exp(np.log(template.typical_wait) + scenario.t_rec_log_sd * rng.normal())
            else:
                raise ValueError(f"Unknown t_rec rule {scenario.t_rec!r}.")
        if template.admissible_fraction < 1.0:
            row["admissible"] = bool(rng.random() < template.admissible_fraction)
        else:
            row["admissible"] = True
        rows.append(row)
        t += T
    seq = pd.DataFrame(rows)
    observed = np.ones(len(seq), bool)
    for start, end in template.gaps:                      # no monitoring, so no glitch is seen
        observed &= ~((seq["mjd"] > start) & (seq["mjd"] < end)).to_numpy()
    extra = [c for c in ("short_window_days", "dnu_microhz") if c in pool.columns]
    return seq[observed].reset_index(drop=True).join(pool[POOL_COLUMNS + extra], on="pool_index")


def _lognormal_noise(true, err, rng, noise=True):
    """Measured = true * exp(e), SD(e) = err/true; reported error = the same relative error.
    noise=False (D14: measurement noise switched off) returns the true values with zero errors."""
    if not noise:
        return true.copy(), np.zeros_like(true)
    rel = err / true
    measured = true * np.exp(rel * rng.normal(size=len(true)))
    return measured, rel * measured


def _catalogue_columns(out, dnu, dnu_err, dnudot, dnudot_err):
    out["dff"], out["dff_err"] = dnu * 1e9, dnu_err * 1e9
    out["df1"], out["df1_err"] = dnudot * 1e3, dnudot_err * 1e3
    out["usable"] = (out["dff"] > 0) & (out["df1"] > 0)       # NaN compares False: unmeasured = not usable
    return out


def _measured_dnu(seq, rng):
    return _lognormal_noise(seq["dnu_frac"].to_numpy(float), seq["dnu_frac_err"].to_numpy(float), rng)


def measure_statistics_layer(seq, rng, noise):
    """S1-S3: Δν/ν and Δν̇/ν̇ each times exp(e) with the resampled glitch's real relative errors
    (noise=False: no measurement noise, for the D14 noise-free check)."""
    dnu, dnu_err = _lognormal_noise(seq["dnu_frac"].to_numpy(float), seq["dnu_frac_err"].to_numpy(float), rng,
                                    noise)
    dnudot, dnudot_err = _lognormal_noise(seq["dnudot_frac"].to_numpy(float),
                                          seq["dnudot_frac_err"].to_numpy(float), rng, noise)
    return _catalogue_columns(seq.copy(), dnu, dnu_err, dnudot, dnudot_err)


def emulation_methods(seq, method_rule, unknown_as):
    """Method class (label) and emulated method per glitch.

    method_rule: "resampled" (D4c: each glitch takes the method of the real glitch it was resampled
    from) or one method name for every glitch (the single-method S4 variants)."""
    if method_rule == "resampled":
        classes = seq["method_class"].to_numpy(object)
    elif method_rule in measure.METHODS:
        classes = np.full(len(seq), method_rule, dtype=object)
    else:
        raise ValueError(f"Unknown method rule {method_rule!r}.")
    emulated = np.where(classes == "unknown", unknown_as, classes)
    bad = set(emulated) - set(measure.METHODS)
    if bad:
        raise ValueError(f"Cannot emulate method classes {sorted(bad)}.")
    return classes, emulated


def window_cadences(seq):
    """D4d: window j >= 1 (after glitch j-1) uses that glitch's source cadence. Window 0 (the pre-glitch
    baseline before the first glitch) has no previous interval and uses the first glitch's cadence."""
    cad = seq["cadence_days"].to_numpy(float)
    return np.concatenate([[cad[0]], cad])


def _emulation_draws(seq, template, method_rule, settings, rng):
    """Everything random about one sequence's emulation except the ν̇ noise level."""
    n = len(seq)
    draws = seq.copy()
    groups, q, tau, present = exponential_ranges(seq, template.name, settings.exponentials)
    draws["tau_d"] = np.exp(rng.uniform(np.log(tau[:, 0]), np.log(tau[:, 1]), n))          # D4o (log-uniform)
    draws["Q"] = np.where(present, np.exp(rng.uniform(np.log(q[:, 0]), np.log(q[:, 1]), n)), 0.0)
    draws["exp_group"], draws["exp_present"] = groups, present
    draws["A"] = measure.exponential_amplitude(draws["Q"], draws["dnu_frac"], draws["tau_d"], template.tau_sd_days)
    draws["method_class"], draws["emulated_method"] = emulation_methods(seq, method_rule, settings.unknown_as)
    if method_rule == "short-window":
        draws["short_window_days"] = settings.short_window_days
    elif "short_window_days" not in draws.columns:
        draws["short_window_days"] = np.nan
    t = measure.sample_times(template.first_mjd - settings.pre_glitch_days, template.obs_end_mjd,
                             draws["mjd"].to_numpy(), window_cadences(seq))
    clean = measure.true_signal(t, draws["mjd"].to_numpy(), draws["dnudot_frac"].to_numpy(),
                                draws["t_rec"].to_numpy(), draws["A"].to_numpy(), draws["tau_d"].to_numpy())
    return draws, t, clean, rng.normal(size=t.size)


def _measure_draws(draws, t, y, settings):
    value = np.full(len(draws), np.nan)
    se = np.full(len(draws), np.nan)
    fitted = np.zeros(len(draws), bool)
    for method in pd.unique(draws["emulated_method"]):
        pick = (draws["emulated_method"] == method).to_numpy()
        w = np.where(pick, draws["short_window_days"].to_numpy(float), 1.0) if method == "short-window" else None
        v, s, f = measure.measure(method, t, y, draws["mjd"].to_numpy(), draws["tau_d"].to_numpy(),
                                  settings.min_samples_after, settings.min_samples_before_window,
                                  settings.min_aic_drop, w)
        pick = (draws["emulated_method"] == method).to_numpy()
        value[pick], se[pick], fitted[pick] = v[pick], s[pick], f[pick]
    return value, se, fitted


def measure_emulated_layer(seq, template, method_rule, settings, noise_sd, rng):
    """S1m and S4: build the synthetic ν̇ series, measure each glitch by its method."""
    draws, t, clean, z = _emulation_draws(seq, template, method_rule, settings, rng)
    value, se, fitted = _measure_draws(draws, t, clean + noise_sd * z, settings)
    draws["exp_fitted"] = fitted
    dnu, dnu_err = _measured_dnu(seq, rng)
    return _catalogue_columns(draws, dnu, dnu_err, value, se)


def _median_relative_error(value, se):
    ok = np.isfinite(value) & (value > 0)
    return float(np.median(se[ok] / value[ok])) if ok.any() else np.nan


def target_relative_errors(templates, min_values):
    """D4c: per pulsar, the median real df1_err/df1 of its usable glitches (imputed errors excluded);
    pulsars with fewer than min_values real error values use the median over all pulsars."""
    def real_rel(tp):
        rel = (tp.pool["dnudot_frac_err"] / tp.pool["dnudot_frac"]).to_numpy(float)
        if "dnudot_frac_err_imputed" in tp.pool.columns:
            rel = rel[~tp.pool["dnudot_frac_err_imputed"].to_numpy(bool)]
        return rel[np.isfinite(rel)]
    global_median = float(np.median(np.concatenate([real_rel(tp) for tp in templates])))
    return {tp.name: (float(np.median(real_rel(tp))) if len(real_rel(tp)) >= min_values else global_median)
            for tp in templates}


def calibrate_noise(template, scenario, sigma, method_rule, settings, target, n_sequences, tolerance, rng):
    """D4c/D4d/D4e: choose the ν̇ sample noise by bisection so the median REPORTED (formal) relative
    error of the emulated Δν̇ (reported SE / measured value, over usable glitches; not the scatter of
    measured against true values) is within tolerance * target of target. D4e: the caller calibrates
    under S1m (no model misfit) at each sigma and uses the result in every emulated scenario.

    The same n_sequences sequences and noise draws are re-measured at every trial noise level (common
    random numbers), so each bisection step uses n_sequences sequences. Returns a dict with the noise,
    the achieved median, whether the tolerance was met, and whether zero noise already exceeds the
    target (the fitted model's misfit alone; noise is then 0).
    """
    trials = []
    for _ in range(n_sequences):
        seq = simulate_sequence(template, scenario, sigma, rng)
        trials.append(_emulation_draws(seq, template, method_rule, settings, rng))

    def median_at(noise_sd):
        values, ses = [], []
        for draws, t, clean, z in trials:
            v, s, _ = _measure_draws(draws, t, clean + noise_sd * z, settings)
            values.append(v)
            ses.append(s)
        return _median_relative_error(np.concatenate(values), np.concatenate(ses))

    def result(noise, achieved, at_floor):
        return {"pulsar": template.name, "noise_sd": noise, "target": target, "achieved": achieved,
                "within_tolerance": bool(abs(achieved - target) <= tolerance * target), "at_floor": at_floor,
                "n_sequences": n_sequences}

    floor = median_at(0.0)
    if not np.isfinite(floor):
        raise RuntimeError(f"{template.name}: no usable emulated glitch to calibrate the noise on.")
    if floor >= target:
        return result(0.0, floor, True)
    scale = float(np.median(template.pool["dnudot_frac"]))
    lo, hi = 0.0, scale
    while median_at(hi) < target:
        lo, hi = hi, 2.0 * hi
        if hi > 1e6 * scale:
            raise RuntimeError(f"{template.name}: noise calibration cannot reach target {target}.")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        achieved = median_at(mid)
        if abs(achieved - target) <= tolerance * target:
            return result(mid, achieved, False)
        lo, hi = (mid, hi) if achieved < target else (lo, mid)
    return result(mid, achieved, False)


def simulate_statistics_dataset(templates, scenario, sigma, noise, rng):
    """One replicate of S1, S2 or S3 -> (intervals, glitches). noise=False switches measurement noise off."""
    if scenario.emulated:
        raise ValueError(f"{scenario.name} is an emulated scenario.")
    tables = [measure_statistics_layer(simulate_sequence(tp, scenario, sigma, rng), rng, noise) for tp in templates]
    return _to_intervals(tables, templates)


def simulate_emulated_dataset(templates, scenario, sigma, method_rule, settings, noise_sd, rng):
    """One replicate of S1m or an S4 variant -> (intervals, glitches). noise_sd: pulsar -> ν̇ noise."""
    if not scenario.emulated:
        raise ValueError(f"{scenario.name} is a statistics-only scenario.")
    tables = [measure_emulated_layer(simulate_sequence(tp, scenario, sigma, rng), tp, method_rule, settings,
                                     noise_sd[tp.name], rng) for tp in templates]
    return _to_intervals(tables, templates)


def _to_intervals(tables, templates):
    """Same interval code as the real data: D11 gaps censor, and D4g drops the censored last interval of
    pulsars whose real obs_end is their last glitch epoch."""
    glitches = pd.concat(tables, ignore_index=True)
    gaps = {tp.name: [start for start, _ in tp.gaps] for tp in templates if tp.gaps}
    intervals, _ = build_intervals(glitches, {tp.name: tp.obs_end_mjd for tp in templates}, gaps)
    drop = [tp.name for tp in templates if tp.drop_censored_last]
    last_censored = (~intervals["event"].astype(bool)) & (~intervals["gap_flag"].astype(bool))
    intervals = intervals[~(intervals["pulsar"].isin(drop) & last_censored)].reset_index(drop=True)
    return add_derived(intervals), glitches


def short_window_admission(entries, max_upper, z_95):
    """D4f admission rule. entries: one dict per (w, sigma) with the S4sw forward-slope mean and its
    standard error across replicates. The upper bound is mean + z_95 * SE; short-window is admitted
    only if every bound is <= max_upper. Returns the verdict and the bounds."""
    rows = []
    for e in entries:
        # No estimate (e.g. no glitch measurable at this w): the bound cannot be shown, so the check fails.
        upper = None if e["mean"] is None else e["mean"] + z_95 * e["se_mean"]
        rows.append(dict(e, upper_95=upper, passes=upper is not None and bool(upper <= max_upper)))
    return {"admitted": all(r["passes"] for r in rows), "max_upper_95": max_upper, "checks": rows}
