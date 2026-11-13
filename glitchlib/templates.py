"""Simulation templates built from the real pulsar structure (real glitch structure only, never real waiting times; D4b).

Per included pulsar (D13): the real first-glitch date, end of monitoring (data/meta/obs_end.csv), number
of glitches, monitoring gaps, spin-down timescale, the share of its glitches that are admissible, and a
pool of its admissible usable glitches (Δν/ν, Δν̇/ν̇, their errors with D4d fills, method class, the
source cadence from the audit or the D4c fallback, the D4f short-window span w, and Δν in µHz from dff
and ν for D4o). Also the D4c noise
target: the median real df1_err/df1 of the pulsar's usable glitches (filled errors excluded), or the
median over all simulated pulsars if it has fewer than D4c.emulation.noise_min_error_values_per_pulsar.
Real waiting times are never read: only glitch epochs, counts and the documented end dates are used.
"""
import numpy as np
import pandas as pd

from glitchlib.cleaning import source_cadences
from glitchlib.simulate import PulsarTemplate


def build_templates(glitches, obs_end, gaps, tau_sd, audit, cfg, pulsars, nu_hz):
    """glitches: data/processed/glitches_clean.csv; obs_end: dict pulsar -> mjd; gaps: monitoring_gaps
    table; tau_sd: dict pulsar -> days; audit: source_audit table; pulsars: names to build; nu_hz: dict
    pulsar -> spin frequency (for Δν in µHz = dff 1e-9 ν 1e6, D4o's J0537−6910 threshold).
    Returns (templates, noise targets dict)."""
    cadences = source_cadences(audit)
    fallback = cfg["D4c.emulation.cadence_until_audit_days"]
    sw = cfg["D4f.short_window.sources"]
    templates, targets, real_rel = [], {}, {}
    for name in pulsars:
        g = glitches[glitches["pulsar"] == name].sort_values("mjd", kind="mergesort")
        admissible = g[g["admissible"]]
        if admissible.empty:
            raise ValueError(f"{name}: no admissible glitch to resample.")
        pool = pd.DataFrame({
            "dnu_frac": admissible["dff"].to_numpy(float) * 1e-9,
            "dnu_frac_err": admissible["dff_err"].to_numpy(float) * 1e-9,
            "dnudot_frac": admissible["df1"].to_numpy(float) * 1e-3,
            "dnudot_frac_err": admissible["df1_err"].to_numpy(float) * 1e-3,
            "method_class": admissible["method_class"].to_numpy(object),
            "cadence_days": [cadences.get(k, fallback) for k in admissible["source_key"]],
            "short_window_days": [float(sw[k]["w_days"]) if k in sw else np.nan for k in admissible["source_key"]],
            "dnu_microhz": admissible["dff"].to_numpy(float) * 1e-9 * float(nu_hz[name]) * 1e6})
        end = float(obs_end[name])
        last = float(g["mjd"].max())
        if end < last:
            raise ValueError(f"{name}: obs_end {end} is before its last glitch {last}.")
        rows = gaps[gaps["name"] == name].dropna(subset=["gap_start_mjd", "gap_end_mjd"])
        templates.append(PulsarTemplate(
            name, float(g["mjd"].min()), end, len(g), pool, float(tau_sd[name]),
            tuple(zip(rows["gap_start_mjd"].astype(float), rows["gap_end_mjd"].astype(float))),
            end == last, len(admissible) / len(g)))
        usable = g[g["usable"] & ~g["df1_err_imputed"].astype(bool)]
        real_rel[name] = (usable["df1_err"] / usable["df1"]).to_numpy(float)
    minimum = cfg["D4c.emulation.noise_min_error_values_per_pulsar"]
    global_median = float(np.median(np.concatenate(list(real_rel.values()))))
    for name, rel in real_rel.items():
        targets[name] = float(np.median(rel)) if len(rel) >= minimum else global_median
    return templates, targets
