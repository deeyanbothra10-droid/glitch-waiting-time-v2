"""Mechanical outcome assignment (D4j with D4l), from numbers only.

T1 (D4l) is directional: it rejects when its signed-root LR statistic exceeds the 97.5th percentile of the
same statistic under S2 at the same sigma. T2 (D4n) rejects when d-hat is below the smaller of the 2.5th
percentiles of d-hat under the S3 and S3m null references; T2's signed root is reported alongside and
never changes the outcome. Reference distributions: S1 (catalogue-level noise), S1m (model true, emulated
measurement) and S4-mix (bias only), at the sigma of the analysis (06: the grid sigma; 08: D4n's matched
sigma per scenario). Powers are the S1m rejection rates of these tests. Outcomes are assigned in order; the first
that applies wins:
  1. D_placebo: placebo b-hat below the 2.5th percentile of the S1m placebo b-hat.
  2. C: b-hat below the 2.5th percentile of S1m b-hat (D4l: no power clause).
  3. A: T1 rejects and b-hat is above the 97.5th percentile of S4-mix b-hat, at or above the 2.5th
     percentile of S1m b-hat and at or below the 97.5th percentile of S1 b-hat (D4l). Within A: A2 if T2
     rejects; B if T2 does not reject and T2 power >= 0.8; otherwise A1, labelled (D4p) "association
     consistent with the scaling, not specific to vortex creep; Δν̇ term not established (T2 power at the
     matched σ reported)", with the T2 power reported beside it (no fixed figure). Only A2
     is evidence for the creep-specific Δν̇ dependence (D4m). A2 carries the report-only label
     "short-window-sensitive" if d-hat refitted without short-window glitches is >= 0 (D4n).
  4. D_bias_range: T1 rejects but b-hat is at or below the 97.5th percentile of S4-mix b-hat.
  5. inconclusive: anything else.
Reported alongside, never changing the outcome: b-hat's percentile within S1m and, above S1m's 97.5th
percentile, the label "less attenuated than emulated"; implied beta; the D15 missingness check.
"""
import numpy as np

OUTCOMES = ("D_placebo", "C", "A2", "B", "A1", "D_bias_range", "inconclusive")


def _pct(values, q):
    return float(np.percentile(np.asarray(values, float), q))


def null_threshold(null_signed_roots, cfg):
    """D4l: the rejection threshold, the 97.5th percentile of the signed root under the null scenario."""
    return _pct(null_signed_roots, cfg["D4l.directional_tests.reject_above_null_percentile"])


def t2_d_cut(null_d, cfg):
    """D4n: T2's cut, the smaller of the 2.5th percentiles of d-hat under the null references (S3, S3m).
    null_d: {scenario: d-hat values}. Returns (cut, {scenario: percentile})."""
    t2 = "D4n.T2"
    if cfg[f"{t2}.combine"] != "smaller_percentile":
        raise ValueError(f"D4n.T2.combine = {cfg[f'{t2}.combine']} is not implemented.")
    per_null = {name: _pct(null_d[name], cfg[f"{t2}.rejects_if_d_hat_below_null_percentile"])
                for name in cfg[f"{t2}.null_scenarios"]}
    return min(per_null.values()), per_null


def short_window_label(outcome, d_without_short_window, cfg):
    """D4n item 3 (report only): the label for an A2 whose d-hat refitted without short-window glitches is
    >= 0. None for other outcomes or when the refit was not available (no data left after dropping them)."""
    sw = "D4n.short_window_label"
    if outcome != cfg[f"{sw}.applies_to"] or d_without_short_window is None:
        return None
    return cfg[f"{sw}.label"] if d_without_short_window >= cfg[f"{sw}.if_d_hat_refitted_without_short_window_at_least"] else None


def references(s1m_b, s1m_placebo_b, s4mix_b, s1_b, cfg):
    """Percentile bounds of D4j/D4l from the reference distributions (D4b: percentile intervals)."""
    a = "D4l.A_slope_condition"
    return {"s1m_b_lower": _pct(s1m_b, cfg[f"{a}.at_or_above_S1m_percentile"]),
            "s1m_b_refutation_bound": _pct(s1m_b, cfg["D4l.C_refutation.bhat_below_S1m_percentile"]),
            "s1m_b_label_bound": _pct(s1m_b, cfg[f"{a}.label_if_above_S1m_percentile"]),
            "s1_b_upper": _pct(s1_b, cfg[f"{a}.at_or_below_S1_percentile"]),
            "s4mix_b_support_bound": _pct(s4mix_b, cfg[f"{a}.above_S4mix_percentile"]),
            "s4mix_b_bias_bound": _pct(s4mix_b, cfg["D4j.D_bias_range.bhat_at_or_below_S4mix_percentile"]),
            "s1m_placebo_bound": _pct(s1m_placebo_b, cfg["D4j.D_placebo.placebo_bhat_below_S1m_placebo_percentile"]),
            "s1m_b": np.sort(np.asarray(s1m_b, float)),
            "n": {"s1m": int(len(s1m_b)), "s1m_placebo": int(len(s1m_placebo_b)), "s4mix": int(len(s4mix_b)),
                  "s1": int(len(s1_b))}}


def classify(t1_rejects, t2_rejects, b, placebo_b, refs, t2_power, cfg, d_without_short_window=None):
    """The single outcome of one analysis (real or simulated), with the report-only labels.
    d_without_short_window: d-hat refitted without short-window glitches (D4n item 3), if computed."""
    if placebo_b < refs["s1m_placebo_bound"]:
        outcome = "D_placebo"
    elif b < refs["s1m_b_refutation_bound"]:
        outcome = "C"
    elif (t1_rejects and b > refs["s4mix_b_support_bound"] and refs["s1m_b_lower"] <= b <= refs["s1_b_upper"]):
        if t2_rejects:
            outcome = "A2"
        elif t2_power >= cfg["D4j.A_scaling_supported.within_A.B_min_T2_power"]:
            outcome = "B"
        else:
            outcome = "A1"
    elif t1_rejects and b <= refs["s4mix_b_bias_bound"]:
        outcome = "D_bias_range"
    else:
        outcome = "inconclusive"
    return {"outcome": outcome,
            "outcome_label": cfg["D4p.A1_label"] if outcome == "A1" else None,              # D4p item 4
            "creep_specific_evidence": outcome == "A2",                       # D4m: only A2
            "T2_power": float(t2_power),            # D4p: the power at the analysis sigma, not a fixed figure
            "b_percentile_in_S1m": float(100 * np.mean(refs["s1m_b"] <= b)),
            "attenuation_label": cfg["D4l.A_slope_condition.label"] if b > refs["s1m_b_label_bound"] else None,
            "short_window_label": short_window_label(outcome, d_without_short_window, cfg)}


def rates(outcomes):
    """Share of analyses with each outcome (exclusive), plus A = A1 + A2 + B and D = both D kinds."""
    n = len(outcomes)
    out = {k: float(np.mean([o == k for o in outcomes])) if n else None for k in OUTCOMES}
    if n:
        out["A_total"] = out["A1"] + out["A2"] + out["B"]
        out["D_total"] = out["D_placebo"] + out["D_bias_range"]
    out["n"] = n
    return out


def percentile_of(value, distribution):
    """Where a value falls in a simulated distribution, as a percentile (D4i: missingness vs S1m)."""
    dist = np.asarray(distribution, float)
    dist = dist[np.isfinite(dist)]
    return float(100 * np.mean(dist <= value)) if dist.size else None


def references_for_json(refs):
    """The reference bounds without the sorted S1m sample (kept only for percentiles)."""
    return {k: v for k, v in refs.items() if k != "s1m_b"}
