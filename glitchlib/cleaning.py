"""Cleaning the catalogue (D13, D9), Δν̇ method classes from the source audit (D4b), spin-down timescale.

No time differences between glitches are computed here: the D13 inclusion count uses only the order
of glitches (glitchlib.intervals.completed_flags).
"""
import re

import numpy as np
import pandas as pd

from glitchlib.intervals import completed_flags

AUDIT_METHODS = ("window", "local-relaxation", "global", "short-window", "unknown")
EXP_FITTED = ("yes", "no", "sometimes", "when-detected")
# Spellings of the catalogue maintainers' own measurements, treated as one source.
OUR_MEASUREMENT = "our measurement"
_OUR_MEASUREMENT_SPELLINGS = {"our measurement", "new - our measurement"}


def primary_reference(ref):
    """A glitch's primary reference: the text of its ref before the first "also" (None if missing)."""
    if ref is None or (isinstance(ref, float) and np.isnan(ref)):
        return None
    return re.split(r"[;.]?\s*\balso\b", ref, maxsplit=1, flags=re.I)[0].strip()


def jodrell_bank_pulsars(glitches, cfg):
    """D4g item 10: pulsars with at least one glitch whose primary source is Jodrell Bank (listed source
    keys, or a primary reference starting with the listed prefix). glitches: pulsar and ref columns."""
    rule = "D4g.robustness_item_10.jodrell_bank_rule"
    keys = glitches["ref"].map(source_key)
    primary = glitches["ref"].map(primary_reference).fillna("")
    jb = keys.isin(cfg[f"{rule}.source_keys"]) | primary.str.startswith(cfg[f"{rule}.primary_reference_prefix"])
    return sorted(glitches.loc[jb, "pulsar"].unique())


def source_key(ref):
    """Key of a glitch's PRIMARY reference: the text before the first "also", matched by its [code]
    when it has one, otherwise by the text itself (trimmed, final full stop removed, lower case).
    Every spelling of "(New -) Our measurement" maps to one key, OUR_MEASUREMENT."""
    primary = primary_reference(ref)
    if primary is None:
        return None
    code = re.search(r"\[(.+?)\]", primary)
    if code:
        return code.group(1)
    text = primary.strip().rstrip(".").lower()
    return OUR_MEASUREMENT if text in _OUR_MEASUREMENT_SPELLINGS else text


def audit_method_class(df1_method, exp_fitted):
    """D4b/D4f method class of one audit row. Blank method -> unknown. Global needs exp_fitted:
    yes or when-detected (D4f) -> global (admissible), no -> global-noexp; sometimes or blank -> unknown
    (per-glitch classification from audit notes is not implemented, so such glitches stay unknown)."""
    method = None if df1_method is None or pd.isna(df1_method) else str(df1_method).strip()
    if not method:
        return "unknown"
    if method not in AUDIT_METHODS:
        raise ValueError(f"source_audit df1_method {method!r} is not one of {AUDIT_METHODS}.")
    if method != "global":
        return method
    exp = None if exp_fitted is None or pd.isna(exp_fitted) else str(exp_fitted).strip().lower()
    if exp is not None and exp not in EXP_FITTED:
        raise ValueError(f"source_audit exp_fitted {exp_fitted!r} is not one of {EXP_FITTED}.")
    return {"yes": "global", "when-detected": "global", "no": "global-noexp"}.get(exp, "unknown")


def assigned_classes(cfg):
    """D4f: classes set by the decisions themselves, overriding a blank audit: the short-window
    sources and "Our measurement" (global, assumed)."""
    out = {key: cfg["D4f.short_window.class"] for key in cfg["D4f.short_window.sources"]}
    out[OUR_MEASUREMENT] = cfg["D4f.our_measurement.class"]
    return out


def admissible_methods(cfg, admission):
    """D9 admissible classes, plus short-window if the D4f admission rule passed (admission: the record
    written by scripts/emulation_fixture.py)."""
    methods = list(cfg["D9.admissible_methods"])
    if admission["admitted"]:
        methods.append(cfg["D4f.short_window.class"])
    return methods


def audit_classes(audit, assigned):
    """Method class per source key: from the audit, with the D4f assignments on top. Audit rows sharing
    a key (spelling variants) must agree, and an audit row must not contradict an assignment."""
    rows = audit.assign(key=audit["ref"].map(source_key),
                        method_class=[audit_method_class(m, e) for m, e in zip(audit["df1_method"], audit["exp_fitted"])])
    classes = {}
    for key, g in rows.groupby("key"):
        known = set(g["method_class"]) - {"unknown"}
        if len(known) > 1:
            raise ValueError(f"source_audit rows {list(g['ref'])} share the key {key!r} but disagree: {sorted(known)}.")
        classes[key] = known.pop() if known else "unknown"
    for key, cls in assigned.items():
        if classes.get(key, "unknown") not in ("unknown", cls):
            raise ValueError(f"source_audit gives {classes[key]!r} for {key!r}, but D4f assigns {cls!r}.")
        classes[key] = cls
    return classes


def fill_relative_errors(table, pairs, rows, min_values):
    """D4d: a row without an error uses its pulsar's median relative error (err/value) over the selected
    rows, or the median over all selected rows if the pulsar has fewer than min_values errors.

    table: DataFrame with a pulsar column; pairs: [(value column, error column), ...]; rows: boolean
    mask of the rows that are filled and that define the medians (e.g. usable glitches).
    Returns a copy with the errors filled and a flag column "<error column>_imputed" per pair.
    """
    out = table.copy()
    rows = np.asarray(rows, bool)
    for value, err in pairs:
        rel = (out[err] / out[value]).to_numpy(float)
        known = rows & np.isfinite(rel)
        global_median = float(np.median(rel[known])) if known.any() else np.nan
        fill = pd.Series(global_median, index=out.index)
        for pulsar, idx in out[known].groupby("pulsar").groups.items():
            if len(idx) >= min_values:
                fill[out["pulsar"] == pulsar] = float(np.median(rel[out.index.get_indexer(idx)]))
        missing = rows & ~np.isfinite(rel)
        out[f"{err}_imputed"] = missing
        out.loc[missing, err] = (fill * out[value])[missing]
    return out


def source_cadences(audit):
    """D4c: emulation cadence per source key, from the audit's cadence_num_days (sources without one
    use the D4c fallback, applied by the caller). Spelling variants of one source must agree."""
    rows = audit.assign(key=audit["ref"].map(source_key)).dropna(subset=["cadence_num_days"])
    cadences = {}
    for key, g in rows.groupby("key"):
        values = set(g["cadence_num_days"].astype(float))
        if len(values) > 1:
            raise ValueError(f"source_audit rows {list(g['ref'])} share the key {key!r} but give cadences {sorted(values)}.")
        cadences[key] = values.pop()
    return cadences


def tau_sd_days(nu_hz, nudot_hz_s):
    """Spin-down timescale nu / (2 |nudot|) in days. Refuses nudot >= 0."""
    nu, nudot = np.asarray(nu_hz, float), np.asarray(nudot_hz_s, float)
    if np.any(nudot[np.isfinite(nudot)] >= 0):
        raise ValueError("nudot must be negative for every pulsar.")
    return nu / (2.0 * np.abs(nudot)) / 86400.0


def clean(raw, cfg, audit, gaps, admissible):
    """D13/D9 cleaning, in order; returns (glitches, log).

    1. remove the non-pulsar objects (D13);
    2. mark usable glitches, dff > 0 and df1 > 0 (unusable rows are kept: they still end intervals);
    2b. fill missing dff_err/df1_err of usable glitches (D4d; flagged);
    3. mark admissible glitches: usable and measured by an admissible method (admissible: D9 classes,
       plus short-window if admitted under D4f; classes from the source audit and D4f);
    4. mark included pulsars: at least D13.min_admissible_completed_intervals_per_pulsar admissible
       glitches that start a completed interval.
    The log has one row per step with rows, pulsars and the counts that step is about.
    """
    key = cfg["D13.pulsar_key"]
    log = [{"step": "0 raw catalogue", "rows": len(raw), "pulsars": raw[key].nunique(),
            "df1_missing": int(raw["df1"].isna().sum()), "df1_le_0": int((raw["df1"] <= 0).sum()),
            "dff_missing": int(raw["dff"].isna().sum()), "dff_le_0": int((raw["dff"] <= 0).sum())}]

    excluded = cfg["D13.excluded_objects"]
    removed = raw[key].isin(excluded)
    g = raw[~removed].copy().rename(columns={key: "pulsar"})
    log.append({"step": "1 remove non-pulsar objects (D13)", "rows": len(g), "pulsars": g["pulsar"].nunique(),
                "rows_removed": int(removed.sum()), "objects_removed": int(raw.loc[removed, key].nunique())})

    g["usable"] = (g["dff"] > cfg["D13.usable_min_dff_exclusive"]) & (g["df1"] > cfg["D13.usable_min_df1_exclusive"])
    log.append({"step": "2 mark usable: dff > 0 and df1 > 0 (D13)", "rows": len(g), "pulsars": g["pulsar"].nunique(),
                "usable": int(g["usable"].sum()), "df1_missing": int(g["df1"].isna().sum()),
                "df1_le_0": int((g["df1"] <= 0).sum()), "dff_missing": int(g["dff"].isna().sum()),
                "dff_le_0": int((g["dff"] <= 0).sum())})

    min_values = cfg["D4d.answers.missing_error_min_values_per_pulsar"]
    g = fill_relative_errors(g, [("dff", "dff_err"), ("df1", "df1_err")], g["usable"], min_values)
    log.append({"step": "2b fill missing errors of usable glitches (D4d)", "rows": len(g),
                "pulsars": g["pulsar"].nunique(), "dff_err_filled": int(g["dff_err_imputed"].sum()),
                "df1_err_filled": int(g["df1_err_imputed"].sum())})

    classes = audit_classes(audit, assigned_classes(cfg))
    g["source_key"] = g["ref"].map(source_key)
    g["source_in_audit"] = g["source_key"].isin(classes)
    g["method_class"] = g["source_key"].map(classes).fillna("unknown")
    g["admissible"] = g["usable"] & g["method_class"].isin(admissible)
    by_class = g[g["usable"]]["method_class"].value_counts()
    all_classes = list(cfg["D8.method_classes_separately"]) + [cfg["D4f.short_window.class"]]
    log.append({"step": "3 mark admissible: " + "+".join(admissible) + " (D9, D4f)",
                "rows": len(g), "pulsars": g["pulsar"].nunique(), "admissible": int(g["admissible"].sum()),
                "usable_source_not_in_audit": int((g["usable"] & ~g["source_in_audit"]).sum()),
                **{f"usable_{c}": int(by_class.get(c, 0)) for c in all_classes}})

    g["starts_completed"] = completed_flags(g, gaps)
    counts = g[g["admissible"] & g["starts_completed"]].groupby("pulsar").size()
    minimum = cfg["D13.min_admissible_completed_intervals_per_pulsar"]
    included = set(counts[counts >= minimum].index)
    g["included_pulsar"] = g["pulsar"].isin(included)
    log.append({"step": f"4 include pulsars with >= {minimum} admissible completed intervals (D13)",
                "rows": int(g["included_pulsar"].sum()), "pulsars": len(included),
                "admissible_completed_in_included": int(counts[counts >= minimum].sum())})
    return g.sort_values(["pulsar", "mjd"], kind="mergesort").reset_index(drop=True), pd.DataFrame(log)
