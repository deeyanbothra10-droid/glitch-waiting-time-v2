"""Glitch counts per Δν̇ method class and per source -> results/method_counts.json. Allowed before
freezing: counts only, no waiting times.

Counts usable glitches by method class (from the source audit via 02_clean.py) for all pulsars, for the
pulsars with at least 5 usable glitches, and for the D13-included pulsars; per pulsar, the admissible
usable glitches and the status of the intervals they start (completed, the D13 count; censored at
obs_end or at a gap; last interval dropped for want of a censoring date), found by comparing epochs
only (glitchlib.intervals.interval_status); and lists the sources of usable glitches whose method is
still unknown, largest first, as the audit to-do list.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.config import load_config  # noqa: E402
from glitchlib.intervals import interval_status  # noqa: E402
from glitchlib.io import gap_starts, git_commit, read_meta  # noqa: E402

IN_PATH = ROOT / "data" / "processed" / "glitches_clean.csv"
OUT_PATH = ROOT / "results" / "method_counts.json"
MIN_USABLE_FOR_META = 5      # the 21 pulsars of data/meta (data/meta/README.md), for reference counts only


def by_class(frame, classes):
    counts = frame["method_class"].value_counts()
    return {c: int(counts.get(c, 0)) for c in classes}


def main():
    cfg = load_config()
    classes = list(cfg["D8.method_classes_separately"]) + [cfg["D4f.short_window.class"]]
    glitches = pd.read_csv(IN_PATH)
    usable = glitches[glitches["usable"]]
    n_usable = usable.groupby("pulsar").size()
    meta_pulsars = n_usable[n_usable >= MIN_USABLE_FOR_META].index
    todo = (usable[usable["method_class"] == "unknown"]
            .groupby(["source_key", "source_in_audit"]).size().sort_values(ascending=False))
    obs = read_meta("obs_end")
    gaps, _ = gap_starts(read_meta("monitoring_gaps"))
    glitches["status"] = interval_status(glitches, dict(zip(obs["name"], obs["obs_end_mjd"].astype(float))), gaps)
    admissible = glitches[glitches["admissible"]]
    status_counts = pd.crosstab(admissible["pulsar"], admissible["status"])
    for col in ["completed", "censored_end", "censored_gap", "dropped_no_end", "dropped_end_at_last_glitch"]:
        if col not in status_counts.columns:
            status_counts[col] = 0
    per_pulsar = pd.DataFrame({
        "admissible_usable": admissible.groupby("pulsar").size(),
        "admissible_completed": status_counts["completed"],
        "censored_at_obs_end": status_counts["censored_end"], "censored_at_gap": status_counts["censored_gap"],
        "last_dropped": status_counts["dropped_no_end"] + status_counts["dropped_end_at_last_glitch"],
        "included": admissible.groupby("pulsar")["included_pulsar"].first()}).sort_values(
        "admissible_completed", ascending=False)
    included = per_pulsar[per_pulsar["included"]]
    totals = {key: int(included[key].sum()) for key in
              ["admissible_completed", "censored_at_obs_end", "censored_at_gap", "last_dropped"]}
    record = {
        "usable_all_pulsars": by_class(usable, classes),
        f"usable_pulsars_with_>=_{MIN_USABLE_FOR_META}_usable": by_class(usable[usable["pulsar"].isin(meta_pulsars)], classes),
        "usable_included_pulsars": by_class(usable[usable["included_pulsar"]], classes),
        "admissible_per_pulsar": per_pulsar.reset_index().to_dict("records"),
        "included_pulsars_totals": dict(totals, pulsars=int(len(included))),
        "unknown_method_sources": [{"source_key": k, "in_audit": bool(a), "usable_glitches": int(n)}
                                   for (k, a), n in todo.items()],
        "git_commit": git_commit(), "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(record, indent=2, default=lambda x: x.item()) + "\n", encoding="utf-8")
    table = pd.DataFrame({k: record[k] for k in list(record)[:3]})
    print(table.to_string())
    print("Admissible counts per pulsar (D13 needs >= "
          f"{cfg['D13.min_admissible_completed_intervals_per_pulsar']} admissible completed intervals):")
    print(included.to_string() if len(included) else "  none")
    print(f"Included pulsars: {len(included)}; totals: " + ", ".join(f"{k} {v}" for k, v in totals.items()))
    print(f"Not included: {int((~per_pulsar['included']).sum())} pulsars with admissible glitches "
          "(full table in the JSON).")
    print(f"{len(todo)} sources with unknown method; largest: "
          + ", ".join(f"{k} ({n})" for (k, _), n in todo.head(8).items()))
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
