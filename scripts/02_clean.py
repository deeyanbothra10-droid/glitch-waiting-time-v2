"""Clean the catalogue -> data/processed/glitches_clean.csv, counts after each step -> logs/exclusions.csv.

Steps, in order (glitchlib.cleaning.clean): remove the D13 non-pulsar objects; mark usable glitches;
mark admissible glitches (D9 and D4f: method classes from the source audit and the D4f assignments;
short-window admissible only if results/short_window_admission.json, written by
scripts/emulation_fixture.py, says the D4f admission rule passed); mark included pulsars (D13).
Before any of that, the known counts of the catalogue are checked and the script stops if one differs.
No time differences between glitches are computed.
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.cleaning import admissible_methods, clean  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import gap_starts, read_meta, read_raw_glitches  # noqa: E402

OUT_PATH = ROOT / "data" / "processed" / "glitches_clean.csv"
LOG_PATH = ROOT / "logs" / "exclusions.csv"
ADMISSION_PATH = ROOT / "results" / "short_window_admission.json"
# Known counts of the catalogue as supplied (Stage 4 brief), checked before any rule is applied.
KNOWN = {"rows": 736, "pulsars": 226, "df1_missing": 118, "df1_le_0": 106,
         "rows_of_excluded_objects": 21, "excluded_objects_present": 8}


def main():
    cfg = load_config()
    raw = read_raw_glitches(cfg)
    key = cfg["D13.pulsar_key"]
    excluded = raw[key].isin(cfg["D13.excluded_objects"])
    found = {"rows": len(raw), "pulsars": raw[key].nunique(), "df1_missing": int(raw["df1"].isna().sum()),
             "df1_le_0": int((raw["df1"] <= 0).sum()), "rows_of_excluded_objects": int(excluded.sum()),
             "excluded_objects_present": int(raw.loc[excluded, key].nunique())}
    differ = {k: (KNOWN[k], found[k]) for k in KNOWN if KNOWN[k] != found[k]}
    if differ:
        sys.exit(f"STOP: known counts differ (expected, found): {differ}")

    if not ADMISSION_PATH.exists():
        sys.exit(f"STOP: {ADMISSION_PATH.relative_to(ROOT)} is missing; run scripts/emulation_fixture.py first "
                 "(D4f admission rule for short-window).")
    admissible = admissible_methods(cfg, json.loads(ADMISSION_PATH.read_text(encoding="utf-8")))
    gaps, undated = gap_starts(read_meta("monitoring_gaps"))
    glitches, log = clean(raw, cfg, read_meta("source_audit"), gaps, admissible)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    glitches.to_csv(OUT_PATH, index=False)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    log.to_csv(LOG_PATH, index=False)

    print("Known counts reproduced:", ", ".join(f"{k} {v}" for k, v in found.items()))
    for _, row in log.iterrows():
        counts = ", ".join(f"{k} {int(v)}" for k, v in row.drop("step").items() if pd.notna(v))
        print(f"  {row['step']}: {counts}")
    if len(undated):
        print(f"Monitoring gaps not applied (no start date yet): {', '.join(undated['name'])}")
    unmatched = glitches.loc[glitches["usable"] & ~glitches["source_in_audit"], "ref"].value_counts()
    if len(unmatched):
        print(f"{int(unmatched.sum())} usable glitches have a primary reference not in the source audit "
              f"({len(unmatched)} references).")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)} and {LOG_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
