"""Robustness analyses on the real data -> results/robustness.json. Same preflight as 07 (with
results/robustness.json absent).

Runs the primary analysis (glitchlib.analysis.primary) for every D8 item, D4f's runs (without
J0537-6910, always reported, exploratory below 40 intervals; without "Our measurement"), D4g item 10
(Jodrell Bank pulsars with obs_end = max(last glitch, 60891)), and D7's exploratory Vela-like / Crab-like
split, and D4o's case without the Crab. Every subset gets the D13 inclusion count of the data being
analysed (confirmed reading); a subset with no pulsar left is reported as "no data". D4p: each case reports
b-hat, d-hat, their calibrated intervals (D4m factors of the primary's null references, read from
results/secondary.json; labelled "approximate (primary factors)", D4p clarification), the uncalibrated Wald intervals, the signed roots, and whether each estimate has
the primary's sign (results/primary.json); no case gets an outcome classification. A fit failure on real
data stops the run (D4k). Prints numbers only.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import models, preflight  # noqa: E402
from glitchlib.analysis import primary_data, robustness_case  # noqa: E402
from glitchlib.cleaning import OUR_MEASUREMENT, jodrell_bank_pulsars  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.intervals import add_derived, build_intervals, epoch_uncertainty_exclusion  # noqa: E402
from glitchlib.io import gap_starts, git_commit, read_meta  # noqa: E402

OUT_REL = "results/robustness.json"


def run(data, cfg, label, exploratory_below=None):
    try:
        return robustness_case(data, cfg, label, exploratory_below, FACTORS, PRIMARY)
    except models.FitError as err:
        sys.exit(f"STOP: a fit failed on the real data in '{label}', so the run stops (D4k): {err}")


FACTORS, PRIMARY = {}, {}


def main():
    preflight.require_all(OUT_REL)
    cfg = load_config()
    if (cfg["D4p.robustness_09.outcome_classification"] != "none"
            or cfg["D4p.robustness_09.interval_factors"] != "primary_factors_from_08"):
        sys.exit("STOP: only D4p's unclassified robustness report is implemented.")
    for name in ("primary", "secondary"):
        if not (ROOT / "results" / f"{name}.json").exists():
            sys.exit(f"STOP: results/{name}.json does not exist; run 07 and 08 first (D4p needs their estimates and factors).")
    prim = json.loads((ROOT / "results" / "primary.json").read_text(encoding="utf-8"))
    sec = json.loads((ROOT / "results" / "secondary.json").read_text(encoding="utf-8"))
    PRIMARY.update({prim[t]["parameter"]: prim[t]["estimate"] for t in ("T1", "T2")})
    FACTORS.update({prim[t]["parameter"]: sec["intervals_b_d"][t]["factor"] for t in ("T1", "T2")})
    checksum = check_unblinding_allowed()
    final = pd.read_csv(ROOT / FINAL_DATASET_REL)
    data = primary_data(final)
    d8 = "D8"
    cases = []
    cases.append(run(data[data["pulsar"] != cfg[f"{d8}.exclude_pulsar"]], cfg, "D8 item 1: without J0537-6910",
                     cfg["D4f.robustness_without_J0537.exploratory_if_fewer_intervals_than"]))
    classes = list(cfg[f"{d8}.method_classes_separately"]) + [cfg["D4f.short_window.class"]]
    for cls in classes:
        cases.append(run(final[final["method_class"] == cls], cfg, f"D8 item 2: method class {cls}"))
    counts = data[data["event"].astype(bool)].groupby("pulsar").size()
    for pulsar in counts[counts >= cfg[f"{d8}.per_pulsar_min_admissible_completed"]].index:
        cases.append(run(data[data["pulsar"] == pulsar], cfg, f"D8 item 3: {pulsar} alone"))
    cases.append(run(data[data["event"].astype(bool)], cfg, "D8 item 4: censored intervals dropped"))
    cases.append(run(data[data["sigma_lnR"] <= cfg[f"{d8}.max_sigma_lnR"]], cfg, "D8 item 5: sigma_lnR <= 0.5"))
    cases.append(run(data[~data["gap_flag"].astype(bool)], cfg, "D8 item 6: gap-flagged intervals excluded"))
    cases.append(run(data[np.exp(data["ln_dnu_frac"]) >= cfg[f"{d8}.min_dnu_frac"]], cfg, "D8 item 7: dnu/nu >= 1e-7"))
    cases.append({"case": "D8 item 8: step dnudot", "status": f"{cfg['D4h.robustness_item_8.status'].replace('_', ' ')}: "
                                                              f"{cfg['D4h.robustness_item_8.reason']} (D4h)"})
    item9 = "D4h.robustness_item_9"
    if cfg[f"{item9}.combination"] != "quadrature" or cfg[f"{item9}.missing_error"] != "exclude":
        sys.exit("STOP: only D4h item 9 with quadrature and missing errors excluded is implemented.")
    excluded = epoch_uncertainty_exclusion(data, cfg[f"{item9}.max_fraction_of_T"])
    cases.append(run(data[~excluded], cfg, "D8 item 9: combined epoch uncertainty <= 0.1 T (D4h)"))

    # D4g item 10: Jodrell Bank pulsars with obs_end = max(last glitch, 60891).
    glitches = pd.read_csv(ROOT / "data" / "processed" / "glitches_clean.csv")
    obs = read_meta("obs_end")
    obs_end = dict(zip(obs["name"], obs["obs_end_mjd"].astype(float)))
    item10 = "D4g.robustness_item_10"
    last = glitches.groupby("pulsar")["mjd"].max()
    for pulsar in jodrell_bank_pulsars(glitches, cfg):
        if pulsar in obs_end:
            obs_end[pulsar] = max(float(last[pulsar]), cfg[f"{item10}.mjd"])
    gaps, _ = gap_starts(read_meta("monitoring_gaps"))
    rebuilt, _ = build_intervals(glitches[glitches["included_pulsar"]], obs_end, gaps)
    rebuilt = add_derived(rebuilt)
    cases.append(run(rebuilt[rebuilt["admissible"].astype(bool)], cfg, "D4g item 10: Jodrell Bank obs_end"))

    cases.append(run(data[data["source_key"] != OUR_MEASUREMENT], cfg, 'D4f: without "Our measurement"'))
    for item in cfg["D4o.D8_robustness_additions"]:                       # D4o item 7
        pulsar = item.removeprefix("exclude_")
        cases.append(run(data[data["pulsar"] != pulsar], cfg, f"D4o: without {pulsar}"))
    # D7 (exploratory): Vela-like if the pulsar's median usable dnu/nu >= threshold.
    usable = glitches[glitches["usable"]]
    median_dnu = (usable["dff"] * 1e-9).groupby(usable["pulsar"]).median()
    vela_like = set(median_dnu[median_dnu >= cfg["D7.vela_like_min_dnu_frac"]].index)
    cases.append(run(data[data["pulsar"].isin(vela_like)], cfg, "D7 exploratory: Vela-like pulsars"))
    cases.append(run(data[~data["pulsar"].isin(vela_like)], cfg, "D7 exploratory: Crab-like pulsars"))

    record = {"cases": cases, "final_dataset_sha256": checksum, "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (ROOT / OUT_REL).write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    for c in cases:
        if c["status"] == "run":
            e = c["estimates_D4p"]
            print(f"{c['case']}: n={c['n_intervals']}, " + ", ".join(
                f"{k}={e[k]['estimate']:.4f} [{e[k]['calibrated_interval'][0]:.4f}, {e[k]['calibrated_interval'][1]:.4f}] "
                f"signed root {e[k]['signed_root']:.3f}, same sign as primary: {e[k]['same_sign_as_primary']}"
                for k in ("b", "d")))
        else:
            print(f"{c['case']}: {c['status']}")
    print(f"Wrote {OUT_REL}")


if __name__ == "__main__":
    main()
