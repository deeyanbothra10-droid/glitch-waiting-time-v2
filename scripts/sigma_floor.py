"""Report only (request of 1 Nov): the floor of the fitted sigma-hat of the primary model (D3, M2) for S1m,
S3, S3m and S4-mix at the lower end of the matched-sigma search -> results/sigma_floor.json. Synthetic data
only.

The generating sigma is the lower end of D4p's search range (0.01; it was 0 under D4n) and the replicate
count is D4n's replicates per step, so this is the first step of the bisection for these scenarios. Seeds: 06's job ids (so 06's seeding)
at their own sigma index (validation.floor_sigma_index), disjoint from every other run. The emulated
scenarios use the D4e noise calibrated at that sigma; S3's size coefficient is D4b's, as in 06.
Reports the mean and the 2.5th and 97.5th percentiles of sigma-hat over the analysed replicates, and
the failed replicates with their reasons (a scenario whose fits all fail is reported as such). Decides
nothing.
"""
import json
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import validation  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_PATH = ROOT / "results" / "sigma_floor.json"
SCENARIOS = ("S1m", "S3", "S3m", "S4-mix")          # as requested; report only


def main():
    cfg = load_config()
    workers = os.cpu_count()
    sigma = float(cfg["D4p.matched_sigma_search.search_range"][0])
    n = cfg["D4n.matched_sigma.replicates_per_step"]
    index = validation.floor_sigma_index(cfg)
    templates, targets = validation.real_templates(cfg)
    jobs = [replace(job, n_replicates=n, secondary=False) for job in validation.validation_jobs(cfg)
            if job.name in SCENARIOS]
    calibration = validation.calibrate(templates, targets, cfg, sigma, index, workers)
    noise = {c["pulsar"]: c["noise_sd"] for c in calibration}
    out = {}
    for job in jobs:
        rec = validation.run_job(templates, job, cfg, sigma, index, noise if job.scenario.emulated else None, workers)
        good = validation.ok(rec)
        failures = validation.failure_check(rec, cfg)
        out[job.name] = {"replicates": int(len(rec)), "analysed": int(len(good)), "discarded": int(rec["discarded"].sum()),
                         "failures": {**failures, "reasons": sorted({r["failure"] for r in failures["reasons"]})}}
        if len(good):                    # none analysed: every fit failed (reported, with the reasons)
            values = good["sigma_hat"].to_numpy(float)
            out[job.name].update({"sigma_hat_mean": float(values.mean()), "sigma_hat_median": float(np.median(values)),
                                  "sigma_hat_p2.5": float(np.percentile(values, 2.5)),
                                  "sigma_hat_p97.5": float(np.percentile(values, 97.5))})
    record = {"description": "Report only: fitted sigma-hat (M2) at the lower end of the D4p search range; synthetic data",
              "generating_sigma": sigma, "replicates_per_scenario": n, "sigma_index": index,
              "seed": cfg["D4b.simulation.seed"], "calibration": calibration,
              "calibration_outside_tolerance": [c["pulsar"] for c in calibration if not c["within_tolerance"]],
              "scenarios": out, "pulsars": [tp.name for tp in templates], "config_sha256": cfg["meta.decisions_md_sha256"],
              "git_commit": git_commit(), "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    for name, v in out.items():
        if not v["analysed"]:
            print(f"{name}: no analysed replicate ({v['failures']['failed']} failed: {'; '.join(v['failures']['reasons'])})")
            continue
        print(f"{name}: sigma-hat mean {v['sigma_hat_mean']:.3f} [{v['sigma_hat_p2.5']:.3f}, {v['sigma_hat_p97.5']:.3f}] "
              f"({v['analysed']} analysed, {v['failures']['failed']} failed)")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
