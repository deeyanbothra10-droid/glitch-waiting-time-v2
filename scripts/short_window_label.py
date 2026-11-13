"""Report only (request of 1 Nov): the D4n item 3 label in 06's simulations -> results/short_window_label.json.
Synthetic data only.

At each sigma of the validation grid, classifies the saved 06 replicates exactly as 06 does (D4l T1, D4n
T2, D4j outcomes with the D4l bounds), then regenerates from their own seeds only the replicates
classified A2, for S1m, S3m, S4-mix and both S4sw variants, with the saved D4e calibration. Each
regenerated replicate must reproduce its saved b-hat (else the run stops). Reports, per scenario and
sigma, the fraction labelled "short-window-sensitive" (d-hat refitted without short-window glitches
>= 0) among the A2 replicates that have a refit, and how many refits had no data left or failed (where no
A2 replicate has a refit, as in S4sw where every glitch is short-window, the label is not computable).

Saved stores are read read-only, searched in order (--reuse, repeatable; default: the three stores 06
was last computed from). Nothing else is simulated: a missing job or calibration stops the run.
"""
import argparse
import json
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import decision, validation  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import git_commit, repo_relative  # noqa: E402

OUT_PATH = ROOT / "results" / "short_window_label.json"
DEFAULT_STORES = [ROOT / "results" / d for d in ("validation_partial", "validation_partial_d4l", "validation_partial_d4m")]
SCENARIOS = ("S1m", "S3m", "S4-mix", "S4sw-espinoza", "S4sw-wang")     # as requested; report only
NEEDED = ("S1", "S1m", "S4-mix")                                       # the D4j references


def load(stores, sigma, name, cfg):
    rec, _ = validation.load_saved(stores, sigma, name)
    if rec is None:
        sys.exit(f"STOP: sigma={sigma} {name} is not in the saved stores; it would need new simulation.")
    return validation.from_saved_records(rec, cfg)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse", type=Path, action="append", default=None,
                        help="a saved partial-store directory (read-only); repeatable, searched in order")
    args = parser.parse_args()
    stores = args.reuse or DEFAULT_STORES
    cfg = load_config()
    workers = os.cpu_count()
    templates, _ = validation.real_templates(cfg)
    jobs = {job.name: job for job in validation.validation_jobs(cfg)}
    out = {}
    for s_idx, sigma in enumerate(validation.validation_sigma_grid(cfg)):
        records = {name: load(stores, sigma, name, cfg) for name in set(NEEDED) | set(SCENARIOS)}
        nulls = {job.name: load(stores, sigma, job.name, cfg)
                 for job in validation.null_jobs(cfg, cfg["D4b.simulation.S3.slope"])}
        records, _, refs, powers = validation.classified(records, nulls, cfg)
        calibration = validation.load_saved_calibration(stores, sigma)
        if calibration is None:
            sys.exit(f"STOP: no saved calibration for sigma={sigma}; it would need new simulation.")
        noise = {c["pulsar"]: c["noise_sd"] for c in calibration}
        out[str(sigma)] = {}
        for name in SCENARIOS:
            good = validation.ok(records[name])
            outcomes = validation.replicate_outcomes(good, refs, powers["T2"], cfg)
            a2 = good["replicate"][[o["outcome"] == "A2" for o in outcomes]].to_numpy()
            job = replace(jobs[name], secondary=False)        # same seeds and data; no permutations or LOO
            regen = validation.run_job(templates, job, cfg, sigma, s_idx, noise, workers, replicates=a2)
            if len(a2):
                saved_b = good.set_index("replicate").loc[a2, "b"].to_numpy(float)
                new = validation.ok(regen).set_index("replicate")
                if len(new) != len(a2) or not np.allclose(new.loc[a2, "b"].to_numpy(float), saved_b, rtol=0, atol=1e-9):
                    sys.exit(f"STOP: sigma={sigma} {name}: the regenerated A2 replicates do not reproduce the saved b-hat.")
            labels = [decision.short_window_label("A2", None if d is None or np.isnan(d) else float(d), cfg)
                      for d in (regen["d_without_short_window"].astype(float) if len(regen) else [])]
            n_label = sum(lab is not None for lab in labels)
            no_data = int(regen["d_without_short_window"].isna().sum()) if len(regen) else 0
            fit_failed = int(regen["d_without_short_window_failure"].notna().sum()) if "d_without_short_window_failure" in regen else 0
            computable = len(a2) - no_data          # the label is undefined where no refit exists
            out[str(sigma)][name] = {"analysed": int(len(good)), "A2_replicates": int(len(a2)),
                                     "A2_with_refit": int(computable),
                                     "labelled_short_window_sensitive": int(n_label),
                                     "fraction_labelled": float(n_label / computable) if computable else None,
                                     "refit_no_data_or_failed": no_data, "refit_failed": fit_failed,
                                     "d_without_short_window_mean": float(regen["d_without_short_window"].astype(float).mean())
                                     if len(a2) else None}
            v = out[str(sigma)][name]
            print(f"sigma={sigma} {name}: A2 {v['A2_replicates']}, with refit {v['A2_with_refit']}, labelled "
                  f"{v['labelled_short_window_sensitive']}"
                  + (f" ({v['fraction_labelled']:.3f})" if computable else " (not computable)"), flush=True)
    record = {"description": "Report only: D4n short-window-sensitive label among 06's A2 replicates (regenerated "
                             "from their seeds); synthetic data",
              "stores": [repo_relative(d) for d in stores],
              "store_meta": {repo_relative(d): json.loads((Path(d) / "meta.json").read_text(encoding="utf-8")) for d in stores},
              "seed": cfg["D4b.simulation.seed"], "per_sigma": out, "config_sha256": cfg["meta.decisions_md_sha256"],
              "git_commit": git_commit(), "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
