"""EXPLORATORY, POST HOC, NOT PRE-REGISTERED (requested 4 Nov, after the frozen run 07-09)
-> results/exploratory/single_pulsar.json. Uses real waiting times for the real fits, so the guard must pass.

For each of B0531+21, B1737-30 and J0537-6910 separately, using only that pulsar's template (its real
glitch properties, source cadences, observation window, gaps and method classes; never its waiting times):
1. M4 (and M3, for the size coefficient) fitted to that pulsar's primary intervals alone (D13 count applied);
   d-hat and M4's sigma-hat recorded.
2. 1,000 one-pulsar replicates each of S3m (size only, measured; slope = that pulsar's own M3 a-hat), the
   pulsar's own-method bias-only variant (S4sw-espinoza for B0531+21 and B1737-30, S4l for J0537-6910) and
   S1m (model true), with 06's scenario and emulation definitions. Each scenario is run at the generating
   sigma at which its median M4 sigma-hat (one-pulsar replicates) matches that pulsar's real M4 sigma-hat:
   D4n's bisection with D4p's range and step limit, D4n's 200 replicates per step and 0.01 tolerance. The
   ν̇ noise is calibrated for that pulsar under S1m at each candidate and at the matched sigma (D4e).
   Seeds: SeedSequence((D4b seed, 900, pulsar index, job id, phase, replicate)), phase 0 for search steps
   (same seeds at every step) and 1 for the full run; calibration seed index 2000 + 10 * pulsar index + phase.
3. The real d-hat's percentile within each scenario's d-hat distribution, and that distribution's 2.5th, 50th
   and 97.5th percentiles.
Numbers only; decides nothing; changes no committed result.
"""
import json
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import decision, exploratory, models, simulate, validation  # noqa: E402
from glitchlib.analysis import primary_data  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.guards import FINAL_DATASET_REL, check_unblinding_allowed  # noqa: E402
from glitchlib.intervals import apply_inclusion  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_DIR = ROOT / "results" / "exploratory"
PULSARS = {"B0531+21": "S4sw-espinoza", "B1737-30": "S4sw-espinoza", "J0537-6910": "S4l"}   # as requested
N_REPLICATES = 1000                                                                          # as requested


def main():
    cfg = load_config()
    checksum = check_unblinding_allowed()
    workers = os.cpu_count()
    seed = cfg["D4b.simulation.seed"]
    data = primary_data(pd.read_csv(ROOT / FINAL_DATASET_REL))
    templates, targets = validation.real_templates(cfg)
    by_name = {tp.name: tp for tp in templates}
    jobs = {job.name: job for job in validation.validation_jobs(cfg)}
    me = cfg["D6.primary_errors_in_lnR"]
    out = {}
    for k, (pulsar, bias_variant) in enumerate(PULSARS.items()):
        own, _ = apply_inclusion(data[data["pulsar"] == pulsar], cfg["D13.min_admissible_completed_intervals_per_pulsar"])
        m4, m3 = models.fit("M4", own, me), models.fit("M3", own, me)
        real = {"n_intervals": int(len(own)), "n_completed": int(own["event"].astype(bool).sum()),
                "d_hat": float(m4.estimates["d"]), "sigma_hat_M4": float(m4.estimates["sigma"]),
                "a_hat_M4": float(m4.estimates["a"]), "a_hat_M3": float(m3.estimates["a"])}
        tp = by_name[pulsar]
        out[pulsar] = {"real": real, "scenarios": {}}
        print(f"{pulsar}: d-hat {real['d_hat']:.4f}, M4 sigma-hat {real['sigma_hat_M4']:.4f}, M3 a-hat {real['a_hat_M3']:.4f}",
              flush=True)
        for name in ("S3m", bias_variant, "S1m"):
            job = jobs[name]
            scenario = replace(job.scenario, slope=real["a_hat_M3"]) if name == "S3m" else job.scenario
            settings = replace(simulate.emulation_from_config(cfg, cfg["D4c.mixtures.unknown_emulated_as"]),
                               short_window_days=job.short_window_days)

            def noise_at(sigma, phase):
                cal = validation.calibrate([tp], {pulsar: targets[pulsar]}, cfg, sigma, 2000 + 10 * k + phase, workers)
                return cal[0]

            def median_at(sigma):
                rec = exploratory.run_single_pulsar(tp, scenario, job.method_rule, settings, sigma,
                                                    noise_at(sigma, 0)["noise_sd"],
                                                    cfg["D4n.matched_sigma.replicates_per_step"],
                                                    (seed, 900, k, job.job_id, 0), cfg, workers)
                good = rec[rec["status"] == "ok"]
                if not len(good):
                    raise validation.SigmaMatchError(f"{pulsar} {name}: no analysed replicate at sigma {sigma}")
                return float(good["sigma_hat"].median())

            match = validation.match_sigma(median_at, real["sigma_hat_M4"], cfg)
            cal = noise_at(match["sigma"], 1)
            rec = exploratory.run_single_pulsar(tp, scenario, job.method_rule, settings, match["sigma"], cal["noise_sd"],
                                                N_REPLICATES, (seed, 900, k, job.job_id, 1), cfg, workers)
            good = rec[rec["status"] == "ok"]
            d = good["d"].to_numpy(float)
            entry = {"method_rule": job.method_rule, "short_window_days": job.short_window_days,
                     "slope": scenario.slope, "covariate": scenario.covariate,
                     "matched_generating_sigma": match["sigma"], "run_at_lower_end": match["at_lower_end"],
                     "bisection_converged": match["converged"], "bisection_steps": match["steps"],
                     "noise_sd": cal["noise_sd"], "calibration_within_tolerance": cal["within_tolerance"],
                     "replicates": int(len(rec)), "analysed": int(len(good)),
                     "discarded": int((rec["status"] == "discarded").sum()),
                     "failed": int(rec["status"].str.startswith("failed").sum()),
                     "sigma_hat_M4_median": float(good["sigma_hat"].median()),
                     "d_hat_p2.5": float(np.percentile(d, 2.5)), "d_hat_p50": float(np.percentile(d, 50)),
                     "d_hat_p97.5": float(np.percentile(d, 97.5)),
                     "real_d_hat_percentile": decision.percentile_of(real["d_hat"], d)}
            out[pulsar]["scenarios"][name] = entry
            print(f"  {name}: matched sigma {entry['matched_generating_sigma']:.4f}"
                  f"{' (lower end)' if entry['run_at_lower_end'] else ''}{'' if entry['bisection_converged'] else ' (not converged)'}"
                  f", median sigma-hat {entry['sigma_hat_M4_median']:.4f}; analysed {entry['analysed']} "
                  f"(discarded {entry['discarded']}, failed {entry['failed']}); d-hat 2.5/50/97.5 = "
                  f"{entry['d_hat_p2.5']:.4f} / {entry['d_hat_p50']:.4f} / {entry['d_hat_p97.5']:.4f}; "
                  f"real d-hat percentile {entry['real_d_hat_percentile']:.1f}", flush=True)
    record = {"label": exploratory.LABEL, "pulsars": out, "replicates": N_REPLICATES, "seed": seed,
              "final_dataset_sha256": checksum, "config_sha256": cfg["meta.decisions_md_sha256"], "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "single_pulsar.json").write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    print(exploratory.LABEL)
    print("Wrote results/exploratory/single_pulsar.json")


if __name__ == "__main__":
    main()
