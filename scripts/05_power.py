"""Power of the primary tests on synthetic data with the real pulsar structure -> results/power.json and
results/power_summary.txt. Synthetic data only: no real waiting times are read.

At each sigma in D14.sigma_grid, with templates from the D13-included pulsars (glitchlib.templates):
  S1 (statistics layer, model true) and S1m (S1 through the measurement emulation, D4c): rejection rates
     of the D4l directional T1 (threshold from the S2 null run) and the D4n T2 (d-hat below the smaller
     2.5th percentile of d-hat under the S3 and S3m null runs), all at the same sigma (power), and the
     D4j outcome rates with the D4l bounds (references: S1, S1m and S4-mix at the same sigma);
  S2 (R irrelevant): rejection rates of T1 and T2 (false positives);
  S4-mix (bias only, real method mix): the D4d reference distribution.
Replicates: D4b (S1, S2: 1,000; S4-mix: 500 per sigma) and D4h (S1m: 1,000 per sigma). The ν̇ noise of the emulated
scenarios is calibrated per pulsar under S1m at each sigma (D4e). Seeds: D4b.simulation.seed.
D4n item 4: every power and outcome rate is reported against the generating sigma and the mean fitted
sigma-hat of the primary model (M2).
"""
import json
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import decision, validation  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import git_commit  # noqa: E402

OUT_PATH = ROOT / "results" / "power.json"
SUMMARY_PATH = ROOT / "results" / "power_summary.txt"


def jobs(cfg):
    """06's jobs S1, S2, S1m and S4-mix (same job ids, so the same seeds), without the secondary extras."""
    return [replace(job, secondary=False) for job in validation.validation_jobs(cfg)
            if job.name in ("S1", "S2", "S1m", "S4-mix")]


def outcome_rates(records, refs, cfg, t2_power):
    return decision.rates([r["outcome"] for r in validation.replicate_outcomes(validation.ok(records), refs, t2_power, cfg)])


def main():
    cfg = load_config()
    workers = os.cpu_count()
    templates, targets = validation.real_templates(cfg)
    results = {"per_sigma": {}}
    lines = [f"{'sigma':>5} {'sigma-hat':>9}  {'scenario':8} {'T1 reject':>9} {'T2 reject':>9} {'A1':>6} {'A2':>6} "
             f"{'replicates':>10}"]
    for s_idx, sigma in enumerate(cfg["D14.sigma_grid"]):
        calibration = validation.calibrate(templates, targets, cfg, sigma, s_idx, workers)
        noise = {c["pulsar"]: c["noise_sd"] for c in calibration}
        records = {job.name: validation.run_job(templates, job, cfg, sigma, s_idx, noise, workers) for job in jobs(cfg)}
        nulls = {job.name: validation.run_job(templates, job, cfg, sigma, s_idx, noise, workers)
                 for job in validation.null_jobs(cfg, cfg["D4b.simulation.S3.slope"])}
        records, thr, refs, power = validation.classified(records, nulls, cfg)
        powers = (power["T1"], power["T2"])
        entry = {"calibration": calibration, "directional_thresholds": thr,
                 "powers_from_S1m": {"T1": powers[0], "T2": powers[1], "generating_sigma": sigma,
                                     "S1m_mean_fitted_sigma_hat": validation.sigma_hat_mean(records["S1m"])},
                 "references": decision.references_for_json(refs), "scenarios": {}}
        for name, rec in records.items():
            good = validation.ok(rec)
            summary = {"replicates": int(len(rec)), "analysed": int(len(good)),
                       "discarded": int(rec["discarded"].sum()), "failed_fits": int(rec["failed"].sum()),
                       "pulsars_dropped_total": int(rec["pulsars_dropped"].sum()),
                       "T1_rejection_rate": validation.rate(rec, "T1_rejects"),
                       "T2_rejection_rate": validation.rate(rec, "T2_rejects"),
                       "T2_rejection_rate_signed_root_D4m_report_only": validation.rate(rec, "T2_rejects_D4m"),
                       "sigma": {"generating": sigma, "mean_fitted_sigma_hat": validation.sigma_hat_mean(rec)},
                       "T1_rejection_rate_chi2_old": validation.rate(rec, "T1_rejects_chi2"),
                       "T2_rejection_rate_chi2_old": validation.rate(rec, "T2_rejects_chi2"),
                       "outcomes": outcome_rates(rec, refs, cfg, powers[1]),
                       "median_completed_intervals": float(good["n_completed"].median()) if len(good) else None,
                       "median_pulsars": float(good["n_pulsars"].median()) if len(good) else None}
            entry["scenarios"][name] = summary
            lines.append(f"{sigma:5.1f} {summary['sigma']['mean_fitted_sigma_hat']:9.3f}  {name:8} {summary['T1_rejection_rate']:9.3f} {summary['T2_rejection_rate']:9.3f} "
                         f"{summary['outcomes']['A1']:6.3f} {summary['outcomes']['A2']:6.3f} {summary['replicates']:10d}")
        results["per_sigma"][str(sigma)] = entry
        print(lines[-4:], flush=True)
    record = {"description": "Power (S1, S1m) and false positives (S2) of T1/T2 with the real pulsar structure; synthetic",
              "per_test_alpha": cfg["D4.per_test_alpha"], "seed": cfg["D4b.simulation.seed"],
              "pulsars": [tp.name for tp in templates], "config_sha256": cfg["meta.decisions_md_sha256"],
              "git_commit": git_commit(), "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              **results}
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    SUMMARY_PATH.write_text("Rejection rates at the D4 per-test alpha; synthetic data with the real pulsar structure\n"
                            + "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"Wrote {OUT_PATH.relative_to(ROOT)} and {SUMMARY_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
