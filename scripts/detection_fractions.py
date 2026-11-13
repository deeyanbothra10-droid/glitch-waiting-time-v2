"""D4o item 5 (report only): the fraction of exponential recoveries the emulation detects, per D4o group,
beside the published detection rates -> results/detection_fractions.json. Synthetic data only.

For S1m and S4-mix at every sigma of the validation grid, the emulation of every replicate is rerun from
06's seeds (job ids and sigma indices) with 06's saved D4e calibration, and the glitches are counted per
D4o group (J0537-6910, crab, others at or above / below the Δν/ν split): all glitches, those with an
exponential, those detected (AIC, D4d), and detections without an exponential. Nothing is fitted.
Beside each D4o count, the same seeds under the superseded D14 placeholder ranges (with the placeholder
record's saved calibration), counted by the same D4o groups. Published: D4o.report_only.
"""
import argparse
import json
import os
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import simulate, validation  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import git_commit, repo_relative  # noqa: E402

OUT_PATH = ROOT / "results" / "detection_fractions.json"
SCENARIOS = ("S1m", "S4-mix")


def fractions(counts):
    return {g: {**c, "detected_of_with_exponential": c["detected"] / c["with_exponential"] if c["with_exponential"] else None,
                "detected_of_glitches": (c["detected"] + c["detected_without_exponential"]) / c["glitches"] if c["glitches"] else None}
            for g, c in sorted(counts.items())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--store", type=Path, default=ROOT / "results" / "validation_partial_d4o",
                        help="06's partial store with the D4o calibrations")
    parser.add_argument("--placeholder-store", type=Path, default=ROOT / "results" / "validation_partial",
                        help="the placeholder record's partial store with its calibrations")
    args = parser.parse_args()
    cfg = load_config()
    workers = os.cpu_count()
    templates, _ = validation.real_templates(cfg)
    d4o = simulate.emulation_from_config(cfg, cfg["D4c.mixtures.unknown_emulated_as"])
    placeholder = replace(d4o, exponentials=simulate.placeholder_exponential_spec(cfg))
    jobs = {job.name: job for job in validation.validation_jobs(cfg) if job.name in SCENARIOS}
    out = {}
    for s_idx, sigma in enumerate(validation.validation_sigma_grid(cfg)):
        out[str(sigma)] = {}
        for label, settings, store in (("D4o", d4o, args.store), ("placeholder_D14", placeholder, args.placeholder_store)):
            cal = validation.load_saved_calibration([store], sigma)
            if cal is None:
                sys.exit(f"STOP: no saved calibration for sigma={sigma} in {store}.")
            noise = {c["pulsar"]: c["noise_sd"] for c in cal}
            for name, job in jobs.items():
                counts = validation.emulated_detection_counts(templates, job, cfg, sigma, s_idx, noise, settings, workers,
                                                              label_spec=d4o.exponentials)
                out[str(sigma)].setdefault(name, {})[label] = fractions(counts)
                print(f"sigma={sigma} {name} {label}: " + "; ".join(
                    f"{g} {v['detected']}/{v['with_exponential']} detected" for g, v in fractions(counts).items()), flush=True)
    record = {"description": "D4o item 5, report only: emulated exponential-recovery detection per D4o group; synthetic",
              "published": cfg["D4o.report_only.published_detection_rates"], "per_sigma": out,
              "stores": {"D4o": repo_relative(args.store), "placeholder_D14": repo_relative(args.placeholder_store)},
              "seed": cfg["D4b.simulation.seed"], "config_sha256": cfg["meta.decisions_md_sha256"], "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
