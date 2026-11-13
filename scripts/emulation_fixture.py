"""Measure S4 bias magnitudes (D4d: a measurement, not a test) on a synthetic fixture.

The fixture pulsars are synthetic stand-ins, not the real pulsars: the real-pulsar S4 runs come later
(scripts/05_power.py, scripts/06_validate.py) once the source audit and the real templates exist.
Every emulation setting comes from config (D4b, D4c, D4d, D4e, D4f, D14); only the synthetic pulsars
below are defined here. The only real input is the source audit's cadences (allowed: real glitch structure, never real waiting times);
no glitch data are read.

For each sigma in D14.sigma_grid and each S4 variant (one method for every glitch, or the method mix
taken from the resampled glitch with "unknown" as window and, for sensitivity, as global), runs
D4b.simulation.replicates.S4_per_variant replicates and records the forward slope (M2, ln T on ln R)
and the backward slope (placebo, previous wait on ln R). D4e: the ν̇ noise is calibrated once per
pulsar per sigma under S1m on the pulsar's method mix, and that noise is used in every variant; each
replicate goes through the D13 inclusion rule (pulsars with too few completed intervals in the
analysed data are dropped; the replicate is discarded if fewer than 2 pulsars remain).
D4f: the short-window variants S4sw run at each source type's main and sensitivity span w, sampled at
that source's audited cadence (D4c: cadence per source from the audit, data/meta/source_audit.csv
cadence_num_days; all other variants keep the fixture's fallback cadence), and the admission rule
(forward-slope 95% upper bound = mean + 1.96 SE <= +0.02 at every sigma and every w) is applied.
Outputs (since D4o, post-freeze): results/emulation_s4_fixture_d4o.json and results/short_window_admission_d4o.json.
The pre-D4o files results/emulation_s4_fixture.json and results/short_window_admission.json (D14 placeholder Q and
tau_d) are the record of the D4f admission that scripts/02_clean.py reads; this script never rewrites them.
"""
import dataclasses
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib import diagnostics, models, simulate  # noqa: E402
from glitchlib.cleaning import source_cadences  # noqa: E402
from glitchlib.config import load_config  # noqa: E402
from glitchlib.intervals import apply_replicate_inclusion  # noqa: E402
from glitchlib.io import git_commit, read_meta  # noqa: E402

OUT_PATH = ROOT / "results" / "emulation_s4_fixture_d4o.json"
ADMISSION_PATH = ROOT / "results" / "short_window_admission_d4o.json"

# ---- synthetic fixture (not real pulsars) ------------------------------------------------------
FIXTURE = {
    "n_pulsars": 10, "pool_size": 15, "n_glitches": 25, "first_mjd": 50000.0, "obs_end_mjd": 75000.0,
    "tau_sd_days": 1e6, "ln_dnu_frac": [np.log(1e-6), 0.8], "ln_dnudot_frac": [np.log(5e-3), 0.6],
    "dnu_rel_err_uniform": [0.005, 0.05], "dnudot_rel_err_uniform": [0.02, 0.2],
    "method_classes": ["window", "local-relaxation", "global", "global-noexp", "unknown"], "fixture_seed": 99,
}


def fixture_templates(cadence_days):
    f = FIXTURE
    rng = np.random.default_rng(f["fixture_seed"])
    templates = []
    for p in range(f["n_pulsars"]):
        n = f["pool_size"]
        dnu = np.exp(rng.normal(*f["ln_dnu_frac"], n))
        dnudot = np.exp(rng.normal(*f["ln_dnudot_frac"], n))
        pool = pd.DataFrame({"dnu_frac": dnu, "dnu_frac_err": rng.uniform(*f["dnu_rel_err_uniform"], n) * dnu,
                             "dnudot_frac": dnudot,
                             "dnudot_frac_err": rng.uniform(*f["dnudot_rel_err_uniform"], n) * dnudot,
                             "method_class": rng.choice(f["method_classes"], n), "cadence_days": cadence_days})
        templates.append(simulate.PulsarTemplate(f"F{p}", f["first_mjd"], f["obs_end_mjd"], f["n_glitches"],
                                                 pool, f["tau_sd_days"], (), False, 1.0))
    return templates



def with_cadence(templates, cadence_days):
    """Copies of the fixture templates whose glitches all come from a source with this cadence."""
    return [dataclasses.replace(tp, pool=tp.pool.assign(cadence_days=cadence_days)) for tp in templates]


def summarise(values):
    v = np.asarray(values, float)
    if v.size < 2:          # nothing (or one value) to summarise, e.g. every replicate discarded
        return {"n": int(v.size), "mean": None, "se_mean": None, "sd": None, "p2.5": None, "p97.5": None}
    return {"n": int(v.size), "mean": float(v.mean()), "se_mean": float(v.std(ddof=1) / np.sqrt(v.size)),
            "sd": float(v.std(ddof=1)), "p2.5": float(np.percentile(v, 2.5)),
            "p97.5": float(np.percentile(v, 97.5))}


def slopes(intervals):
    return (models.fit("M2", intervals, "ignore").estimates["b"],
            diagnostics.placebo_fit(intervals, "ignore").estimates["b"])


class Tally:
    """Slopes of one analysed dataset across replicates, with D4e inclusion counts."""

    def __init__(self):
        self.forward, self.backward, self.dropped, self.discarded, self.failed = [], [], [], 0, 0

    def add(self, intervals, min_completed, min_pulsars, pulsars):
        kept, dropped, discarded = apply_replicate_inclusion(intervals, min_completed, min_pulsars, pulsars)
        self.dropped.append(len(dropped))
        if discarded:
            self.discarded += 1
            return
        try:
            f, b = slopes(kept)
        except models.FitError:
            self.failed += 1
            return
        self.forward.append(f)
        self.backward.append(b)

    def summary(self):
        return {"forward": summarise(self.forward), "backward": summarise(self.backward),
                "pulsars_dropped_total": int(np.sum(self.dropped)),
                "pulsars_dropped_mean_per_replicate": float(np.mean(self.dropped)),
                "replicates_discarded": self.discarded, "failed_fits": self.failed}


def main():
    cfg = load_config()
    seed = cfg["D14.seed"]
    n_rep = cfg["D4b.simulation.replicates.S4_per_variant"]
    admissible = cfg["D9.admissible_methods"]
    min_completed = cfg["D13.min_admissible_completed_intervals_per_pulsar"]
    min_pulsars = cfg["D4e.replicate_inclusion.discard_replicate_if_fewer_pulsars_than"]
    tolerance = cfg["D4d.answers.calibration_tolerance_fraction_of_target"]
    templates = fixture_templates(cfg["D4c.emulation.cadence_until_audit_days"])
    targets = simulate.target_relative_errors(templates, cfg["D4c.emulation.noise_min_error_values_per_pulsar"])
    s4 = simulate.scenario_from_config(cfg, "S4")
    calibration_scenario = simulate.scenario_from_config(cfg, cfg["D4e.calibration.scenario"])
    mix_settings = simulate.emulation_from_config(cfg, cfg["D4c.mixtures.unknown_emulated_as"])
    variants = {"S4w": ("window", "window", None, None), "S4l": ("local-relaxation", "window", None, None),
                "S4g": ("global", "window", None, None), "S4g0": ("global-noexp", "window", None, None),
                "S4-mix": ("resampled", cfg["D4c.mixtures.unknown_emulated_as"], None, None),
                "S4-mix-unknown-as-global": ("resampled", cfg["D4c.mixtures.unknown_sensitivity_emulated_as"], None, None)}
    cadences = source_cadences(read_meta("source_audit"))
    fallback = cfg["D4c.emulation.cadence_until_audit_days"]
    for key, source in cfg["D4f.short_window.sources"].items():
        cadence = cadences.get(key, fallback)
        for w in (source["w_days"], source["w_sensitivity_days"]):
            variants[f"S4sw-{source['source_type']}-w{w:g}"] = ("short-window", "window", float(w), cadence)
    admission_entries = []
    seeds = np.random.SeedSequence(seed)
    results = {"calibration": {}, "variants": {}}
    for sigma in cfg["D14.sigma_grid"]:
        sigma_seq = seeds.spawn(1)[0]
        cal_rng = np.random.default_rng(sigma_seq.spawn(1)[0])
        noise, cal_rows = {}, []
        for tp in templates:
            out = simulate.calibrate_noise(tp, calibration_scenario, sigma, "resampled", mix_settings, targets[tp.name],
                                           cfg["D4d.answers.calibration_sequences_per_step"], tolerance, cal_rng)
            noise[tp.name] = out["noise_sd"]
            cal_rows.append(out)
        results["calibration"][str(sigma)] = {
            "scenario": calibration_scenario.name, "per_pulsar": cal_rows,
            "outside_tolerance": [r["pulsar"] for r in cal_rows if not r["within_tolerance"]]}
        for name, (rule, unknown_as, w, cadence) in variants.items():
            settings = dataclasses.replace(simulate.emulation_from_config(cfg, unknown_as), short_window_days=w)
            run_templates = templates if cadence is None else with_cadence(templates, cadence)
            rng = np.random.default_rng(sigma_seq.spawn(1)[0])
            all_data, admissible_only = Tally(), Tally()
            names = [tp.name for tp in run_templates]
            for _ in range(n_rep):
                intervals, _ = simulate.simulate_emulated_dataset(run_templates, s4, sigma, rule, settings, noise, rng)
                all_data.add(intervals, min_completed, min_pulsars, names)
                if rule == "resampled":
                    admissible_only.add(intervals[intervals["method_class"].isin(admissible)], min_completed,
                                        min_pulsars, names)
            entry = {"sigma": sigma, "method_rule": rule, "unknown_as": unknown_as, "short_window_days": w,
                     "cadence_days": cadence if cadence is not None else fallback,
                     "all_data": all_data.summary()}
            if rule == "short-window":
                forward = entry["all_data"]["forward"]
                admission_entries.append({"variant": name, "w": w, "cadence_days": cadence, "sigma": sigma,
                                          "mean": forward["mean"],
                                          "se_mean": forward["se_mean"]})
            if rule == "resampled":
                entry["admissible_only"] = admissible_only.summary()
            results["variants"][f"{name}|sigma={sigma}"] = entry
            a = entry["all_data"]
            def fmt(x):
                return "   n/a        " if x["mean"] is None else f"{x['mean']:+.3f} ± {x['se_mean']:.3f}"
            print(f"sigma={sigma} {name:26s} forward {fmt(a['forward'])}   backward {fmt(a['backward'])}   "
                  f"dropped {a['pulsars_dropped_total']} discarded {a['replicates_discarded']} failed {a['failed_fits']}")
    record = {"description": "S4 bias magnitudes on a SYNTHETIC fixture (not real pulsars); D4d measurement, D4e rules",
              "seed": seed, "replicates_per_variant_per_sigma": n_rep, "fixture": FIXTURE,
              "config_sha256": cfg["meta.decisions_md_sha256"], "git_commit": git_commit(),
              "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), **results}
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(record, indent=2, default=float) + "\n", encoding="utf-8")
    verdict = simulate.short_window_admission(admission_entries, cfg["D4f.admission.max_upper_95"],
                                              float(stats.norm.ppf(0.975)))
    verdict.update({"rule": "D4f: S4sw forward-slope upper 95% bound (mean + 1.96 SE) <= max at every sigma and w",
                    "fixture": "synthetic", "seed": seed, "git_commit": record["git_commit"],
                    "generated_utc": record["generated_utc"]})
    ADMISSION_PATH.write_text(json.dumps(verdict, indent=2, default=float) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)} and {ADMISSION_PATH.relative_to(ROOT)}: "
          f"short-window {'ADMITTED' if verdict['admitted'] else 'NOT admitted'} to the primary analysis.")


if __name__ == "__main__":
    main()
