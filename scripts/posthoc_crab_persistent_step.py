"""Post hoc, NOT pre-registered: Crab implied beta with a persistent-step Delta-nudot.

Reproduces the numbers quoted in Section 9.2 of the paper. Inputs are committed files only:
  data/processed/final_dataset.csv  (the frozen primary data; Crab admissible intervals)
  data/meta/pulsar_params.csv        (ATNF nu, nudot)
  results/secondary.json             (frozen M1 intercept c_p for B0531+21)
plus one literature constant: |Delta nudot_p| ~= 7e-8 * Delta nu (Lyne et al. 2015, MNRAS 446, 857),
and the fraction 0.54 of the persistent step that appears promptly (same paper: 46 per cent grows
on a ~320-d time-scale).

Method: replacing the catalogue df1 by the persistent-step value multiplies each R by
r = df1_cat / df1_persistent, so the M1 intercept shifts by -mean(ln r). The implied beta is
beta = 2e-3 * tau_sd / exp(c_p) - 1/2, with tau_sd = nu / (2 |nudot|) in days.
Writes results/exploratory/crab_persistent_step.json.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LYNE_COEFF = 7e-8          # |Delta nudot_p| / Delta nu  [s^-1]
PROMPT_FRACTION = 0.54     # share of the persistent step that is instantaneous


def main():
    fd = pd.read_csv(ROOT / "data/processed/final_dataset.csv")
    crab = fd[(fd.pulsar == "B0531+21") & (fd.included_pulsar == True) & (fd.admissible == True)]
    pp = pd.read_csv(ROOT / "data/meta/pulsar_params.csv").set_index("name").loc["B0531+21"]
    nu, nudot = float(pp.nu_hz), float(pp.nudot_hz_s)
    tau_sd_days = nu / (2 * abs(nudot)) / 86400.0
    sec = json.loads((ROOT / "results/secondary.json").read_text())
    c_p = next(b["c_p"] for b in sec["implied_beta"] if b["pulsar"] == "B0531+21")

    dnu_hz = crab.dff * 1e-9 * nu
    df1_persistent = LYNE_COEFF * dnu_hz / abs(nudot) * 1e3      # in units of 1e-3
    ratio = crab.df1 / df1_persistent

    def beta(shift):
        return 2e-3 * tau_sd_days / np.exp(c_p - shift) - 0.5

    out = {
        "label": "EXPLORATORY, POST HOC, NOT PRE-REGISTERED",
        "pulsar": "B0531+21",
        "n_intervals": int(len(crab)),
        "inputs": {"c_p_M1": c_p, "tau_sd_days": tau_sd_days, "nu_hz": nu, "nudot_hz_s": nudot,
                   "lyne_coeff_per_s": LYNE_COEFF, "prompt_fraction": PROMPT_FRACTION},
        "ratio_catalogue_to_persistent": {"median": float(ratio.median()),
                                          "min": float(ratio.min()), "max": float(ratio.max())},
        "implied_beta": {
            "catalogue": float(beta(0.0)),
            "persistent_step_total": float(beta(np.mean(np.log(ratio)))),
            "persistent_step_prompt_only": float(beta(np.mean(np.log(ratio / PROMPT_FRACTION)))),
        },
    }
    dest = ROOT / "results/exploratory/crab_persistent_step.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
