"""Post hoc Crab persistent-step calculation (scripts/posthoc_crab_persistent_step.py): a known answer on
synthetic inputs, and exact reproduction of the committed results/exploratory/crab_persistent_step.json.
Both run the script against a temporary project root, so no committed file is read for writing or changed."""
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "posthoc_crab_persistent_step.py"


def _load(root):
    spec = importlib.util.spec_from_file_location("posthoc_crab_persistent_step", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT = Path(root)                      # read inputs from, and write the output to, the temporary root
    return module


def _run(root):
    _load(root).main()
    return (Path(root) / "results/exploratory/crab_persistent_step.json").read_text()


def test_known_answer_on_synthetic_inputs(tmp_path):
    # nu = 30 Hz, nudot = -3e-10 Hz/s: persistent df1 (1e-3) = 7e-8 * dff 1e-9 * 30 / 3e-10 * 1e3 = 0.007 * dff,
    # so dff = 1000, df1 = 28 gives a catalogue-to-persistent ratio of exactly 4 for every Crab interval.
    nu, nudot = 30.0, -3e-10
    tau = nu / (2 * abs(nudot)) / 86400.0
    c_p = np.log(2e-3 * tau)                      # implied beta with the catalogue values = 1 - 1/2 = 0.5
    (tmp_path / "data/processed").mkdir(parents=True)
    (tmp_path / "data/meta").mkdir(parents=True)
    (tmp_path / "results/exploratory").mkdir(parents=True)
    pd.DataFrame({"pulsar": ["B0531+21"] * 3 + ["B0531+21", "J0537-6910"],
                  "included_pulsar": [True, True, True, True, True], "admissible": [True, True, True, False, True],
                  "dff": [1000.0, 1000.0, 1000.0, 5.0, 1000.0], "df1": [28.0, 28.0, 28.0, 99.0, 1.0]}
                 ).to_csv(tmp_path / "data/processed/final_dataset.csv", index=False)
    pd.DataFrame({"name": ["B0531+21"], "nu_hz": [nu], "nudot_hz_s": [nudot]}).to_csv(tmp_path / "data/meta/pulsar_params.csv", index=False)
    (tmp_path / "results/secondary.json").write_text(json.dumps({"implied_beta": [{"pulsar": "B0531+21", "c_p": c_p}]}))
    out = json.loads(_run(tmp_path))
    assert out["label"] == "EXPLORATORY, POST HOC, NOT PRE-REGISTERED"
    assert out["n_intervals"] == 3                                    # the inadmissible and other-pulsar rows are excluded
    assert out["ratio_catalogue_to_persistent"] == pytest.approx({"median": 4.0, "min": 4.0, "max": 4.0})
    assert out["implied_beta"]["catalogue"] == pytest.approx(0.5)
    assert out["implied_beta"]["persistent_step_total"] == pytest.approx(4 * 1.0 - 0.5)          # intercept shift ln 4
    assert out["implied_beta"]["persistent_step_prompt_only"] == pytest.approx(4 / 0.54 - 0.5)   # shift ln(4 / 0.54)


def test_reproduces_the_committed_output_exactly(tmp_path):
    for rel in ("data/processed/final_dataset.csv", "data/meta/pulsar_params.csv", "results/secondary.json"):
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, tmp_path / rel)
    (tmp_path / "results/exploratory").mkdir(parents=True, exist_ok=True)
    committed = (REPO / "results/exploratory/crab_persistent_step.json").read_text()
    assert _run(tmp_path) == committed
