"""Cleaning, source matching, method classes and spin-down timescale on synthetic inputs."""
import numpy as np
import pandas as pd
import pytest

from glitchlib import cleaning
from glitchlib.config import load_config
from glitchlib.io import RAW_COLUMNS


@pytest.mark.parametrize("ref,key", [
    ("Espinoza+2011 [elsk11]", "elsk11"),
    ("Antonopoulou+2018 [aeka18] Also in Ferdman+2018 [fagk18]", "aeka18"),
    ("Basu+2021 [bsa+21]; also in Lower+2021 [ljd+21]", "bsa+21"),
    ("Shaw+2018; also Basu+2019 [bkj+19]", "shaw+2018"),
    ("New - Our measurement. Also in Keith+2024 [kjk+24]", "our measurement"),
    ("New - Our Measurement", "our measurement"),
    ("New - Our measurement.", "our measurement"),
    ("Our measurement", "our measurement"),
    (None, None),
])
def test_source_key_uses_the_primary_reference(ref, key):
    assert cleaning.source_key(ref) == key


@pytest.mark.parametrize("method,exp,expected", [
    ("window", None, "window"), ("local-relaxation", "no", "local-relaxation"),
    ("global", "yes", "global"), ("global", "no", "global-noexp"), ("global", "sometimes", "unknown"),
    ("global", None, "unknown"), (None, None, "unknown"), ("unknown", "yes", "unknown"),
    ("global", "when-detected", "global"), ("short-window", None, "short-window"),
])
def test_audit_method_class_rules(method, exp, expected):
    assert cleaning.audit_method_class(method, exp) == expected


def test_audit_method_class_refuses_unknown_labels():
    with pytest.raises(ValueError, match="df1_method"):
        cleaning.audit_method_class("timing", None)
    with pytest.raises(ValueError, match="exp_fitted"):
        cleaning.audit_method_class("global", "maybe")


def test_audit_spelling_variants_must_agree():
    audit = pd.DataFrame({"ref": ["New - Our measurement", "New - Our Measurement."],
                          "df1_method": ["window", "local-relaxation"], "exp_fitted": [None, None]})
    with pytest.raises(ValueError, match="disagree"):
        cleaning.audit_classes(audit, {})
    audit["df1_method"] = ["window", None]
    assert cleaning.audit_classes(audit, {}) == {"our measurement": "window"}


def test_d4f_assignments_override_a_blank_audit_and_refuse_a_contradiction():
    cfg = load_config()
    assigned = cleaning.assigned_classes(cfg)
    assert assigned == {"elsk11": "short-window", "wmp+00": "short-window", "our measurement": "global"}
    audit = pd.DataFrame({"ref": ["Espinoza+2011 [elsk11]", "New - Our measurement", "Basu+2021 [bsa+21]"],
                          "df1_method": [None, None, "local-relaxation"], "exp_fitted": [None, None, None]})
    assert cleaning.audit_classes(audit, assigned) == {"elsk11": "short-window", "our measurement": "global",
                                                       "bsa+21": "local-relaxation", "wmp+00": "short-window"}
    audit.loc[0, "df1_method"] = "window"
    with pytest.raises(ValueError, match="D4f"):
        cleaning.audit_classes(audit, assigned)


def test_admissible_methods_follow_the_admission_record():
    cfg = load_config()
    assert cleaning.admissible_methods(cfg, {"admitted": False}) == ["local-relaxation", "global"]
    assert cleaning.admissible_methods(cfg, {"admitted": True}) == ["local-relaxation", "global", "short-window"]


def test_tau_sd_known_answers():
    year = 365.25
    assert cleaning.tau_sd_days(11.19, -1.56e-11) / year == pytest.approx(11_400, rel=0.01)    # Vela
    assert cleaning.tau_sd_days(29.95, -3.77e-10) / year == pytest.approx(1_259, rel=0.01)     # Crab
    with pytest.raises(ValueError, match="negative"):
        cleaning.tau_sd_days([11.19], [1e-11])


def _catalogue():
    # Pulsar A: 6 glitches, 5 usable, all from an admissible (local-relaxation) source.
    # Pulsar B: 4 usable glitches from a window source. Object M82-X2 is on the D13 exclusion list.
    rows = []
    for i, df1 in enumerate([5.0, 6.0, None, 4.0, 7.0, 3.0]):
        rows.append([i, "A", "A", i + 1, 50000 + 500 * i, 1.0, 1000.0, 10.0, df1, 0.1, "Paper+2020 [loc20]"])
    for i, df1 in enumerate([5.0, -1.0, 4.0, 4.0, 2.0]):
        rows.append([10 + i, "B", "B", i + 1, 50000 + 700 * i, 1.0, 900.0, 9.0, df1, 0.1, "Win+2010 [win10]"])
    rows.append([20, "M82-X2", "X", 1, 55000, 1.0, 5.0, 1.0, 1.0, 0.1, "Other [oth]"])
    df = pd.DataFrame(rows, columns=RAW_COLUMNS)
    df.loc[df["df1"].isna(), "df1"] = np.nan
    return df


def test_clean_steps_and_log():
    cfg = load_config()
    audit = pd.DataFrame({"ref": ["Paper+2020 [loc20]", "Win+2010 [win10]"],
                          "df1_method": ["local-relaxation", "window"], "exp_fitted": [None, None]})
    glitches, log = cleaning.clean(_catalogue(), cfg, audit, {}, ["local-relaxation", "global"])
    steps = log.set_index("step")
    assert list(log["rows"])[:2] == [12, 11]
    assert log.iloc[0]["df1_missing"] == 1 and log.iloc[0]["df1_le_0"] == 1
    assert log.iloc[1]["objects_removed"] == 1
    assert log.iloc[2]["usable"] == 9
    assert log.iloc[3]["df1_err_filled"] == 0
    assert log.iloc[4]["admissible"] == 5 and log.iloc[4]["usable_window"] == 4
    # A: admissible glitches starting a completed interval = 4 (the last usable one is last overall).
    assert log.iloc[5]["pulsars"] == 1 and log.iloc[5]["admissible_completed_in_included"] == 4
    assert set(glitches.loc[glitches["included_pulsar"], "pulsar"]) == {"A"}
    assert len(glitches) == 11                       # unusable rows are kept: they still end intervals


def test_fill_relative_errors_known_answer():
    # Pulsar A: 4 usable rows with relative errors 0.1, 0.2, 0.3 and one missing -> A's median 0.2.
    # Pulsar B: 1 known (0.5) and 1 missing -> fewer than 3 values -> median of all known usable rows
    # (0.1, 0.2, 0.3, 0.5) = 0.25. The unusable row is neither filled nor used for the medians.
    table = pd.DataFrame({"pulsar": ["A"] * 5 + ["B"] * 2,
                          "df1": [1.0, 2.0, 1.0, 4.0, 1.0, 2.0, 2.0],
                          "df1_err": [0.1, 0.4, 0.3, np.nan, 9.0, 1.0, np.nan],
                          "usable": [True, True, True, True, False, True, True]})
    out = cleaning.fill_relative_errors(table, [("df1", "df1_err")], table["usable"], 3)
    assert out.loc[3, "df1_err"] == pytest.approx(0.2 * 4.0)
    assert out.loc[6, "df1_err"] == pytest.approx(0.25 * 2.0)
    assert list(out["df1_err_imputed"]) == [False, False, False, True, False, False, True]
    assert out.loc[4, "df1_err"] == 9.0 and np.isnan(table.loc[3, "df1_err"])


def test_source_cadences_from_the_audit():
    audit = pd.DataFrame({"ref": ["Espinoza+2011 [elsk11]", "New - Our measurement", "Our measurement", "X [x1]"],
                          "cadence_num_days": [6.0, 12.0, 12.0, np.nan]})
    assert cleaning.source_cadences(audit) == {"elsk11": 6.0, "our measurement": 12.0}
    audit.loc[2, "cadence_num_days"] = 14.0
    with pytest.raises(ValueError, match="cadences"):
        cleaning.source_cadences(audit)


def test_real_audit_classes_and_cadences():
    # The audited sources map to the classes and cadences written in data/meta/source_audit.csv.
    from glitchlib.io import read_meta
    cfg = load_config()
    audit = read_meta("source_audit")
    classes = cleaning.audit_classes(audit, cleaning.assigned_classes(cfg))
    expected = {"elsk11": "short-window", "aeka18": "local-relaxation", "bsa+21": "global",
                "hke+26": "local-relaxation", "ljd+21": "global", "ymh+12": "local-relaxation",
                "wmp+00": "short-window", "ggyz22": "global", "our measurement": "global"}
    assert {k: classes[k] for k in expected} == expected
    assert {k: v for k, v in classes.items() if k not in expected} == {k: "unknown" for k in classes if k not in expected}
    cadences = cleaning.source_cadences(audit)
    assert cadences == {"elsk11": 6.0, "aeka18": 14.0, "bsa+21": 12.0, "hke+26": 14.0, "ljd+21": 30.0,
                        "ymh+12": 21.0, "wmp+00": 28.0, "ggyz22": 15.0, "our measurement": 12.0}


def test_jodrell_bank_rule_known_answer():
    cfg = load_config()
    glitches = pd.DataFrame({"pulsar": ["A", "A", "B", "C", "D", "E", "F"],
                             "ref": ["Yu+2012 [ymh+12]", "Espinoza+2011 [elsk11]",
                                     "Shaw+2018 (ATel #11625); also Basu+2019 [bkj+19]",
                                     "Lower+2021 [ljd+21]; also in Basu+2021 [bsa+21]",   # JB only as "also"
                                     "Our measurement", "Basu+2019 [bkj+19]", None]})
    assert cleaning.jodrell_bank_pulsars(glitches, cfg) == ["A", "B", "D", "E"]


def test_jodrell_bank_rule_reproduces_the_recorded_list():
    # D4g records the included pulsars the rule gave on 31 Oct; the code must reproduce them.
    import pathlib
    cfg = load_config()
    path = pathlib.Path(__file__).resolve().parents[1] / "data" / "processed" / "glitches_clean.csv"
    glitches = pd.read_csv(path)
    included = set(glitches.loc[glitches["included_pulsar"], "pulsar"])
    jb = [p for p in cleaning.jodrell_bank_pulsars(glitches, cfg) if p in included]
    assert jb == cfg["D4g.robustness_item_10.jodrell_bank_rule.included_pulsars_computed_31_oct"]
