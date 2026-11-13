"""The validator accepts the real transcription and refuses every kind of broken config."""
import copy

import pytest
import yaml

from glitchlib.config import (CONFIG_PATH, DECISIONS_MD_PATH, SCHEMA_PATH, ConfigError, load_config,
                              validate)


def _real():
    with open(CONFIG_PATH, encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    with open(SCHEMA_PATH, encoding="utf-8") as fh:
        schema = yaml.safe_load(fh)
    return config, schema


def _refused(config, schema, decision_id):
    with pytest.raises(ConfigError) as err:
        validate(config, schema)
    assert f"Decision {decision_id}:" in str(err.value)
    return str(err.value)


def test_real_config_is_valid():
    cfg = load_config()
    # Spot checks against sentences in DECISIONS.md.
    assert cfg["D4.per_test_alpha"] == 0.025
    assert cfg["D14.seed"] == 20261031
    assert cfg["D14.sigma_grid"] == [0.2, 0.4, 0.6]
    assert cfg["D12.gate.min_admissible_completed_intervals"] == 40
    assert len(cfg["D13.excluded_objects"]) == 8
    assert cfg["D4a.permutation.n_permutations"] == 10000
    assert cfg["D4a.spearman.per_pulsar_min_completed_admissible"] == 4
    assert cfg["D4b.simulation.validation_permutation.n_permutations"] == 999
    assert cfg["D4b.diagnostics.missingness.method"] == "conditional_logit"


def test_no_provisional_items_remain():
    # D4p item 5 removed the D9, D10 and D11 tags; the D14.exponentials tag was removed because D4o replaces it.
    assert load_config().provisional_items() == []


def test_missing_leaf_names_its_decision():
    config, schema = _real()
    del config["D14"]["seed"]
    message = _refused(config, schema, "D14")
    assert "D14.seed" in message and "missing" in message


def test_missing_nested_leaf_names_its_decision():
    config, schema = _real()
    del config["D12"]["outcomes"]["C_refutation"]["power_sigma"]
    assert "D12.outcomes.C_refutation.power_sigma" in _refused(config, schema, "D12")


def test_null_value_counts_as_missing():
    config, schema = _real()
    config["D4"]["per_test_alpha"] = None
    assert "missing" in _refused(config, schema, "D4")


def test_missing_whole_decision():
    config, schema = _real()
    del config["D7"]
    _refused(config, schema, "D7")


def test_unknown_key_is_refused():
    config, schema = _real()
    config["D8"]["extra_cut"] = 0.3
    assert "not in the schema" in _refused(config, schema, "D8")


def test_unknown_top_level_is_refused():
    config, schema = _real()
    config["D16"] = {"provisional": False}
    with pytest.raises(ConfigError, match="D16"):
        validate(config, schema)


def test_wrong_type_is_refused():
    config, schema = _real()
    config["D14"]["seed"] = "20261031"
    assert "integer" in _refused(config, schema, "D14")


def test_boolean_is_not_a_number():
    config, schema = _real()
    config["D4"]["familywise_alpha"] = True
    _refused(config, schema, "D4")


def test_disallowed_value_is_refused():
    config, schema = _real()
    config["D6"]["primary_errors_in_lnR"] = "deming"
    assert "allowed values" in _refused(config, schema, "D6")


def test_disallowed_list_item_is_refused():
    config, schema = _real()
    config["D9"]["admissible_methods"] = ["local-relaxation", "window-ish"]
    _refused(config, schema, "D9")


def test_out_of_range_is_refused():
    config, schema = _real()
    config["D4"]["per_test_alpha"] = 1.5
    assert "maximum" in _refused(config, schema, "D4")


def test_empty_list_is_refused():
    config, schema = _real()
    config["D14"]["sigma_grid"] = []
    _refused(config, schema, "D14")


def test_changed_decisions_md_is_refused(tmp_path):
    config, schema = _real()
    edited = tmp_path / "DECISIONS.md"
    edited.write_text(DECISIONS_MD_PATH.read_text(encoding="utf-8") + "\nchanged\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="changed since it was transcribed"):
        validate(config, schema, edited)


def test_accessor_refuses_unknown_path():
    cfg = load_config()
    with pytest.raises(ConfigError, match="Decision D14"):
        cfg["D14.number_of_permutations"]


def test_validation_does_not_modify_input():
    config, schema = _real()
    before = copy.deepcopy(config)
    validate(config, schema)
    assert config == before
