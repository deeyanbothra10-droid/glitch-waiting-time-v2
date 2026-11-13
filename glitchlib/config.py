"""Load config/decisions.yaml and check it against config/decisions_schema.yaml.

Every script gets its decisions through load_config(). Nothing here supplies a default: a missing,
extra, mistyped or disallowed value stops the run with an error naming the decision ID.
"""
from pathlib import Path

import yaml

from glitchlib.io import sha256_file

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "decisions.yaml"
SCHEMA_PATH = PROJECT_ROOT / "config" / "decisions_schema.yaml"
DECISIONS_MD_PATH = PROJECT_ROOT / "DECISIONS.md"

DECISION_IDS = ["D1", "D2", "D3", "D4", "D4a", "D4b", "D4c", "D4d", "D4e", "D4f", "D4g", "D4h", "D4i", "D4j", "D4k", "D4l", "D4m", "D4n", "D4o", "D4p"] + [f"D{i}" for i in range(5, 16)]
LEAF_TYPES = {"string", "integer", "number", "boolean", "list"}


class ConfigError(Exception):
    """A decision value is missing, unknown, or not of the type and values the schema allows."""


def _decision_of(path):
    return path.split(".")[0]


def _fail(path, problem):
    raise ConfigError(
        f"Decision {_decision_of(path)}: '{path}' {problem}. Every value must come from DECISIONS.md "
        f"(transcribed into config/decisions.yaml); fix the transcription or ask, never add a default."
    )


def _is_leaf(spec):
    # A schema leaf is a mapping whose "type" is one of LEAF_TYPES; a group must not use "type" as a key.
    return isinstance(spec, dict) and isinstance(spec.get("type"), str) and spec["type"] in LEAF_TYPES


def _check_leaf(value, spec, path):
    kind = spec["type"]
    # bool is a subclass of int in Python, so it is excluded explicitly from the numeric types.
    ok = {
        "string": isinstance(value, str),
        "integer": isinstance(value, int) and not isinstance(value, bool),
        "number": isinstance(value, (int, float)) and not isinstance(value, bool),
        "boolean": isinstance(value, bool),
        "list": isinstance(value, list),
    }[kind]
    if not ok:
        _fail(path, f"must be of type {kind}, got {type(value).__name__} ({value!r})")
    if "allowed" in spec and value not in spec["allowed"]:
        _fail(path, f"has value {value!r}; allowed values are {spec['allowed']}")
    if "min" in spec and value < spec["min"]:
        _fail(path, f"has value {value!r}, below the minimum {spec['min']}")
    if "max" in spec and value > spec["max"]:
        _fail(path, f"has value {value!r}, above the maximum {spec['max']}")
    if kind == "list":
        if len(value) == 0:
            _fail(path, "is an empty list")
        for i, item in enumerate(value):
            _check_leaf(item, spec["items"], f"{path}[{i}]")


def _check_group(values, schema, path):
    if not isinstance(values, dict):
        _fail(path, f"must be a group of keys, got {type(values).__name__}")
    for key in values:
        if key not in schema:
            _fail(f"{path}.{key}", "is not in the schema (DECISIONS.md does not define it)")
    for key, spec in schema.items():
        sub = f"{path}.{key}"
        if key not in values or values[key] is None:
            _fail(sub, "is missing")
        if _is_leaf(spec):
            _check_leaf(values[key], spec, sub)
        else:
            _check_group(values[key], spec, sub)


def validate(config, schema, decisions_md_path=DECISIONS_MD_PATH):
    """Raise ConfigError unless config matches schema exactly and DECISIONS.md is the transcribed version."""
    if not isinstance(config, dict):
        raise ConfigError("config/decisions.yaml is empty or not a mapping of decision IDs.")
    for top in ["meta"] + DECISION_IDS:
        if top not in schema:
            raise ConfigError(f"Decision {top}: missing from config/decisions_schema.yaml.")
        if top not in config:
            raise ConfigError(f"Decision {top}: missing from config/decisions.yaml.")
    for top in config:
        if top not in schema:
            raise ConfigError(f"'{top}' in config/decisions.yaml is not a decision ID ({', '.join(DECISION_IDS)}) or meta.")
    for top in ["meta"] + DECISION_IDS:
        _check_group(config[top], schema[top], top)
    actual = sha256_file(decisions_md_path)
    recorded = config["meta"]["decisions_md_sha256"]
    if actual != recorded:
        raise ConfigError(
            f"DECISIONS.md has changed since it was transcribed (SHA-256 now {actual}, recorded {recorded}). "
            "Re-transcribe config/decisions.yaml from it, then update meta.decisions_md_sha256."
        )


class Config:
    """Validated decisions. Read values with cfg["D14.seed"]; unknown paths raise ConfigError."""

    def __init__(self, data):
        self._data = data

    def __getitem__(self, path):
        node = self._data
        walked = []
        for part in path.split("."):
            walked.append(part)
            if not isinstance(node, dict) or part not in node:
                _fail(".".join(walked), f"is not in config/decisions.yaml (needed as '{path}')")
            node = node[part]
        return node

    def provisional_items(self):
        """Paths of every group tagged provisional: true (items still to be confirmed after the audit)."""
        found = []

        def walk(node, path):
            if isinstance(node, dict):
                if node.get("provisional") is True:
                    found.append(path)
                for key, value in node.items():
                    walk(value, f"{path}.{key}")

        for top in DECISION_IDS:
            walk(self._data[top], top)
        return found

    def as_dict(self):
        return self._data


def load_config(config_path=CONFIG_PATH, schema_path=SCHEMA_PATH, decisions_md_path=DECISIONS_MD_PATH):
    with open(config_path, encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    with open(schema_path, encoding="utf-8") as fh:
        schema = yaml.safe_load(fh)
    validate(config, schema, decisions_md_path)
    return Config(config)
