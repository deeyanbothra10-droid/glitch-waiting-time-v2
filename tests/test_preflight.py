"""Preflight of the real-run scripts refuses when any condition fails."""
import hashlib
import json
import subprocess

import pytest

from glitchlib import preflight


@pytest.fixture(autouse=True)
def _normal_mode(monkeypatch):
    # These tests check the normal preflight even when the suite itself runs inside run_all.py --reproduce
    # (GLITCH_REPRODUCE=1 inherited from the environment); the reproduction-mode test sets it itself.
    monkeypatch.delenv("GLITCH_REPRODUCE", raising=False)


def _project(tmp_path, frozen=True, validation="pass", output=False):
    (tmp_path / "prereg").mkdir()
    (tmp_path / "data" / "processed").mkdir(parents=True)
    (tmp_path / "results").mkdir()
    data = b"pulsar,t\nX,1\n"
    (tmp_path / "data" / "processed" / "final_dataset.csv").write_bytes(data)
    (tmp_path / "prereg" / "data_freeze.json").write_text(json.dumps({"sha256": hashlib.sha256(data).hexdigest()}))
    if frozen:
        (tmp_path / "prereg" / "FROZEN").write_text("frozen\n")
    if validation:
        result = "PASS" if validation == "pass" else "FAIL"
        (tmp_path / "results" / "validation.json").write_text(
            json.dumps({"criteria": [{"name": "c1", "result": "PASS"}, {"name": "c2", "result": result}]}))
    if output:
        (tmp_path / "results" / "primary.json").write_text("{}")
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-A"], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x"], cwd=tmp_path,
                   check=True)
    return tmp_path


def test_all_conditions_hold(tmp_path):
    preflight.require_all("results/primary.json", _project(tmp_path), run_tests=False)


@pytest.mark.parametrize("kwargs,reason", [({"frozen": False}, "FROZEN"), ({"validation": "fail"}, "not all PASS"),
                                           ({"validation": None}, "validation.json does not exist"),
                                           ({"output": True}, "already exists")])
def test_refuses(tmp_path, kwargs, reason):
    root = _project(tmp_path, **kwargs)
    with pytest.raises(SystemExit, match=reason):
        preflight.require_all("results/primary.json", root, run_tests=False)


class _Cfg:
    def __init__(self, tagged):
        self.tagged = tagged

    def provisional_items(self):
        return self.tagged


def test_d4p_refuses_provisional_decisions(tmp_path):
    root = _project(tmp_path)
    preflight.require_all("results/primary.json", root, run_tests=False, refuse_provisional_cfg=_Cfg([]))
    with pytest.raises(SystemExit, match="still tagged provisional: D14.exponentials"):
        preflight.require_all("results/primary.json", root, run_tests=False, refuse_provisional_cfg=_Cfg(["D14.exponentials"]))


def test_refuses_dirty_tree_and_changed_data(tmp_path):
    root = _project(tmp_path)
    (root / "data" / "processed" / "final_dataset.csv").write_bytes(b"pulsar,t\nX,2\n")
    with pytest.raises(SystemExit) as err:
        preflight.require_all("results/primary.json", root, run_tests=False)
    assert "not clean" in str(err.value) and "SHA-256" in str(err.value)


def test_reproduction_mode_skips_only_clean_tree_and_run_once(tmp_path, monkeypatch):
    root = _project(tmp_path, output=True)
    (root / "untracked.txt").write_text("dirty\n")
    with pytest.raises(SystemExit) as err:
        preflight.require_all("results/primary.json", root, run_tests=False)
    assert "not clean" in str(err.value) and "already exists" in str(err.value)
    monkeypatch.setenv("GLITCH_REPRODUCE", "1")
    preflight.require_all("results/primary.json", root, run_tests=False)            # both skipped
    (root / "prereg" / "FROZEN").unlink()
    with pytest.raises(SystemExit, match="FROZEN"):                                    # the guard still applies
        preflight.require_all("results/primary.json", root, run_tests=False)
