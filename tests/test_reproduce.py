"""Reproduction comparison (post-freeze tooling) on small known cases."""
import json
import subprocess
import time

from glitchlib import reproduce


def test_json_differences_ignore_commit_timestamp_and_timing_fields():
    old = {"a": 1.0, "git_commit": "x", "generated_utc": "t1", "nested": {"seconds": 3, "b": [1, 2]}}
    new = {"a": 1.0, "git_commit": "y", "generated_utc": "t2", "nested": {"seconds": 9, "b": [1, 3]}}
    assert reproduce.json_differences(old, new) == [(".nested.b[1]", 2, 3)]
    assert reproduce.json_differences({"a": 1}, {"a": 1, "c": 2}) == [(".c", "<absent>", 2)]


def test_pdf_creation_dates_are_ignored():
    a = b"%PDF /CreationDate (D:20261101120000Z) stream"
    b = b"%PDF /CreationDate (D:20261102130000Z) stream"
    assert reproduce.compare_file("f.pdf", a, b) == []
    assert reproduce.compare_file("f.pdf", a, b.replace(b"stream", b"other")) != []
    assert reproduce.compare_file("f.csv", b"x\n1\n", b"x\n1\n") == []
    assert reproduce.compare_file("f.csv", b"x\n1\n", b"x\n2\n") != []


def test_compare_tree_classifies_files(tmp_path):
    (tmp_path / "results").mkdir()
    for name, content in (("same.json", {"v": 1, "generated_utc": "a"}), ("diff.json", {"v": 1}), ("old.json", {"v": 1})):
        (tmp_path / "results" / name).write_text(json.dumps(content))
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-A"], cwd=tmp_path, check=True)
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x"], cwd=tmp_path, check=True)
    cloned_at = time.time()
    time.sleep(1.1)
    (tmp_path / "results" / "same.json").write_text(json.dumps({"v": 1, "generated_utc": "b"}))
    (tmp_path / "results" / "diff.json").write_text(json.dumps({"v": 2}))
    (tmp_path / "results" / "extra.json").write_text("{}")
    report = reproduce.compare_tree(tmp_path, cloned_at, dirs=("results",))
    assert report["reproduced"] == ["results/same.json"]
    assert list(report["differs"]) == ["results/diff.json"]
    assert report["not_regenerated"] == ["results/old.json"] and report["new_files"] == ["results/extra.json"]
