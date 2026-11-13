"""The no-peeking guard refuses unless FROZEN exists and the dataset checksum matches the freeze record."""
import hashlib
import json

import pytest

from glitchlib.guards import PeekingError, check_unblinding_allowed, unblinding_allowed

DATASET = b"name,T\nsynthetic,100\n"


def _project(tmp_path, frozen=True, record=True, dataset=True, recorded_sha=None):
    (tmp_path / "prereg").mkdir()
    (tmp_path / "data" / "processed").mkdir(parents=True)
    if frozen:
        (tmp_path / "prereg" / "FROZEN").write_text("frozen\n")
    if dataset:
        (tmp_path / "data" / "processed" / "final_dataset.csv").write_bytes(DATASET)
    if record:
        sha = recorded_sha or hashlib.sha256(DATASET).hexdigest()
        (tmp_path / "prereg" / "data_freeze.json").write_text(json.dumps({"sha256": sha}))
    return tmp_path


def test_allows_when_frozen_and_checksum_matches(tmp_path):
    root = _project(tmp_path)
    assert check_unblinding_allowed(root) == hashlib.sha256(DATASET).hexdigest()
    assert unblinding_allowed(root)


def test_refuses_without_frozen(tmp_path):
    root = _project(tmp_path, frozen=False)
    with pytest.raises(PeekingError, match="FROZEN does not exist"):
        check_unblinding_allowed(root)
    assert not unblinding_allowed(root)


def test_refuses_when_checksum_differs(tmp_path):
    root = _project(tmp_path, recorded_sha="0" * 64)
    with pytest.raises(PeekingError, match="SHA-256"):
        check_unblinding_allowed(root)


def test_refuses_when_dataset_edited_after_freeze(tmp_path):
    root = _project(tmp_path)
    (root / "data" / "processed" / "final_dataset.csv").write_bytes(DATASET + b"synthetic,200\n")
    with pytest.raises(PeekingError, match="SHA-256"):
        check_unblinding_allowed(root)


def test_refuses_without_freeze_record(tmp_path):
    with pytest.raises(PeekingError, match="data_freeze.json does not exist"):
        check_unblinding_allowed(_project(tmp_path, record=False))


def test_refuses_with_unreadable_freeze_record(tmp_path):
    root = _project(tmp_path)
    (root / "prereg" / "data_freeze.json").write_text('{"date": "2027-01-10"}')
    with pytest.raises(PeekingError, match="sha256"):
        check_unblinding_allowed(root)


def test_refuses_without_dataset(tmp_path):
    with pytest.raises(PeekingError, match="final_dataset.csv does not exist"):
        check_unblinding_allowed(_project(tmp_path, dataset=False))


def test_real_project_is_frozen_with_the_recorded_checksums():
    # Since the freeze (prereg-v2, fe07d83) and the data freeze (0b68712), the guard must allow the real
    # project, return the recorded data-freeze checksum, and FROZEN must record the frozen DECISIONS.md.
    from glitchlib.guards import FREEZE_RECORD_REL, FROZEN_REL, PROJECT_ROOT
    assert unblinding_allowed()
    assert check_unblinding_allowed() == "27a63d37ab877d79e28340bf113e768fe612e290bba1a8f4f2a1ae7ff8bdbe3b"
    from glitchlib.reproduce import reproducing
    if not reproducing():             # a reproduction's 04 records its own clone's commit (run_all.py --reproduce)
        assert json.loads((PROJECT_ROOT / FREEZE_RECORD_REL).read_text(encoding="utf-8"))["git_commit"].startswith("fe07d83")
    frozen = (PROJECT_ROOT / FROZEN_REL).read_text(encoding="utf-8")
    assert "DECISIONS.md SHA-256: 2fd950febe3648b20098cca04a4ebe2742c905c43218d7977fa4894198f00bed" in frozen


def test_freeze_record_written_by_the_pipeline_is_accepted(tmp_path):
    from glitchlib.guards import write_freeze_record
    root = _project(tmp_path, record=False)
    record = write_freeze_record("2027-01-10", "abc123", root)
    assert record["sha256"] == hashlib.sha256(DATASET).hexdigest()
    assert check_unblinding_allowed(root) == record["sha256"]
    (root / "data" / "processed" / "final_dataset.csv").write_bytes(DATASET + b"edited,1\n")
    with pytest.raises(PeekingError, match="SHA-256"):
        check_unblinding_allowed(root)
