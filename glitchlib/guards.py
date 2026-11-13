"""The no-peeking guard.

Nothing may relate R, dnu or dnudot (or their logs) to REAL waiting times unless the
pre-registration is frozen (prereg/FROZEN exists) and the analysis dataset is exactly the one
whose checksum was recorded at the data freeze (prereg/data_freeze.json). This is the only place
that check is written; every script that touches real waiting times calls check_unblinding_allowed()
first. Synthetic data never needs it.
"""
import json
from pathlib import Path

from glitchlib.io import sha256_file

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FROZEN_REL = Path("prereg") / "FROZEN"
FREEZE_RECORD_REL = Path("prereg") / "data_freeze.json"
FINAL_DATASET_REL = Path("data") / "processed" / "final_dataset.csv"


class PeekingError(RuntimeError):
    """Raised when real waiting times are requested before the freeze conditions hold."""


def check_unblinding_allowed(root=PROJECT_ROOT):
    """Return the verified checksum if real waiting times may be analysed; otherwise raise PeekingError."""
    root = Path(root)
    frozen = root / FROZEN_REL
    record = root / FREEZE_RECORD_REL
    dataset = root / FINAL_DATASET_REL
    refuse = "Refusing to analyse real waiting times: "
    if not frozen.is_file():
        raise PeekingError(refuse + f"{FROZEN_REL} does not exist (the pre-registration is not frozen).")
    if not record.is_file():
        raise PeekingError(refuse + f"{FREEZE_RECORD_REL} does not exist (the data are not frozen).")
    try:
        recorded = json.loads(record.read_text(encoding="utf-8"))["sha256"]
    except (json.JSONDecodeError, KeyError, TypeError) as err:
        raise PeekingError(refuse + f"{FREEZE_RECORD_REL} has no readable 'sha256' entry ({err}).")
    if not dataset.is_file():
        raise PeekingError(refuse + f"{FINAL_DATASET_REL} does not exist.")
    actual = sha256_file(dataset)
    if actual != recorded:
        raise PeekingError(
            refuse + f"SHA-256 of {FINAL_DATASET_REL} is {actual}, but {FREEZE_RECORD_REL} records {recorded}."
        )
    return actual


def write_freeze_record(date, git_commit, root=PROJECT_ROOT):
    """Record the SHA-256 of data/processed/final_dataset.csv in prereg/data_freeze.json (the data
    freeze). check_unblinding_allowed later compares the dataset against this record."""
    root = Path(root)
    record = {"file": str(FINAL_DATASET_REL), "sha256": sha256_file(root / FINAL_DATASET_REL),
              "date": date, "git_commit": git_commit}
    (root / FREEZE_RECORD_REL).parent.mkdir(parents=True, exist_ok=True)
    (root / FREEZE_RECORD_REL).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def unblinding_allowed(root=PROJECT_ROOT):
    """True/False form of check_unblinding_allowed, for run_all.py to decide what to skip."""
    try:
        check_unblinding_allowed(root)
    except PeekingError:
        return False
    return True
