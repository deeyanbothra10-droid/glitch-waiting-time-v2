"""Post-freeze tooling: compare a reproduction (run_all.py --reproduce DIR, a fresh clone in which the whole
pipeline was rerun) with the committed files of the same commit.

Rules: JSON files must match exactly except fields named in JSON_IGNORE (commit, timestamp and timing fields);
CSV and other text files byte for byte; PDFs byte for byte once their /CreationDate entries are removed.
A tracked file the run did not rewrite (modification time not after the clone) is listed as "not regenerated".
"""
import json
import os
import re
import subprocess
from pathlib import Path

REPRODUCE_ENV = "GLITCH_REPRODUCE"          # set to "1" by run_all.py --reproduce for the pipeline it runs
JSON_IGNORE = {"git_commit", "generated_utc", "seconds", "seconds_max", "date"}
COMPARED_DIRS = ("results", "data/processed", "figures", "prereg", "logs")


def reproducing():
    return os.environ.get(REPRODUCE_ENV) == "1"


def json_differences(old, new, path="", ignore=JSON_IGNORE):
    """Every differing leaf as (path, old, new), skipping keys in ignore at any depth."""
    out = []
    if isinstance(old, dict) and isinstance(new, dict):
        for key in sorted(set(old) | set(new), key=str):
            if key in ignore:
                continue
            if key not in old or key not in new:
                out.append((f"{path}.{key}", old.get(key, "<absent>"), new.get(key, "<absent>")))
            else:
                out.extend(json_differences(old[key], new[key], f"{path}.{key}", ignore))
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            out.append((path, f"<{len(old)} items>", f"<{len(new)} items>"))
        for i, (a, b) in enumerate(zip(old, new)):
            out.extend(json_differences(a, b, f"{path}[{i}]", ignore))
    elif old != new and not (isinstance(old, float) and isinstance(new, float) and old != old and new != new):
        out.append((path, old, new))
    return out


def strip_pdf_dates(data):
    return re.sub(rb"/CreationDate *\(D:[^)]*\)", b"", data)


def compare_file(relpath, old_bytes, new_bytes):
    """Differences for one file: a list (empty when it reproduces)."""
    if relpath.endswith(".json"):
        return json_differences(json.loads(old_bytes), json.loads(new_bytes))
    if relpath.endswith(".pdf"):
        return [] if strip_pdf_dates(old_bytes) == strip_pdf_dates(new_bytes) else [("", "pdf content differs", "")]
    return [] if old_bytes == new_bytes else [("", "bytes differ", "")]


def compare_tree(directory, cloned_at, dirs=COMPARED_DIRS):
    """Compare every tracked file under dirs in the clone `directory` with its committed version (HEAD)."""
    directory = Path(directory)
    tracked = subprocess.run(["git", "ls-files", "--", *dirs], cwd=directory, capture_output=True, text=True,
                             check=True).stdout.split()
    untracked = subprocess.run(["git", "ls-files", "--others", "--exclude-standard", "--", *dirs], cwd=directory,
                               capture_output=True, text=True, check=True).stdout.split()
    report = {"reproduced": [], "differs": {}, "not_regenerated": [], "new_files": untracked}
    for rel in tracked:
        if rel.endswith(".gitkeep"):
            continue
        path = directory / rel
        if path.stat().st_mtime <= cloned_at:
            report["not_regenerated"].append(rel)
            continue
        old = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=directory, capture_output=True, check=True).stdout
        diffs = compare_file(rel, old, path.read_bytes())
        if diffs:
            report["differs"][rel] = [[str(p), repr(a), repr(b)] for p, a, b in diffs]
        else:
            report["reproduced"].append(rel)
    return report
