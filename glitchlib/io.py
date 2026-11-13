"""Readers for the raw catalogue and the data/meta files, with column checks.

Files are only ever read here, never written. Missing-value codes come from the config (D13).
"""
import hashlib
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_GLITCHES_PATH = PROJECT_ROOT / "data" / "raw" / "glitches.csv"
META_DIR = PROJECT_ROOT / "data" / "meta"

# Checksum of data/raw/glitches.csv as supplied (data/raw/README.md).
RAW_GLITCHES_SHA256 = "5ba59d0160fe1f00eb70d464aa048c8b9968980c251118dca6deed825f849de7"

RAW_COLUMNS = ["idx", "name", "jname", "no", "mjd", "mjd_err", "dff", "dff_err", "df1", "df1_err", "ref"]
RAW_NUMERIC = ["idx", "no", "mjd", "mjd_err", "dff", "dff_err", "df1", "df1_err"]

META_FILES = {
    "obs_end": {
        "file": "obs_end.csv",
        "columns": ["name", "jname", "obs_end_mjd", "source", "notes"],
        "numeric": ["obs_end_mjd"],
    },
    "monitoring_gaps": {
        "file": "monitoring_gaps.csv",
        "columns": ["name", "jname", "gap_start_mjd", "gap_end_mjd", "source", "notes"],
        "numeric": ["gap_start_mjd", "gap_end_mjd"],
    },
    "pulsar_params": {
        "file": "pulsar_params.csv",
        "columns": ["name", "jname", "nu_hz", "nudot_hz_s", "source", "catalogue_version"],
        "numeric": ["nu_hz", "nudot_hz_s"],
    },
    "source_audit": {
        "file": "source_audit.csv",
        "columns": ["ref", "n_glitches", "telescope", "cadence_days", "cadence_num_days", "exp_fitted",
                    "df1_definition", "df1_method", "tau_d_days", "Q", "span", "quote", "notes"],
        # cadence_days is the text from the paper (e.g. "2-10"); cadence_num_days is the number the emulation uses.
        "numeric": ["n_glitches", "cadence_num_days"],
    },
}


class DataFormatError(Exception):
    """A file does not have the expected columns or contains an unparseable value."""


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _check_columns(df, expected, path):
    found = list(df.columns)
    if found != expected:
        missing = [c for c in expected if c not in found]
        extra = [c for c in found if c not in expected]
        raise DataFormatError(
            f"{path}: columns {found} do not match the expected {expected} "
            f"(missing {missing}, unexpected {extra}, or a different order)."
        )


def _to_numeric(df, columns, path):
    for col in columns:
        parsed = pd.to_numeric(df[col], errors="coerce")
        bad = df[col][parsed.isna() & df[col].notna()]
        if len(bad):
            raise DataFormatError(f"{path}: column '{col}' has non-numeric values {sorted(set(bad))[:5]}.")
        df[col] = parsed.astype(float)
    return df


def read_raw_glitches(cfg, path=RAW_GLITCHES_PATH):
    """Read the catalogue: columns checked, D13 missing codes -> NaN, numeric columns as float.

    Rows keep the file's order; sorting by mjd (D13) is done by the cleaning step, not here.
    """
    # keep_default_na=False: only the D13 codes mean "missing", so no other string is silently dropped.
    df = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=cfg["D13.missing_codes"])
    _check_columns(df, RAW_COLUMNS, path)
    df = df.astype(object).where(df.notna(), None)
    return _to_numeric(df, RAW_NUMERIC, path)


def read_meta(kind, meta_dir=META_DIR):
    """Read one data/meta file ('obs_end', 'monitoring_gaps', 'pulsar_params', 'source_audit').

    Blank cells become NaN (these files are filled by hand and use blanks, not X, for missing values).
    """
    if kind not in META_FILES:
        raise DataFormatError(f"Unknown meta file '{kind}'; expected one of {sorted(META_FILES)}.")
    spec = META_FILES[kind]
    path = Path(meta_dir) / spec["file"]
    df = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[""])
    _check_columns(df, spec["columns"], path)
    df = df.astype(object).where(df.notna(), None)
    return _to_numeric(df, spec["numeric"], path)


def describe_catalogue(df, key):
    """Structure-only summary for provenance: no statistic involving dff, df1 or waits between glitches."""
    return {
        "n_rows": int(len(df)),
        "n_pulsars": int(df[key].nunique()),
        "pulsar_key": key,
        "columns": list(df.columns),
        "mjd_min": float(df["mjd"].min()),
        "mjd_max": float(df["mjd"].max()),
    }


def git_commit(root=PROJECT_ROOT):
    """Current commit hash, or None outside a git checkout (recorded in every output)."""
    import subprocess
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def repo_relative(path, root=PROJECT_ROOT):
    """A path as recorded in outputs: relative to the repository root, with '/' separators, so that no
    machine-specific absolute path is written. A path outside the repository is written relative to the
    root as well (with '..'). Only the recorded text changes; the path that is read is unchanged."""
    import os
    return Path(os.path.relpath(Path(path).resolve(), Path(root).resolve())).as_posix()


def gap_starts(gaps_table):
    """D11 gap list from data/meta/monitoring_gaps.csv: dict pulsar -> list of gap-start mjds, and the
    rows that cannot be applied because their start date is still missing."""
    dated = gaps_table[gaps_table["gap_start_mjd"].notna()]
    starts = {p: sorted(g["gap_start_mjd"].astype(float)) for p, g in dated.groupby("name")}
    return starts, gaps_table[gaps_table["gap_start_mjd"].isna()]
