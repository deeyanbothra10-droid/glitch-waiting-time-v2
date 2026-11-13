"""Readers check columns, map only the D13 missing codes to missing, and parse numbers."""
import hashlib
import math

import pytest

from glitchlib.config import load_config
from glitchlib.io import (META_FILES, RAW_COLUMNS, DataFormatError, describe_catalogue, read_meta,
                          read_raw_glitches, repo_relative, sha256_file)

HEADER = ",".join(RAW_COLUMNS)
ROWS = [
    "0,B0000+00,0000+0000,1,50000.5,0.1,2100,10,6.0,0.5,ref1",
    "1,B0000+00,0000+0000,2,51000,X,30,x,X,X,ref1",
    "2,J1111-1111,X,1,49000,2,1800,5,4.0,0.2,ref2",
]


@pytest.fixture
def cfg():
    return load_config()


def _write(tmp_path, lines, name="glitches.csv"):
    path = tmp_path / name
    path.write_text("\r\n".join(lines) + "\r\n")
    return path


def test_sha256_known_answer(tmp_path):
    path = tmp_path / "abc.txt"
    path.write_bytes(b"abc")
    assert sha256_file(path) == hashlib.sha256(b"abc").hexdigest()
    assert sha256_file(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_repo_relative_known_answers(tmp_path):
    # Known answers on a synthetic layout: inside the root -> plain relative path; the root itself -> ".";
    # outside the root -> relative with "..". No absolute path is ever returned.
    root = tmp_path / "repo"
    (root / "results" / "store_a").mkdir(parents=True)
    (tmp_path / "elsewhere").mkdir()
    assert repo_relative(root / "results" / "store_a", root) == "results/store_a"
    assert repo_relative(str(root / "results" / "store_a"), root) == "results/store_a"
    assert repo_relative(root, root) == "."
    assert repo_relative(tmp_path / "elsewhere", root) == "../elsewhere"
    assert not repo_relative(root / "results" / "store_a", root).startswith("/")


def test_raw_reader_parses_synthetic_file(tmp_path, cfg):
    df = read_raw_glitches(cfg, _write(tmp_path, [HEADER] + ROWS))
    assert list(df.columns) == RAW_COLUMNS
    assert df.loc[0, "dff"] == 2100.0 and df.loc[0, "df1"] == 6.0
    assert math.isnan(df.loc[1, "mjd_err"]) and math.isnan(df.loc[1, "df1"])
    assert math.isnan(df.loc[1, "dff_err"])          # lower-case x is also missing
    assert df.loc[2, "jname"] is None                 # X in a text column is missing too
    assert df.loc[0, "name"] == "B0000+00"


def test_raw_reader_keeps_file_order(tmp_path, cfg):
    df = read_raw_glitches(cfg, _write(tmp_path, [HEADER] + ROWS))
    assert list(df["mjd"]) == [50000.5, 51000.0, 49000.0]


def test_raw_reader_rejects_wrong_columns(tmp_path, cfg):
    bad_header = HEADER.replace("df1_err", "df1_error")
    with pytest.raises(DataFormatError, match="df1_err"):
        read_raw_glitches(cfg, _write(tmp_path, [bad_header] + ROWS))


def test_raw_reader_rejects_unknown_code(tmp_path, cfg):
    rows = ROWS + ["3,B0000+00,0000+0000,3,52000,1,?,1,1.0,0.1,ref1"]
    with pytest.raises(DataFormatError, match="'dff'"):
        read_raw_glitches(cfg, _write(tmp_path, [HEADER] + rows))


def test_describe_catalogue_known_answer(tmp_path, cfg):
    df = read_raw_glitches(cfg, _write(tmp_path, [HEADER] + ROWS))
    summary = describe_catalogue(df, "name")
    assert summary["n_rows"] == 3 and summary["n_pulsars"] == 2
    assert summary["mjd_min"] == 49000.0 and summary["mjd_max"] == 51000.0


def test_meta_reader_blank_is_missing(tmp_path):
    lines = ["name,jname,obs_end_mjd,source,notes", "B0000+00,0000+0000,55000,paper,", "J1111-1111,X,,,"]
    _write(tmp_path, lines, META_FILES["obs_end"]["file"])
    df = read_meta("obs_end", tmp_path)
    assert df.loc[0, "obs_end_mjd"] == 55000.0
    assert math.isnan(df.loc[1, "obs_end_mjd"])


def test_meta_reader_rejects_wrong_columns(tmp_path):
    _write(tmp_path, ["name,obs_end_mjd", "B0000+00,55000"], META_FILES["obs_end"]["file"])
    with pytest.raises(DataFormatError):
        read_meta("obs_end", tmp_path)


def test_real_files_have_expected_structure(cfg):
    # Structure only: counts stated in the brief, no statistic on dff, df1 or waits.
    df = read_raw_glitches(cfg)
    assert len(df) == 736 and df["name"].nunique() == 226
    for kind in META_FILES:
        read_meta(kind)
