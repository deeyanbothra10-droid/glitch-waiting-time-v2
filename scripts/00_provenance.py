"""Record where the raw catalogue came from and what it contains -> data/processed/provenance.json.

Structure only: checksum, row and pulsar counts, columns, mjd range. No statistic involving dff,
df1 or the time between glitches is computed. The raw file is only read.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from glitchlib.config import load_config  # noqa: E402
from glitchlib.io import (RAW_GLITCHES_PATH, RAW_GLITCHES_SHA256, describe_catalogue,  # noqa: E402
                          git_commit, read_raw_glitches, sha256_file)

DOWNLOAD_DATE = "2026-10-30"
SOURCE_URL = "https://www.jb.man.ac.uk/pulsar/glitches/gTable.html"
OUT_PATH = ROOT / "data" / "processed" / "provenance.json"


def main():
    cfg = load_config()
    checksum = sha256_file(RAW_GLITCHES_PATH)
    if checksum != RAW_GLITCHES_SHA256:
        sys.exit(f"STOP: SHA-256 of {RAW_GLITCHES_PATH} is {checksum}, "
                 f"not the supplied {RAW_GLITCHES_SHA256}. The raw file has changed.")
    df = read_raw_glitches(cfg)
    record = {
        "file": str(RAW_GLITCHES_PATH.relative_to(ROOT)),
        "sha256": checksum,
        "sha256_matches_supplied": True,
        **describe_catalogue(df, cfg["D13.pulsar_key"]),
        "download_date": DOWNLOAD_DATE,
        "source_url": SOURCE_URL,
        "random_seed": None,  # this script uses no random numbers
        "git_commit": git_commit(),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}: {record['n_rows']} rows, {record['n_pulsars']} pulsars, "
          f"mjd {record['mjd_min']:.0f}-{record['mjd_max']:.0f}, SHA-256 matches the supplied file.")


if __name__ == "__main__":
    main()
