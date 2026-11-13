"""Templates from a synthetic cleaned-glitch table (no waiting times are read)."""
import numpy as np
import pandas as pd
import pytest

from glitchlib.config import load_config
from glitchlib.templates import build_templates


def _glitches():
    # A: 5 glitches, 4 admissible (3 local, 1 short-window source elsk11), one unusable.
    # B: 3 glitches, 2 admissible, one of them with a filled (imputed) df1_err.
    rows = [("A", 50000, True, True, "aeka18", "local-relaxation", 0.10, False),
            ("A", 50400, True, True, "aeka18", "local-relaxation", 0.20, False),
            ("A", 50900, False, False, "aeka18", "unknown", np.nan, False),
            ("A", 51500, True, True, "elsk11", "short-window", 0.30, False),
            ("A", 52000, True, True, "zzz99", "local-relaxation", 0.40, False),
            ("B", 51000, True, True, "bsa+21", "global", 0.50, False),
            ("B", 52000, True, True, "bsa+21", "global", 0.25, True),
            ("B", 53000, True, False, "wjm+10", "unknown", 0.60, False)]
    df = pd.DataFrame(rows, columns=["pulsar", "mjd", "usable", "admissible", "source_key", "method_class",
                                     "rel", "df1_err_imputed"])
    df["dff"], df["dff_err"] = 1000.0, 10.0
    df["df1"] = 5.0
    df["df1_err"] = df["rel"] * df["df1"]
    return df.drop(columns="rel")


def _audit():
    return pd.DataFrame({"ref": ["Antonopoulou+2018 [aeka18]", "Espinoza+2011 [elsk11]", "Basu+2021 [bsa+21]"],
                         "cadence_num_days": [14.0, 6.0, 12.0]})


def _gaps(rows=()):
    return pd.DataFrame(list(rows), columns=["name", "gap_start_mjd", "gap_end_mjd"])


def test_templates_known_answer():
    cfg = load_config()
    tps, targets = build_templates(_glitches(), {"A": 53000.0, "B": 53000.0}, _gaps([("A", 50500.0, 50800.0)]),
                                   {"A": 1e6, "B": 2e6}, _audit(), cfg, ["A", "B"], {"A": 10.0, "B": 20.0})
    a, b = tps
    assert (a.first_mjd, a.obs_end_mjd, a.n_glitches) == (50000.0, 53000.0, 5)
    assert a.typical_wait == pytest.approx(600.0)
    assert a.admissible_fraction == pytest.approx(4 / 5) and len(a.pool) == 4
    assert list(a.pool["cadence_days"]) == [14.0, 14.0, 6.0, 14.0]          # zzz99 not audited: D4c fallback
    assert a.pool["short_window_days"].iloc[2] == 40.0 and np.isnan(a.pool["short_window_days"].iloc[0])
    assert a.gaps == ((50500.0, 50800.0),) and not a.drop_censored_last
    assert b.drop_censored_last                                             # obs_end = its last glitch
    assert b.pool["dnudot_frac"].iloc[0] == pytest.approx(5e-3)
    # A has 4 real usable error values: its median 0.25. B has 2 non-imputed usable values (0.5, 0.6;
    # 0.25 is imputed): fewer than 3, so the global median of 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 = 0.35.
    assert targets["A"] == pytest.approx(0.25) and targets["B"] == pytest.approx(0.35)
    # D4o: Δν in µHz = Δν/ν × ν × 1e6
    assert a.pool["dnu_microhz"].tolist() == pytest.approx((a.pool["dnu_frac"] * 10.0 * 1e6).tolist())


def test_obs_end_before_last_glitch_is_refused():
    with pytest.raises(ValueError, match="before its last glitch"):
        build_templates(_glitches(), {"A": 51000.0, "B": 53000.0}, _gaps(), {"A": 1e6, "B": 1e6}, _audit(),
                        load_config(), ["A", "B"], {"A": 10.0, "B": 10.0})
