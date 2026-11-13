# Forward forecast v2 (post-freeze; NOT part of the pre-registration)

Written by `scripts/forecast_v2.py`; the forecasts are in `forecast_v2.json`. Forecast v1 (`forecast_v1.json`, `README.md`, tag `forecast-v1`) is kept unchanged.

## Why v1 was superseded

v1 conditioned every pulsar on no glitch before the 2026-10-30 snapshot. That is valid only where monitoring is known to continue up to the snapshot, which is documented only for the pulsars of the ongoing Jodrell Bank timing programme (the D4g Jodrell Bank rule). For the others, a missing glitch after the documented end of monitoring may mean no one was looking. v2 was made the same day as v1, before any outcome was known.

## Scored set

11 pulsars are scored (forecastable and covered by the D4g Jodrell Bank rule). With this few pulsars and windows this wide, the score is a weak test; it is reported whatever it shows.

## Rule (fixed before any outcome)

1. Label: post-freeze forward forecast v2; NOT part of the pre-registration (prereg-v2). Supersedes forecast v1 (tag forecast-v1, kept unchanged) before any outcome was known.
2. Snapshot: the Jodrell Bank glitch catalogue downloaded 2026-10-30 (MJD 61343); for each of the 21 included pulsars, its last catalogued glitch g. Pulsars whose g lacks dff > 0 or df1 > 0 are not forecastable.
3. Forecast: ln T ~ Normal(mu, sigma), T = wait after g in days, from the frozen fits in results/primary.json: M2 mu = c_p + b ln R_g with R_g = (dff/1000)/df1; M0 (baseline) mu = c_p; each model's own sigma-hat. Median and central 50% (0.25-0.75) and 90% (0.05-0.95) windows of the conditional distribution.
4. Scored set: only the forecastable pulsars covered by the D4g Jodrell Bank rule (ongoing Jodrell Bank timing programme), conditioned on T > t0 = 61343 - epoch(g).
5. Reported, not scored: every other forecastable pulsar, conditioned on T > documented obs_end - epoch(g), or unconditioned if obs_end equals epoch(g); labelled 'reported, not scored: monitoring after obs_end not documented'.
6. Scoring, on 2028-11-01 (MJD 62076): download the catalogue again. For each scored pulsar the outcome is the first glitch after g: the first catalogued glitch of that pulsar whose epoch exceeds g's epoch by more than g's catalogued epoch uncertainty (0 if none is given).
7. Void: if that glitch's epoch is before MJD 61343, the pulsar is void and not scored.
8. Score per scored pulsar: log predictive density of the observed wait (density of T in days) under the conditional forecast, or, if no glitch by MJD 62076, log of the conditional survival at the wait to that date; M2 minus M0. Total score = sum over scored pulsars; positive means M2 forecasts better.
9. Coverage: for the 50% and 90% windows, the fraction of scored pulsars whose observed wait lies inside the window (ends included). A pulsar with no glitch by the scoring date counts as outside if the window ends before the scoring date, and as not determinable (excluded from that coverage) otherwise; the counts are reported.
10. Functions: glitchlib.forecast (conditional_quantiles, log_predictive, coverage, outcome_for, forecast_window).

## Scored forecasts

| Pulsar | Last glitch g | Conditioning | Model | Median | 50% window | 90% window |
|---|---|---|---|---|---|---|
| B0531+21 | 2025-08-06 | no glitch before MJD 61343 (snapshot) | M2 | 2027-05-26 | 2027-01-23 to 2028-01-04 | 2026-11-14 to 2029-08-31 |
| B0531+21 | 2025-08-06 | no glitch before MJD 61343 (snapshot) | M0 | 2027-07-18 | 2027-02-12 to 2028-05-05 | 2026-11-17 to 2030-08-18 |
| B1757-24 | 2019-03-22 | no glitch before MJD 61343 (snapshot) | M2 | 2030-01-14 | 2028-02-19 to 2033-07-06 | 2027-01-22 to 2042-12-03 |
| B1757-24 | 2019-03-22 | no glitch before MJD 61343 (snapshot) | M0 | 2030-01-06 | 2028-02-08 to 2033-08-15 | 2027-01-19 to 2043-12-01 |
| B1800-21 | 2020-04-19 | no glitch before MJD 61343 (snapshot) | M2 | 2030-11-23 | 2028-07-08 to 2035-02-19 | 2027-02-18 to 2046-07-06 |
| B1800-21 | 2020-04-19 | no glitch before MJD 61343 (snapshot) | M0 | 2030-06-21 | 2028-04-18 to 2034-07-20 | 2027-02-01 to 2046-03-18 |
| B1823-13 | 2020-01-31 | no glitch before MJD 61343 (snapshot) | M2 | 2031-05-17 | 2028-09-22 to 2036-01-26 | 2027-03-04 to 2048-08-17 |
| B1823-13 | 2020-01-31 | no glitch before MJD 61343 (snapshot) | M0 | 2032-04-24 | 2029-01-29 to 2038-04-13 | 2027-03-26 to 2055-04-11 |
| J0205+6449 | 2018-07-21 | no glitch before MJD 61343 (snapshot) | M2 | 2028-04-22 | 2027-06-05 to 2029-12-04 | 2026-12-07 to 2034-04-10 |
| J0205+6449 | 2018-07-21 | no glitch before MJD 61343 (snapshot) | M0 | 2028-07-29 | 2027-07-13 to 2030-07-12 | 2026-12-13 to 2035-12-15 |
| J0631+1036 | 2018-08-22 | no glitch before MJD 61343 (snapshot) | M2 | 2028-04-05 | 2027-05-29 to 2029-10-30 | 2026-12-06 to 2034-01-13 |
| J0631+1036 | 2018-08-22 | no glitch before MJD 61343 (snapshot) | M0 | 2028-06-10 | 2027-06-23 to 2030-03-29 | 2026-12-10 to 2035-03-21 |
| J1737-3137 | 2022-09-26 | no glitch before MJD 61343 (snapshot) | M2 | 2028-12-08 | 2027-09-11 to 2031-03-08 | 2026-12-25 to 2037-03-29 |
| J1737-3137 | 2022-09-26 | no glitch before MJD 61343 (snapshot) | M0 | 2029-04-20 | 2027-10-30 to 2032-01-23 | 2027-01-02 to 2039-12-11 |
| J1814-1744 | 2011-03-06 | no glitch before MJD 61343 (snapshot) | M2 | 2030-07-25 | 2028-05-02 to 2034-09-01 | 2027-02-03 to 2045-10-16 |
| J1814-1744 | 2011-03-06 | no glitch before MJD 61343 (snapshot) | M0 | 2030-04-05 | 2028-03-15 to 2034-02-06 | 2027-01-25 to 2044-10-21 |
| J1841-0524 | 2020-09-19 | no glitch before MJD 61343 (snapshot) | M2 | 2028-11-26 | 2027-09-02 to 2031-03-05 | 2026-12-23 to 2037-05-06 |
| J1841-0524 | 2020-09-19 | no glitch before MJD 61343 (snapshot) | M0 | 2028-09-18 | 2027-08-01 to 2030-11-07 | 2026-12-17 to 2036-11-24 |
| J2021+3651 | 2023-12-11 | no glitch before MJD 61343 (snapshot) | M2 | 2029-10-04 | 2028-02-08 to 2032-08-06 | 2027-01-27 to 2040-02-02 |
| J2021+3651 | 2023-12-11 | no glitch before MJD 61343 (snapshot) | M0 | 2029-08-09 | 2027-12-25 to 2032-07-28 | 2027-01-14 to 2040-12-29 |
| J2229+6114 | 2026-01-31 | no glitch before MJD 61343 (snapshot) | M2 | 2028-07-31 | 2027-09-13 to 2030-01-08 | 2027-01-16 to 2033-09-30 |
| J2229+6114 | 2026-01-31 | no glitch before MJD 61343 (snapshot) | M0 | 2028-06-04 | 2027-07-23 to 2029-12-09 | 2026-12-25 to 2034-03-09 |

## Reported, not scored (reported, not scored: monitoring after obs_end not documented)

| Pulsar | Last glitch g | Conditioning | Model | Median | 50% window | 90% window |
|---|---|---|---|---|---|---|
| B1046-58 | 2024-04-06 | unconditioned (documented obs_end equals g) | M2 | 2026-08-18 | 2025-09-28 to 2028-01-17 | 2025-01-06 to 2031-09-11 |
| B1046-58 | 2024-04-06 | unconditioned (documented obs_end equals g) | M0 | 2026-06-01 | 2025-07-12 to 2027-12-03 | 2024-11-07 to 2032-02-12 |
| J0537-6910 | 2025-03-28 | no glitch before documented obs_end MJD 60827.0 | M2 | 2025-08-02 | 2025-06-27 to 2025-10-01 | 2025-06-06 to 2026-03-11 |
| J0537-6910 | 2025-03-28 | no glitch before documented obs_end MJD 60827.0 | M0 | 2025-08-11 | 2025-07-01 to 2025-10-25 | 2025-06-06 to 2026-05-25 |
| J1023-5746 | 2019-12-06 | no glitch before documented obs_end MJD 59269.0 | M2 | 2022-03-19 | 2021-08-09 to 2023-04-10 | 2021-03-27 to 2026-01-30 |
| J1023-5746 | 2019-12-06 | no glitch before documented obs_end MJD 59269.0 | M0 | 2022-06-10 | 2021-09-09 to 2023-10-21 | 2021-04-01 to 2027-09-02 |
| J1420-6048 | 2019-03-13 | unconditioned (documented obs_end equals g) | M2 | 2020-12-20 | 2020-04-21 to 2022-01-12 | 2019-10-05 to 2024-10-09 |
| J1420-6048 | 2019-03-13 | unconditioned (documented obs_end equals g) | M0 | 2022-02-01 | 2020-11-23 to 2024-02-11 | 2019-12-27 to 2029-10-03 |
| J1617-5055 | 2012-12-06 | unconditioned (documented obs_end equals g) | M2 | 2014-05-26 | 2013-11-06 to 2015-04-12 | 2013-05-25 to 2017-07-17 |
| J1617-5055 | 2012-12-06 | unconditioned (documented obs_end equals g) | M0 | 2013-11-14 | 2013-06-25 to 2014-07-12 | 2013-03-10 to 2016-05-12 |

Note, J0537-6910: its glitch rate (one glitch every 100-200 d) makes a 580 d glitch-free stretch to the 2026-10-30 snapshot implausible, so the snapshot is not used for conditioning.

## Not forecastable

- B1338-62: last glitch lacks dff and df1 > 0 (g at MJD 58214.0)
- B1727-47: last glitch lacks df1 > 0 (g at MJD 57984.0)
- B1737-30: last glitch lacks df1 > 0 (g at MJD 60980.0)
- B1758-23: last glitch lacks df1 > 0 (g at MJD 58993.5)
- J0729-1448: last glitch lacks df1 > 0 (g at MJD 58266.6)
