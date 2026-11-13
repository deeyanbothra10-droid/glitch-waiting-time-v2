# Forward forecast v1 (post-freeze; NOT part of the pre-registration)

Written by `scripts/forecast_v1.py`; the forecasts are in `forecast_v1.json`.

## Rule (fixed before any outcome)

1. Label: post-freeze forward forecast v1; NOT part of the pre-registration (prereg-v2).
2. Snapshot: the Jodrell Bank glitch catalogue downloaded 2026-10-30 (MJD 61343); for each of the 21 included pulsars, its last catalogued glitch g. Pulsars whose g lacks dff > 0 or df1 > 0 are not forecastable.
3. Forecast: ln T ~ Normal(mu, sigma), T = wait after g in days, from the frozen fits in results/primary.json: M2 mu = c_p + b ln R_g with R_g = (dff/1000)/df1; M0 (baseline) mu = c_p; each model's own sigma-hat. Conditioned on T > t0 = 61343 - epoch(g). Median and central 50% (0.25-0.75) and 90% (0.05-0.95) windows of the conditional distribution.
4. Scoring, on 2028-11-01 (MJD 62076): download the catalogue again. For each forecast pulsar the outcome is the first glitch after g: the first catalogued glitch of that pulsar whose epoch exceeds g's epoch by more than g's catalogued epoch uncertainty (0 if none is given).
5. Void: if that glitch's epoch is before MJD 61343, the pulsar is void and not scored.
6. Score per pulsar: log predictive density of the observed wait (density of T in days) under the conditional forecast, or, if no glitch by MJD 62076, log of the conditional survival at the wait to that date; M2 minus M0. Total score = sum over scored pulsars; positive means M2 forecasts better.
7. Coverage: for the 50% and 90% windows, the fraction of scored pulsars whose observed wait lies inside the window (ends included). A pulsar with no glitch by the scoring date counts as outside if the window ends before the scoring date, and as not determinable (excluded from that coverage) otherwise; the counts are reported.
8. Functions: glitchlib.forecast (conditional_quantiles, log_predictive, coverage, outcome_for).

## Forecasts (conditional on no glitch before 2026-10-30)

| Pulsar | Last glitch g | Model | Median | 50% window | 90% window |
|---|---|---|---|---|---|
| B0531+21 | 2025-08-06 | M2 | 2027-05-26 | 2027-01-23 to 2028-01-04 | 2026-11-14 to 2029-08-31 |
| B0531+21 | 2025-08-06 | M0 | 2027-07-18 | 2027-02-12 to 2028-05-05 | 2026-11-17 to 2030-08-18 |
| B1046-58 | 2024-04-06 | M2 | 2028-03-31 | 2027-05-31 to 2029-09-30 | 2026-12-07 to 2033-10-14 |
| B1046-58 | 2024-04-06 | M0 | 2028-05-14 | 2027-06-14 to 2030-01-31 | 2026-12-09 to 2034-12-29 |
| B1757-24 | 2019-03-22 | M2 | 2030-01-14 | 2028-02-19 to 2033-07-06 | 2027-01-22 to 2042-12-03 |
| B1757-24 | 2019-03-22 | M0 | 2030-01-06 | 2028-02-08 to 2033-08-15 | 2027-01-19 to 2043-12-01 |
| B1800-21 | 2020-04-19 | M2 | 2030-11-23 | 2028-07-08 to 2035-02-19 | 2027-02-18 to 2046-07-06 |
| B1800-21 | 2020-04-19 | M0 | 2030-06-21 | 2028-04-18 to 2034-07-20 | 2027-02-01 to 2046-03-18 |
| B1823-13 | 2020-01-31 | M2 | 2031-05-17 | 2028-09-22 to 2036-01-26 | 2027-03-04 to 2048-08-17 |
| B1823-13 | 2020-01-31 | M0 | 2032-04-24 | 2029-01-29 to 2038-04-13 | 2027-03-26 to 2055-04-11 |
| J0205+6449 | 2018-07-21 | M2 | 2028-04-22 | 2027-06-05 to 2029-12-04 | 2026-12-07 to 2034-04-10 |
| J0205+6449 | 2018-07-21 | M0 | 2028-07-29 | 2027-07-13 to 2030-07-12 | 2026-12-13 to 2035-12-15 |
| J0537-6910 | 2025-03-28 | M2 | 2027-02-10 | 2026-12-10 to 2027-06-02 | 2026-11-06 to 2028-03-31 |
| J0537-6910 | 2025-03-28 | M0 | 2027-03-10 | 2026-12-21 to 2027-08-05 | 2026-11-08 to 2028-09-21 |
| J0631+1036 | 2018-08-22 | M2 | 2028-04-05 | 2027-05-29 to 2029-10-30 | 2026-12-06 to 2034-01-13 |
| J0631+1036 | 2018-08-22 | M0 | 2028-06-10 | 2027-06-23 to 2030-03-29 | 2026-12-10 to 2035-03-21 |
| J1023-5746 | 2019-12-06 | M2 | 2028-04-23 | 2027-06-05 to 2029-12-08 | 2026-12-07 to 2034-05-02 |
| J1023-5746 | 2019-12-06 | M0 | 2028-10-03 | 2027-08-07 to 2030-12-06 | 2026-12-18 to 2037-01-25 |
| J1420-6048 | 2019-03-13 | M2 | 2028-05-25 | 2027-06-18 to 2030-02-14 | 2026-12-09 to 2034-10-10 |
| J1420-6048 | 2019-03-13 | M0 | 2029-05-31 | 2027-11-11 to 2032-05-05 | 2027-01-03 to 2040-08-30 |
| J1617-5055 | 2012-12-06 | M2 | 2028-10-29 | 2027-08-21 to 2031-01-01 | 2026-12-20 to 2036-10-10 |
| J1617-5055 | 2012-12-06 | M0 | 2028-12-29 | 2027-09-14 to 2031-05-21 | 2026-12-24 to 2037-11-18 |
| J1737-3137 | 2022-09-26 | M2 | 2028-12-08 | 2027-09-11 to 2031-03-08 | 2026-12-25 to 2037-03-29 |
| J1737-3137 | 2022-09-26 | M0 | 2029-04-20 | 2027-10-30 to 2032-01-23 | 2027-01-02 to 2039-12-11 |
| J1814-1744 | 2011-03-06 | M2 | 2030-07-25 | 2028-05-02 to 2034-09-01 | 2027-02-03 to 2045-10-16 |
| J1814-1744 | 2011-03-06 | M0 | 2030-04-05 | 2028-03-15 to 2034-02-06 | 2027-01-25 to 2044-10-21 |
| J1841-0524 | 2020-09-19 | M2 | 2028-11-26 | 2027-09-02 to 2031-03-05 | 2026-12-23 to 2037-05-06 |
| J1841-0524 | 2020-09-19 | M0 | 2028-09-18 | 2027-08-01 to 2030-11-07 | 2026-12-17 to 2036-11-24 |
| J2021+3651 | 2023-12-11 | M2 | 2029-10-04 | 2028-02-08 to 2032-08-06 | 2027-01-27 to 2040-02-02 |
| J2021+3651 | 2023-12-11 | M0 | 2029-08-09 | 2027-12-25 to 2032-07-28 | 2027-01-14 to 2040-12-29 |
| J2229+6114 | 2026-01-31 | M2 | 2028-07-31 | 2027-09-13 to 2030-01-08 | 2027-01-16 to 2033-09-30 |
| J2229+6114 | 2026-01-31 | M0 | 2028-06-04 | 2027-07-23 to 2029-12-09 | 2026-12-25 to 2034-03-09 |

## Not forecastable

- B1338-62: last glitch lacks dff and df1 > 0 (g at MJD 58214.0)
- B1727-47: last glitch lacks df1 > 0 (g at MJD 57984.0)
- B1737-30: last glitch lacks df1 > 0 (g at MJD 60980.0)
- B1758-23: last glitch lacks df1 > 0 (g at MJD 58993.5)
- J0729-1448: last glitch lacks df1 > 0 (g at MJD 58266.6)
