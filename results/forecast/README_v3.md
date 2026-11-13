# Forward forecast v3 (post-freeze; NOT part of the pre-registration)

Written by `scripts/forecast_v3.py`; the forecasts are in `forecast_v3.json`. Forecasts v1 and v2 (`forecast_v1.json`, `README.md`, `forecast_v2.json`, `README_v2.md`; tags `forecast-v1`, `forecast-v2`) are kept unchanged.

## Why v2 was superseded

v2 scored every Jodrell Bank pulsar conditioned on no glitch before the 2026-10-30 snapshot. For some of them that assumption is not supported: J1814-1744's last glitch is 5,717 d before the snapshot, 4.9x its longest catalogued wait (1,166 d); and J0205+6449, J0631+1036 and J1814-1744 have a documented obs_end in Dec 2018 (end of the Basu et al. 2022 Jodrell Bank sample). v3 adds rules (b) and (c) below. It was made before any outcome was known.

## Scored set

6 pulsars are scored: B0531+21, B1800-21, B1823-13, J1737-3137, J2021+3651, J2229+6114. With so few pulsars and windows this wide, the score is a weak test; it is reported whatever it shows.

## Rule (fixed before any outcome)

1. Label: post-freeze forward forecast v3; NOT part of the pre-registration (prereg-v2). Supersedes forecast v2 (tag forecast-v2) and v1 (tag forecast-v1), both kept unchanged, before any outcome was known.
2. Snapshot: the Jodrell Bank glitch catalogue downloaded 2026-10-30 (MJD 61343); for each of the 21 included pulsars, its last catalogued glitch g. Pulsars whose g lacks dff > 0 or df1 > 0 are not forecastable.
3. Forecast: ln T ~ Normal(mu, sigma), T = wait after g in days, from the frozen fits in results/primary.json: M2 mu = c_p + b ln R_g with R_g = (dff/1000)/df1; M0 (baseline) mu = c_p; each model's own sigma-hat. Median and central 50% (0.25-0.75) and 90% (0.05-0.95) windows of the conditional distribution.
4. Scored set: a forecastable pulsar is scored only if (a) it is covered by the D4g Jodrell Bank rule; (b) it has no documented obs_end before MJD 61343 in data/meta/obs_end.csv (an obs_end equal to the last catalogued glitch epoch is the D4g placeholder, not a documented end of monitoring); and (c) 61343 - epoch(g) is no longer than its longest catalogued wait between consecutive glitches. Scored pulsars are conditioned on T > 61343 - epoch(g).
5. Reported, not scored: every other forecastable pulsar, with the rules it fails and their numbers; conditioned on T > documented obs_end - epoch(g), or unconditioned if obs_end equals epoch(g).
6. Scoring, on 2028-11-01 (MJD 62076): download the catalogue again. For each scored pulsar the outcome is the first glitch after g: the first catalogued glitch of that pulsar whose epoch exceeds g's epoch by more than g's catalogued epoch uncertainty (0 if none is given).
7. Void: if that glitch's epoch is before MJD 61343, the pulsar is void and not scored.
8. Score per scored pulsar: log predictive density of the observed wait (density of T in days) under the conditional forecast, or, if no glitch by MJD 62076, log of the conditional survival at the wait to that date; M2 minus M0. Total score = sum over scored pulsars; positive means M2 forecasts better.
9. Coverage: for the 50% and 90% windows, the fraction of scored pulsars whose observed wait lies inside the window (ends included). A pulsar with no glitch by the scoring date counts as outside if the window ends before the scoring date, and as not determinable (excluded from that coverage) otherwise; the counts are reported.
10. Functions: glitchlib.forecast (conditional_quantiles, log_predictive, coverage, outcome_for, forecast_window, scoring_eligibility).

## Scored forecasts

| Pulsar | Last glitch g | Conditioning | Model | Median | 50% window | 90% window |
|---|---|---|---|---|---|---|
| B0531+21 | 2025-08-06 | no glitch before MJD 61343 (snapshot) | M2 | 2027-05-26 | 2027-01-23 to 2028-01-04 | 2026-11-14 to 2029-08-31 |
| B0531+21 | 2025-08-06 | no glitch before MJD 61343 (snapshot) | M0 | 2027-07-18 | 2027-02-12 to 2028-05-05 | 2026-11-17 to 2030-08-18 |
| B1800-21 | 2020-04-19 | no glitch before MJD 61343 (snapshot) | M2 | 2030-11-23 | 2028-07-08 to 2035-02-19 | 2027-02-18 to 2046-07-06 |
| B1800-21 | 2020-04-19 | no glitch before MJD 61343 (snapshot) | M0 | 2030-06-21 | 2028-04-18 to 2034-07-20 | 2027-02-01 to 2046-03-18 |
| B1823-13 | 2020-01-31 | no glitch before MJD 61343 (snapshot) | M2 | 2031-05-17 | 2028-09-22 to 2036-01-26 | 2027-03-04 to 2048-08-17 |
| B1823-13 | 2020-01-31 | no glitch before MJD 61343 (snapshot) | M0 | 2032-04-24 | 2029-01-29 to 2038-04-13 | 2027-03-26 to 2055-04-11 |
| J1737-3137 | 2022-09-26 | no glitch before MJD 61343 (snapshot) | M2 | 2028-12-08 | 2027-09-11 to 2031-03-08 | 2026-12-25 to 2037-03-29 |
| J1737-3137 | 2022-09-26 | no glitch before MJD 61343 (snapshot) | M0 | 2029-04-20 | 2027-10-30 to 2032-01-23 | 2027-01-02 to 2039-12-11 |
| J2021+3651 | 2023-12-11 | no glitch before MJD 61343 (snapshot) | M2 | 2029-10-04 | 2028-02-08 to 2032-08-06 | 2027-01-27 to 2040-02-02 |
| J2021+3651 | 2023-12-11 | no glitch before MJD 61343 (snapshot) | M0 | 2029-08-09 | 2027-12-25 to 2032-07-28 | 2027-01-14 to 2040-12-29 |
| J2229+6114 | 2026-01-31 | no glitch before MJD 61343 (snapshot) | M2 | 2028-07-31 | 2027-09-13 to 2030-01-08 | 2027-01-16 to 2033-09-30 |
| J2229+6114 | 2026-01-31 | no glitch before MJD 61343 (snapshot) | M0 | 2028-06-04 | 2027-07-23 to 2029-12-09 | 2026-12-25 to 2034-03-09 |

## Reported, not scored

| Pulsar | Last glitch g | Conditioning | Model | Median | 50% window | 90% window | Fails |
|---|---|---|---|---|---|---|---|
| B1046-58 | 2024-04-06 | unconditioned (obs_end equals g) | M2 | 2026-08-18 | 2025-09-28 to 2028-01-17 | 2025-01-06 to 2031-09-11 | (a) not covered by the D4g Jodrell Bank rule |
| B1046-58 | 2024-04-06 | unconditioned (obs_end equals g) | M0 | 2026-06-01 | 2025-07-12 to 2027-12-03 | 2024-11-07 to 2032-02-12 | |
| B1757-24 | 2019-03-22 | unconditioned (obs_end equals g) | M2 | 2024-05-30 | 2022-06-20 to 2027-07-11 | 2020-11-15 to 2035-07-15 | (c) snapshot - g = 2778.2 d exceeds its longest catalogued wait 2318.0 d (1.20x) |
| B1757-24 | 2019-03-22 | unconditioned (obs_end equals g) | M0 | 2023-03-21 | 2021-07-27 to 2026-01-06 | 2020-04-25 to 2033-10-21 | |
| J0205+6449 | 2018-07-21 | no glitch before documented obs_end MJD 58483.0 | M2 | 2020-01-22 | 2019-07-11 to 2020-12-05 | 2019-02-16 to 2023-03-07 | (b) documented obs_end MJD 58483.0 is before MJD 61343; (c) snapshot - g = 3023.0 d exceeds its longest catalogued wait 1987.0 d (1.52x) |
| J0205+6449 | 2018-07-21 | no glitch before documented obs_end MJD 58483.0 | M0 | 2019-12-02 | 2019-06-02 to 2020-10-19 | 2019-02-01 to 2023-04-09 | |
| J0537-6910 | 2025-03-28 | no glitch before documented obs_end MJD 60827.0 | M2 | 2025-08-02 | 2025-06-27 to 2025-10-01 | 2025-06-06 to 2026-03-11 | (a) not covered by the D4g Jodrell Bank rule; (b) documented obs_end MJD 60827.0 is before MJD 61343 |
| J0537-6910 | 2025-03-28 | no glitch before documented obs_end MJD 60827.0 | M0 | 2025-08-11 | 2025-07-01 to 2025-10-25 | 2025-06-06 to 2026-05-25 | |
| J0631+1036 | 2018-08-22 | no glitch before documented obs_end MJD 58483.0 | M2 | 2020-01-17 | 2019-07-15 to 2020-11-12 | 2019-02-22 to 2022-12-27 | (b) documented obs_end MJD 58483.0 is before MJD 61343; (c) snapshot - g = 2990.9 d exceeds its longest catalogued wait 2639.7 d (1.13x) |
| J0631+1036 | 2018-08-22 | no glitch before documented obs_end MJD 58483.0 | M0 | 2019-10-12 | 2019-05-10 to 2020-07-07 | 2019-01-27 to 2022-08-03 | |
| J1023-5746 | 2019-12-06 | no glitch before documented obs_end MJD 59269.0 | M2 | 2022-03-19 | 2021-08-09 to 2023-04-10 | 2021-03-27 to 2026-01-30 | (a) not covered by the D4g Jodrell Bank rule; (b) documented obs_end MJD 59269.0 is before MJD 61343; (c) snapshot - g = 2520.0 d exceeds its longest catalogued wait 961.0 d (2.62x) |
| J1023-5746 | 2019-12-06 | no glitch before documented obs_end MJD 59269.0 | M0 | 2022-06-10 | 2021-09-09 to 2023-10-21 | 2021-04-01 to 2027-09-02 | |
| J1420-6048 | 2019-03-13 | unconditioned (obs_end equals g) | M2 | 2020-12-20 | 2020-04-21 to 2022-01-12 | 2019-10-05 to 2024-10-09 | (a) not covered by the D4g Jodrell Bank rule; (c) snapshot - g = 2788.0 d exceeds its longest catalogued wait 1345.0 d (2.07x) |
| J1420-6048 | 2019-03-13 | unconditioned (obs_end equals g) | M0 | 2022-02-01 | 2020-11-23 to 2024-02-11 | 2019-12-27 to 2029-10-03 | |
| J1617-5055 | 2012-12-06 | unconditioned (obs_end equals g) | M2 | 2014-05-26 | 2013-11-06 to 2015-04-12 | 2013-05-25 to 2017-07-17 | (a) not covered by the D4g Jodrell Bank rule |
| J1617-5055 | 2012-12-06 | unconditioned (obs_end equals g) | M0 | 2013-11-14 | 2013-06-25 to 2014-07-12 | 2013-03-10 to 2016-05-12 | |
| J1814-1744 | 2011-03-06 | no glitch before documented obs_end MJD 58483.0 | M2 | 2021-12-31 | 2020-03-20 to 2025-04-07 | 2019-03-19 to 2034-02-17 | (b) documented obs_end MJD 58483.0 is before MJD 61343; (c) snapshot - g = 5717.0 d exceeds its longest catalogued wait 1166.0 d (4.90x) |
| J1814-1744 | 2011-03-06 | no glitch before documented obs_end MJD 58483.0 | M0 | 2021-06-06 | 2019-12-21 to 2024-03-08 | 2019-03-03 to 2031-12-26 | |
| J1841-0524 | 2020-09-19 | unconditioned (obs_end equals g) | M2 | 2023-11-22 | 2022-09-14 to 2025-10-16 | 2021-09-23 to 2030-09-08 | (c) snapshot - g = 2231.4 d exceeds its longest catalogued wait 1637.1 d (1.36x) |
| J1841-0524 | 2020-09-19 | unconditioned (obs_end equals g) | M0 | 2022-09-13 | 2021-11-19 to 2024-02-01 | 2021-04-05 to 2027-12-13 | |

Note, J0537-6910: its glitch rate (one glitch every 100-200 d) makes a 580 d glitch-free stretch to the 2026-10-30 snapshot implausible, so the snapshot is not used for conditioning.

## Not forecastable

- B1338-62: last glitch lacks dff and df1 > 0 (g at MJD 58214.0)
- B1727-47: last glitch lacks df1 > 0 (g at MJD 57984.0)
- B1737-30: last glitch lacks df1 > 0 (g at MJD 60980.0)
- B1758-23: last glitch lacks df1 > 0 (g at MJD 58993.5)
- J0729-1448: last glitch lacks df1 > 0 (g at MJD 58266.6)
