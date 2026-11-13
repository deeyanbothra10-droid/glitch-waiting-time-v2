# Do neutron-star glitches keep time?

A pre-registered, measurement-aware test of the vortex-creep prediction of the time to the next pulsar glitch.
This is the research repository for the paper. It contains the manuscript ([PAPER.pdf](PAPER.pdf), built from
[paper/glitch_paper_revised.tex](paper/glitch_paper_revised.tex); [PAPER_FINAL_REVIEW.pdf](PAPER_FINAL_REVIEW.pdf) is
the approved review build and is identical to it), the analysis code, the catalogue snapshot and the data derived
from it, the frozen pre-registration records, the results, figures and forecasts, and their documentation. Where the
paper and the repository differ, the committed files in this repository are the record.

## Research question

The vortex-creep model predicts the wait T to the next glitch from the fractional steps in spin frequency
(Δν/ν) and spin-down rate (Δν̇/ν̇) of the last one (Zhou et al. 2022, equation 48). Writing
R = (Δν/ν)/(Δν̇/ν̇), the prediction is ln T = ln R + c_p within each pulsar p. The test asks whether, within pulsars,
waiting times scale with R (T1), and whether the spin-down step carries the predicted negative coefficient at fixed
glitch size (T2), the part of the prediction that a dependence on glitch size alone cannot produce.

## Data sources

- Jodrell Bank Glitch Catalogue (https://www.jb.man.ac.uk/pulsar/glitches/gTable.html), downloaded on
  2026-10-30: 736 glitches in 226 objects, stored unchanged as `data/raw/glitches.csv`.
- ATNF Pulsar Catalogue, version 2.9.0: spin frequencies and their derivatives in `data/meta/pulsar_params.csv`.
- A source audit of how each Δν̇ was measured: `data/meta/source_audit.csv` (its provenance is described in
  `data/meta/README.md`).
- Documented ends of monitoring and monitoring gaps: `data/meta/obs_end.csv`, `data/meta/monitoring_gaps.csv`.

## Sample

225 admissible inter-glitch intervals (219 completed, 6 right-censored) in 21 pulsars, each with at least four
admissible completed intervals. By Δν̇ measurement method: 70 local relaxation, 55 global fits with exponentials and
100 short window. Counts per pulsar and method are in `data/processed/final_dataset.csv`; the exclusion steps are in
`logs/exclusions.csv` and `logs/dropped_last_intervals.csv`.

## Models

Log-normal accelerated-failure-time models with right censoring, one intercept c_p per pulsar and one scatter σ
shared by all pulsars:

| Model | Linear predictor |
|---|---|
| M0 | c_p |
| M1 | ln R + c_p |
| M2 | b ln R + c_p |
| M3 | a ln(Δν/ν) + c_p |
| M4 | a ln(Δν/ν) + d ln(Δν̇/ν̇) + c_p |

The prediction is b > 0 in M2 and d < 0 in M4 (b = 1 and d = −1 without measurement error).

## Simulated populations

Every decision refers a real statistic to its distribution over simulated catalogues that share the real pulsars,
glitch counts, observing windows, method mix and uncertainties (1,000 replicates per population, 500 for the
bias-only S4 populations; seed 20261031):

- **S1**: the prediction is true, with catalogue-level measurement noise only.
- **S1m**: as S1, but Δν̇ is measured through an emulation of each fitting method, with the recovery time t_rec = T.
- **S2**: waits independent of glitch properties.
- **S3**: waits depend on glitch size alone, ln T = a ln(Δν/ν) + c_p + σε.
- **S3m**: as S3, but measured through the emulation.
- **S4 variants**: waits as in S2, with Δν̇ measured through one method and t_rec independent of T, so that any
  association is a measurement artefact: S4w (window), S4l (local relaxation), S4g (global), S4g0 (global without
  exponentials), S4sw-espinoza and S4sw-wang (short window, w = 40 d and 125 d).
- **S4-mix**: as the S4 variants, with the real method mix.

## Tests

- **T1 (scaling)**: the signed root of the likelihood-ratio statistic for M2 against M0 must exceed its 97.5th
  percentile under S2.
- **T2 (the Δν̇ term)**: d̂ in M4 must fall below the smaller of its 2.5th percentiles under S3 and S3m, whose size
  coefficient is set to the real M3 estimate.

In the real run, each reference population is simulated at the generating σ that reproduces the real fitted σ̂.
The outcome is the first of D (artefact signature), C (refutation), A (scaling supported; within A: A2, B or A1),
D (signal within bias range) and inconclusive whose rule applies (DECISIONS.md, D4j as amended by D4l–D4p).

## Result

The pre-registered outcome is **A2**: scaling supported, with the Δν̇ term established.

- b̂ = 0.268, calibrated 95 per cent interval [0.200, 0.336]; T1 signed root 7.398 against 2.037 (S2 97.5th
  percentile).
- d̂ = −0.257, calibrated 95 per cent interval [−0.384, −0.130]; below −0.130 (S3m 2.5th percentile; S3: −0.094).
- σ̂ = 0.696 (M2). Power at the matched σ: T1 0.997, T2 0.769.

The full record is in `results/primary.json` (07), `results/secondary.json` (08) and `results/robustness.json` (09).

## How to reproduce

Requirements: Python 3.13 with numpy, pandas, scipy, matplotlib, pyyaml and pytest; lifelines is needed only for
one cross-check test. The versions used were Python 3.13.1, numpy 2.2.1, pandas 2.2.3, scipy 1.15.1,
matplotlib 3.10.0, pyyaml 6.0.3, pytest 9.2.0 and lifelines 0.30.3; `run_all.py` writes the versions it finds to
`results/environment.txt`.

```bash
python3 -m pytest
```

runs every test (synthetic data, except structure checks of the inputs).

```bash
python3 run_all.py --reproduce /path/to/empty/folder
```

clones this repository's HEAD into the empty folder and reruns every step of `run_all.py` there (`00_provenance`,
`01_pulsar_params`, `emulation_fixture`, `02_clean`, `method_counts`, `collinearity`, `05_power`, `06_validate`,
`03_intervals`, `04_uncertainty`, `07_primary`, `08_secondary`, `09_robustness`, `10_figures`). It then compares every file it writes with the committed version (JSON exactly, apart from commit, timestamp and timing
fields; CSV and text byte for byte; PDF without creation dates) and writes `reproduction_report.json` in that folder.
Everything is written under that folder. In this mode only, the 07–09 preflight skips its run-once and clean-tree
checks; the no-peeking guard and its checksum check stay in force. A full run takes several hours on 8 cores.
`logs/reproduction/reproduction_report_5605ad4.json` is the report of a full run made during the project's earlier
development (identifier 5605ad4, a commit of that history; see "Git history of this repository"). No full run has
been made from this repository's release commit.

`python3 run_all.py` without `--reproduce` runs the same steps in place and rewrites committed outputs; in this
repository 07–09 then stop at their preflight, because their outputs already exist (each runs once).

The post-freeze scripts are not part of `run_all.py`: `scripts/exploratory_posthoc.py`,
`scripts/exploratory_influence.py`, `scripts/exploratory_single_pulsar.py`, `scripts/exploratory_heterogeneity_a.py`,
`scripts/sigma_floor.py`, `scripts/forecast_v1.py`, `scripts/forecast_v2.py` and `scripts/forecast_v3.py` run from
committed files and rewrite their committed outputs. `scripts/detection_fractions.py` and
`scripts/short_window_label.py` read the saved partial results of 06 runs (`results/validation_partial*`), which are
not in the repository.

`results/short_window_label.json`, written by `scripts/short_window_label.py`, is a report-only record made on
2026-11-01 from 06's saved replicates under the settings in force before the freeze. It is kept in the repository
with its numbers as generated; only its store paths were rewritten as repository-relative paths. The script cannot
currently regenerate it: run under the frozen configuration, it stops at its first check because the regenerated A2
replicates do not reproduce the saved b-hat, so this result has not been independently reproduced. The script as it
was before the change to repository-relative paths stops in the same way, so that change did not cause the failure.

```bash
python3 scripts/posthoc_crab_persistent_step.py
```

is an exploratory, post hoc calculation, not part of the pre-registered primary analysis. It recomputes the Crab
pulsar's implied β with a persistent-step Δν̇ in place of the catalogue values (Section 9.2 of the paper), from
committed files (`data/processed/final_dataset.csv`, `data/meta/pulsar_params.csv` and the frozen M1 intercept in
`results/secondary.json`) and one literature scaling (Lyne et al. 2015). It rewrites
`results/exploratory/crab_persistent_step.json`; `tests/test_posthoc_crab.py` checks that it reproduces that file
exactly.

Every decision is read from `config/decisions.yaml`, transcribed from `DECISIONS.md` (see `config/TRANSCRIPTION.md`);
seeds come from the config and are recorded in every output.

## Building the manuscript

`paper/glitch_paper_revised.tex` builds `PAPER.pdf` with pdfTeX (TeX Live 2026, pdfTeX 1.40.29). With `mnras.cls` and
`newtx` not installed, its preamble falls back to the plain article class with `mathptmx`, which is the layout of
`PAPER.pdf` (14 pages); with `mnras.cls` installed the same source builds the MNRAS layout instead. From the
repository root:

```bash
mkdir -p build
SOURCE_DATE_EPOCH=1794553842 FORCE_SOURCE_DATE=1 pdflatex -interaction=nonstopmode -output-directory=build paper/glitch_paper_revised.tex
```

Run the second command three times, then copy `build/glitch_paper_revised.pdf` to `PAPER.pdf` and
`PAPER_FINAL_REVIEW.pdf`. The fixed date (2026-11-13 07:10:42 UTC) and the fixed output path make the build
byte-identical. `build/` is not part of the repository.

## Where each part of the paper lives

| Paper | Repository |
|---|---|
| Sec. 2, the prediction | `DECISIONS.md` (D1–D3); `glitchlib/models.py` |
| Sec. 3.1, catalogue and intervals | `data/raw/glitches.csv`; `data/meta/`; `scripts/00_provenance.py`–`04_uncertainty.py`; `glitchlib/cleaning.py`, `glitchlib/intervals.py`; `data/processed/`; `logs/exclusions.csv`, `logs/dropped_last_intervals.csv` |
| Sec. 3.2 and Appendix B, source audit | `data/meta/source_audit.csv`; `glitchlib/cleaning.py`; `results/short_window_admission.json` |
| Sec. 3.3 and Table 1, sample | `data/processed/final_dataset.csv`, `data/processed/pulsar_params.csv` |
| Sec. 4, emulation of Δν̇ measurement | `glitchlib/measure.py`, `glitchlib/simulate.py`, `glitchlib/templates.py`; `scripts/emulation_fixture.py`; `results/emulation_s4_fixture*.json`, `results/short_window_admission*.json`, `results/detection_fractions.json`; `logs/emulation_findings.md` |
| Sec. 5.1–5.3, models, populations, tests, outcomes | `glitchlib/models.py`, `glitchlib/simulate.py`, `glitchlib/validation.py`, `glitchlib/decision.py`, `glitchlib/analysis.py`, `glitchlib/diagnostics.py` |
| Sec. 5.4 and Appendix A, pre-registration and amendments | `DECISIONS.md`; `config/`; `prereg/FROZEN`, `prereg/data_freeze.json`; `glitchlib/guards.py`, `glitchlib/preflight.py` |
| Sec. 6, Table 2 and Fig. 1, validation | `scripts/05_power.py`, `scripts/06_validate.py`; `results/validation.json`, `results/power.json`, `results/sigma_floor.json`; `figures/F5_validation.pdf` |
| Sec. 7 and Table 3, pre-registered result | `scripts/07_primary.py`, `scripts/08_secondary.py`; `results/primary.json`, `results/secondary.json` |
| Table 4 and Fig. 4, implied β | `results/secondary.json`; `figures/F6_implied_beta.pdf` |
| Fig. 2, Fig. 3, Fig. 5 | `figures/F2_lnT_vs_lnR.pdf`, `figures/F4_model_comparison.pdf`, `figures/F3_kaplan_meier.pdf` (all from `scripts/10_figures.py`) |
| Sec. 7.1 and Table 5, robustness | `scripts/09_robustness.py`; `results/robustness.json` |
| Sec. 8 and Table 6, post hoc analyses | `glitchlib/exploratory.py`; `scripts/exploratory_*.py`; `results/exploratory/` |
| Sec. 9.5 and Appendix C, forward forecast | `glitchlib/forecast.py`; `scripts/forecast_v1.py`–`forecast_v3.py`; `results/forecast/` (versions v1–v3; v3 is the final version) |
| Reproduction | `run_all.py`; `glitchlib/reproduce.py`; `logs/reproduction/` |

`figures/F1_exclusion_flow.pdf` is produced by `scripts/10_figures.py` but not used in the paper. Tests are in
`tests/`.

## Limitations

As listed in Section 9.4 of the paper:

(i) Shared scatter. The primary model uses one σ for all pulsars, which the single-pulsar fits contradict. The
simulations share the assumption, so calibration does not correct for it. A hierarchical model with pulsar-specific
scatter and partially pooled slopes would be the natural next step.

(ii) Emulation. See Section 4.3: no phase-connected fitting, no red noise, two uncited settings, and four to six
pulsars, depending on σ, whose emulated uncertainties exceed the catalogue's.

(iii) Sample. PSR J0537−6910 contributes 28 per cent of the intervals, though the result survives its removal. The
method behind the catalogue maintainers' own measurements is assumed (robustness case included).

(iv) Constant β. Equation (3) assumes β_p is constant within a pulsar.

(v) Self-held pre-registration. The freeze is attested by the author's version-control history and an e-mailed record
of the checksums, not by a public registry. The rules were amended sixteen times before the freeze, partly in response to simulated power.

## Licence, citation and AI assistance

See [LICENSE](LICENSE) (MIT for code, CC BY 4.0 for results, copyright reserved for the manuscript; no licence is
granted for the catalogue snapshot or the input tables in `data/raw/` and `data/meta/`), [CITATION.cff](CITATION.cff)
and [AI_ASSISTANCE.md](AI_ASSISTANCE.md).

## How to verify the pre-registration

The analysis decisions were frozen before any statistic relating glitch size to waiting time was computed on the
real data. The pre-registration and the research timeline are documented by the dated, frozen records in this
repository; every check below uses only this repository.

### Frozen records and checksums (SHA-256)

| File | SHA-256 |
|---|---|
| `data/raw/glitches.csv` (catalogue downloaded 2026-10-30) | `5ba59d0160fe1f00eb70d464aa048c8b9968980c251118dca6deed825f849de7` |
| `DECISIONS.md` (the frozen rules) | `2fd950febe3648b20098cca04a4ebe2742c905c43218d7977fa4894198f00bed` |
| `config/decisions.yaml` (their frozen transcription) | `3cd2b2f9254697aabaf60a6264d6c7856f8918d3577b6877dfb108fdfb1c671f` |
| `data/processed/final_dataset.csv` (data freeze) | `27a63d37ab877d79e28340bf113e768fe612e290bba1a8f4f2a1ae7ff8bdbe3b` |
| `results/forecast/forecast_v1.json` | `1aebf098db140138b701735fea8e5a6bd33333f198924d3eb1ab762aa91a52c6` |
| `results/forecast/forecast_v2.json` | `44b7c64fa8c0bc4d31b1caa3d4aa959913568971bb2d165e2ed19b922ada3b8c` |
| `results/forecast/forecast_v3.json` | `9f46b89465ce63096ef68822f618c8c6c08b87ed40df8c81f05fb3878a3f43f0` |

```bash
shasum -a 256 data/raw/glitches.csv DECISIONS.md config/decisions.yaml data/processed/final_dataset.csv \
  results/forecast/forecast_v1.json results/forecast/forecast_v2.json results/forecast/forecast_v3.json
cat prereg/FROZEN prereg/data_freeze.json
```

`prereg/FROZEN` records the freeze date (2026-11-01) and the checksums of `DECISIONS.md` and `config/decisions.yaml`;
`prereg/data_freeze.json` records the checksum of `data/processed/final_dataset.csv`. The digests of the rules, the
configuration, the analysis data set and the catalogue file (Section 5.4 of the paper) were also sent to an external
e-mail account on 2026-11-01; that e-mail is held by the author.

### Documented timeline

| Date | Event | Record in this repository |
|---|---|---|
| 2026-10-30 | catalogue downloaded | `data/raw/README.md` |
| 2026-10-31 | rules D1–D15 written; amendments D4a–D4k (D4a's heading is undated) | `DECISIONS.md` |
| 2026-11-01 | amendments D4l–D4p | `DECISIONS.md` |
| 2026-11-01 07:21 UTC | validation on simulated data (105 of 105 criteria pass) | `results/validation.json` |
| 2026-11-01 | pre-registration frozen; data freeze | `prereg/FROZEN`, `prereg/data_freeze.json` |
| 2026-11-01 07:58 UTC | 07 primary analysis, run once on the real data | `results/primary.json` |
| 2026-11-01 08:09 UTC | 08 secondary analyses and 09 robustness cases | `results/secondary.json`, `results/robustness.json` |
| 2026-11-01 to 2026-11-04 | post hoc exploratory analyses | `results/exploratory/` |
| 2026-11-12 15:08–15:16 UTC | forward forecasts v1–v3 (v3 final; scored on 2028-11-01) | `results/forecast/` |

The times are the `generated_utc` fields written by the analysis code. They come from the author's computer clock,
so they are self-recorded, not independent timestamps. The order of the freeze and the real-data analysis is also
enforced in code: `glitchlib/guards.py` refuses any statistic relating glitch properties to real waiting times
unless `prereg/FROZEN` exists and the analysis data set matches the checksum in `prereg/data_freeze.json`, and the
primary script refuses to run once its output exists. Everything after the 09 robustness run (figures, the
reproduction tooling, `results/exploratory/`, `results/forecast/`) is post-freeze and labelled as such; none of it
changes a pre-registered result.

### Git history of this repository

This repository's Git history begins with a single consolidated release commit made on 2026-11-13. It does not
establish the chronology of the research, which is documented by the dated and frozen records above. The project's
earlier development history, in which the freeze was made, is not included. Identifiers of that history are kept
unchanged where the records contain them (the `git_commit` fields of the results, the commit named in
`prereg/FROZEN` and `prereg/data_freeze.json`, and the identifier 5605ad4 of the reproduction report); they cannot be
looked up in this repository.
