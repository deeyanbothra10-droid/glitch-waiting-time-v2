# Measurement-emulation findings

Findings from the Δν̇ measurement emulation (glitchlib/measure.py) on synthetic data only.
Each entry: date, what was run, numbers, explanation, and what was decided.

---
## 2026-10-31 · S4 local-relaxation shows a small negative forward bias

**Run.** tests/test_simulate.py fixture, before D4c: 10 pulsars, W_p = 1000 d, span 25,000 d,
cadence 10 d, ν̇ noise 1e-3, 1,000 d pre-glitch monitoring, ≥ 5 samples per window, σ = 0.4,
t_rec ~ lognormal(1000 d, 0.3) independent of T, Q ~ log-uniform(0.005, 0.05),
τ_d ~ log-uniform(1, 100) d, 12 replicates (seed 20261031). Exponentials fitted at every τ_d.

**Numbers** (mean slope across replicates ± its standard error):

| Method | Forward b̂ (ln T on ln R) | Backward b̂ (placebo) | Expected | Agrees |
|---|---|---|---|---|
| window | +0.051 ± 0.006 | −0.026 ± 0.006 | positive forward | yes |
| global-noexp | +0.004 ± 0.007 | −0.119 ± 0.006 | negative backward | yes |
| global | −0.001 ± 0.005 | −0.004 ± 0.004 | none | yes |
| local-relaxation | −0.036 ± 0.006 | −0.006 ± 0.005 | none | **no** |

**Explanation.** When t_rec < T, the step finishes relaxing partway through the window, so ν̇ falls
and then stays flat. That shape resembles an exponential decay, and the local fit's exponential term
absorbs the bend, so the reported total Δν̇ is too large. The effect grows with T/t_rec: with noise
switched off, ln(measured/true) rises with ln T_k (coefficient +0.30). Long waits therefore get
too-small R, which gives a negative forward slope. With t_rec = T (S1) there is no bend and the
method is exact. A second effect: exponentials faster than the cadence cannot be resolved, and
fitting their amplitudes gives wild values (34% of glitches with τ_d < 3 d at a 10-d cadence were
off by more than 2 J); D4c now leaves τ_d < 2 × cadence out of every fit.

**Decision (D4c).** Code not changed to remove the bias. Accepted as a finding: it works against
T1 (b > 0). Revised S4l forward expectation: not positive, upper end of the 95% interval ≤ +0.02.

---
## 2026-10-31 · Under the D4c settings, local-relaxation and global show backward biases

**Run.** Same fixture and seed as above, with the D4c emulation settings: cadence 14 d, 365 d of
pre-glitch monitoring, ≥ 4 samples after a glitch (window: also ≥ 4 before), exponentials with
τ_d < 28 d left out of every fit but kept in the data, c_p = ln W_p − σ²/2. ν̇ noise fixed at 1e-3.

| Method | Forward b̂ | Backward b̂ | Expectation | Agrees |
|---|---|---|---|---|
| window | +0.051 ± 0.006 | −0.017 ± 0.007 | positive forward | yes |
| global-noexp | −0.006 ± 0.008 | −0.102 ± 0.007 | negative backward | yes |
| global | −0.001 ± 0.006 | **−0.075 ± 0.008** | none | **no** (backward) |
| local-relaxation | **+0.019 ± 0.004** | **−0.046 ± 0.008** | forward 95% upper ≤ +0.02; backward none | **no** (both) |

**Explanation.** With τ_d log-uniform on 1–100 d, about 72% of exponentials (τ_d < 28 d) are now
unresolved: left out of the fit but present in the data. Noise-free, for glitches whose exponential is
not fitted, ln(measured/true) rises with the previous wait ln T_(k−1) (global +0.74, local +0.34)
and, for local, falls with the following wait ln T_k (−0.33).
- Backward: the previous glitch's unfitted exponential sits at the start of the previous window and
  distorts that window's fitted level, which is the baseline for the next step. How much depends on
  the previous window's length, so R_k depends on T_(k−1): a negative backward slope.
- Local forward: a glitch's own unfitted exponential raises the fitted line at the start of its
  window, more so for a short window, so Δν̇ is overstated for short waits: a small positive forward
  slope (+0.019, 95% upper end +0.028).
The earlier negative forward bias of local-relaxation came from resolved exponentials absorbing the
end of the linear recovery; with most exponentials now unresolved, that effect is outweighed.
These biases depend on the provisional τ_d and Q ranges (D14) and on the cadence.

**Decision.** Pending; code not tuned.

**Decision (D4d, 31 Oct).** S4 is a measurement, not a test: the pass/fail S4 tests were removed and
replaced by noise-free correctness tests of the emulator. The D4c resolvability rule was replaced by
AIC-based exponential detection.

---
## 2026-10-31 · S4 magnitudes on the synthetic fixture under D4d (results/emulation_s4_fixture.json)

**Run.** `python3 scripts/emulation_fixture.py`. Synthetic fixture (not real pulsars): 10 pulsars,
W_p = 1000 d, 25,000 d span, 15-glitch pools, method classes drawn at random. All emulation settings
from config: cadence 14 d (D4c fallback), 365 d pre-glitch, ≥ 4 samples after (window: and before),
exponentials fitted only if they lower the window's AIC by ≥ 10 (D4d), Q and τ_d from D14
(provisional), Δν/ν and Δν̇ noise from the real relative errors. ν̇ noise calibrated per pulsar on
the S4 mix (50 sequences per bisection step, tolerance 5%). 500 replicates per variant per σ,
seed 20261031. Mean slope across replicates ± its standard error:

| Variant | σ = 0.2 fwd / bwd | σ = 0.4 fwd / bwd | σ = 0.6 fwd / bwd |
|---|---|---|---|
| S4w (window) | +0.017 / −0.004 | +0.055 / −0.011 | +0.100 / −0.018 |
| S4l (local-relaxation) | +0.005 / −0.009 | +0.006 / −0.026 | −0.002 / −0.037 |
| S4g (global) | −0.005 / −0.020 | −0.005 / −0.051 | −0.005 / −0.061 |
| S4g0 (global-noexp) | −0.004 / −0.038 | −0.007 / −0.112 | +0.024 / −0.179 |
| S4-mix, all data (unknown as window) | −0.002 / −0.006 | −0.007 / −0.021 | −0.007 / −0.029 |
| S4-mix, all data (unknown as global) | −0.003 / −0.014 | −0.003 / −0.041 | −0.000 / −0.054 |
| S4-mix, admissible only (unknown as window) | +0.000 / −0.017 | −0.001 / −0.043 | +0.000 / −0.048 |

Standard errors of these means are 0.001–0.003. Central 95% of the admissible-only S4-mix forward
b̂: [−0.041, +0.048] at σ = 0.2, [−0.082, +0.078] at σ = 0.4, [−0.126, +0.127] at σ = 0.6.

**Observations.**
- Window gives a positive forward bias that grows with σ. Every emulated method gives a negative
  backward bias that grows with σ; global-noexp is largest.
- The admissible-only mix (the primary analysis's data) shows no forward bias but a backward bias
  of about −0.02 to −0.05. Under D4d the placebo is judged against the S1m placebo distribution,
  which is simulated with the same measurement mix, so this bias is part of that reference.
- Failed fits: 1–4 of 500 replicates in the mixtures (a pulsar with no completed interval in that
  replicate); 0 in the single-method variants. They are counted in the JSON and left out of the
  summaries.
- Noise calibration hit its 5% tolerance for 10/10 pulsars at σ = 0.2, 8/10 at 0.4, 5/10 at 0.6.
  The others are "at floor": with zero noise the fitted model's misfit alone already gives a median
  relative error above the target, because in S4 the recovery (t_rec) does not end at the next
  glitch, and this mismatch grows with σ. Those pulsars get zero added noise and a larger relative
  error than the real data.

---
## 2026-10-31 · Marginal exponential detection adds scatter beyond the reported error bars

**Run.** tests/test_measure.py: four glitches, permanent steps, exponentials A = 1–5 × 10⁻³ with
τ_d = 20 d, ν̇ noise 10⁻³, cadence 5 d, 300 noise realisations, local-relaxation.

**Numbers.** Detection rates per glitch 0.81, 0.03, 0.32, 1.00; ratio of actual scatter to mean
reported error 2.24, 2.12, 3.66, 1.09. With exponentials 10 × larger (always detected) the ratios are
1.02, 0.96, 1.04, 1.09.

**Explanation.** When an exponential is near the detection threshold, whether it is fitted depends
on the noise, and the measured Δν̇ jumps between two models (step only, or step + exponential). The
reported error is the least-squares error of whichever model was chosen, so it leaves out this
model-selection scatter. Effects: σ_lnR from emulated errors is too small for such glitches, and the
noise calibration (which matches reported errors) slightly understates the true scatter. Real papers
that choose their models on the data have the same property.

**Decision (D4e).** Stated as a limitation; it affects robustness item 5 (σ_lnR ≤ 0.5). Code not changed.

---
## 2026-10-31 · S4 magnitudes on the synthetic fixture under D4e (results/emulation_s4_fixture.json)

**Run.** As entry 3, with the D4e rules: ν̇ noise calibrated once per pulsar per σ under S1m on the
pulsar's own method mix (unknown as window) and used in every variant; the D13 inclusion rule
(≥ 4 completed intervals in the analysed data) applied to every replicate; a replicate discarded if
fewer than 2 pulsars remain. 500 replicates per variant per σ, seed 20261031.

| Variant | σ = 0.2 fwd / bwd | σ = 0.4 fwd / bwd | σ = 0.6 fwd / bwd |
|---|---|---|---|
| S4w (window) | +0.017 / −0.004 | +0.055 / −0.011 | +0.101 / −0.018 |
| S4l (local-relaxation) | +0.001 / −0.007 | −0.003 / −0.022 | −0.009 / −0.035 |
| S4g (global) | +0.001 / −0.014 | +0.003 / −0.041 | −0.002 / −0.054 |
| S4g0 (global-noexp) | −0.005 / −0.038 | −0.007 / −0.112 | +0.024 / −0.180 |
| S4-mix, all data (unknown as window) | −0.002 / −0.003 | −0.007 / −0.015 | −0.008 / −0.026 |
| S4-mix, all data (unknown as global) | −0.001 / −0.010 | −0.001 / −0.033 | −0.000 / −0.048 |
| S4-mix, admissible only (unknown as window) | +0.003 / −0.012 | +0.003 / −0.034 | −0.000 / −0.042 |

Standard errors of the means 0.001–0.003. Central 95% of the admissible-only S4-mix forward b̂:
[−0.041, +0.054] at σ = 0.2, [−0.093, +0.096] at σ = 0.4, [−0.125, +0.124] at σ = 0.6.

**Inclusion (D4e).** No pulsar dropped and no replicate discarded in the single-method and all-data
variants. Admissible-only: about 0.73 of 10 pulsars dropped per replicate (358–382 per 500 replicates);
no replicate discarded. No fit failed.

**Calibration under S1m (D4e) does not reach the target for most pulsars.** Outside the 5% tolerance:
6 of 10 pulsars at every σ (σ = 0.2: F1, F3, F5, F6, F7, F9; σ = 0.4 and 0.6: F1, F2, F3, F5, F6, F9),
all "at floor" (zero added noise already gives a median reported relative error above the target;
e.g. F6: 0.13–0.16 against 0.076). Reason: S1m removes the misfit only for local-relaxation. With
t_rec = T and noise switched off, the median reported relative error is about 1e-13 for local-relaxation
but 0.08–0.12 for global (one ν̈ cannot follow recovery slopes J_k/T_k that differ from window to
window), 0.10–0.18 for global-noexp and 0.25–0.38 for window (no recovery term). A pulsar whose mix
contains these methods therefore sits at the floor, and its emulated errors are larger than the real ones.

**Decision.** Pending.

---
## 2026-10-31 · D4f admission rule for short-window: NOT admitted (results/short_window_admission.json)

**Run.** `python3 scripts/emulation_fixture.py`, same synthetic fixture, seed and D4e rules as entry 5
(the S4 variants of entry 5 reproduce exactly). New: S4sw, a constant ν̇ level over w days before and
after each glitch, each span cut at the neighbouring glitch, no exponentials (D4f), at the fixture
cadence of 14 d (D4c fallback, since the audit has no cadences yet), 500 replicates per w per σ.
Rule (D4f): admitted only if mean + 1.96 SE of the forward slope ≤ +0.02 at every σ and every w.

| Variant | σ = 0.2 fwd / bwd | σ = 0.4 fwd / bwd | σ = 0.6 fwd / bwd | Upper 95% bound (σ = 0.2 / 0.4 / 0.6) |
|---|---|---|---|---|
| Espinoza-type, w = 40 d | n/a | n/a | n/a | none: all 500 replicates discarded |
| Espinoza-type, w = 150 d | +0.000 / −0.012 | +0.004 / −0.040 | +0.006 / −0.070 | +0.001 / +0.006 / +0.009 |
| Wang-type, w = 125 d | −0.000 / −0.008 | +0.000 / −0.031 | +0.006 / −0.056 | +0.001 / +0.003 / +0.010 |
| Wang-type, w = 300 d | +0.001 / −0.021 | +0.004 / −0.061 | +0.013 / −0.110 | +0.002 / +0.006 / +0.017 |

**Verdict: not admitted.** The only failing check is w = 40 d: at a 14-d cadence a 40-d span holds 3
samples (days 0, 14, 28), fewer than the 4 that D4c requires after a glitch, so no glitch is measurable,
every replicate has no pulsar left and is discarded (D4e), and no bound exists. A check with no estimate
counts as failed, because the rule admits short-window "only if" the bound is shown to be ≤ +0.02.
Every check that could be computed passes, the largest bound being +0.017 (w = 300 d, σ = 0.6).
Short-window is therefore analysed in robustness item 2 only (D4f).

**Observations.** Short-window behaves like window over a short span: little forward bias (the spans
rarely reach the next glitch), and a negative backward bias that grows with σ and with w (the "before"
level still holds part of the previous step's recovery, which depends on the previous wait).

**Note.** The verdict rests on the provisional fixture cadence. Espinoza+2011 used Jodrell Bank timing,
whose cadence is typically much shorter than 14 d; at a cadence of 13 d or less a 40-d span holds 4
samples. The D4f rule is applied as written; whether the fixture should use the audited source cadence
for S4sw is for the user to decide.

**Superseded by entry 7** (kept here unchanged). That run used the 14-d fallback cadence because the
source audit did not exist yet. Entry 7 reruns S4sw at the audited cadences, which is D4c's rule
("Cadence: per source from the audit. Until the audit exists, 14 days"), now that the audit exists.
The rerun applies that existing rule; it is not a change made because of this result.
Saved copies of this run: results/short_window_admission_cadence14.json,
results/emulation_s4_fixture_sw_cadence14.json.

---
## 2026-10-31 · D4f admission rule rerun at the audited cadences: ADMITTED (results/short_window_admission.json)

**Why rerun.** D4c: "Cadence: per source from the audit. Until the audit exists, 14 days." The audit
now exists (data/meta/source_audit.csv, cadence_num_days): Espinoza+2011 6 d, Wang+2000 28 d. This
rerun applies that rule to the S4sw variants. It is not a change made because of entry 6's result,
which is kept above.

**Run.** `python3 scripts/emulation_fixture.py`. Everything as in entry 6 (synthetic fixture, seed
20261031, ν̇ noise calibrated under S1m at the 14-d fixture cadence, D4e inclusion, 500 replicates per
w per σ), except that the S4sw variants are sampled at their source's audited cadence. The other 18
variants reproduce entry 5 exactly.

| Variant (cadence) | σ = 0.2 fwd / bwd | σ = 0.4 fwd / bwd | σ = 0.6 fwd / bwd | Upper 95% bound (σ = 0.2 / 0.4 / 0.6) |
|---|---|---|---|---|
| Espinoza-type, w = 40 d (6 d) | +0.001 / −0.003 | +0.002 / −0.011 | −0.000 / −0.020 | +0.002 / +0.004 / +0.004 |
| Espinoza-type, w = 150 d (6 d) | +0.000 / −0.012 | +0.004 / −0.037 | +0.006 / −0.070 | +0.002 / +0.007 / +0.010 |
| Wang-type, w = 125 d (28 d) | +0.000 / −0.009 | +0.000 / −0.028 | +0.002 / −0.054 | +0.001 / +0.003 / +0.006 |
| Wang-type, w = 300 d (28 d) | +0.002 / −0.019 | +0.005 / −0.060 | +0.013 / −0.109 | +0.003 / +0.007 / +0.016 |

No replicate discarded, no pulsar dropped, no fit failed.

**Verdict: admitted.** All 12 checks pass (mean + 1.96 SE ≤ +0.02 at every σ and every w, both source
types; largest bound +0.016 at w = 300 d, σ = 0.6). Under D4f, short-window enters the primary analysis.
The backward slopes are negative and grow with σ and w, as for window (entry 6); they enter the S1m
placebo reference through the measurement mix (D4d).

---
## 2026-11-01 · 06_validate (D4j/D4k): T2 false positives under S3 and a slightly liberal LR test

**Run.** `caffeinate -i python3 scripts/06_validate.py` at commit ab790d0; real structure of the 21 D13
pulsars; sigma grid 0.2, 0.4, 0.6, 0.8, 1.2; replicates per D4b/D4h; seed 20261031. No replicate failed
in any of the 60 jobs. results/validation.json.

**S3 T2 false-positive rate fails the D4h bound (0.040) at sigma 0.2-0.8:** 0.109, 0.069, 0.053, 0.044;
0.030 at sigma 1.2. Diagnostic (not part of 06): S3 at sigma 0.2, 400 replicates, T2 rejection 0.105 with
measurement noise and 0.037 without it; mean a(M4) 0.990 vs 0.999, mean d +0.0050 vs -0.0006.
Explanation: the lognormal noise on Δν/ν (real dff_err) slightly attenuates the size coefficient a, and
ln(Δν̇/ν̇), which is correlated with ln(Δν/ν), takes up the remainder (errors in variables; D6 ignores
errors in ln R in the primary fit). The bias in d is tiny, but at small sigma its standard error is small
enough for T2 to reject it. The effect shrinks as sigma grows.

**S2 false-positive rates are slightly above the nominal 0.025 everywhere** (T1 and T2 at five sigmas:
0.029-0.041, mean 0.0345, about 6 Monte-Carlo SE above 0.025); two cells exceed the bound marginally
(T1 at sigma 0.4: 0.040; T2 at sigma 0.8: 0.041). The chi-square approximation of the LR test is
slightly liberal with about 20 pulsar intercepts and about 200 intervals.

**Other observations.** S1 (statistics layer) is "inconclusive" in 94-100% of replicates because its
b-hat (about 0.79) lies above the central 95% of the S1m reference (about 0.04-0.40), which the
measurement emulation attenuates. Under S1m at sigma 0.2, 13.3% of replicates are classified C because
T1 power (0.864) is just above 0.8, so every non-rejection counts as a refutation. S4w (window bias) gives
A2 in 30-67% of replicates at sigma >= 0.4; the real admissible data contain no window-class glitches.
S1m missingness: estimable in 1000 of 1000 replicates at every sigma; mean log odds ratio about -0.93.

**Decision.** Pending.

---
## 2026-11-01 · 06 under D4l: all criteria pass; A1 is reached under the size-only world (S3)

**Run.** `caffeinate -i python3 scripts/06_validate.py --reuse-evaluation results/validation_partial` at
596adca: the evaluation replicates of ab790d0 read-only; new runs only for the D4l nulls (S2, S3; 1,000
per sigma, job ids 41-42) and the S1/S3/S1m diagnostic reruns (same seeds; b-hat reproduced exactly).
No replicate failed in any of the 85 jobs. results/validation.json.

**Criteria.** 105 of 105 pass. Directional false-positive rates: S2 T1 0.017-0.031, S2 T2 0.008-0.036,
S3 T2 0.020-0.030 (bound 0.040). Thresholds (97.5th percentiles of the signed root under the nulls):
T1 1.95-2.20, T2 1.95-2.49. S1m powers: T1 0.894, 0.791, 0.752, 0.615, 0.190 and T2 0.397, 0.366,
0.259, 0.221, 0.063 at sigma 0.2, 0.4, 0.6, 0.8, 1.2.

**Outcome rates (selected).** S1m: A1 + A2 = 0.871, 0.769, 0.734, 0.597, 0.185; C = 0.023-0.025.
S4-mix: A1 + A2 = 0.000, 0.006, 0.024, 0.026, 0.012. **S3 (size only): A1 = 0.834, 0.811, 0.869, 0.872,
0.863**, A2 = 0.012-0.022: under the rival "size does the work", T1 rejects and b-hat (about 0.80, from
the correlation of R with Δν/ν) falls inside the A bounds, so A1 is reached most of the time. Only A2
(T2) separates the model from size-only, and T2 power is 0.06-0.40. S1: A2 = 0.94-0.97, with the label
"less attenuated than emulated" in every replicate.

**Report-only diagnostics (D4l item 4).**
- S3 rejections under the old chi-square T2: 27-36% had d-hat < 0 (most leaked in the positive direction).
- SD of b-hat across replicates / mean model SE: S1 2.77 to 1.57 (sigma 0.2 to 1.2), S3 1.96 to 1.52,
  S1m 1.27 to 1.12; for d-hat: S1 3.09 to 1.61, S3 1.45 to 1.09, S1m 1.21 to 1.10. With measurement error
  in ln R the model-based (Wald) SEs understate the replicate-to-replicate spread, most at small sigma.
- S1m, slope of measured on true ln(Δν̇/ν̇): global 0.61-0.68, local-relaxation 1.02-1.10, short-window
  1.19-1.23.

**Decision.** Pending.

---
## 2026-11-01 · 06 under D4m: measured size-only null (S3m); b-hat does not separate creep from size-only

**Run.** `caffeinate -i python3 scripts/06_validate.py --reuse results/validation_partial --reuse
results/validation_partial_d4l --store results/validation_partial_d4m` at 004dffc. New jobs only: S3m as an
evaluation scenario (job 7, 1,000 per sigma) and as a T2 null (job 43, 1,000 per sigma). No failed
replicates in any of the 95 jobs. 120 of 120 criteria pass; A2 rate under S3m 0.009, 0.020, 0.015, 0.021,
0.019 at sigma 0.2, 0.4, 0.6, 0.8, 1.2 (bound 0.040).

**T2 threshold (D4m: the larger of S3 and S3m).** S3 / S3m 97.5th percentiles: 2.49/1.90, 2.15/2.09,
2.17/2.02, 1.95/2.13, 1.97/1.82; S3m is the larger only at sigma 0.8, where S1m T2 power falls from 0.221
(D4l) to 0.161. T1 power unchanged: 0.894, 0.791, 0.752, 0.615, 0.190.

**S3m outcome rates (report only).** A1 = 0.796, 0.667, 0.625, 0.521, 0.212; A2 = 0.009-0.021;
inconclusive 0.13-0.73; C 0.02-0.05; D_placebo 0.02-0.04. Mean b-hat under S3m (0.229, 0.208, 0.194,
0.173, 0.129) is almost the same as under S1m (0.230, 0.217, 0.204, 0.182, 0.110): once measured, the size-only
world and the creep world give the same slope on ln R.

**S1m against S3m (report only).** Fraction of S1m b-hat above S3m's 97.5th percentile: 0.037, 0.038,
0.033, 0.027, 0.013. Fraction of S1m d-hat below S3m's 2.5th percentile: 0.905, 0.666, 0.447, 0.234,
0.069. The separation lies in d, not in b.

**Interval factors (D4m; SD / mean SE under the nulls).** b (S2): 1.05, 1.09, 1.10, 1.04, 1.07;
d (S3m): 1.03, 1.07, 1.10, 1.08, 1.06. These are much smaller than the SD/SE ratios under the model-true
scenarios (S1: 1.6-2.8 for b; S1m: 1.1-1.3), because measurement error inflates the spread of b-hat and d-hat
mainly when there is a real slope to attenuate.

**Decision.** Pending.

---
## 2026-11-01 · 06 under D4n: T2 on d-hat against the S3/S3m null references (from saved stores only)

**Run.** `python3 scripts/06_validate.py --reuse results/validation_partial --reuse results/validation_partial_d4l
--reuse results/validation_partial_d4m --store results/validation_partial_d4n --saved-only` at 0ecf505. No
simulation: every job and calibration came from the three saved stores (the new store holds only meta.json).
120 of 120 criteria pass.

**T2 (D4n).** Cut on d-hat (the smaller 2.5th percentile under S3-null / S3m-null): -0.0408, -0.0632, -0.0913,
-0.1250, -0.1605 at sigma 0.2, 0.4, 0.6, 0.8, 1.2 (S3's is the smaller at 0.2 and 1.2). S1m T2 power: 0.866,
0.666, 0.447, 0.234, 0.062 (D4m signed-root T2: 0.397, 0.366, 0.259, 0.161, 0.063). T2 false-positive rates:
S2 0.001-0.029, S3 0.014-0.021, S3m 0.009-0.025, S4-mix 0.004-0.018.

**Outcome changes.** S1m A2: 0.799, 0.607, 0.424, 0.225, 0.051 (D4m: 0.383, 0.347, 0.252, 0.156, 0.057). At
sigma 0.2 the T2 power is >= 0.8, so B becomes reachable: B = 0.072 under S1m, 0.836 under S3, 0.797 under S3m.
A2 under S3m: 0.008, 0.022, 0.023, 0.019, 0.022 (bound 0.040). S4w A2: 0.002, 0.334, 0.594, 0.706, 0.794
(D4m: 0.002, 0.314, 0.592, 0.704, 0.804).

**Mean fitted sigma-hat of M2 (D4n item 4)** at generating sigma 0.2 / 0.4 / 0.6 / 0.8 / 1.2: S1m 0.824 / 0.855
/ 0.909 / 0.972 / 1.100; S3m 0.876 / 0.909 / 0.953 / 0.997 / 1.094; S2 0.189 / 0.378 / 0.563 / 0.748 / 1.112;
S4-mix 0.189 / 0.375 / 0.544 / 0.695 / 0.961; S1 0.624 / 0.704 / 0.815 / 0.958 / 1.264; S3 1.102 / 1.149 /
1.224 / 1.320 / 1.563.

**Short-window label (D4n item 3).** Not available from the saved stores: the replicates predate D4n and the
datasets are not saved, so the refit without short-window glitches needs new simulation.

**Decision.** Pending.

---
## 2026-11-01 · Placeholder-settings record: sigma-hat floors, short-window label, 05 rerun (D14 provisional Q, tau_d)

**sigma-hat floors (report only; scripts/sigma_floor.py at 84998ac; generating sigma = 0, 200 replicates,
seed index 7).** Mean fitted sigma-hat of M2 [2.5th, 97.5th percentile]: S1m 0.799 [0.636, 0.962]; S3 1.081
[0.898, 1.261]; S3m 0.890 [0.760, 1.039]. S4-mix: all 200 replicates fail (M0 does not converge): with sigma = 0
and no covariate in the waiting times, each pulsar's completed waits are identical, so sigma-hat -> 0. The
same holds for S2. Consequence for 08: D4n's bisection starts at sigma = 0, and for S2 and S4-mix that step has
no analysed replicate, so the matching stops (SigmaMatchError); D4n has no rule for this.

**Short-window label (D4n item 3; scripts/short_window_label.py at 84998ac).** A2 replicates regenerated from
their seeds; every one reproduced its saved b-hat. Fraction labelled "short-window-sensitive" among A2
replicates with a refit, sigma 0.2 / 0.4 / 0.6 / 0.8 / 1.2: S1m 0.220 (176/799) / 0.211 (128/607) / 0.146
(62/424) / 0.089 (20/225) / 0.020 (1/51); S3m 0.000 (0/8) / 0.000 (0/22) / 0.000 (0/23) / 0.000 (0/19) /
0.045 (1/22); S4-mix none in A2 / none / 0/4 / 0/5 / 0/2. S4sw-espinoza and S4sw-wang: not computable (every
glitch is short-window, so no data remain).

**05 rerun (5a19496).** Reproduces 06 at sigma 0.2, 0.4, 0.6: 756 values compared (thresholds, references,
calibration, powers, rates, outcomes, sigma-hat), 0 differences.

---
## 2026-11-01 · Full rerun under D4o (literature exponential recoveries) and D4p, against the placeholder record 7e5b134

**Runs (0a59496).** `python3 scripts/06_validate.py --store results/validation_partial_d4o` (nothing reused), then
05, sigma_floor.py, detection_fractions.py. No failed replicates. 06: all 105 criteria pass (25 D14/D4h/D4m
checks + 80 failed-replicate checks; the placeholder record's 120 included 15 diagnostic-rerun jobs). 05
reproduces 06 at sigma 0.2, 0.4, 0.6. Statistics-layer jobs (S1, S2, S3, S1-noise-off, S2-null, S3-null) are
unchanged replicate for replicate; their outcomes move only through the new references.

**S1m powers (new / placeholder), sigma 0.2 / 0.4 / 0.6 / 0.8 / 1.2.** T1 0.990/0.894, 0.980/0.791, 0.952/0.752,
0.902/0.615, 0.492/0.190. T2 (D4n) 0.954/0.866, 0.828/0.666, 0.597/0.447, 0.371/0.234, 0.109/0.062.

**S1m b-hat 95% percentile interval (new | placeholder).** [0.216, 0.668] | [0.077, 0.402]; [0.202, 0.621] |
[0.061, 0.386]; [0.134, 0.586] | [0.043, 0.375]; [0.093, 0.572] | [-0.000, 0.369]; [-0.031, 0.467] | [-0.083, 0.305].
Mean 0.442, 0.407, 0.370, 0.325, 0.222 (placeholder 0.230, 0.217, 0.204, 0.182, 0.110).

**Outcomes.** S1m A2 0.916, 0.799, 0.581, 0.357, 0.104 (placeholder 0.799, 0.607, 0.424, 0.225, 0.051); B 0.034,
0.149 at sigma 0.2, 0.4. S3m: B 0.895, 0.858 at 0.2, 0.4; A1 0.812, 0.784, 0.455 at 0.6-1.2; A2 0.016-0.024
(bound 0.040). S3: B 0.871, 0.876 at 0.2, 0.4. S4w A2 0.000, 0.000, 0.094, 0.510, 0.596 (placeholder 0.002, 0.334,
0.594, 0.706, 0.794). S4-mix A2 <= 0.010.

**Mean fitted sigma-hat (M2), new / placeholder.** S1m 0.760/0.824, 0.804/0.855, 0.867/0.909, 0.943/0.972,
1.087/1.100; S3m 0.836/0.876, 0.876/0.909, 0.929/0.953, 0.979/0.997, 1.081/1.094; S2 and S4-mix unchanged to
3 decimals (S2 identical; S4-mix within 0.003).

**Floors (D4p lower end 0.01; placeholder at 0).** S1m 0.736 [0.608, 0.847] (0.799 [0.636, 0.962]); S3 1.080
[0.884, 1.246] (1.081); S3m 0.831 [0.721, 0.939] (0.890); S4-mix 0.009 [0.009, 0.010], all 200 analysed
(placeholder at 0: all 200 failed).

**Detection (D4o item 5; detected / all glitches, S1m at sigma 0.2; placeholder in brackets).** J0537-6910 0.273
(0.415); others 0.176 (0.271); Crab 0.041 (0.008). Published: J0537 12/66 = 0.18; Yu et al. 27/107 = 0.25.
J0537 has false detections (Q = 0 below 20 µHz): 1289 at sigma 0.2. 23 of 24 Crab pool glitches are
short-window, which never fits exponentials.

**Short-window label (06 records).** S1m 296/916, 246/799, 149/581, 66/357, 12/104 of A2 labelled.

**Decision.** Pending.

---
## 2026-11-01 · Post-freeze: S4 fixture and short-window admission rerun under D4o (new files; originals kept)

**Run.** `python3 scripts/emulation_fixture.py` at 44f7b14, which now writes results/emulation_s4_fixture_d4o.json
and results/short_window_admission_d4o.json. The originals (results/emulation_s4_fixture.json,
results/short_window_admission.json, committed at 05e112e under the D14 placeholder Q and tau_d) are untouched:
they are the record of the D4f admission, and 02_clean.py reads the original admission file. The fixture's
synthetic pulsars all fall in D4o's "others" group.

**Admission (D4f rule, upper 95% bound <= +0.02).** Still admitted. Largest upper bound 0.0129 (S4sw-wang w300,
sigma 0.6), against 0.0158 in the original; all 12 bounds within 0.0029 of the original.

**S4 forward slopes, original -> D4o (sigma 0.2 / 0.4 / 0.6).** S4w +0.017/+0.055/+0.101 -> +0.016/+0.051/+0.091;
S4l +0.001/-0.003/-0.009 -> -0.000/-0.010/-0.024; S4g +0.001/+0.003/-0.002 -> +0.001/+0.006/+0.005;
S4g0 -0.005/-0.007/+0.024 -> -0.003/-0.005/+0.024; S4-mix -0.002/-0.007/-0.008 -> -0.001/-0.004/-0.002.
