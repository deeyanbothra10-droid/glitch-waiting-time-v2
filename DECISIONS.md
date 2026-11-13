# Analysis decisions (D1–D15)

Status: PROPOSED by Claude, 31 Oct 2026, for a personal (non-competition) project. Not yet frozen.
Checked against the catalogue file (structure and counts only, no waiting-time statistics), but not
against the source audit. Items marked **[after audit]**
must be confirmed once the source audit and end-of-monitoring dates exist. Revising them then is fine:
nothing here depends on waiting times. Freeze (git tag + prereg/FROZEN) before any statistic relating R to T.

Notation: admissible data = intervals whose Δν̇ was measured by a method that models post-glitch
relaxation (`local-relaxation`, or `global` with exponentials included).

---

## D1 · Research question
**Primary:** Within individual glitching pulsars, does the waiting time to the next glitch scale with
R = (Δν/ν)/(Δν̇/ν̇) as the vortex-creep model predicts, and does the Δν̇ term add predictive power beyond
glitch size alone, when Δν̇ is measured by relaxation-modelled fits?
**Secondary (exploratory):** Is the relation stronger in Vela-like than in Crab-like pulsars?
*Why:* this asks the part only vortex creep claims (the Δν̇ term), and restricts the confirmatory test
to data that the bias emulation says can discriminate.

## D2 · Hypotheses (tested on admissible data)
- **T1 (scaling):** in M2, b > 0. Null: b = 0.
- **T2 (the Δν̇ term):** in M4, d < 0. Null: d = 0.
- The model predicts b = 1 and d = −1, but noise attenuates estimates, so "consistent with the model" is
  judged against the S1 simulation (see D6, D12), not against the raw value 1.

## D3 · Primary model
Lognormal accelerated-failure-time model, right-censored, one intercept per pulsar, one shared σ.
T1 uses M2 (b ln R + c_p); T2 uses M4 (a ln Δν/ν + d ln Δν̇/ν̇ + c_p).
*Why:* handles censoring; gives effect sizes; M4 asks the sharpest question. Shared σ because most
pulsars have only 4–10 intervals.

## D4 · Criterion and significance
Likelihood-ratio tests: T1 (M2 vs M0, H0 b = 0) and T2 (M4 vs M3, H0 d = 0), two-sided.
Family-wise α = 0.05 with Bonferroni: each test at α = 0.025.
Everything else (Spearman, permutation, LOO-CV, AIC, implied β, class comparison) is secondary:
reported with estimates and intervals, no significance claims.

## D4a · Secondary-analysis details
- Permutation test: within each pulsar, permute ln R across ALL of that pulsar's intervals, completed and
  censored together; each interval keeps its own T and censoring flag, only ln R moves. Statistic: the M2
  slope estimate b̂ from the full censored likelihood. Two-sided p = (1 + #{|b̂_perm| ≥ |b̂_obs|}) /
  (1 + N), with N = 10,000 permutations and seed 20261031. The placebo permutation (D15) uses the same
  scheme with the previous wait as outcome.
  Why: under the null, R is exchangeable within a pulsar whatever the censoring, so permuting R (not T)
  keeps every censored interval valid and uses all the data.
- Spearman: completed intervals only, computed per pulsar for pulsars with ≥ 4 completed admissible
  intervals; pooled as a Fisher-z weighted mean with weights (n_p − 3). Descriptive only: no p-values,
  and always labelled "completed intervals only".
  Why: rank correlation cannot handle censoring; dropping censored intervals biases towards short waits,
  so it is reported as a description, never as a test.
- LOO-CV: leave out one interval at a time, refit, score the held-out interval by log f(T) if completed
  or log S(c) if censored; report the mean over all intervals for each model. Not leave-one-pulsar-out,
  because a held-out pulsar's intercept cannot be predicted.

## D4b · Clarifications (31 Oct)
### Interval estimates (D4)
- All reported intervals are 95%: Wald intervals from the inverse Hessian for model parameters; percentile intervals for simulated quantities; the pooled Spearman uses Fisher z ± 1.96/sqrt(Σ(n_p − 3)).
- Spearman edge case (D4a): a pulsar with |r_s| = 1 is reported and flagged, and excluded from the pooled Fisher-z mean. The excluded pulsars are listed.

### Simulation (D14)
- Typical synthetic wait per pulsar: W_p = (obs_end − first-glitch mjd) / (number of glitches). This is an explicit, permitted exception to rule 3: it carries the pulsar's glitch RATE only and no information about any R–T relation.
- Centring: in every scenario, covariates are centred within pulsar, and c_p = ln W_p.
  S1: ln T = (ln R − mean_p ln R) + c_p + σε.
  S2: ln T = c_p + σε.
  S3: ln T = a·(ln dnu_frac − mean_p ln dnu_frac) + c_p + σε, with a = 1.
- S4 recovery time: t_rec ~ lognormal(median W_p, log-sd 0.3), drawn independently of T. It never uses real or simulated waits.
- Replicates: 1,000 per scenario per σ value for S1–S3, and 500 per S4 variant per σ value, with seed 20261031.
- Validation cost limits: the LR tests run on every replicate. Permutation tests (primary and placebo) run with 999 permutations on the first 200 replicates of each scenario × σ. LOO-CV runs on the first 100. The real run uses 10,000 permutations, as in D4a.
- Method class: an S4g0 (global without exponentials) case is df1_method = global with exp_fitted = no. If exp_fitted = sometimes, classify per glitch when the source audit notes say which glitches had exponentials; otherwise the class is unknown. For admissibility (D9), global requires exp_fitted = yes.

### Diagnostics (D15)
- Collinearity: Pearson correlation of ln(Δν/ν) and ln(Δν̇/ν̇) over usable glitches, for pulsars with ≥ 4 admissible usable glitches. The D15 rule uses the unweighted median of |r| across those pulsars. Also report the pooled within-pulsar (demeaned) correlation.
- Placebo: outcome = the previous wait T_(k−1), for every interval whose starting glitch has a predecessor (that wait is always completed). Test: the D4a permutation scheme, two-sided. "Significantly negative" means p < 0.05 and b̂ < 0; "significantly positive" means p < 0.05 and b̂ > 0.
- Missingness: conditional (pulsar-stratified) logistic regression of "df1 missing" on ln(following wait), over glitches in the included pulsars whose following wait is completed. Glitches whose following wait is censored are excluded and counted. Two-sided α = 0.05. Report the odds ratio with a 95% interval. Also report descriptively the same model with "df1 ≤ 0" as the outcome.
  CODE CHANGE: replace the per-pulsar-intercept logistic fit with the conditional logit likelihood, written with numpy/scipy, and test it against a hand-worked example.

## D4c · Simulation design (31 Oct)
### S4 expectations, revised before any real data
- S4l (local-relaxation), forward: not positive (upper end of the 95% interval ≤ +0.02). Reason: when t_rec < T, the step's recovery stops partway through the window. The local fit's exponential term absorbs that bend and overstates Δν̇ for long waits, which gives a small negative bias. The bias works against T1 and is reported as a finding. S4w, S4g and S4g0: expectations unchanged.

### Two layers of simulation
- S1, S2 and S3 test the statistics code only. They do not use the measurement emulation, and t_rec is not used. Measured Δν̇ = true Δν̇ × exp(e), where e is normal with SD = the real relative error (df1_err/df1) of the glitch that was resampled.
- S1m = S1 passed through the measurement emulation, with t_rec = T. Its b̂ distribution is the reference for D6 and D12 ("the S1-simulated b̂ distribution" means S1m from now on).
- S4 variants: as before, through the emulation, with t_rec independent of T.
- Mixtures (S1m and the S4 mix): each simulated glitch takes the method of the real glitch it was resampled from. For the primary analysis, only admissible glitches are used. In all-data mixtures, "unknown" is emulated as window (worst case), with a sensitivity run that emulates it as global.

### Emulation settings
- Cadence: per source from the audit. Until the audit exists, 14 days.
- Noise: per pulsar, choose the ν̇ sample noise (by bisection) so that the median relative error of the emulated Δν̇ matches the median real df1_err/df1 of that pulsar's usable glitches. Pulsars with fewer than 3 error values use the global median.
- Pre-glitch monitoring: 365 days before the first glitch.
- Measurability: a glitch needs at least 4 samples after it (and, for window, at least 4 before it). Otherwise its Δν̇ is missing: it is unusable but still ends the previous interval.
- Resolvability: exponentials with τ_d < 2 × cadence are not included in any fit. Their effect stays in the data, as it would in a real paper.

### Glitch counts
- Use c_p = ln W_p − σ²/2, so the MEAN synthetic wait equals W_p and simulated counts match real counts on average.

## D4d · S4 as measurement; calibrated decision rules (31 Oct)
### S4 is a measurement, not a test (restores D14)
- Remove every pass/fail test on S4 bias directions or sizes. S4 magnitudes depend on provisional settings. They are computed, saved to results/ and logs/emulation_findings.md, and reported.
- Replace them with correctness tests of the emulator on noise-free toy cases with known answers:
  (a) window with no exponentials gives J(1 − T/(2·t_rec)) for T ≤ t_rec;
  (b) window with one exponential and constant step gives Δν̇_p + (Δν_d/T)(1 − e^(−T/τ_d));
  (c) global and local with every exponential fitted and t_rec = T are exact;
  (d) an exponential below the detection rule is never fitted.

### Exponential detection (replaces the D4c resolvability rule)
- In each fitted window, include an exponential term only if including it lowers AIC by at least 10 (τ_d known to the fit). This replaces "τ_d < 2 × cadence is never fitted". Analysts fit what is significant; what is not significant is, by construction, small compared with the noise.

### Decision rules calibrated against simulation (amends D12 and D15)
- Reference distributions come from the same pulsars, counts, measurement mix (admissible only for the primary analysis) and noise as the real data: S1m (model true) and S4-mix (bias only, t_rec independent of T), with 1,000 replicates each at the σ estimated from the real primary fit.
- Support (A) additionally requires: observed b̂ above the 97.5th percentile of the S4-mix b̂ distribution (beats bias alone), and observed placebo b̂ not below the 2.5th percentile of the S1m placebo distribution.
- Artefact (D): observed b̂ inside the central 95% of the S4-mix distribution, or observed placebo b̂ below the 2.5th percentile of the S1m placebo distribution.
- The raw placebo permutation p-value is still reported, but no longer decides anything on its own.

### Answers to the four open questions
1. Δν/ν noise: yes, in every layer: multiply by exp(e), with e normal and SD equal to the real relative error dff_err/dff of the resampled glitch.
2. Calibration: 50 simulated sequences per pulsar per bisection step; tolerance 5% of the target; record the achieved value.
3. Cadence within a pulsar: each interval uses the cadence of the source of the glitch that starts it; the pre-glitch baseline uses the previous interval's cadence.
4. Missing errors: a usable glitch without df1_err (or dff_err) uses its pulsar's median relative error (the global median if the pulsar has fewer than 3 values), and is flagged.

## D4e · Noise calibration and replicate rules (31 Oct)
- Readings confirmed: global uses the same per-window AIC decision to choose which exponentials enter the global fit; the 365-day baseline before the first glitch uses the first glitch's cadence.
- Calibration target: match the median of the fit's REPORTED (formal) relative error of Δν̇ to the median real df1_err/df1, not the scatter of measured against true values. Real df1_err values are formal errors from noise and do not include model misfit.
- Calibration scenario: calibrate once per pulsar under S1m (t_rec = T, so no misfit), at each σ, and use that same ν̇ noise level in every scenario (S1m and all S4 variants). Noise is a property of the data, not of the physics. Report, per pulsar, the achieved value and any pulsar outside the 5% tolerance.
- Replicates with too few intervals: apply the D13 inclusion rule (≥ 4 admissible completed intervals) to each simulated replicate, exactly as for real data. Pulsars that fail it are dropped from that replicate's fit; the replicate itself is kept. Report how many pulsars were dropped per scenario. A replicate is discarded only if fewer than 2 pulsars remain, and those are counted.
- Limitation recorded: when detection is marginal, the reported Δν̇ errors understate the true scatter by 2–4× (emulation_findings entry 4). Real catalogue errors probably share this. It affects robustness item 5 (σ_lnR ≤ 0.5) and is stated as a limitation; no code change.

## D4f · Audit-driven method classes (31 Oct)
- New class `short-window`: Δν̇ from separate pre-/post-glitch local solutions over spans fixed by the analyst (not ending at the next glitch), differenced at the glitch epoch, transients not fitted. Applies to Espinoza+2011 (post-glitch span ideally 30–50 d) and Wang+2000 (~100–300 d).
- Emulation S4sw: constant ν̇ level fitted over w days before and w days after the glitch (each truncated at the neighbouring glitch), no exponentials; Δν̇ = post level − pre level. w = 40 d (sensitivity 150 d) for Espinoza-type; w = 125 d (sensitivity 300 d) for Wang-type.
- Admission rule, fixed before any real data: short-window enters the primary analysis only if, in the synthetic fixture, the S4sw forward slope's 95% upper bound is ≤ +0.02 at every σ, for both the main and the sensitivity w. Otherwise it is analysed in robustness item 2 only. The result is recorded in emulation_findings.md either way.
- exp_fitted = "when-detected" counts as "yes" for admissibility (it is what D4d's detection rule models).
- "Our measurement" rows are global (assumed, following Basu et al. 2022); admitted; a robustness run excludes them; to be confirmed with the catalogue maintainers.
- Robustness: the primary analysis without J0537−6910 is always reported, labelled exploratory if it has fewer than 40 intervals.
- Sources not yet audited remain `unknown` (excluded from the primary analysis).

## D4g · Censoring dates and the J0537 gap (31 Oct)
- J0537−6910 gap: RXTE monitoring ended MJD 55926.9 (31 Dec 2011; Antonopoulou et al. 2018); NICER monitoring began MJD 57982 (17 Aug 2017; Ho et al. 2026, Table 2). The interval starting at MJD 55819 is censored at 55926.9 (D11).
- obs_end, documented spans only:
  - J0537−6910: 60827 (NICER data to 2025 Jun 1, Ho et al. 2026).
  - J1023−5746: 59269 (Fermi data end, Gügercinoğlu et al. 2022).
  - Jodrell Bank pulsars whose last glitch is before 58483: 58483 (31 Dec 2018, end of the Basu et al. 2022 sample): J0631+1036, J0205+6449, J0729−1448, J1814−1744.
  - All other pulsars: obs_end = their last glitch epoch, i.e. no censoring information; that last interval is dropped and logged, never entered with c = 0. To be improved if the maintainers or the source papers give later coverage dates.
- New robustness item 10: Jodrell Bank pulsars use obs_end = max(last glitch, 60891), i.e. the latest catalogue measurement (MJD 61071) minus 180 days for reporting lag.
- Jodrell Bank pulsars (for item 10): a pulsar counts as Jodrell Bank-monitored if at least one of its catalogue glitches has a Jodrell Bank primary source (Espinoza+2011, Basu+2021, Basu+2019, any Shaw+ entry, or any "Our measurement" spelling). Resulting included (D13) pulsars, computed 31 Oct: B0531+21, B1727−47, B1737−30, B1757−24, B1758−23, B1800−21, B1823−13, J0205+6449, J0631+1036, J0729−1448, J1737−3137, J1814−1744, J1841−0524, J2021+3651, J2229+6114 (15 of 21; 144 of 218 catalogue pulsars).

## D4h · Intercept, false-positive bound and open items (31 Oct)
- Intercept (replaces D4c's formula, keeps its aim): c_p = ln W_p − σ²/2 − ln m_p, where m_p is the mean over the pulsar's resampled pool of exp(slope × centred covariate) for that scenario (S1/S1m: covariate ln R, slope 1; S3: covariate ln dnu_frac, slope a; S2/S4: m_p = 1). The mean synthetic wait then equals W_p in every scenario.
- False-positive bound (corrects D14): at the per-test α = 0.025, a test passes if its false-positive rate ≤ 0.025 + 3·sqrt(0.025·0.975/n), where n is the number of replicates (≈ 0.040 for n = 1,000). Applies to T1 and T2 under S2, and to T2 under S3.
- Robustness item 9: combined epoch uncertainty = sqrt(mjd_err_start² + mjd_err_end²). An interval is excluded if that exceeds 0.1·T, or if either error is missing (X).
- Robustness item 8: not run, because no step-only Δν̇ values are recorded. Reported as "not run" with that reason.
- Confirmed readings: S1m uses 1,000 replicates per σ; Support (A) uses D4d's percentile placebo condition, with the raw permutation p reported only; the implied-β check counts all admissible intervals, completed and censored.

## D4i · Tiered decision rules after the power analysis (31 Oct)
Reason (recorded before freezing): 05 showed realistic power of about 0.78 for T1 and about 0.34 for T2 at σ = 0.4. It also showed that the implied-β condition depends on pulsar glitch rates, sizes and τ_sd, all known before the freeze, rather than on the R–T relation under test. This revision was made after seeing those simulated values, and no real R–T statistic.
- A1, scaling supported: T1 rejects; observed b̂ above the 97.5th percentile of S4-mix b̂; observed b̂ within the central 95% of S1m b̂; placebo condition as in D4d.
- A2, Δν̇ term supported: A1 holds and T2 rejects (d < 0).
- If T2 does not reject, the Δν̇ term is "not established (T2 power ≈ 0.34)". That is NOT outcome B, because B requires T2 power ≥ 0.8, which this design does not have.
- C, refutation: observed b̂ below the 2.5th percentile of S1m b̂; or T1 does not reject AND the T1 power at the σ estimated from the real primary fit is ≥ 0.8.
- D, artefact: as in D4d.
- Implied β: reported per pulsar with 95% intervals, as a physical-consistency check discussed by class (Vela-like vs Crab-like). It no longer decides any outcome.
- In simulations, outcome rates are computed with these rules.
- D15 missingness: judged against its S1m distribution, like the placebo. Under a true model, short following waits make Δν̇ unmeasurable, so some association is expected.

## D4j · Exclusive outcome order; corrected artefact definition (31 Oct)
Reason: under D4d, D fired whenever b̂ lay inside the central 95% of S4-mix, which also captures a genuine null (b̂ ≈ 0) and overlapped with C. An artefact means a signal that bias could explain, not the absence of a signal. Recorded before freezing; no real R–T statistic has been computed.
Outcomes are assigned in this order; the first that applies wins:
1. D (artefact signature): observed placebo b̂ below the 2.5th percentile of the S1m placebo distribution.
2. C (refutation): observed b̂ below the 2.5th percentile of S1m b̂; or T1 does not reject AND T1 power at the fitted σ is ≥ 0.8.
3. A (scaling supported): T1 rejects, b̂ above the 97.5th percentile of S4-mix b̂, and b̂ within the central 95% of S1m b̂. Within A:
   - A2 if T2 also rejects (d < 0);
   - B if T2 does not reject and T2 power at the fitted σ is ≥ 0.8;
   - otherwise A1 (Δν̇ term not established).
4. D (signal within bias range): T1 rejects but b̂ is at or below the 97.5th percentile of S4-mix b̂.
5. Inconclusive: anything else.
Implied β and the D15 missingness check are reported alongside and never change the outcome.

## D4k · Failed replicates, σ grid, reference σ (31 Oct)
1. A failed replicate (fit failure, timeout, unexpected error) is excluded from its job's rates and listed with its reason. If more than 1% of replicates in any job fail, that job is FAIL until the cause is explained in logs/emulation_findings.md.
2. Add σ = 0.8 and 1.2 to the validation grid (keep existing values). Reason: ln of a memoryless wait has SD 1.28, so real within-pulsar scatter may exceed 0.6, where false-positive control is untested.
3. For the real classification, every percentile and power in D4j comes from simulations run at the real fitted σ̂ of the primary model (D3), never interpolated from the grid; replicate count and seeding as in 06.
4. A missingness check that is not estimable (separation) is reported as such; the S1m reference uses estimable replicates only, with the count reported.
5. In 07–09, any fit failure on real data stops the run; it is never skipped.

## D4l · Calibrated directional tests; corrected C and A bounds (1 Nov)
Reason: 06 (3921bdf) showed (i) T2 false positives under S3 driven by measurement noise, with the leaked bias in d̂ positive, opposite to the predicted sign; (ii) the LR test slightly liberal under S2 (mean 0.0345); (iii) D4j's power clause in C gave 13.3% false refutations under S1m at σ = 0.2; (iv) S1, which has only catalogue-level noise, was "inconclusive" 94–100% of the time because b̂ exceeded S1m's range. No real R–T statistic has been computed.
1. T1 and T2 are directional. Statistic: the signed root of the LR statistic, signed so that positive means the predicted direction (b > 0 for T1, d < 0 for T2). Each rejects when its statistic exceeds the 97.5th percentile of the same statistic under its null scenario (S2 for T1, S3 for T2), simulated at the same σ with 1,000 replicates and seeds disjoint from the evaluation replicates. In 08 the null references are run at the real σ̂, with S3's size coefficient set to the real M3 estimate. Powers in D4j/D4k use these tests.
2. C is: observed b̂ below the 2.5th percentile of S1m b̂. The power clause is removed.
3. A's slope condition becomes: b̂ above the 97.5th percentile of S4-mix b̂, at or above the 2.5th percentile of S1m b̂, and at or below the 97.5th percentile of S1 b̂. Report b̂'s percentile within S1m; if it is above S1m's 97.5th percentile, add the label "less attenuated than emulated" (report only).
4. Report only, never outcome-changing: for S3 T2 rejections under the old test, the fraction with d̂ < 0; for S1, S3 and S1m, the SD of b̂ and d̂ across replicates against the mean model SE; for S1m, the slope of measured on true ln(Δν̇/ν̇) per method class.

## D4m · Measured size-only null; A1 relabelled; calibrated intervals (1 Nov)
Reason: 9c36d64 showed S3 (size only) reaches A1 in 81–87% of replicates because ln R is correlated with ln(Δν/ν), so A1 cannot separate vortex creep from a size-only world. S3 also skips the measurement emulation, so T2's null omits measurement-induced bias in d (S4w reaches A2). Model SEs understate replicate spread by up to 3.1×. No real R–T statistic has been computed.
1. New scenario S3m: S3 passed through the measurement emulation with the source-audit method mix, as for S4-mix. Run in 06 as an evaluation scenario and as a null reference with disjoint seeds; in 08, at the real σ̂ with the real M3 size coefficient.
2. T2 rejects only if its signed root exceeds the 97.5th percentile under both S3 and S3m.
3. A1 is renamed "association consistent with the scaling, not specific to vortex creep (a size-only world produces it)". Only A2 is reported as evidence for the creep-specific Δν̇ dependence.
4. D4h's false-positive bound also applies to the A2 rate under S3m at every σ.
5. Reported intervals for b and d multiply the model SE by (SD of the estimate ÷ mean model SE) from the matching null reference at the same σ (S2 for b, S3m for d). Plain Wald intervals are shown alongside, labelled uncalibrated. Decisions are unaffected.
6. Report only: S3m outcome rates; the fraction of S1m b̂ and d̂ beyond S3m's 97.5th/2.5th percentiles.
Clarification (1 Nov): S3m uses all S4-mix emulation settings, including t_rec independent of T; S3m evaluation uses 1,000 replicates per σ, and its null reference uses D4l's 1,000.

## D4n · σ matched on fitted scatter; T2 on d̂ against null references (1 Nov)
Reason: the d̂ report after e53a4c3 showed (i) the primary model's fitted σ̂ far exceeds the generating σ (about 0.82–0.88 at σ = 0.2), so D4k item 3, which uses the real σ̂ as the generating σ, would build references for a much noisier world; (ii) the D4l signed root scales d̂ by each fit's own SE, which under S1m is 2.8× the null spread of d̂. A test on d̂ against the null distribution of d̂ is calibrated under S2, S3, S3m and S4-mix (rates ≤ 0.029) and treats d̂ as A and C already treat b̂. Chosen after seeing power on simulated data only; no real R–T statistic has been computed.
1. Replaces D4k item 3. In 08, each reference scenario (S1, S1m, S2, S3, S3m, S4-mix) is simulated at the generating σ for which that scenario's median fitted σ̂ of the primary model equals the real σ̂, found by bisection on [0, 2] with 200 replicates per step to within 0.01 in σ̂, then run with the full replicate count. If a scenario's fitted σ̂ at σ = 0 already exceeds the real σ̂, it is run at σ = 0 and this is reported.
2. T2 rejects when d̂ is below the smaller of the 2.5th percentiles of d̂ under the S3 and S3m null references at the same σ (06) or matched σ (08). The D4l signed root is reported alongside and never changes the outcome.
3. Report only: A2 is labelled "short-window-sensitive" if d̂ refitted without short-window glitches is ≥ 0.
4. 05 and 06 report every power and outcome rate against both the generating σ and the mean fitted σ̂ of the primary model.

## D4o · Literature-based exponential recoveries; replaces D14's provisional Q and τ_d (1 Nov)
Reason: D14's Q ~ log-uniform(0.005, 0.05), τ_d ~ log-uniform(1, 100) d was provisional. Published fits are pulsar-dependent: most young pulsars have τ_d of tens to hundreds of days, far above the placeholder's 10-day median, while J0537−6910 and the Crab have short τ_d and the Crab a large Q. Decided before any real R–T statistic.
Per glitch, Q and τ_d are drawn log-uniform within its group:
Each range is the envelope (minimum to maximum) of the single-component fits tabulated in the cited source; multi-component fits and the J1617−5055 event are excluded.
1. J0537−6910: exponential only if Δν ≥ 20 µHz (from dff and ν); Q ∈ [0.002, 0.04]; τ_d ∈ [4, 37] d. Source: arXiv:2512.13426 (A&A 706, A61, 2026), Table 2.
2. B0531+21: Q ∈ [0.58, 1.00]; τ_d ∈ [2.0, 21.1] d. Sources: Wong, Backer & Lyne 2001 (astro-ph/0010010), Table 3; Wang et al. 2012 (arXiv:1203.4291), Table 2.
3. All other included pulsars, from Lower et al. 2021 (arXiv:2109.07612), Table 2: if Δν/ν ≥ 1e-7, Q ∈ [0.002, 0.06] and τ_d ∈ [5.8, 488] d; if Δν/ν < 1e-7, Q ∈ [0.10, 0.15] and τ_d ∈ [54, 2340] d. Yu et al. 2013 (arXiv:1211.2035) is cited for the detection rate only.
4. Detection as D4d: AIC ≥ 10 only. The 2 × cadence rule is not reinstated.
5. Report only: the fraction of exponentials the emulation detects, per group, beside the published detection rates (J0537: 12 of 66; Yu et al.: 27 of 107).
6. Not emulated, recorded as limitations: the Crab's slow growth of the persistent Δν̇ step (≈46%, τ ≈ 320 d; Lyne et al. 2015, arXiv:1410.0886); multi-component recoveries; the J1617−5055 Q ≈ 0.98 event.
7. Add to D8 robustness: the Crab left out (J0537−6910 is already D8 item 1).
Recorded limitations (1 Nov, report only): the short-window label fires on 11.5–32% of model-true A2 replicates (S1m, σ = 1.2 to 0.2), so it signals lost power, not an artefact; the emulation detects exponentials in J0537 at 0.19–0.27 of glitches across σ (0.27 at σ = 0.2) against 12/66 published, including on glitches given none.

## D4p · Pre-freeze gaps (1 Nov)
Reason: the open-items list after 7e5b134 found cases 07–09 cannot handle. No real R–T statistic has been computed.
1. D4n's search runs on [0.01, 2]. If a scenario's fitted σ̂ at 0.01 exceeds the real σ̂, it runs at 0.01 and this is reported. If bisection has not converged after 30 steps, it uses the midpoint of the final bracket and reports it.
2. Confirmed readings for 08: emulated scenarios are noise-calibrated at each candidate and at the matched σ; interval factors use each null at its own matched σ; the search uses its own seed index with the same seeds at every step.
3. 09 robustness cases report b̂, d̂, their calibrated intervals and signed roots, and whether each estimate has the primary's sign. They receive no outcome classification.
4. A1 label text: "association consistent with the scaling, not specific to vortex creep; Δν̇ term not established (T2 power at the matched σ reported)". No fixed power figure.
5. Provisional tags removed from D10 and D11 (values from D4g) and from D9 (admissible set fixed; unaudited sources are excluded, not assumed). 07's preflight refuses if any decision is still tagged provisional.
6. Recorded as limitations, not changed: the S3m/S4 t_rec log-SD 0.3 and the short-window spans (40 and 125 d) are not taken from tabulated values.
Clarification (1 Nov), confirmed reading for 09: robustness intervals reuse the primary's interval factors from 08 and are labelled "approximate (primary factors)".

## D5 · Baselines
M0 (pulsar-specific typical wait), M3 (size only), and a per-pulsar homogeneous Poisson model
(exponential AFT) reported descriptively. *Why:* M3 is the rival that separates the model's
claim from the known size–wait link.

## D6 · Measurement uncertainty
Primary fit ignores errors in ln R. Attenuation is handled by **calibration**: S1 (model true) is
simulated with the real error bars and real measurement methods, giving the distribution of b̂ the
model would produce on data like yours. Observed b̂ is compared with that distribution (D12).
Robustness: refit with intervals where σ_lnR ≤ 0.5. *Why:* errors-in-variables for censored AFT models is
complex and fragile; simulation calibration is transparent and uses the same code.

## D7 · Vela-like / Crab-like rule
Vela-like if the pulsar's median usable Δν/ν ≥ 1×10⁻⁶; otherwise Crab-like.
**Exploratory only.** Report T1 and T2 estimates per class with intervals; no test.
*Why:* glitch size is the distinguishing property in the literature, easy to compute, and fixed
before outcomes; too few Crab-like pulsars for a confirmatory test.

## D8 · Robustness set (all reported, whatever they show)
1. Excluding J0537−6910.
2. Each Δν̇ method class separately (window, local-relaxation, global, global without exponentials, unknown).
3. Each pulsar with ≥ 10 admissible completed intervals, separately.
4. Censored intervals dropped entirely.
5. σ_lnR ≤ 0.5 only.
6. Intervals flagged for monitoring gaps excluded.
7. Δν/ν ≥ 1×10⁻⁷ only.
8. Step (permanent) Δν̇ instead of the catalogue total, where the source paper reports it.
9. Intervals whose combined start and end epoch uncertainty exceeds 10% of T excluded.

## D9 · Which Δν̇
Primary: catalogue Δν̇/ν̇ (`df1`) for admissible data only. Window, unknown and global-without-exponentials
measurements are excluded from the primary test and analysed in robustness item 2.
Secondary: step values from source papers where reported (robustness item 8).
*Why:* the catalogue value is the only one available consistently; restricting to relaxation-modelled
methods removes the artefact the emulation found. Known risk: totals that include short exponentials
dilute the signal (bias toward null); item 8 checks it. **[after audit]**

## D10 · Censoring dates
Per pulsar, in order of preference: last observation date from the catalogue maintainers; else the end of
the data span in the paper reporting the most recent glitch; when sources disagree, the earliest
defensible date. If none exists, that pulsar's last interval is dropped (logged).
*Why:* a too-early date only loses information; a too-late one is wrong. **[after audit]**

## D11 · Monitoring gaps
An interval that spans a listed gap is **censored at the gap start**. Gap list from source papers,
including J0537−6910's RXTE→NICER gap (verify dates). Robustness item 6 excludes such intervals instead.
*Why:* uses exactly what is known: no glitch from the last glitch until monitoring stopped. **[after audit]**

## D12 · Decision rules
- **Gate:** at least 40 admissible completed intervals across at least 5 pulsars. If not met: the
  method-stratified bias analysis (D14 S4 scenarios + D8 item 2) becomes the main result, and T1/T2 on
  all data are reported as exploratory.
- **Support (outcome A):** T1 and T2 both reject; observed b̂ within the central 95% of the S1-simulated
  b̂ distribution; placebo not significantly negative (D15); implied β > 0 for every pulsar with
  ≥ 5 admissible intervals.
- **Size does the work (B):** T1 rejects, T2 does not, and T2 had power ≥ 0.8 in S1.
- **Refutation (C):** T1 does not reject and T1 power in S1 ≥ 0.8 at σ = 0.4; or b̂ lies outside the
  central 95% of the S1 distribution.
- **Artefact (D):** placebo significantly negative, or support appears only in window-type data.
- **Inconclusive:** anything else, including T2 when collinearity makes it underpowered (D15).

## D13 · Inclusion
- Pulsar key: the `name` column (one excluded object has jname "X"). Missing values are coded `X`/`x`.
- Exclude these 8 objects (exact `name` strings): 4U_0142+61, 1E_1048.1-5937, CXO_J164710.2-455216,
  1RXS_J1708-4009, AX_1838.0-0655, 1E_1841-045, 1E_2259+586, M82-X2.
- Usable glitch: dff > 0 and df1 > 0. Glitches are ordered by mjd, never by `no` (J1837-0604's
  numbering is out of time order).
- Pulsars enter if they have ≥ 4 admissible completed intervals.
- Epoch uncertainty: no glitch is ever removed, because removing one would merge two intervals into a
  false long one. All glitches keep their role in ending intervals; 18 glitches in the 21 pulsars have
  mjd_err > 30 days and 6 have it unknown (`X`). Intervals touching them are flagged; robustness item 9
  excludes intervals whose combined epoch uncertainty exceeds 10% of T.
- No fractional-error cut in the primary analysis (robustness item 5 applies one).

## D14 · Simulation settings
- σ ∈ {0.2, 0.4, 0.6}; 1,000 replicates per S1–S3 scenario, 500 per S4 variant; seed 20261031.
- S1: t_rec = T (recovery completes at the next glitch). S4: t_rec independent of T, drawn from a
  lognormal with the pulsar's median wait and log-sd 0.3.
- Exponentials: Q ~ log-uniform(0.005, 0.05), τ_d ~ log-uniform(1, 100) days.
  **Provisional: replace with ranges from Gügercinoğlu et al. (2022) and the source audit.**
- Cadence and Δν̇ noise per source from the source audit; noise scaled so simulated df1 relative errors
  match the real df1_err distribution.
- Method mix in S4 = proportions in the source audit.
- **Pass criteria:** S2 false-positive rate for T1 ≤ 0.07 and for T2 ≤ 0.07 (0.05 plus Monte-Carlo
  tolerance); S3 false-positive rate for T2 ≤ 0.07; S1 b̂ median within 0.05 of the noise-free
  expectation when measurement noise is switched off. S4 results are reported, not pass/fail.

## D15 · Diagnostics
- **Placebo:** same model as M2 but with the previous wait T_(k−1) as outcome, admissible data.
  Significantly negative (p < 0.05): artefact flag, support cannot be declared. Significantly positive:
  reported as possible "longer wait, bigger glitch" physics; does not block support. Otherwise clean.
- **Missingness:** within-pulsar logistic regression of "Δν̇ missing" on ln(following wait). Report the
  odds ratio; if significant, state the selection caveat.
- **Collinearity:** within-pulsar correlation of ln Δν/ν and ln Δν̇/ν̇ (allowed before freezing). If the
  median |r| > 0.9, T2 is declared underpowered in advance and treated as exploratory.
