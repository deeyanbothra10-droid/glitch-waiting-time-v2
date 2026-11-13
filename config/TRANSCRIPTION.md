# Transcription check: config/decisions.yaml ← DECISIONS.md

Each key in config/decisions.yaml, its value, and the sentence of DECISIONS.md it comes from.
Labels in `code font` are names for the code; the value is the sentence. Every group also has
`provisional` (false everywhere: D4p item 5 removed the D9, D10 and D11 tags; the D14.exponentials tag was removed on 1 Nov because D4o replaces it).

| Key | Value | Source sentence in DECISIONS.md |
|---|---|---|
| meta.status | proposed | "Status: PROPOSED by Claude, 31 Oct 2026 … Not yet frozen." |
| meta.decisions_md_sha256 | 2fd9…0bed | (checksum of DECISIONS.md at transcription; not a decision) |
| D1.primary_question | verbatim | D1 "**Primary:** Within individual glitching pulsars, …" |
| D1.secondary_question | verbatim | D1 "**Secondary (exploratory):** Is the relation stronger …" |
| D1.secondary_exploratory | true | D1 "Secondary (exploratory)" |
| D2.data | admissible | D2 heading "Hypotheses (tested on admissible data)" |
| D2.T1.model / parameter / expected_sign | M2 / b / positive | D2 "T1 (scaling): in M2, b > 0." |
| D2.T1.null_value | 0 | D2 "Null: b = 0." |
| D2.T2.model / parameter / expected_sign | M4 / d / negative | D2 "T2 (the Δν̇ term): in M4, d < 0." |
| D2.T2.null_value | 0 | D2 "Null: d = 0." |
| D2.model_prediction.b / d | 1 / −1 | D2 "The model predicts b = 1 and d = −1" |
| D2.consistency_judged_against | S1_simulation | D2 "judged against the S1 simulation (see D6, D12), not against the raw value 1" |
| D3.family / censoring | lognormal_aft / right | D3 "Lognormal accelerated-failure-time model, right-censored" |
| D3.intercepts / sigma | per_pulsar / shared | D3 "one intercept per pulsar, one shared σ" |
| D3.T1_model / T2_model | M2 / M4 | D3 "T1 uses M2 (b ln R + c_p); T2 uses M4 (a ln Δν/ν + d ln Δν̇/ν̇ + c_p)." |
| D4.test | likelihood_ratio | D4 "Likelihood-ratio tests" |
| D4.T1_comparison | full M2, reduced M0 | D4 "T1 (M2 vs M0, H0 b = 0)" |
| D4.T2_comparison | full M4, reduced M3 | D4 "T2 (M4 vs M3, H0 d = 0)" |
| D4.sidedness | two-sided | D4 "two-sided" |
| D4.familywise_alpha / correction | 0.05 / bonferroni | D4 "Family-wise α = 0.05 with Bonferroni" |
| D4.per_test_alpha | 0.025 | D4 "each test at α = 0.025" |
| D4.secondary_analyses | spearman, permutation, loo_cv, aic, implied_beta, class_comparison | D4 "Everything else (Spearman, permutation, LOO-CV, AIC, implied β, class comparison) is secondary" |
| D4.secondary_reporting | estimates_and_intervals_no_significance | D4 "reported with estimates and intervals, no significance claims" |
| D4a.permutation.strata / permuted_variable | pulsar / ln_R | D4a "within each pulsar, permute ln R" |
| D4a.permutation.intervals | all | D4a "across ALL of that pulsar's intervals, completed and censored together; each interval keeps its own T and censoring flag, only ln R moves" |
| D4a.permutation.statistic | M2_bhat_full_censored_likelihood | D4a "Statistic: the M2 slope estimate b̂ from the full censored likelihood." |
| D4a.permutation.sidedness / p_value | two-sided / one_plus_count_abs_ge_over_one_plus_N | D4a "Two-sided p = (1 + #{\|b̂_perm\| ≥ \|b̂_obs\|}) / (1 + N)" |
| D4a.permutation.n_permutations / seed | 10000 / 20261031 | D4a "with N = 10,000 permutations and seed 20261031" |
| D4a.permutation.placebo_uses_same_scheme | true | D4a "The placebo permutation (D15) uses the same scheme with the previous wait as outcome." |
| D4a.spearman.intervals | completed_only | D4a "Spearman: completed intervals only" |
| D4a.spearman.per_pulsar_min_completed_admissible | 4 | D4a "computed per pulsar for pulsars with ≥ 4 completed admissible intervals" |
| D4a.spearman.pooling / pooling_weight | fisher_z_weighted_mean / n_minus_3 | D4a "pooled as a Fisher-z weighted mean with weights (n_p − 3)" |
| D4a.spearman.p_values / label | none / completed intervals only | D4a "Descriptive only: no p-values, and always labelled "completed intervals only"" |
| D4a.loo_cv.unit | interval | D4a "leave out one interval at a time, refit" |
| D4a.loo_cv.score_completed / score_censored | log_f_T / log_S_c | D4a "score the held-out interval by log f(T) if completed or log S(c) if censored" |
| D4a.loo_cv.summary | mean_over_all_intervals_per_model | D4a "report the mean over all intervals for each model" |
| D4b.interval_estimates.level | 0.95 | D4b "All reported intervals are 95%" |
| D4b.interval_estimates.model_parameters | wald_inverse_hessian | D4b "Wald intervals from the inverse Hessian for model parameters" |
| D4b.interval_estimates.simulated_quantities | percentile | D4b "percentile intervals for simulated quantities" |
| D4b.interval_estimates.pooled_spearman | fisher_z_pm_z_over_sqrt_sum_n_minus_3 | D4b "the pooled Spearman uses Fisher z ± 1.96/sqrt(Σ(n_p − 3))" |
| D4b.interval_estimates.spearman_abs_r_one | flag_and_exclude_from_pooled_and_list | D4b "a pulsar with \|r_s\| = 1 is reported and flagged, and excluded from the pooled Fisher-z mean. The excluded pulsars are listed." |
| D4b.simulation.typical_wait | obs_span_over_n_glitches | D4b "W_p = (obs_end − first-glitch mjd) / (number of glitches)" |
| D4b.simulation.covariate_centring / intercept | within_pulsar / ln_W_p | D4b "covariates are centred within pulsar, and c_p = ln W_p" |
| D4b.simulation.S1 | ln_R, slope 1 | D4b "S1: ln T = (ln R − mean_p ln R) + c_p + σε." |
| D4b.simulation.S2 | none | D4b "S2: ln T = c_p + σε." |
| D4b.simulation.S3 | ln_dnu_frac, slope 1 | D4b "S3: ln T = a·(ln dnu_frac − mean_p ln dnu_frac) + c_p + σε, with a = 1." |
| D4b.simulation.S4_t_rec | lognormal, median W_p, log_sd 0.3, independent of T, uses_waits false | D4b "t_rec ~ lognormal(median W_p, log-sd 0.3), drawn independently of T. It never uses real or simulated waits." |
| D4b.simulation.replicates_per / replicates / seed | scenario_and_sigma / 1000, 1000, 1000, 500 / 20261031 | D4b "1,000 per scenario per σ value for S1–S3, and 500 per S4 variant per σ value, with seed 20261031" |
| D4b.simulation.validation_lr_tests | all_replicates | D4b "the LR tests run on every replicate" |
| D4b.simulation.validation_permutation | primary, placebo / 999 / 200 | D4b "Permutation tests (primary and placebo) run with 999 permutations on the first 200 replicates of each scenario × σ" |
| D4b.simulation.validation_loo_cv_first_n_replicates | 100 | D4b "LOO-CV runs on the first 100" |
| D4b.simulation.real_run_n_permutations | 10000 | D4b "The real run uses 10,000 permutations, as in D4a." |
| D4b.simulation.method_class.global_exp_fitted_no | global-noexp | D4b "an S4g0 (global without exponentials) case is df1_method = global with exp_fitted = no" |
| D4b.simulation.method_class.global_exp_fitted_sometimes | per_glitch_from_audit_notes_else_unknown | D4b "If exp_fitted = sometimes, classify per glitch when the source audit notes say which glitches had exponentials; otherwise the class is unknown." |
| D4b.simulation.method_class.global_exp_fitted_yes / admissible_global_requires_exp_fitted | global / yes | D4b "For admissibility (D9), global requires exp_fitted = yes." |
| D4b.diagnostics.collinearity.correlation / over | pearson / usable_glitches | D4b "Pearson correlation of ln(Δν/ν) and ln(Δν̇/ν̇) over usable glitches" |
| D4b.diagnostics.collinearity.min_admissible_usable_glitches_per_pulsar | 4 | D4b "for pulsars with ≥ 4 admissible usable glitches" |
| D4b.diagnostics.collinearity.rule_summary | unweighted_median_abs_r | D4b "The D15 rule uses the unweighted median of \|r\| across those pulsars." |
| D4b.diagnostics.collinearity.also_report | pooled_within_pulsar_demeaned_r | D4b "Also report the pooled within-pulsar (demeaned) correlation." |
| D4b.diagnostics.placebo.outcome / rows / previous_wait_always_completed | previous_wait / intervals_whose_start_glitch_has_predecessor / true | D4b "outcome = the previous wait T_(k−1), for every interval whose starting glitch has a predecessor (that wait is always completed)" |
| D4b.diagnostics.placebo.test / sidedness | D4a_permutation / two-sided | D4b "Test: the D4a permutation scheme, two-sided." |
| D4b.diagnostics.placebo.alpha / significantly_negative / significantly_positive | 0.05 / p < α and b̂ < 0 / p < α and b̂ > 0 | D4b ""Significantly negative" means p < 0.05 and b̂ < 0; "significantly positive" means p < 0.05 and b̂ > 0." |
| D4b.diagnostics.missingness.method / strata | conditional_logit / pulsar | D4b "conditional (pulsar-stratified) logistic regression" |
| D4b.diagnostics.missingness.response / predictor | df1_missing / ln_following_wait | D4b "of "df1 missing" on ln(following wait)" |
| D4b.diagnostics.missingness.rows / censored_following_wait | included_pulsars_glitches_with_completed_following_wait / exclude_and_count | D4b "over glitches in the included pulsars whose following wait is completed. Glitches whose following wait is censored are excluded and counted." |
| D4b.diagnostics.missingness.sidedness / alpha | two-sided / 0.05 | D4b "Two-sided α = 0.05." |
| D4b.diagnostics.missingness.report | odds_ratio, interval | D4b "Report the odds ratio with a 95% interval." |
| D4b.diagnostics.missingness.descriptive_also | df1_le_0 | D4b "Also report descriptively the same model with "df1 ≤ 0" as the outcome." |
| D4b.simulation.intercept | ln_W_p (superseded) | D4b "c_p = ln W_p"; replaced by D4c.glitch_counts.intercept |
| D4c.s4_expectations.S4l_forward / S4l_forward_max_upper_95 | not_positive / 0.02 | D4c "S4l (local-relaxation), forward: not positive (upper end of the 95% interval ≤ +0.02)." |
| D4c.s4_expectations.S4l_bias_reporting | finding_works_against_T1 | D4c "The bias works against T1 and is reported as a finding." |
| D4c.s4_expectations.S4w / S4g / S4g0 | unchanged | D4c "S4w, S4g and S4g0: expectations unchanged." |
| D4c.layers.statistics_only | S1, S2, S3 | D4c "S1, S2 and S3 test the statistics code only. They do not use the measurement emulation, and t_rec is not used." |
| D4c.layers.statistics_measured_dnudot | true_times_exp_normal_sd_real_relative_error | D4c "Measured Δν̇ = true Δν̇ × exp(e), where e is normal with SD = the real relative error (df1_err/df1) of the glitch that was resampled." |
| D4c.layers.emulated / S1m_t_rec | S1m, S4 / equals_T | D4c "S1m = S1 passed through the measurement emulation, with t_rec = T." ; "S4 variants: as before, through the emulation" |
| D4c.layers.S1_reference_for_D6_D12 | S1m | D4c "Its b̂ distribution is the reference for D6 and D12 ("the S1-simulated b̂ distribution" means S1m from now on)." |
| D4c.layers.S4_t_rec | independent_of_T | D4c "with t_rec independent of T" |
| D4c.mixtures.method_assignment | resampled_glitch_method | D4c "each simulated glitch takes the method of the real glitch it was resampled from" |
| D4c.mixtures.primary_analysis_glitches | admissible_only | D4c "For the primary analysis, only admissible glitches are used." |
| D4c.mixtures.unknown_emulated_as / unknown_sensitivity_emulated_as | window / global | D4c ""unknown" is emulated as window (worst case), with a sensitivity run that emulates it as global" |
| D4c.emulation.cadence / cadence_until_audit_days | per_source_from_audit / 14 | D4c "Cadence: per source from the audit. Until the audit exists, 14 days." |
| D4c.emulation.noise_calibration | per_pulsar_bisection_on_median_relative_error | D4c "per pulsar, choose the ν̇ sample noise (by bisection) so that the median relative error of the emulated Δν̇ matches the median real df1_err/df1 of that pulsar's usable glitches" |
| D4c.emulation.noise_min_error_values_per_pulsar / noise_fallback | 3 / global_median | D4c "Pulsars with fewer than 3 error values use the global median." |
| D4c.emulation.pre_glitch_monitoring_days | 365 | D4c "Pre-glitch monitoring: 365 days before the first glitch." |
| D4c.emulation.min_samples_after / min_samples_before_window | 4 / 4 | D4c "a glitch needs at least 4 samples after it (and, for window, at least 4 before it)" |
| D4c.emulation.unmeasurable_glitch | missing_unusable_still_ends_previous_interval | D4c "Otherwise its Δν̇ is missing: it is unusable but still ends the previous interval." |
| D4c.emulation.exp_min_tau_over_cadence | 2 | D4c "exponentials with τ_d < 2 × cadence are not included in any fit. Their effect stays in the data" |
| D4c.glitch_counts.intercept | ln_W_p_minus_half_sigma_squared | D4c "Use c_p = ln W_p − σ²/2, so the MEAN synthetic wait equals W_p" |
| D4d.s4_role.pass_fail_tests | none | D4d "Remove every pass/fail test on S4 bias directions or sizes." (restores D14 "S4 results are reported, not pass/fail") |
| D4d.s4_role.outputs | results, emulation_findings_log | D4d "They are computed, saved to results/ and logs/emulation_findings.md, and reported." |
| D4d.s4_role.emulator_correctness_tests | four tests (a)–(d) | D4d "(a) window with no exponentials gives J(1 − T/(2·t_rec)) for T ≤ t_rec; (b) window with one exponential and constant step gives Δν̇_p + (Δν_d/T)(1 − e^(−T/τ_d)); (c) global and local with every exponential fitted and t_rec = T are exact; (d) an exponential below the detection rule is never fitted." |
| D4d.exponential_detection.unit / criterion / min_aic_drop | each_fitted_window / aic_drop / 10 | D4d "In each fitted window, include an exponential term only if including it lowers AIC by at least 10" |
| D4d.exponential_detection.tau_d_known_to_fit | true | D4d "(τ_d known to the fit)" |
| D4d.calibrated_rules.reference_scenarios | S1m, S4-mix | D4d "S1m (model true) and S4-mix (bias only, t_rec independent of T)" |
| D4d.calibrated_rules.reference_inputs / primary_measurement_mix | same_pulsars_counts_measurement_mix_noise_as_real / admissible_only | D4d "from the same pulsars, counts, measurement mix (admissible only for the primary analysis) and noise as the real data" |
| D4d.calibrated_rules.replicates / sigma | 1000 / estimated_from_real_primary_fit | D4d "with 1,000 replicates each at the σ estimated from the real primary fit" |
| D4d.calibrated_rules.support_requires | b̂ above S4-mix 97.5th percentile; placebo b̂ not below S1m placebo 2.5th percentile | D4d "Support (A) additionally requires: observed b̂ above the 97.5th percentile of the S4-mix b̂ distribution …, and observed placebo b̂ not below the 2.5th percentile of the S1m placebo distribution." |
| D4d.calibrated_rules.artefact_if | b̂ inside S4-mix central 95%; or placebo b̂ below S1m placebo 2.5th percentile | D4d "Artefact (D): observed b̂ inside the central 95% of the S4-mix distribution, or observed placebo b̂ below the 2.5th percentile of the S1m placebo distribution." |
| D4d.calibrated_rules.raw_placebo_permutation_p | reported_not_deciding | D4d "The raw placebo permutation p-value is still reported, but no longer decides anything on its own." |
| D4d.answers.dnu_noise | times_exp_normal_sd_real_relative_error_all_layers | D4d "Δν/ν noise: yes, in every layer: multiply by exp(e), with e normal and SD equal to the real relative error dff_err/dff of the resampled glitch." |
| D4d.answers.calibration_sequences_per_step / calibration_tolerance_fraction_of_target / calibration_record_achieved | 50 / 0.05 / true | D4d "Calibration: 50 simulated sequences per pulsar per bisection step; tolerance 5% of the target; record the achieved value." |
| D4d.answers.cadence_within_pulsar / pre_glitch_baseline_cadence | interval_uses_starting_glitch_source_cadence / previous_interval | D4d "each interval uses the cadence of the source of the glitch that starts it; the pre-glitch baseline uses the previous interval's cadence" |
| D4d.answers.missing_error_fill / missing_error_min_values_per_pulsar / missing_error_fallback / missing_error_flagged | pulsar_median_relative_error / 3 / global_median / true | D4d "a usable glitch without df1_err (or dff_err) uses its pulsar's median relative error (the global median if the pulsar has fewer than 3 values), and is flagged" |
| D4e.readings_confirmed.global_exponential_detection | same_per_window_aic_decision | D4e "global uses the same per-window AIC decision to choose which exponentials enter the global fit" |
| D4e.readings_confirmed.pre_glitch_baseline_cadence | first_glitch_cadence | D4e "the 365-day baseline before the first glitch uses the first glitch's cadence" |
| D4e.calibration.target_statistic / target_reference | median_reported_formal_relative_error / median_real_df1_err_over_df1 | D4e "match the median of the fit's REPORTED (formal) relative error of Δν̇ to the median real df1_err/df1, not the scatter of measured against true values" |
| D4e.calibration.scenario / per | S1m / pulsar_and_sigma | D4e "calibrate once per pulsar under S1m (t_rec = T, so no misfit), at each σ" |
| D4e.calibration.noise_used_in | all_emulated_scenarios | D4e "use that same ν̇ noise level in every scenario (S1m and all S4 variants)" |
| D4e.calibration.report | achieved_per_pulsar, pulsars_outside_tolerance | D4e "Report, per pulsar, the achieved value and any pulsar outside the 5% tolerance." |
| D4e.replicate_inclusion.rule | D13_min_admissible_completed_intervals_per_pulsar | D4e "apply the D13 inclusion rule (≥ 4 admissible completed intervals) to each simulated replicate, exactly as for real data" |
| D4e.replicate_inclusion.failing_pulsars / replicate / report | dropped_from_replicate_fit / kept / pulsars_dropped_per_scenario | D4e "Pulsars that fail it are dropped from that replicate's fit; the replicate itself is kept. Report how many pulsars were dropped per scenario." |
| D4e.replicate_inclusion.discard_replicate_if_fewer_pulsars_than / report_discarded | 2 / true | D4e "A replicate is discarded only if fewer than 2 pulsars remain, and those are counted." |
| D4e.limitations.marginal_detection_errors_understated | stated_limitation_no_code_change | D4e "the reported Δν̇ errors understate the true scatter by 2–4× ... It affects robustness item 5 (σ_lnR ≤ 0.5) and is stated as a limitation; no code change." |
| D4f.short_window.class / definition | short-window / pre_post_local_solutions_fixed_spans_differenced_no_transients | D4f "New class `short-window`: Δν̇ from separate pre-/post-glitch local solutions over spans fixed by the analyst (not ending at the next glitch), differenced at the glitch epoch, transients not fitted." |
| D4f.short_window.sources.elsk11 (source_type, w_days, w_sensitivity_days) | espinoza, 40, 150 | D4f "Applies to Espinoza+2011 (post-glitch span ideally 30–50 d)" ; "w = 40 d (sensitivity 150 d) for Espinoza-type" |
| D4f.short_window.sources.wmp+00 (source_type, w_days, w_sensitivity_days) | wang, 125, 300 | D4f "and Wang+2000 (~100–300 d)" ; "w = 125 d (sensitivity 300 d) for Wang-type" |
| D4f.short_window.emulation | constant_level_w_before_and_after_truncated_at_neighbouring_glitch_no_exponentials | D4f "Emulation S4sw: constant ν̇ level fitted over w days before and w days after the glitch (each truncated at the neighbouring glitch), no exponentials; Δν̇ = post level − pre level." |
| D4f.admission.statistic / max_upper_95 / every_sigma / main_and_sensitivity_w / fixture | S4sw_forward_slope_upper_95 / 0.02 / true / true / synthetic | D4f "short-window enters the primary analysis only if, in the synthetic fixture, the S4sw forward slope's 95% upper bound is ≤ +0.02 at every σ, for both the main and the sensitivity w" |
| D4f.admission.if_passed / otherwise / record | primary_analysis / robustness_item_2_only / emulation_findings_log | D4f "Otherwise it is analysed in robustness item 2 only. The result is recorded in emulation_findings.md either way." |
| D4f.exp_fitted_when_detected_counts_as | yes | D4f "exp_fitted = "when-detected" counts as "yes" for admissibility" |
| D4f.our_measurement | global, assumed, admitted, robustness run excluding | D4f ""Our measurement" rows are global (assumed, following Basu et al. 2022); admitted; a robustness run excludes them; to be confirmed with the catalogue maintainers." |
| D4f.robustness_without_J0537 | always reported; exploratory if < 40 intervals | D4f "the primary analysis without J0537−6910 is always reported, labelled exploratory if it has fewer than 40 intervals" |
| D4f.unaudited_sources | unknown | D4f "Sources not yet audited remain `unknown` (excluded from the primary analysis)." |
| D4g.gaps.J0537-6910.gap_start_mjd | 55926.9 | D4g "RXTE monitoring ended MJD 55926.9 (31 Dec 2011; Antonopoulou et al. 2018)" |
| D4g.gaps.J0537-6910.gap_end_mjd | 57982 | D4g "NICER monitoring began MJD 57982 (17 Aug 2017; Ho et al. 2026, Table 2)" |
| D4g.gaps.J0537-6910.censored_interval_start_mjd | 55819 | D4g "The interval starting at MJD 55819 is censored at 55926.9 (D11)." |
| D4g.obs_end_documented.J0537-6910 | 60827 | D4g "J0537−6910: 60827 (NICER data to 2025 Jun 1, Ho et al. 2026)." |
| D4g.obs_end_documented.J1023-5746 | 59269 | D4g "J1023−5746: 59269 (Fermi data end, Gügercinoğlu et al. 2022)." |
| D4g.jodrell_bank_basu2022 | 58483 if last glitch before 58483; J0631+1036, J0205+6449, J0729-1448, J1814-1744 | D4g "Jodrell Bank pulsars whose last glitch is before 58483: 58483 (31 Dec 2018, end of the Basu et al. 2022 sample): J0631+1036, J0205+6449, J0729−1448, J1814−1744." |
| D4g.all_other_pulsars | last_glitch_epoch; dropped_and_logged_never_zero_censoring | D4g "All other pulsars: obs_end = their last glitch epoch, i.e. no censoring information; that last interval is dropped and logged, never entered with c = 0." |
| D4g.robustness_item_10 | Jodrell Bank pulsars: obs_end = max(last glitch, 60891) | D4g "New robustness item 10: Jodrell Bank pulsars use obs_end = max(last glitch, 60891), i.e. the latest catalogue measurement (MJD 61071) minus 180 days for reporting lag." |
| D4g.robustness_item_10.jodrell_bank_rule.source_keys / primary_reference_prefix | elsk11, bsa+21, bkj+19, our measurement / Shaw+ | D4g "a pulsar counts as Jodrell Bank-monitored if at least one of its catalogue glitches has a Jodrell Bank primary source (Espinoza+2011, Basu+2021, Basu+2019, any Shaw+ entry, or any "Our measurement" spelling)" |
| D4g.robustness_item_10.jodrell_bank_rule.included_pulsars_computed_31_oct | 15 pulsars | D4g "Resulting included (D13) pulsars, computed 31 Oct: B0531+21, … J2229+6114 (15 of 21; 144 of 218 catalogue pulsars)." (a check on the rule, not a separate decision) |
| D4h.intercept | ln W_p − σ²/2 − ln m_p; m_p = pool mean of exp(slope × centred covariate); 1 without covariate | D4h "c_p = ln W_p − σ²/2 − ln m_p, where m_p is the mean over the pulsar's resampled pool of exp(slope × centred covariate) for that scenario (S1/S1m: covariate ln R, slope 1; S3: covariate ln dnu_frac, slope a; S2/S4: m_p = 1)" |
| D4h.false_positive_bound | α 0.025, multiplier 3, n = replicates; S2 T1, S2 T2, S3 T2 | D4h "a test passes if its false-positive rate ≤ 0.025 + 3·sqrt(0.025·0.975/n), where n is the number of replicates ... Applies to T1 and T2 under S2, and to T2 under S3." |
| D4h.robustness_item_9 | quadrature; 0.1; missing → exclude | D4h "combined epoch uncertainty = sqrt(mjd_err_start² + mjd_err_end²). An interval is excluded if that exceeds 0.1·T, or if either error is missing (X)." |
| D4h.robustness_item_8 | not_run; no step-only values | D4h "Robustness item 8: not run, because no step-only Δν̇ values are recorded. Reported as "not run" with that reason." |
| D4h.confirmed | S1m 1000 per σ; D4d percentile placebo; raw permutation p reported only; implied β over all admissible intervals | D4h "S1m uses 1,000 replicates per σ; Support (A) uses D4d's percentile placebo condition, with the raw permutation p reported only; the implied-β check counts all admissible intervals, completed and censored." |
| D4i.reason_recorded_before_freeze | true | D4i "Reason (recorded before freezing): 05 showed realistic power of about 0.78 for T1 and about 0.34 for T2 at σ = 0.4 ... This revision was made after seeing those simulated values, and no real R–T statistic." |
| D4i.A1_scaling_supported | T1 rejects; b̂ > S4-mix 97.5th; b̂ within S1m central 95%; placebo as D4d | D4i "A1, scaling supported: T1 rejects; observed b̂ above the 97.5th percentile of S4-mix b̂; observed b̂ within the central 95% of S1m b̂; placebo condition as in D4d." |
| D4i.A2_dnudot_supported | A1 and T2 rejects with d < 0 | D4i "A2, Δν̇ term supported: A1 holds and T2 rejects (d < 0)." |
| D4i.T2_not_rejected | label "not established (T2 power ≈ 0.34)"; not outcome B | D4i "If T2 does not reject, the Δν̇ term is "not established (T2 power ≈ 0.34)". That is NOT outcome B" |
| D4i.B_requires_T2_power_at_least | 0.8 | D4i "because B requires T2 power ≥ 0.8, which this design does not have" |
| D4i.C_refutation | b̂ < S1m 2.5th; or T1 not rejected and T1 power ≥ 0.8 at the real σ̂ | D4i "C, refutation: observed b̂ below the 2.5th percentile of S1m b̂; or T1 does not reject AND the T1 power at the σ estimated from the real primary fit is ≥ 0.8." |
| D4i.D_artefact | as_D4d | D4i "D, artefact: as in D4d." |
| D4i.implied_beta | reported, not deciding; 95% intervals; by D7 class | D4i "Implied β: reported per pulsar with 95% intervals, as a physical-consistency check discussed by class (Vela-like vs Crab-like). It no longer decides any outcome." |
| D4i.simulation_outcome_rates | D4i_rules | D4i "In simulations, outcome rates are computed with these rules." |
| D4i.missingness.judged_against | S1m_distribution | D4i "D15 missingness: judged against its S1m distribution, like the placebo." |
| D4j.reason_recorded_before_freeze | true | D4j "Reason: under D4d, D fired whenever b̂ lay inside the central 95% of S4-mix ... Recorded before freezing; no real R–T statistic has been computed." |
| D4j.order | D_placebo, C, A, D_bias_range, inconclusive | D4j "Outcomes are assigned in this order; the first that applies wins" |
| D4j.D_placebo | placebo b̂ < S1m placebo 2.5th | D4j "1. D (artefact signature): observed placebo b̂ below the 2.5th percentile of the S1m placebo distribution." |
| D4j.C_refutation | b̂ < S1m 2.5th; or T1 not rejected and T1 power ≥ 0.8 | D4j "2. C (refutation): observed b̂ below the 2.5th percentile of S1m b̂; or T1 does not reject AND T1 power at the fitted σ is ≥ 0.8." |
| D4j.A_scaling_supported | T1 rejects; b̂ > S4-mix 97.5th; b̂ within S1m central 95% | D4j "3. A (scaling supported): T1 rejects, b̂ above the 97.5th percentile of S4-mix b̂, and b̂ within the central 95% of S1m b̂." |
| D4j.A_scaling_supported.within_A | A2: T2 rejects (d < 0); B: T2 not rejected and T2 power ≥ 0.8; A1 otherwise | D4j "A2 if T2 also rejects (d < 0); B if T2 does not reject and T2 power at the fitted σ is ≥ 0.8; otherwise A1 (Δν̇ term not established)." |
| D4j.D_bias_range | T1 rejects and b̂ ≤ S4-mix 97.5th | D4j "4. D (signal within bias range): T1 rejects but b̂ is at or below the 97.5th percentile of S4-mix b̂." |
| D4j.inconclusive | anything_else | D4j "5. Inconclusive: anything else." |
| D4j.never_change_outcome | implied_beta, D15_missingness | D4j "Implied β and the D15 missingness check are reported alongside and never change the outcome." |
| D4k.failed_replicates | excluded from rates; listed with reason; > 1% per job → FAIL until explained in logs/emulation_findings.md | D4k "1. A failed replicate (fit failure, timeout, unexpected error) is excluded from its job's rates and listed with its reason. If more than 1% of replicates in any job fail, that job is FAIL until the cause is explained in logs/emulation_findings.md." |
| D4k.validation_sigma_grid_additional | 0.8, 1.2 | D4k "2. Add σ = 0.8 and 1.2 to the validation grid (keep existing values)." |
| D4k.real_reference_sigma / real_reference_replicates_and_seeds | fitted σ̂ of the primary model / as in 06 | D4k "3. For the real classification, every percentile and power in D4j comes from simulations run at the real fitted σ̂ of the primary model (D3), never interpolated from the grid; replicate count and seeding as in 06." |
| D4k.missingness_not_estimable | reported as not estimable; S1m reference from estimable replicates; count reported | D4k "4. A missingness check that is not estimable (separation) is reported as such; the S1m reference uses estimable replicates only, with the count reported." |
| D4k.real_data_fit_failure | stops_run | D4k "5. In 07–09, any fit failure on real data stops the run; it is never skipped." |
| D4l.reason_recorded_before_freeze | true | D4l "Reason: 06 (3921bdf) showed (i) ... (iv) ... No real R–T statistic has been computed." |
| D4l.directional_tests (statistic, directions) | signed root of LR; b > 0 for T1, d < 0 for T2 | D4l "T1 and T2 are directional. Statistic: the signed root of the LR statistic, signed so that positive means the predicted direction (b > 0 for T1, d < 0 for T2)." |
| D4l.directional_tests (reject, nulls, replicates, seeds) | > 97.5th percentile under S2 (T1) / S3 (T2); 1,000 replicates; disjoint seeds | D4l "Each rejects when its statistic exceeds the 97.5th percentile of the same statistic under its null scenario (S2 for T1, S3 for T2), simulated at the same σ with 1,000 replicates and seeds disjoint from the evaluation replicates." |
| D4l.directional_tests (real run) | nulls at σ̂ in 08; S3 size coefficient = real M3 estimate | D4l "In 08 the null references are run at the real σ̂, with S3's size coefficient set to the real M3 estimate." |
| D4l.directional_tests.powers_use_these_tests | true | D4l "Powers in D4j/D4k use these tests." |
| D4l.C_refutation | b̂ < S1m 2.5th; power clause removed | D4l "2. C is: observed b̂ below the 2.5th percentile of S1m b̂. The power clause is removed." |
| D4l.A_slope_condition | b̂ > S4-mix 97.5th, ≥ S1m 2.5th, ≤ S1 97.5th; label above S1m 97.5th | D4l "3. A's slope condition becomes: b̂ above the 97.5th percentile of S4-mix b̂, at or above the 2.5th percentile of S1m b̂, and at or below the 97.5th percentile of S1 b̂. Report b̂'s percentile within S1m; if it is above S1m's 97.5th percentile, add the label "less attenuated than emulated" (report only)." |
| D4l.report_only | three diagnostics | D4l "4. Report only, never outcome-changing: for S3 T2 rejections under the old test, the fraction with d̂ < 0; for S1, S3 and S1m, the SD of b̂ and d̂ across replicates against the mean model SE; for S1m, the slope of measured on true ln(Δν̇/ν̇) per method class." |
| D4m.reason_recorded_before_freeze | true | D4m "Reason: 9c36d64 showed S3 (size only) reaches A1 in 81–87% ... No real R–T statistic has been computed." |
| D4m.S3m (base, measurement, method_mix) | S3; emulated; source-audit mix | D4m "1. New scenario S3m: S3 passed through the measurement emulation with the source-audit method mix, as for S4-mix." |
| D4m.S3m.t_rec | as_S4_mix | D4m clarification "S3m uses all S4-mix emulation settings, including t_rec independent of T" |
| D4m.S3m.evaluation_replicates | as_S3 | D4m clarification "S3m evaluation uses 1,000 replicates per σ, and its null reference uses D4l's 1,000." |
| D4m.S3m (roles, seeds, real_run) | evaluation and null reference; disjoint seeds; 08 at σ̂ with the real M3 coefficient | D4m "Run in 06 as an evaluation scenario and as a null reference with disjoint seeds; in 08, at the real σ̂ with the real M3 size coefficient." |
| D4m.T2_rejects_if_exceeds_null_percentile_under | S3, S3m | D4m "2. T2 rejects only if its signed root exceeds the 97.5th percentile under both S3 and S3m." |
| D4m.A1_label / creep_specific_evidence | renamed label / A2 only | D4m "3. A1 is renamed "association consistent with the scaling, not specific to vortex creep (a size-only world produces it)". Only A2 is reported as evidence for the creep-specific Δν̇ dependence." |
| D4m.false_positive_bound_also_applies_to | A2 rate under S3m, every σ | D4m "4. D4h's false-positive bound also applies to the A2 rate under S3m at every σ." |
| D4m.calibrated_intervals | SE × (SD ÷ mean SE) from S2 (b) and S3m (d) at the same σ; plain Wald shown, labelled uncalibrated; decisions unaffected | D4m "5. Reported intervals for b and d multiply the model SE by (SD of the estimate ÷ mean model SE) from the matching null reference at the same σ (S2 for b, S3m for d). Plain Wald intervals are shown alongside, labelled uncalibrated. Decisions are unaffected." |
| D4m.report_only | S3m outcome rates; S1m b̂ above S3m 97.5th; S1m d̂ below S3m 2.5th | D4m "6. Report only: S3m outcome rates; the fraction of S1m b̂ and d̂ beyond S3m's 97.5th/2.5th percentiles." |
| D4n.reason_recorded_before_freeze / chosen_after | true / simulated_power_only | D4n "Reason: the d̂ report after e53a4c3 showed … Chosen after seeing power on simulated data only; no real R–T statistic has been computed." |
| D4n.matched_sigma (replaces, scenarios, statistic, target) | D4k item 3; S1, S1m, S2, S3, S3m, S4-mix; median fitted σ̂ of the primary model = real σ̂ | D4n "1. Replaces D4k item 3. In 08, each reference scenario (S1, S1m, S2, S3, S3m, S4-mix) is simulated at the generating σ for which that scenario's median fitted σ̂ of the primary model equals the real σ̂," |
| D4n.matched_sigma (search, search_range, replicates_per_step, tolerance_in_sigma_hat, then_replicates) | bisection; [0, 2]; 200; 0.01; full count | D4n "found by bisection on [0, 2] with 200 replicates per step to within 0.01 in σ̂, then run with the full replicate count." |
| D4n.matched_sigma.if_fitted_sigma_hat_at_zero_exceeds_target | run_at_zero_and_report | D4n "If a scenario's fitted σ̂ at σ = 0 already exceeds the real σ̂, it is run at σ = 0 and this is reported." |
| D4n.T2 (statistic, percentile, null_scenarios, combine, null σ) | d̂ below the smaller 2.5th percentile under S3 and S3m; same σ in 06, matched σ in 08 | D4n "2. T2 rejects when d̂ is below the smaller of the 2.5th percentiles of d̂ under the S3 and S3m null references at the same σ (06) or matched σ (08)." |
| D4n.T2.signed_root | reported, never changes the outcome | D4n "The D4l signed root is reported alongside and never changes the outcome." |
| D4n.short_window_label | A2; "short-window-sensitive"; d̂ without short-window ≥ 0; report only | D4n "3. Report only: A2 is labelled "short-window-sensitive" if d̂ refitted without short-window glitches is ≥ 0." |
| D4n.report_against_sigma / report_against_sigma_in | generating σ and mean fitted σ̂ of the primary model / 05, 06 | D4n "4. 05 and 06 report every power and outcome rate against both the generating σ and the mean fitted σ̂ of the primary model." |
| D4o.reason_recorded_before_freeze / replaces / draw | true / D14.exponentials / log-uniform within group | D4o "Reason: D14's Q ~ log-uniform(0.005, 0.05), τ_d ~ log-uniform(1, 100) d was provisional. … Decided before any real R–T statistic." and "Per glitch, Q and τ_d are drawn log-uniform within its group:" |
| D4o.ranges / excluded_from_envelopes | envelope of single-component tabulated fits / multi-component fits, J1617−5055 event | D4o "Each range is the envelope (minimum to maximum) of the single-component fits tabulated in the cited source; multi-component fits and the J1617−5055 event are excluded." |
| D4o.groups.J0537-6910 | Δν ≥ 20 µHz; Q [0.002, 0.04]; τ_d [4, 37] d | D4o "1. J0537−6910: exponential only if Δν ≥ 20 µHz (from dff and ν); Q ∈ [0.002, 0.04]; τ_d ∈ [4, 37] d. Source: arXiv:2512.13426 (A&A 706, A61, 2026), Table 2." |
| D4o.groups.crab | Q [0.58, 1.00]; τ_d [2.0, 21.1] d | D4o "2. B0531+21: Q ∈ [0.58, 1.00]; τ_d ∈ [2.0, 21.1] d. Sources: Wong, Backer & Lyne 2001 (astro-ph/0010010), Table 3; Wang et al. 2012 (arXiv:1203.4291), Table 2." |
| D4o.groups.others | split 1e-7; ≥: Q [0.002, 0.06], τ_d [5.8, 488] d; <: Q [0.10, 0.15], τ_d [54, 2340] d | D4o "3. All other included pulsars, from Lower et al. 2021 (arXiv:2109.07612), Table 2: if Δν/ν ≥ 1e-7, Q ∈ [0.002, 0.06] and τ_d ∈ [5.8, 488] d; if Δν/ν < 1e-7, Q ∈ [0.10, 0.15] and τ_d ∈ [54, 2340] d." |
| D4o.detection | D4d AIC only | D4o "4. Detection as D4d: AIC ≥ 10 only. The 2 × cadence rule is not reinstated." |
| D4o.report_only | detection fraction per group; 12/66; Yu 27/107 (detection rate only) | D4o "5. Report only: the fraction of exponentials the emulation detects, per group, beside the published detection rates (J0537: 12 of 66; Yu et al.: 27 of 107)." and item 3 "Yu et al. 2013 (arXiv:1211.2035) is cited for the detection rate only." |
| D4o.limitations_not_emulated | three limitations | D4o "6. Not emulated, recorded as limitations: the Crab's slow growth of the persistent Δν̇ step (≈46%, τ ≈ 320 d; Lyne et al. 2015, arXiv:1410.0886); multi-component recoveries; the J1617−5055 Q ≈ 0.98 event." |
| D4o.D8_robustness_additions | the Crab left out | D4o "7. Add to D8 robustness: the Crab left out (J0537−6910 is already D8 item 1)." |
| D4o.limitations_recorded_report_only | label 0.115–0.32 of S1m A2; J0537 detection 0.19–0.27 of glitches | D4o "Recorded limitations (1 Nov, report only): the short-window label fires on 11.5–32% of model-true A2 replicates (S1m, σ = 1.2 to 0.2), so it signals lost power, not an artefact; the emulation detects exponentials in J0537 at 0.19–0.27 of glitches across σ (0.27 at σ = 0.2) against 12/66 published, including on glitches given none." |
| D4p.reason_recorded_before_freeze | true | D4p "Reason: the open-items list after 7e5b134 found cases 07–09 cannot handle. No real R–T statistic has been computed." |
| D4p.matched_sigma_search | [0.01, 2]; run at 0.01 if exceeded, reported; 30 steps, then the midpoint of the final bracket, reported | D4p "1. D4n's search runs on [0.01, 2]. If a scenario's fitted σ̂ at 0.01 exceeds the real σ̂, it runs at 0.01 and this is reported. If bisection has not converged after 30 steps, it uses the midpoint of the final bracket and reports it." |
| D4p.confirmed_readings_08 | calibration at each candidate and matched σ; factors from each null at its matched σ; own seed index, same seeds every step | D4p "2. Confirmed readings for 08: emulated scenarios are noise-calibrated at each candidate and at the matched σ; interval factors use each null at its own matched σ; the search uses its own seed index with the same seeds at every step." |
| D4p.robustness_09 | b̂, d̂, calibrated intervals, signed roots, same sign as primary; no classification | D4p "3. 09 robustness cases report b̂, d̂, their calibrated intervals and signed roots, and whether each estimate has the primary's sign. They receive no outcome classification." |
| D4p.robustness_09.interval_factors / interval_label | primary's factors from 08 / "approximate (primary factors)" | D4p "Clarification (1 Nov), confirmed reading for 09: robustness intervals reuse the primary's interval factors from 08 and are labelled "approximate (primary factors)"." |
| D4p.A1_label / A1_label_fixed_power_figure | new text / false | D4p "4. A1 label text: "association consistent with the scaling, not specific to vortex creep; Δν̇ term not established (T2 power at the matched σ reported)". No fixed power figure." |
| D4p.provisional_tags_removed / D9_unaudited_sources / preflight_07_refuses_if_any_provisional | D9, D10, D11 / excluded, not assumed / true | D4p "5. Provisional tags removed from D10 and D11 (values from D4g) and from D9 (admissible set fixed; unaudited sources are excluded, not assumed). 07's preflight refuses if any decision is still tagged provisional." |
| D4p.limitations_recorded | t_rec log-SD 0.3; short-window spans 40 and 125 d | D4p "6. Recorded as limitations, not changed: the S3m/S4 t_rec log-SD 0.3 and the short-window spans (40 and 125 d) are not taken from tabulated values." |
| D5.baselines | M0, M3, poisson_exponential_aft | D5 "M0 …, M3 (size only), and a per-pulsar homogeneous Poisson model (exponential AFT)" |
| D5.poisson_scope / poisson_reporting | per_pulsar / descriptive | D5 "a per-pulsar homogeneous Poisson model … reported descriptively" |
| D6.primary_errors_in_lnR | ignore | D6 "Primary fit ignores errors in ln R." |
| D6.attenuation_handling / calibration_scenario | simulation_calibration / S1 | D6 "Attenuation is handled by calibration: S1 (model true) is simulated …" |
| D6.calibration_inputs | real_error_bars, real_measurement_methods | D6 "with the real error bars and real measurement methods" |
| D6.robustness_max_sigma_lnR | 0.5 | D6 "Robustness: refit with intervals where σ_lnR ≤ 0.5." |
| D7.statistic / vela_like_min_dnu_frac | median_usable_dnu_frac / 1.0e-6 | D7 "Vela-like if the pulsar's median usable Δν/ν ≥ 1×10⁻⁶; otherwise Crab-like." |
| D7.exploratory_only | true | D7 "Exploratory only." |
| D7.report_per_class / test | T1_estimate, T2_estimate / none | D7 "Report T1 and T2 estimates per class with intervals; no test." |
| D8.report_all | true | D8 heading "(all reported, whatever they show)" |
| D8.exclude_pulsar | J0537-6910 | D8 item 1 "Excluding J0537−6910." (written with an ASCII hyphen, as in the `name` column) |
| D8.method_classes_separately | window, local-relaxation, global, global-noexp, unknown | D8 item 2 "Each Δν̇ method class separately (window, local-relaxation, global, global without exponentials, unknown)." |
| D8.per_pulsar_min_admissible_completed | 10 | D8 item 3 "Each pulsar with ≥ 10 admissible completed intervals, separately." |
| D8.drop_censored | true | D8 item 4 "Censored intervals dropped entirely." |
| D8.max_sigma_lnR | 0.5 | D8 item 5 "σ_lnR ≤ 0.5 only." |
| D8.exclude_gap_flagged | true | D8 item 6 "Intervals flagged for monitoring gaps excluded." |
| D8.min_dnu_frac | 1.0e-7 | D8 item 7 "Δν/ν ≥ 1×10⁻⁷ only." |
| D8.step_dnudot_where_reported | true | D8 item 8 "Step (permanent) Δν̇ instead of the catalogue total, where the source paper reports it." |
| D8.max_epoch_uncertainty_fraction_of_T | 0.10 | D8 item 9 "Intervals whose combined start and end epoch uncertainty exceeds 10% of T excluded." |
| D9.primary_dnudot | catalogue_df1 | D9 "Primary: catalogue Δν̇/ν̇ (`df1`) for admissible data only." |
| D9.admissible_methods | local-relaxation, global | Notation "admissible data = … (`local-relaxation`, or `global` with exponentials included)" |
| D9.excluded_from_primary | window, unknown, global-noexp | D9 "Window, unknown and global-without-exponentials measurements are excluded from the primary test" |
| D9.secondary_dnudot | source_paper_step | D9 "Secondary: step values from source papers where reported (robustness item 8)." |
| D9.provisional | false | D9 "**[after audit]**"; removed by D4p "5. Provisional tags removed from D10 and D11 (values from D4g) and from D9 …" |
| D10.preference_order | maintainer_last_observation, most_recent_glitch_paper_data_span_end | D10 "last observation date from the catalogue maintainers; else the end of the data span in the paper reporting the most recent glitch" |
| D10.on_disagreement | earliest_defensible | D10 "when sources disagree, the earliest defensible date" |
| D10.if_none | drop_last_interval_and_log | D10 "If none exists, that pulsar's last interval is dropped (logged)." |
| D10.provisional | false | D10 "**[after audit]**"; removed by D4p item 5 |
| D11.interval_spanning_gap | censor_at_gap_start | D11 "An interval that spans a listed gap is censored at the gap start." |
| D11.gap_list_source / gaps_to_verify | source_papers / J0537-6910 RXTE->NICER | D11 "Gap list from source papers, including J0537−6910's RXTE→NICER gap (verify dates)." |
| D11.robustness_alternative | exclude | D11 "Robustness item 6 excludes such intervals instead." |
| D11.provisional | false | D11 "**[after audit]**"; removed by D4p item 5 |
| D12.gate.min_admissible_completed_intervals / min_pulsars | 40 / 5 | D12 "Gate: at least 40 admissible completed intervals across at least 5 pulsars." |
| D12.gate.if_not_met | verbatim | D12 "If not met: the method-stratified bias analysis … reported as exploratory." |
| D12.outcomes.A_support.rule | verbatim | D12 "Support (outcome A): …" |
| D12.outcomes.A_support.s1_central_fraction | 0.95 | A: "within the central 95% of the S1-simulated b̂ distribution" |
| D12.outcomes.A_support.implied_beta_min / …_min_admissible_intervals | 0 / 5 | A: "implied β > 0 for every pulsar with ≥ 5 admissible intervals" |
| D12.outcomes.B_size_does_the_work.rule / T2_min_power_S1 | verbatim / 0.8 | D12 "Size does the work (B): … T2 had power ≥ 0.8 in S1." |
| D12.outcomes.C_refutation.rule | verbatim | D12 "Refutation (C): …" |
| D12.outcomes.C_refutation.T1_min_power_S1 / power_sigma | 0.8 / 0.4 | C: "T1 power in S1 ≥ 0.8 at σ = 0.4" |
| D12.outcomes.C_refutation.s1_central_fraction | 0.95 | C: "b̂ lies outside the central 95% of the S1 distribution" |
| D12.outcomes.D_artefact.rule | verbatim | D12 "Artefact (D): placebo significantly negative, or support appears only in window-type data." |
| D12.outcomes.inconclusive.rule | verbatim | D12 "Inconclusive: anything else, including T2 when collinearity makes it underpowered (D15)." |
| D13.pulsar_key | name | D13 "Pulsar key: the `name` column" |
| D13.missing_codes | X, x | D13 "Missing values are coded `X`/`x`." |
| D13.excluded_objects | 8 names | D13 "Exclude these 8 objects (exact `name` strings): …" |
| D13.usable_min_dff_exclusive / usable_min_df1_exclusive | 0 / 0 | D13 "Usable glitch: dff > 0 and df1 > 0." |
| D13.order_by | mjd | D13 "Glitches are ordered by mjd, never by `no`" |
| D13.min_admissible_completed_intervals_per_pulsar | 4 | D13 "Pulsars enter if they have ≥ 4 admissible completed intervals." |
| D13.remove_glitches_for_epoch_uncertainty | false | D13 "no glitch is ever removed" |
| D13.flag_mjd_err_above_days / flag_mjd_err_unknown | 30 / true | D13 "have mjd_err > 30 days and 6 have it unknown (`X`). Intervals touching them are flagged" |
| D13.expected_counts_21_pulsars | 18 / 6 | D13 "18 glitches in the 21 pulsars have mjd_err > 30 days and 6 have it unknown" |
| D13.primary_fractional_error_cut | none | D13 "No fractional-error cut in the primary analysis" |
| D14.sigma_grid | 0.2, 0.4, 0.6 | D14 "σ ∈ {0.2, 0.4, 0.6}" |
| D14.replicates.S1/S2/S3 / S4_per_variant | 1000 / 500 | D14 "1,000 replicates per S1–S3 scenario, 500 per S4 variant" |
| D14.seed | 20261031 | D14 "seed 20261031" |
| D14.S1_t_rec | equals_T | D14 "S1: t_rec = T (recovery completes at the next glitch)." |
| D14.S4_t_rec | independent, lognormal, median = pulsar_median_wait, log_sd 0.3 | D14 "S4: t_rec independent of T, drawn from a lognormal with the pulsar's median wait and log-sd 0.3." |
| D14.exponentials.Q | log-uniform 0.005–0.05 | D14 "Q ~ log-uniform(0.005, 0.05)" |
| D14.exponentials.tau_d_days | log-uniform 1–100 | D14 "τ_d ~ log-uniform(1, 100) days" |
| D14.exponentials.provisional | false | D14 "Provisional: replace with ranges from Gügercinoğlu et al. (2022) and the source audit."; replaced by D4o "Literature-based exponential recoveries; replaces D14's provisional Q and τ_d" (0269544), so the tag was removed (1 Nov) |
| D14.cadence_and_noise / noise_scaling | per_source_from_source_audit / match_real_df1_relative_errors | D14 "Cadence and Δν̇ noise per source from the source audit; noise scaled so simulated df1 relative errors match the real df1_err distribution." |
| D14.S4_method_mix | source_audit_proportions | D14 "Method mix in S4 = proportions in the source audit." |
| D14.pass_criteria.S2_T1_… / S2_T2_… | 0.07 / 0.07 | D14 "S2 false-positive rate for T1 ≤ 0.07 and for T2 ≤ 0.07" |
| D14.pass_criteria.S3_T2_max_false_positive_rate | 0.07 | D14 "S3 false-positive rate for T2 ≤ 0.07" |
| D14.pass_criteria.S1_bhat_median_max_abs_deviation_noise_free | 0.05 | D14 "S1 b̂ median within 0.05 of the noise-free expectation when measurement noise is switched off" |
| D14.pass_criteria.S4 | reported_not_pass_fail | D14 "S4 results are reported, not pass/fail." |
| D15.placebo.model / outcome / data | M2 / previous_wait / admissible | D15 "same model as M2 but with the previous wait T_(k−1) as outcome, admissible data" |
| D15.placebo.p_threshold | 0.05 | D15 "Significantly negative (p < 0.05)" |
| D15.placebo.if_significantly_negative | artefact_flag_support_blocked | D15 "artefact flag, support cannot be declared" |
| D15.placebo.if_significantly_positive | report_does_not_block_support | D15 "reported as possible … physics; does not block support" |
| D15.placebo.otherwise | clean | D15 "Otherwise clean." |
| D15.missingness.method / response / predictor | within_pulsar_logistic_regression / dnudot_missing / ln_following_wait | D15 "within-pulsar logistic regression of "Δν̇ missing" on ln(following wait)" |
| D15.missingness.report / if_significant | odds_ratio / state_selection_caveat | D15 "Report the odds ratio; if significant, state the selection caveat." |
| D15.collinearity.measure / variables | within_pulsar_correlation / ln_dnu_frac, ln_dnudot_frac | D15 "within-pulsar correlation of ln Δν/ν and ln Δν̇/ν̇" |
| D15.collinearity.allowed_before_freeze | true | D15 "(allowed before freezing)" |
| D15.collinearity.summary / underpowered_if_above / if_underpowered | median_abs_r / 0.9 / T2_exploratory | D15 "If the median \|r\| > 0.9, T2 is declared underpowered in advance and treated as exploratory." |
