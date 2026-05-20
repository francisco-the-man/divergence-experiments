
Changes vs parent code (per Case-B follow-up discipline):
- KEPT IDENTICAL: _nk_fitness_table, _nk_fitness, _make_rmf, _rmf_fitness, _origin_fixation_step, _hill_estimator, _ks_to_exponential, _summary_stats, _parse_cell, _neighbor_context, the saddle-aware stasis-detector logic in run(), trailing-epoch dropping policy, eps_jump_rel default 0.005.
- CHANGED: default max_substitutions 400→2000 (parent's 400 default replaced; task params override anyway).
- ADDED: _ccdf_features (20-bin empirical CCDF on fixed log10 edges -2..6), _jump_log_stats, _stasis_jump_rankcorr, _lag1_autocorr_log. These are pure feature-extraction helpers; they do not alter simulation dynamics.
- ADDED to run(): tracking of jump_mags, waiting_times, entry_contexts during the epoch loop (parent only tracked stasis lengths). Returned in result so the aggregator can compute stratified CCDFs and the jump-vs-stasis boxplot.
- ADDED to result dict: 'features' now contains 28 keys (4 original + 20 CCDF + 2 jump + 1 rankcorr + 1 autocorr); 'feature_groups' tags each for the permutation-importance plot; 'stasis_lengths_stratified' for the saddle-vs-plateau confound plot; 'jump_magnitudes', 'waiting_times', 'entry_contexts' as raw arrays for the cross-task plots in analysis_plan.

Seed offset (+1000) is applied by the task generator (seed_start=1000 in the plan); the code consumes whatever seed is passed.

The aggregator (downstream) is expected to: stack `features` dicts across all 320 tasks into a DataFrame, build the 8-way `param_bin_label` target, drop rows with NaN features, run the three feature-set ablation (4-scalar / +CCDF / full) with stratified 5-fold CV RandomForest, and compute balanced-accuracy CIs + permutation importance. Per-task observables here are exactly what those tests reference.

Observed during profiling: trajectories often terminate at local peaks well before max_subs=2000 (NK_K2 with seed 1000 produced 5 epochs in 0.34s; the parent reported similar small epoch counts on rugged cells). This is inherent to the origin-fixation dynamics + saddle-aware threshold and matches parent behavior — so the 5× max_subs bump is a ceiling, not a guarantee. If the aggregator finds the null still holds, the diagnostic interpretation in the plan ("within-family ruggedness not recoverable from single-trajectory stasis distributions") is well-supported.

PROFILE:
  sample_params: {"seed": 1000, "landscape_family_and_param": "NK_K2", "N": 20, "pop_size": 1000, "max_substitutions": 2000, "stasis_detector": "saddle_aware"}
  measured_elapsed_sec: 0.344
  scale_factor: 1   # already ran the full max_substitutions=2000; the loop terminates naturally at peaks. The worst-case cell (a smooth landscape where trajectories DON'T peak-trap) would scale linearly with steps actually taken. Even at the full 2000 step budget the per-step cost is ~N+ε fitness evaluations + occasional context check (also N), so worst case ≈ 2000 × 40 fitness evals ≈ negligible. NK fitness is O(N) per call, RMF is O(N) + dict hit. A maximally-active 2000-step trajectory would take roughly 2000 × 0.344 / 5 ≈ 140 s (linear extrapolation from 5 effective steps to 2000).
  projected_per_task_sec: ~150 (loose upper bound; typical task <1s based on probe)
  verdict: fits
  notes: 320 parallel tasks, each well under the 3600s per-task ceiling. Result dict is small (<10KB per task even for 2000-step trajectories: 2000 ints + 2000 floats + 2000 floats + 2000 short strings ≈ 80KB worst case, still under the 100KB guideline).
