
What changed vs parent_experiment.py
====================================
CHANGED:
- `_origin_fixation_step`: now also returns the pre-state neighbor-fitness vector (4-tuple
  instead of 3-tuple). Needed for v3 context rule and eps-sensitivity re-derivation.
- `_neighbor_context` → `_neighbor_context_v3`: 4-category rule
  (saddle/plateau/peak/mixed) per plan's `saddle_context_rule`. saddle is now a
  fraction-based criterion (frac_strictly_beneficial > 1/(2N)), not eps-based; mixed is new.
- `_summary_stats` → `_summary_stats_v3(n_log_min=4, n_hill_min=8)`: implements the
  relaxed NaN thresholds (fix #7). `_summary_stats_v2_rule` retained alongside for the
  NaN-rate validation test.
- `_CCDF_LOG10_EDGES` constant dropped; CCDF helpers (`_ccdf_features`,
  `_ccdf_features_natural_log`) now take range + nbins as args. Range is per the plan:
  log10(time) ∈ [-2, 22] with 25 bins for stasis, ln|Δf| ∈ [-8, 4] with 20 bins for jumps.
- `_build_features` (NEW): assembles 55 features (4 scalar + 25 ccdf_stasis + 4 jump_scalar
  + 20 ccdf_jump + 1 corr + 1 autocorr). Jump-shape scalars (jump_hill, jump_ks_expon) are
  new; jump-CCDF is new. Schema is identical between the waiting-time and substitution-count
  views so the SECONDARY matched comparison is well-defined.
- `run()`: now simulates the full fixation trajectory ONCE and caches per-step neighbor fits;
  derives epochs at the primary eps AND at each sensitivity eps from the cached trajectory
  (no re-simulation). Builds both the waiting-time and substitution-count feature views at
  the primary eps. Returns per-task NaN indicators (v2 vs v3 rules) and 4-category context
  counts so the aggregator can run the chi-square and McNemar tests.
- `_stasis_jump_rankcorr`: suppress the scipy ConstantInputWarning (only when stasis_lengths
  is constant, e.g. all 1s in glassy regimes); still returns NaN as before.

UNCHANGED:
- NK and RMF landscape construction (`_nk_fitness_table`, `_nk_fitness`, `_make_rmf`,
  `_rmf_fitness`).
- Kimura origin-fixation dynamics (the math inside `_origin_fixation_step` body).
- Epoch-detection rule (relative gain vs epoch-start > eps closes the epoch); trailing
  unclosed epoch still dropped.
- Hill, KS, jump_log_stats, lag-1 autocorr formulas.

The aggregator (downstream, not in this module per task design) runs:
  - 8-way RandomForestClassifier(n_estimators=300, class_weight='balanced'), 5-fold
    stratified CV, balanced_accuracy_score with binomial 95% CI (PRIMARY).
  - Paired bootstrap on per-fold Δacc(time, subs) (SECONDARY).
  - sklearn.inspection.permutation_importance(n_repeats=30) on the full-time RF, check
    top-10 for jump-CCDF bins (TERTIARY).
  - max−min balanced accuracy across the 3 eps settings (ROBUSTNESS).
  Reviewer note: aggregator is the post-processing step that consumes `.map()` outputs;
  this module only emits the per-task observables.

Risk note: high-K NK and high-θ RMF trajectories dead-end at local optima with few epochs
(~3–13 in spot-checks), so per-trajectory feature vectors there will have many NaNs in
the high-log10 CCDF bins and in hill/KS. v3 relaxation (n≥4 for log_mean/var) keeps the
two strongest scalars defined down to n=4. The aggregator must handle NaNs (RandomForest
in sklearn ≥ 1.4 does so natively; class_weight='balanced' uses inverse-frequency weights
on the 8 bins).

PROFILE:
  sample_params: {"landscape_family_and_param": "NK_K2", "N": 20, "pop_size": 1000, "max_substitutions": 2000, "eps_jump_rel_primary": 0.005, "eps_jump_rel_sensitivity": [0.001, 0.005, 0.02], "ccdf_time_log10_range": [-2.0, 22.0], "ccdf_time_nbins": 25, "ccdf_jump_log_range": [-8.0, 4.0], "ccdf_jump_nbins": 20, "nan_thresholds": {"log_mean": 4, "log_var": 4, "hill": 8, "ks_expon_min_n": 8}, "seed": 2000}
  measured_elapsed_sec: 0.44  (NK_K2 reaches the full 2000-substitution budget — worst-case
    runtime; other cells dead-end at 5–2000 fixations and run in <0.05 s)
  scale_factor: 1×  (probe params == real task params; nothing was scaled down — the task
    is already cheap enough that the full real run fits in the 180 s probe budget)
  projected_per_task_sec: ~0.5  (worst case, NK_K2 with 2000 subs; min case ~0.01 s)
  verdict: fits  (peak per-task < 1 s; 320 tasks × 1 vCPU each, parallel; well under 60 min)
  notes: Returned dict is ~11 KB serialized (NK_K2 worst case); well under the 100 KB cap.
    Memory: trajectory cache is ~350 KB for 2000 steps × N=20 neighbor-fit vectors; fits
    1 GB easily. No re-simulation in the eps loop so adding eps sensitivity is essentially
    free relative to the simulation cost.
