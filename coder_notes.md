
Changes vs parent code (explicit diff):
  CHANGED:
    * Dropped _safe_div and removed "/" from BINOPS entirely (PI: "drop _safe_div").
    * primitive_set() rewritten: now three conditions M / S / C instead of
      {minimal, rich}. S adds {neg, id} + duplicate "1.0" terminal
      (size-only). C adds {sin, cos, square, sqrtp} (genuine new portals).
    * UNOPS split into UNOPS_S and UNOPS_C; ALL_UN merges for evaluator.
    * run_gp(): tournament k=2 (was 4) for true neutral drift; pop_size and
      max_gen are now params-driven; added visited_samples collection
      (log-spaced sub-sampling of current-best trees, serialized as strings).
    * Renamed plateau_durations -> plateau_events; now emits
      (start_gen, end_gen, duration, fitness_level, observed_flag) so
      downstream code can do censored MLE log-normal σ and CSN xmin search,
      and band by fitness level (fitness_banded_post_τ definition of late phase).
    * Added egraph_banded_measure(): the new banded e-graph proxy task,
      computing per-band mean class size (semantic-fingerprint membership-
      weighted) AND out-degree portal density (P(1-node mutation -> better
      band)). This is the dissociated static landscape probe.
    * run() now dispatches on kind in {"gp_run","egraph_banded"} per the
      plan's two-task-kind design; returns the full per-gen best trace +
      plateau event log + visited trees, leaving τ_target calibration,
      σ-fitting, Wilcoxon/Spearman/CSN tests to post-processing across tasks.
  KEPT (verbatim):
    * gen_tree, tree_size, tree_nodes, replace_node, evaluate (extended to
      new unaries via ALL_UN dispatch), fitness, crossover, mutate.
    * Noise injection rationale (additive y-noise creates MSE floor).
    * Elitism: best-of-population preserved each generation
      ("preserve_only" matches plan).
    * px = 0.85 crossover rate, 0.15 secondary mutation, max_depth=4 init.

Null-criteria observables this task emits (post-processing handles aggregation):
  * plateau_events with observed/censored flag and fitness_level
    -> pooled censored MLE log-normal σ per (target, condition);
    -> CSN powerlaw.Fit with xmin search + distribution_compare on pooled
       late-phase durations per cell.
  * full best_hist (subsampled) -> per-target τ_target calibration across
    conditions, fitness banding for late-phase selection.
  * egraph_banded results (size, portal) -> Δsize_S and Δconnectivity_C
    for the two Spearman tests.

PROFILE:
  sample_params: {"kind":"gp_run","target":{"id":"T4_deep_mul","expr":"(x0 + x1) * (x2 + x0) * (x1 + x2)","n_vars":3},"condition":"C","seed":7,"pop_size":200,"max_generations":150,"size_cap":30,"noise_sigma":0.05}
  measured_elapsed_sec: 2.27
  scale_factor: 10  # max_gen 150 -> 1500 is 10x; pop_size already full at 200; T4/C is heaviest target/cond combo (largest trees, transcendentals)
  projected_per_task_sec: ~23
  verdict: fits
  notes: Initial small-pop probe (pop=50, max_gen=75, T1/C) was 0.30s; scaling
    to full (pop=200, max_gen=1500) by ×4×20 = ×80 projects ~24s. Heavy-target
    probe (pop=200, max_gen=150, T4/C) at 2.27s × 10 = ~23s confirms. With 120
    reps × 6 targets × 3 conds = 2160 GP tasks at ~25s each = ~54k task-seconds,
    each task is a separate parallel container so per-container wall ≈ 25s << 60min
    cap. The 18 egraph_banded tasks (n_samples=4000) are similar or smaller cost.
    Comfortable margin; no need to reduce max_generations or n_reps.
