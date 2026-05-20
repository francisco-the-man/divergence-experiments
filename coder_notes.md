
PROFILE:
  sample_params: {"seed": 0, "landscape_family_and_param": "NK_K2", "N": 20, "pop_size": 1000, "max_substitutions": 400, "stasis_detector": "saddle_aware"}
  measured_elapsed_sec: 0.081
  scale_factor: 1  # profile ran a FULL task (max_substitutions=400, n_replicates=1 — that's already one task)
  projected_per_task_sec: ~0.1 (NK_K2). NK_K16 builds a 20×2^17 ≈ 21 MB contribs table once; per-step cost is the same O(N) neighbor enumeration so still well under 1s. RMF caches eta lazily; trajectory visits ≤401 unique genotypes so cache stays small.
  verdict: fits
  notes:
    - Each task = one trajectory on one landscape; plan specifies 320 such tasks (8 cells × 40 reps) for parallel map. The aggregator/report-writer collects the per-task summary_stats and runs the classifier + Spearman tests.
    - "stasis_lengths_subs" mostly equal 1 because beneficial fixations dominate at SSWM with Ne=1000 (i.e. once a beneficial neighbor exists, it fixes before drift through many neutrals). The richer observable is "stasis_lengths_time" (the waiting time in dimensionless units of 1/(mu*Ne)) — this is what carries the landscape-topology signature per Sibani/Koonin. I return both; the analysis layer can use either, but I recommend "summary_stats_time" as the primary feature set.
    - eps_jump_rel = 0.005 (relative to |f|) defines a "beneficial fixation". I made it scale-relative so it transfers between NK (f ∈ [0,1]) and RMF (f roughly in [-N*theta, +noise]).
    - Saddle-aware stratification implemented per lit-reviewer's concern: `stasis_contexts` and `summary_stats_saddle*` give the saddle-only feature set for the stratification_drop test.
    - Right-censored trailing epoch is dropped (standard handling for waiting-time data).
    - mu_per_genome from the plan does not appear in dynamics: under origin-fixation Gillespie, only relative fixation rates matter for which neighbor fixes, and the time unit is 1/(mu*Ne) which is set by external normalization. The classifier sees only distribution shape, which is mu-invariant. Documenting this for the reviewer.
    - Memory for NK_K=16: 20*2^17*8 bytes ≈ 21 MB, comfortable in 1 GB.
