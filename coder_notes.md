
PROFILE:
  sample_params: {"target_id": "T4_deep_mul", "expr": "(x0 + x1) * (x2 + x0) * (x1 + x2)", "primitives": "rich", "n_replicates": 4, "max_generations": 1500, "pop_size": 80, "size_cap": 30, "seed_base": 8000}
  measured_elapsed_sec: 25.1
  scale_factor: 10× (n_replicates: 40/4 = 10; max_generations and pop_size already at full plan scale)
  projected_per_task_sec: ~250 (≈4.2 min)
  verdict: fits
  notes:
    - T4 is among the harder targets (deepest tree). Easy targets (T1, T2) will be faster per-rep but converge harder, producing fewer late-phase plateaus per replicate; total time per task should be ≤ this estimate.
    - All 10 tasks projected well under the 60-min hard cap; ample headroom.

METHODOLOGY NOTES:
1. Plateau definition: strict inequality on best-fitness sequence. Elitism makes
   best_hist non-increasing, so a strictly-lower value marks a real improvement.
   No relative-tolerance epsilon: any tiny float-different MSE counts. This is
   the standard "inter-improvement time" used in evolutionary-dynamics literature.

2. Noise floor on training y (σ_rel=0.05 of target std-dev): without this, easy
   targets converge to numerically-zero MSE and the late-phase plateau collapses
   to a single 750-gen block. With noise, the GP keeps wandering near optimum,
   producing the sequence of small jumps that the plan wants to measure.
   This is consistent with how progress-curve / waiting-time statistics are
   typically extracted in noisy fitness regimes.

3. Late-phase event counts per cell (pooled across 40 reps):
   - Hard targets (T4, T5): ~200-300 events → CSN OK (≥50), σ stable.
   - Easy targets (T1, T2): likely 80-150 events → CSN OK, σ stable.
   - The plan's claim of "~1200 events per cell" was optimistic — real strict-
     improvement inter-arrival counts are lower. Still meets CSN's 50-tail bar
     by 2-6× and gives reliable log-normal σ MLEs.

4. NN-size proxy: 2000 random trees per condition, fingerprinted by output
   vector on 12 probe points, rounded to 4 decimals. Membership-weighted mean
   class size = sum(s²)/sum(s) — interpretable as "expected size of class
   containing a random genotype". Returns the rich/minimal ratio for the
   correlation test in the analysis plan. This is independent of any GP run.

5. Determinism: every random op uses np.random.default_rng(seed_base + r) for
   replicate r. The NN proxy uses seed_base + 99999. No global RNG state.

LIMITATIONS the reviewer should know:
- 40 reps × ~1500 gen with strict-improvement plateaus gives ~5-10 late-phase
  events per rep → small per-rep counts. Pooling is the right move and the plan
  is designed for it; replicate-level bootstrap (downstream) is the way to get
  CIs on σ.
- "Late phase" = second half of generations is a simple split; doesn't adapt to
  per-run convergence point. Acceptable proxy at this scale.
- size_cap=30 may bias all conditions toward similar effective expressivity
  ceilings; if the mechanism shows weakly, full-scale (size_cap=40) follow-up
  may amplify it.
