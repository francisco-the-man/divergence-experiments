Revisions per reviewer:
(1) Training data (X,y) seeded by TARGET-only seed (TARGET_DATA_SEED) so rich and minimal cells for the same target see identical data.
(2) Fitness checkpoints are var_y * fixed_fractions [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001] — target-determined, condition-independent. Both cells of a target probe the neutral-network proxy at identical MSE thresholds.
(3) var_y, checkpoint_fractions, checkpoints all returned for downstream sanity-checking.
GP-run RNG still uses seed_base+r (per-condition trajectories differ, as intended).
Neutral-proxy sampling RNG uses seed_base+99999 (sampling differs across conditions, but the grid it's evaluated against is identical) — this seemed acceptable since the observable is class-count at matched thresholds, not the specific trees sampled.