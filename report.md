## TL;DR

The mechanism is **not demonstrated**. Although the per-target shift in plateau heavy-tailedness (Δσ_late) tracks the e-graph neutral-network-size ratio in rank order (Spearman ρ = 0.80, 95% bootstrap CI [0.11, 1.00]), the *direction* of the shift is mixed — only 2 of 5 targets move the predicted way (rich > minimal), and the paired one-sided Wilcoxon fails decisively (V = 7, p = 0.5938). One of three preregistered null criteria triggered ⇒ mechanism not demonstrated.

## Headline figure

![Paired σ shifts](sigma_paired_by_target.png)

## Results table

| Test | Observed | Null criterion (fails if…) | Pass? |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, n=5 targets) | V = 7, p = 0.5938; 2/5 targets in predicted direction | p > 0.05 *or* direction mixed | ✗ |
| Spearman ρ(Δσ_late, NN-size ratio) with bootstrap CI | ρ = 0.80, 95% CI [0.11, 1.00] | 95% CI crosses 0 | ✓ |
| CSN log-normal vs exponential (rich condition, per target) | 3/5 targets significant at α=0.05 (T1, T4, T5) | none significant in any target | ✓ |
| **Overall mechanism demonstration** | criterion (i) tripped | any single criterion fails | ✗ |

## Per-test figures and interpretation

### 1. Paired σ across conditions

![Paired σ](sigma_paired_by_target.png)

The hypothesis predicts every line to slope upward (rich primitives → larger σ_late). Instead three of five targets slope *down* — T1_poly_sep (Δσ = -0.3209), T3_rational (Δσ = -0.0504), T5_transcend (Δσ = -0.3322) — and only T2_coupled (+0.2662) and T4_deep_mul (+0.4541) move as predicted. The one-sided paired Wilcoxon gives V = 7, p = 0.5938. There is no aggregate shift in heavy-tailedness when the primitive set is enriched.

### 2. Δσ vs. neutral-network ratio

![Sigma shift vs NN ratio](sigma_shift_vs_nn_ratio.png)

The e-graph NN-size ratio (rich / minimal) ranges from 9.13 (T5) to 14.14 (T2). Within that range, the *rank* of Δσ tracks the rank of NN-ratio surprisingly well: Spearman ρ = 0.80, bootstrap 95% CI [0.11, 1.00] — the CI just clears zero. This is consistent with the mechanism: where the primitive enrichment inflates neutral-network sizes the most, σ does shift most positively. But on absolute scale the effect lives near zero — the regression goes from "moderately negative Δσ" at small ratios to "moderately positive Δσ" at large ratios, crossing zero in the middle. So the *correlation* claim of the hypothesis is consistent with the data; the *level* claim (rich > minimal overall) is not.

### 3. CCDF of late-phase plateau durations by target

![CCDF by target](tail_ccdf_by_target.png)

Visually, the rich and minimal CCDFs largely overlap; the censoring shoulder at duration ≈ 750 (max generations cap) is prominent on every panel. The CSN log-likelihood-ratio test (log-normal vs exponential) is significant in 3/5 targets, but only 2 of those 3 prefer log-normal (T1 R = -3.92, T5 R = -2.95); T4 actually prefers *exponential* (R = +3.21, p = 0.0013). So criterion (iii) ("not distinguishable in any target") passes — the distributions are distinguishable from exponential — but the direction of distinguishability is not uniformly toward heavy-tailed log-normal.

## Discussion

Holding the target fixed and enriching the primitive set with identity-equivalent operators inflates e-graph-measured neutral-network sizes by 6×–17× across the five targets, yet does not produce a corresponding inflation in plateau-duration σ_late. The rank-order correlation (ρ = 0.80) is suggestive — it's exactly what the mechanism predicts conditional on a real effect existing — but the magnitudes that produce that rank are tiny (|Δσ| ≤ 0.45) and the mean shift is near zero. The most direct read is **M2 from the plan: heavy plateau tails are largely generic GP dynamics, not a readout of landscape neutrality at this scale.** Possible alternative reads: (a) the chosen NN-proxy (e-graph mean class size on 2000 random samples) measures the wrong slice of the landscape — actual GP search occupies a non-uniform region; (b) censoring at gen=1500 truncates the very tail that distinguishes the distributions; (c) the scale-down (pop=80, 40 reps) lost discriminating power in σ even though tail-event counts cleared the CSN bar. Follow-up would be: scale up to tk_9bdd20b5's full pop and generations on just T2/T4 (the two targets that moved as predicted) to see whether their positive Δσ holds.

## Caveats

- **Heavy censoring** at duration = 750 generations (= half of max_generations, the late-phase window). On several panels >25% of pooled late durations are right-censored at 750; σ estimates treat these as observed and therefore *underestimate* the true heavy-tailedness of both conditions equally — but the bias may not cancel in the difference Δσ.
- **n = 5 targets** is the unit of analysis for both Wilcoxon and Spearman. The Spearman CI [0.11, 1.00] is wide and the lower bound is fragile to a single target swap; do not over-read.
- **Neutral-proxy** is a uniform-random-sample estimate of mean equivalence-class size, not the actual measure on search-occupied trees. The 6×–17× ratio it reports is plausibly the right *direction* but unknown calibration.
- **CSN test is rich-condition-only** per plan; minimal-condition tail behavior was not tested for distinguishability.
- T4_deep_mul prefers exponential over log-normal in the rich condition (R = +3.21), which is at odds with the mechanism's qualitative prediction even though it passes "distinguishable from exponential."
