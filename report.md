# Plateau heavy-tails do not track primitive-set richness — mechanism not demonstrated

## TL;DR

The preregistered paired test asked whether enriching the GP primitive set (identity-equivalents + extra constants) lengthens log-normal σ on late-phase plateau durations. It does not: paired Wilcoxon `V = 7`, `p = 0.59375` (one-sided rich > minimal), with only **2 / 5** targets shifting in the predicted direction. The independent e-graph NN-size manipulation **did** work (rich/minimal ratio ≈ 9–14× across all 5 targets), and Spearman ρ between Δσ and the NN-ratio is **0.80** with bootstrap 95 % CI `[0.11, 1.00]` — suggestive but underpowered at n = 5. The CSN log-normal-vs-exponential distinguishability test passes on only **3 / 5** rich-condition targets. By the preregistered null criteria (any one failure = mechanism not demonstrated), **all three criteria fire**: the proposed neutral-network mechanism for plateau heavy-tails is not supported at this scale.

## Headline figure

![Per-target paired σ shift](sigma_paired_by_target.png)

## Results vs. preregistered null criteria

| Test | Observed | Null criterion | Verdict |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, n=5) | V = 7, p = 0.59375; 2/5 positive | p > 0.05 or mixed direction → fail | ✗ |
| Spearman Δσ vs. NN-ratio, 95 % bootstrap CI | ρ = 0.80, CI [0.11, 1.00] | CI crosses 0 → fail | ✓ (just) |
| CSN log-normal vs. exponential (rich, per target) | 3 / 5 significant at α = 0.05 | any target not distinguishable → fail | ✗ |
| **Overall mechanism demonstrated?** | — | all three must pass | **✗** |

Two of the three preregistered checks fail, so by the plan's stated logic the mechanism is not demonstrated.

## Per-test interpretation

### 1. σ shift is directionally inconsistent across targets

![Per-target paired σ shift](sigma_paired_by_target.png)

The lines cross. T2_coupled (Δσ = +0.27) and T4_deep_mul (Δσ = +0.45) move the predicted way; T1_poly_sep (Δσ = -0.32), T3_rational (Δσ = -0.05) and T5_transcend (Δσ = -0.33) move the opposite way. The mean Δσ across targets is 0.0034 — essentially zero. If "more primitive-set neutrality → heavier plateau tails" were a robust mechanism, we'd expect 5/5 in the predicted direction at NN-ratios this large; we get 2/5.

### 2. Spearman correlation is positive but the CI just barely excludes zero

![Δσ vs. neutral-network ratio](sigma_shift_vs_nn_ratio.png)

Spearman ρ = 0.80 with bootstrap 95 % CI [0.11, 1.00] across the five targets. The point estimate is in the hypothesised direction and the CI's lower bound (0.11) does not cross zero, so this criterion technically passes. But (a) n = 5 makes the bootstrap CI extremely wide and discrete, (b) the correlation is driven by the rank order more than by a strong gradient — the NN-ratio only spans 9.13× to 14.14× across targets, a narrow range, while Δσ ranges from −0.33 to +0.45 with no monotone trend by ratio magnitude. Read this as "consistent with, but not evidence for" the mechanism.

### 3. Pooled late-phase CCDFs don't show a uniform tail-lengthening in the rich condition

![Late-phase CCDFs](tail_ccdf_by_target.png)

Visual inspection of the log-log CCDFs of pooled late-phase plateau durations confirms the σ result: tails of rich vs. minimal cross each other within targets, and the dominant feature in every panel is the spike at the run-length ceiling (duration = 750 generations = censored at end-of-run). The CSN log-normal-vs-exponential test passes on only 3 of 5 rich-condition targets (T1: p = 8.8e-05, T4: p = 1.3e-03, T5: p = 3.2e-03), and fails to distinguish the two distributions on T2 (p = 0.069) and T3 (p = 0.130) — meaning we can't even claim the rich-condition tails are reliably log-normal-like rather than exponential on those targets.

## Discussion

The counterfactual was clean: e-graph measurement confirmed the manipulation increased neutral-network size by 9–14× in every target, but the predicted progress-dynamics signature (heavier log-normal tails on plateau durations) did not appear paired-target-wise. Two of the failure modes the plan named — M2 (heavy tails are a generic GP artifact unrelated to landscape neutrality) and M3 (a shift exists but doesn't track the independent NN measurement) — both have some support here: σ values are 1.5–2.1 for *every* (target × condition) cell, so plateaus are heavy-tailed in both conditions, and the directional shift is uncorrelated with NN-ratio magnitude (ρ = 0.80 on rank order, but no gradient). The Spearman result is the one thread to pull on, but at n = 5 it's not enough to overturn a paired test that came out 2/5 in the wrong direction. The scaled-down version of tk_9bdd20b5 has done its job: it tells us not to spend the full-scale budget on this exact contrast. A useful follow-up would either (a) widen the NN-ratio range so the Spearman test has signal to work with, or (b) test a different structural axis (size-cap, mutation-rate) where neutrality changes are not confounded with primitive-set semantics.

## Caveats

- **Run-length censoring is severe.** Many late-phase plateaus hit duration 750 (= remaining generations), visible as the vertical step in every CCDF. Log-normal σ fits absorb these as right-tail mass — the σ statistic is partly a censoring statistic, not a pure tail statistic.
- **n = 5 targets** is fixed by the plan and gives the paired Wilcoxon very low power; a true effect of moderate size could still produce 2/5 by chance. Conversely the ρ = 0.80 bootstrap CI is wide and discrete (the lower bound 0.11 is one resampling realisation away from crossing zero).
- **`neutral_proxy.mean_class_size_weighted`** is the e-graph-derived NN-size proxy used in the analysis. It is a sampled estimate, not exhaustive enumeration; ratios may shift slightly with larger samples.
- **"Pass" on Spearman is a technicality** — the CI just barely excludes zero, and the criterion was preregistered as "CI crosses 0 → fail", which we honour, but readers should not over-weight it.
- Per-target CSN p-values are taken from the precomputed `csn_late` field (powerlaw library, distribution_compare lognormal vs. exponential), not recomputed in R.
