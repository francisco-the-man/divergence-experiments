## TL;DR

The preregistered mechanism — *richer primitive set → larger neutral networks → heavier plateau-duration tails* — **was not demonstrated** at this scale. The paired Wilcoxon test on log-normal σ (rich vs. minimal, 5 targets) gave V = 7, p = 0.606 with the direction mixed (only 2/5 targets shifted positive). The correlation between Δσ and the e-graph NN-size ratio was high (Spearman ρ = 0.80) and its bootstrap CI just cleared zero, but the underlying signal it correlates with is itself absent. The CSN log-normal-vs-exponential comparison was significant in only 3/5 targets in the rich condition, and in 2 of those 3 (T2, T4) it pointed *away* from log-normal. Two of three null-criteria fired ⇒ mechanism not demonstrated.

## Headline figure

![Paired σ shift vs. e-graph NN ratio](sigma_shift_vs_nn_ratio.png)

## Results table

| Test | Observed | Null criterion | Verdict |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich vs. minimal, one-sided greater) | V = 7, p = 0.606; 2/5 targets positive, direction mixed | p > 0.05 or direction mixed ⇒ null | ✗ FAIL |
| Spearman ρ(Δσ, NN-ratio) with bootstrap 95% CI | ρ = 0.80, 95% CI [0.11, 1.00] | CI crosses 0 ⇒ null | ✓ PASS (CI excludes 0) |
| CSN log-normal vs exponential (rich, per target) | 3/5 significant at α=0.05; only 2/5 (T1, T5) favor log-normal | LR not significant in any target ⇒ null | ✗ FAIL (strict reading: not all 5) |
| **Overall mechanism demonstrated** | — | All three must pass | **✗ NO** |

## Per-test figures + interpretation

### 1. Paired σ_late by target (the Wilcoxon)

![Paired sigma by target](sigma_paired_by_target.png)

The hypothesis predicts every line slopes *up* from minimal → rich. Instead we see three targets sloping *down* (T1_poly_sep Δσ = -0.32, T3_rational Δσ = -0.05, T5_transcend Δσ = -0.33) and two sloping up (T2_coupled Δσ = +0.27, T4_deep_mul Δσ = +0.45). Mean σ is essentially unchanged (1.875 minimal vs. 1.878 rich). The paired Wilcoxon (V = 7, p = 0.606, one-sided rich > minimal) cannot reject the null, and the *direction-mixed* clause of the preregistered null fires regardless of the p-value.

### 2. Δσ vs. e-graph neutral-network-size ratio (the Spearman)

The headline scatter shows a strong rank correlation (ρ = 0.80, bootstrap 95% CI [0.11, 1.00], p = 0.104). Taken at face value this *passes* the preregistered null-criterion (CI excludes zero). But the y-axis quantity it correlates with — Δσ — has no consistent sign across targets, so what we're seeing is a rank ordering between two small-n variables (n = 5) where the NN-ratio itself spans only ~9.1× to ~14.1×. Read cautiously: with five points and a CI lower bound at 0.11, this is suggestive but not strong evidence; the bootstrap distribution is heavily ceiling-bounded (upper CI = 1.00).

### 3. CSN log-normal vs. exponential, rich condition, per target

![Tail CCDF by target](tail_ccdf_by_target.png)

Across the five rich-condition late-phase tails: T1 strongly favors log-normal (R = -3.92, p = 8.8e-05), T5 favors log-normal (R = -2.95, p = 3.2e-03), but T2 (R = +1.82, p = 0.069), T3 (R = +1.51, p = 0.130) and T4 (R = +3.21, p = 1.3e-03) point in the opposite direction — T4 significantly so. So the "log-normal is the right heavy-tailed model" assumption baked into the σ-as-summary-statistic isn't uniformly supported: in three of five targets the *exponential* fits at least as well. That undermines the σ-comparison itself: when the distribution isn't log-normal, "σ of the log-normal fit" doesn't cleanly mean "tail heaviness."

## Discussion

The mechanism predicted a coherent story: enrich the primitive set → e-graph NN-size ratio jumps ~10–14× → plateau-duration tails get heavier in a way that tracks the NN-size jump per target. We see the middle step (NN ratios are large and roughly comparable across targets, 9.1–14.1) but neither the predicted endpoint (no consistent σ inflation) nor a clean distributional signature (3/5 targets don't even prefer log-normal in the rich condition). The high Spearman ρ is the one positive finding, but with n = 5 targets and a CI that essentially spans [0.11, 1.00], it's better read as "consistent with the mechanism" than "supportive of it." A cleaner reading: at this scale, primitive-set richness doesn't translate into the predicted progress-dynamics signature, even though it clearly inflates the independently-measured neutral-network proxy.

## Caveats

- **n = 5 targets** for the paired test and the Spearman. Power on the Wilcoxon is fixed by n_targets, not n_replicates — scaling up replicates won't help this analysis.
- **σ as a tail-heaviness proxy is only valid when log-normal is the right model.** In 3/5 rich targets it isn't (T2, T3, T4 favor exponential by the CSN LR). The preregistered statistic is partially mis-specified by its own diagnostic test.
- **Plateau-truncation at max_generations=1500.** Many tail durations hit 750 (the late-phase window length) — visible as the right-edge spike in the CCDFs. Censoring will compress σ estimates, possibly differentially across conditions.
- **NN-size proxy is a sampled e-graph count** (`mean_class_size_weighted` over ~2000 random trees), not the full induced fitness landscape's neutral structure. Treat it as ordinal, not metric.
- **Spearman bootstrap CI is ceiling-bound** at 1.00 because n = 5 admits a max ρ of exactly 1; the lower bound (0.11) is the more informative number.
- **Scale-down rationale (pop=80, gens=1500, reps=40)** was preregistered as "if the mechanism shows here we can scale up." It didn't show — scaling up is not automatically warranted; the right next step is fixing the σ-mis-specification issue (e.g. comparing tail exponents under whichever distribution actually fits per cell).
