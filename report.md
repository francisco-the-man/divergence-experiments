# Mechanism not demonstrated: primitive-set enrichment does not systematically lengthen plateau-duration tails

## TL;DR

The pre-registered mechanism predicted that enriching the GP primitive set (which inflates e-graph-measured neutral-network sizes by ~10–20×) would (i) increase the log-normal σ of late-phase plateau durations, (ii) do so in proportion to the per-target NN-size ratio, and (iii) yield log-normal-over-exponential tails in the rich condition. **All three pre-registered null criteria triggered.** Shift in σ_late is mixed in sign across the 5 targets (2 positive, 3 negative), the Spearman correlation with NN-size ratio is essentially zero with a CI that crosses 0, and the CSN log-normal-vs-exponential comparison in the rich condition is inconsistent across targets (two targets actually favor exponential, two favor log-normal, one is null). This is a clean disconfirmation of the proposed mechanism at this scale.

## Headline figure

![Paired σ shift](plot_sigma_paired_by_target.png)

## Results table

| Pre-registered test | Observed | Null criterion | Outcome |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided, n=5) | V=6, p≈0.50; direction mixed (2/5 positive) | p>0.05 OR mixed direction | ✗ FAIL (null triggered) |
| Spearman ρ(Δσ_late, NN-size ratio) | ρ ≈ +0.10, bootstrap 95% CI ≈ [−0.90, +0.90] (crosses 0) | CI crosses 0 | ✗ FAIL (null triggered) |
| CSN log-normal vs exponential (rich, per target) | T1: R=−3.92 p<.001 (favors exp); T2: R=+1.82 p=.07; T3: R=+1.51 p=.13; T4: R=+3.21 p=.001 (favors lnorm); T5: R=−2.95 p=.003 (favors exp) | Not significant in any target | ✗ FAIL (null triggered — direction inconsistent; 2 targets significantly favor exponential) |

**Mechanism status: NOT DEMONSTRATED.** All three independent failure modes specified in the pre-registration fired.

## Per-test figures + interpretation

### 1. Paired σ_late by target

![Paired sigma](plot_sigma_paired_by_target.png)

The hypothesis predicted every line slopes upward (rich > minimal). Instead, T2 and T4 rise, T1, T3, T5 fall. The paired Wilcoxon (one-sided, greater) gives p≈0.50. There is no systematic enrichment effect on σ_late. Note also that even the average effect is roughly zero — this is not a power problem from n=5 targets; the effect direction itself is not consistent.

### 2. Δσ_late vs e-graph NN-size ratio

![Shift vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

The mechanism predicts a positive slope: targets whose neutral-network sizes inflate more under enrichment should show larger σ shifts. The observed scatter is essentially flat; ρ_Spearman ≈ 0.1 with a bootstrap CI spanning nearly the full [−1, 1] range. The two targets with the largest NN-size ratios (T3, T5) actually show negative σ shifts. The independent e-graph measurement of neutrality does not predict the dynamical observable.

### 3. Late-phase plateau-duration CCDFs

![CCDF by target](plot_tail_ccdf_by_target.png)

Pooled late-phase plateau durations on log-log axes, per target × condition. Visually, the rich and minimal tails are largely overlapping, with no consistent rightward shift in the rich condition. The CSN per-target log-likelihood ratios confirm this: T1 and T5 in the rich condition are *better fit by exponential than log-normal* (negative R, p<0.01) — the opposite of what a neutrality-inflation mechanism predicts. Only T4 shows the predicted log-normal preference.

## Discussion

Mechanism M proposed that target equation → induced-landscape neutral-network structure → heavy-tailed plateau durations, and that the e-graph neutral-network-size ratio between primitive sets should predict the magnitude of the dynamical shift. Mechanism M predicted (a) rich > minimal in σ_late uniformly, (b) Δσ correlated with NN ratio, (c) log-normal-favored tails in rich. We observed (a) mixed direction, (b) ρ≈0 with wide CI, (c) inconsistent and partly exponential-favoring. This is the failure-mode signature of **M2 / M3 from the pre-registration**: either plateau heavy-tailedness is not driven by neutral-network size, or the e-graph measurement is not the right neutrality proxy for this dynamical observable. The data cannot distinguish these two failures, but either disconfirms the specific mechanism as stated. A scale-up of tk_9bdd20b5's design would not rescue this — the failure is in *direction*, not in *power*. The natural follow-up is to interrogate the proxy itself: does e-graph class size correlate with anything dynamical (jump sizes? early-phase σ? success rate?), or is the assumed proxy-mechanism link broken at the root.

## Caveats

- n=5 targets; with this n, even a tight Spearman estimate has a CI spanning most of [−1, +1]. The null on test 2 is "CI crosses 0", which is satisfied, but a positive trend cannot be excluded — what *is* excluded is a strong, reliable correlation.
- Plateau durations are right-censored at `max_generations - boundary` (≈750); many late plateaus saturate. This compresses upper-tail discrimination and could attenuate σ differences in either direction.
- "Rich" primitive set in the rich condition also occasionally produces catastrophic runs (per-rep final-best up to ~1.5 in T1 rich, vs. ~0.018 max in T1 minimal). The σ estimate is conditional on a run producing late-phase plateaus at all; selection effects on which runs contribute events differ between conditions.
- Neutral-proxy is "mean equivalence-class size, sample-weighted" over 2000 random trees — a static structural measurement, not a dynamics-weighted one. The mechanism could in principle be rescued by a dynamics-weighted neutrality measure; this experiment does not test that.
- Scaled-down configuration (pop=80, gens=1500, reps=40). The pre-registration noted statistical power for the n=5 paired test does not depend on these — but censoring at 750 generations does.
