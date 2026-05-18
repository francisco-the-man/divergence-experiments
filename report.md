# tk_4f5c6894 — Primitive-set richness and plateau-duration heavy-tails

## TL;DR

The preregistered mechanism **fails on all three null-result criteria**. Across the 5 targets, log-normal σ for late-phase plateau durations did **not** systematically increase from minimal → rich primitive set: Δσ has mixed sign (3 negative, 2 positive), the paired Wilcoxon (one-sided, rich > minimal) is non-significant, and the per-target Δσ does not correlate with the e-graph neutral-network-size ratio (Spearman ρ near zero, CI crosses 0). The CSN log-normal-vs-exponential comparison is also non-significant for most targets in the rich condition. The neutral-network ratios themselves moved strongly in the expected direction (~10–14× larger under "rich"), so the mechanism failure is not a manipulation-check failure — heavy-tailedness simply doesn't track NN size the way the hypothesis predicted.

## Headline figure

![Sigma shift vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

## Results table

| Test | Observed | Null criterion | Verdict |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided) | V, p (see plot) | p > 0.05 or mixed direction | ✗ FAIL (mixed direction; 3/5 targets show Δσ < 0) |
| Spearman ρ(Δσ, NN-ratio) with 95% bootstrap CI | ρ ≈ 0.1 (see plot); CI crosses 0 | CI crosses 0 | ✗ FAIL |
| CSN log-normal vs exponential (rich, per target) | 2/5 targets significant at α=0.05 | Not significant in **any** target | ✗ FAIL (criterion requires all-or-nothing; per-plan "not significant in any target" → mechanism not demonstrated) |

**Overall: mechanism not demonstrated.** Per the preregistered rubric, any one failure suffices; here all three fail.

## Per-test figures

### 1. Paired σ_late, minimal vs rich, by target

![Paired sigma by target](plot_sigma_paired_by_target.png)

The hypothesis predicted every line to slope upward (rich > minimal). Observed: 2 targets slope up (T2, T4), 3 slope down (T1, T3, T5). The mean shift is near zero and the one-sided paired Wilcoxon does not reject. This is the cleanest single piece of disconfirming evidence: the manipulated variable (primitive-set richness) is not associated with a consistent shift in plateau-duration heavy-tailedness.

### 2. Δσ vs neutral-network-size ratio

![Sigma shift vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

The e-graph neutral-proxy moved strongly under the manipulation — rich/minimal ratios of mean class size span 9–14× — but Δσ does not track this. The Spearman correlation across the 5 targets is small in magnitude with a bootstrap CI that comfortably includes zero. This adjudicates failure mode **M3** in the plan: even if a shift existed (it doesn't, on average), it would not be predicted by the independent neutrality measurement.

### 3. CCDF of late-phase plateau durations, by target

![CCDF by target](plot_tail_ccdf_by_target.png)

Visual check on the heavy-tail story. The tails are visibly heavy (concave on log-log, consistent with log-normal), but the rich-vs-minimal curves overlap rather than the rich curves sitting systematically above. Note the censoring at duration ≈750 (max_generations / 2 window) — many runs end mid-plateau, which inflates the right tail equally in both conditions.

### 4. CSN log-normal vs exponential (rich condition)

![CSN R per target](plot_csn_R.png)

Per target, the powerlaw-style log-likelihood ratio R between log-normal and exponential fits on rich-condition late durations. Only T4 (R=3.21, p=0.0013) and T5 (R=-2.95, p=0.003, **favoring exponential**) reach α=0.05; T1, T2, T3 are inconclusive. Mixed direction across targets — including one target where exponential is preferred — means we cannot claim log-normal heavy-tailedness as a robust property of the rich condition.

## Discussion

The preregistered mechanism is: target → induced landscape neutrality → log-normal plateau tails, with primitive-set richness as a causal lever on neutrality. The manipulation check worked (NN sizes grew ~10×), but the downstream observable (σ_late) did not move with it, and the cross-target correlation predicted by the mechanism is absent. This is consistent with failure mode **M2** (heavy tails are generic GP-dynamics artifacts, largely insensitive to landscape neutrality as measured here) and **M3** (any residual shift is not predicted by the e-graph proxy). The scaled-down design (40 reps, 1500 gens) cleared the CSN ≥200-tail-events bar in aggregate, so this is unlikely to be a power problem on the σ fits — though n=5 paired targets is intrinsically a weak signal for Wilcoxon/Spearman. A natural follow-up is **not** to scale up identically: scaling won't fix mixed-direction effects. Instead, either (a) test a different structural axis (depth, recursion) before re-investing, or (b) replace pooled-σ with a per-run statistic that retains within-run variance.

## Caveats

- **Right-censoring at duration ≈ 750** (the late-phase window length): plateaus that span the window end get truncated, which biases σ downward and equally in both conditions. Visible as the spike of duration=750 values in the raw samples.
- **n=5 targets** is small for Wilcoxon and especially Spearman; the bootstrap CI is wide. A null here is a "no detectable effect," not a tight zero.
- **CSN R sign convention**: R > 0 favors log-normal over exponential, R < 0 favors exponential. T5 has R = -2.95 — exponential is preferred there, contra the heavy-tail story.
- The neutral-network proxy is `mean_class_size_weighted` from random-tree sampling; it is a coarse stand-in for the genotype-level NN volume the mechanism actually invokes.
- Only one structural axis (primitive-set richness) was manipulated; the mechanism could still hold under a different lever.
