# GA plateau heavy-tails vs primitive-set richness — NULL on the headline test

## TL;DR

The mechanism predicted that enriching the GP primitive set with identity-equivalent operators would systematically **increase** late-phase plateau-duration heavy-tailedness (log-normal σ) on every target, with the per-target shift tracking an independent e-graph neutral-network-size ratio. We see neither: the paired Wilcoxon is null (direction mixed, **2/5 targets shift in the predicted direction**), and Spearman ρ between Δσ and NN-ratio is essentially zero. The e-graph manipulation **did** robustly enlarge neutral-network proxies by ~6–23×, so the manipulation worked — the proposed dynamical consequence didn't follow. **Mechanism not demonstrated at this scale.**

## Headline figure

![Paired σ shift by target](plot_sigma_paired_by_target.png)

Lines cross. If the mechanism held, every line would slope up (minimal → rich). Three of five slope **down**.

## Results table

| Pre-registered test | Observed | Null criterion | Pass? |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided) | V=7, p≈0.50; 2/5 targets positive | p>0.05 OR mixed direction → fail | ✗ FAIL |
| Spearman ρ(Δσ, NN-ratio), bootstrap 95% CI | ρ≈0.10, 95% CI crosses 0 (n=5) | CI crosses 0 → fail | ✗ FAIL |
| CSN log-normal vs exponential (rich, per target) | T1 p<0.001 (log-normal preferred), T4 p=0.001, T5 p=0.003; T2 p=0.07, T3 p=0.13 | not significant in ANY target → fail | ✓ PASS (3/5 targets) |

Per the pre-registered `null_result_criteria` ("any one of i/ii/iii = mechanism not demonstrated"), **the mechanism is not demonstrated**: criteria (i) and (ii) both fail. Criterion (iii) passes for 3 of 5 targets — heavy tails are real where they appear, just not modulated by primitive richness in the predicted way.

## Per-test figures

### Δσ vs neutral-network-size ratio (the correlation test)

![Sigma shift vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

The e-graph proxy shows the manipulation worked: mean equivalence-class size jumps from ~3–9 (minimal) to ~50–86 (rich), a 6–23× inflation. But Δσ scatters around zero with no monotone relationship to the ratio. The hypothesis required a positive, ordered slope from this plot; we get noise. This is the cleanest disconfirmation of the proposed mechanism — the independent landscape measurement moves a lot, the dynamical observable doesn't track it.

### Pooled late-phase CCDFs

![Tail CCDFs by target](plot_tail_ccdf_by_target.png)

Heavy tails are visually present in every cell (curves decay slower than exponential on log-log), consistent with criterion (iii) passing in most targets. But the rich-vs-minimal CCDF orderings are inconsistent across panels: T2 and T4 show the predicted "rich has fatter tail," T1 and T5 show the opposite, T3 is overlapping. This is the same story as the σ plot, just visualised on the raw durations.

### Early-phase σ (exploratory, not pre-registered)

![Early phase sigma](plot_sigma_early.png)

*Flagged as exploratory.* The early-phase σ pattern is similarly mixed and slightly compressed in the rich condition (4/5 targets have σ_early_rich < σ_early_minimal). If anything, richer primitives produce slightly **less** heavy-tailed plateau distributions early on — opposite to the hypothesis. Not in the pre-registered plan; included only to show the null isn't an artifact of the late-phase windowing.

## Discussion

The pre-registered mechanism — *primitive-set neutrality → larger neutral networks → heavier plateau tails* — fails its directional prediction at this scale. The manipulation check (e-graph NN-size ratios) confirms the rich condition really does induce vastly more genotype-phenotype redundancy, so the null isn't "the knob didn't turn." What we can say: heavy-tailed plateau durations are a robust phenomenon (log-normal beats exponential on the pooled tails for 3/5 targets in the rich condition), but their σ is not a simple readout of static neutral-network size. The most likely culprits are (a) other landscape properties (basin geometry, mutation operator coupling) dominate plateau dynamics, (b) the "neutral_proxy" mean-class-size statistic is not the right summary of the relevant landscape feature, or (c) richer primitives also increase exploration variance, which can *shorten* some plateaus by accident — partially offsetting the predicted lengthening. Follow-up: don't scale this experiment up; redesign the landscape observable.

## Caveats

- n=5 targets gives weak power on both the Wilcoxon and Spearman; CIs are wide. A positive result here would have been suggestive, not definitive — but a directional split of 2/5 cannot be rescued by power alone.
- "Late phase" is defined by the experiment code's windowing; we did not re-derive it from raw curves.
- `neutral_proxy.mean_class_size_weighted` is a sampled e-graph estimate (n_valid=2000), not exact; rich condition has fewer classes covering the same sample, so the ratio is well-identified at least directionally.
- Plateau durations are right-censored at max_generations=1500; the saturating "750"s in pooled samples are this censoring. Log-normal fits ignore the censoring — this is consistent across conditions and shouldn't bias the rich-vs-minimal comparison, but absolute σ values should not be taken as estimates of an uncensored tail.
- Scaled-down run (pop=80, gen=1500, reps=40). The pre-registered tail-event bar (≥200) is **not** cleared per cell (most cells have 70–150 late events), contrary to the plan's projection. Still well above 30, but worth noting.
