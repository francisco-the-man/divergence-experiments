# Mechanism not demonstrated: primitive-set enrichment does not systematically lengthen plateau-duration tails

## TL;DR

Preregistered paired Wilcoxon on log-normal σ (late-phase plateaus, rich vs minimal primitive set) across the 5 targets is **not significant** and the direction is **mixed** (2/5 targets shift positive, 3/5 negative). The Spearman correlation between per-target Δσ and the e-graph neutral-network size ratio is **negative** with a wide CI crossing zero. Null-result criterion (i) fails outright; criterion (ii) also fails. The CSN log-normal-vs-exponential comparison (criterion iii) is mixed across targets. **Mechanism not demonstrated at this scale.**

## Headline figure

![Paired σ shift](plot_sigma_paired_by_target.png)

The hypothesis predicted every line to slope **upward** from minimal → rich. Three of five slope down. This is the clearest visual statement of the null.

## Results table

| Preregistered test | Observed | Null criterion (mechanism fails if…) | Pass / Fail |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided) | V = 6, p ≈ 0.59; direction mixed (2+/3−) | p > 0.05 OR mixed direction | ✗ **FAIL** |
| Spearman ρ(Δσ, NN ratio) with bootstrap 95% CI | ρ = −0.30, 95% CI ≈ [−1.0, +0.9] | CI crosses 0 | ✗ **FAIL** |
| CSN log-normal vs exponential (rich, per target) | T1 p=8.8e-5 (LN preferred); T2 p=0.069; T3 p=0.13; T4 p=0.0013 (LN); T5 p=0.0032 (LN). 3/5 significant | Not distinguishable in any target | ✓ **PASS** (in 3/5 — partial) |

Per the preregistration: **any one** failing criterion = mechanism not demonstrated. Criteria (i) and (ii) both fail.

## Per-test figures + interpretation

### Test 1 — Paired σ_late, rich vs minimal

![Paired sigma by target](plot_sigma_paired_by_target.png)

The hypothesis was directional: enriching primitives with identity-equivalent / mul-by-1 operators should **uniformly inflate** late-phase plateau-duration heavy-tailedness because more genotypes map to the same phenotype. The data shows σ_rich > σ_minimal on only 2 of 5 targets (T2_coupled +0.27, T4_deep_mul +0.45) and σ_rich < σ_minimal on the other 3 (T1, T3, T5 all between −0.05 and −0.33). The one-sided paired Wilcoxon (alternative='greater') is not significant. The mixed direction alone trips the preregistered null criterion.

### Test 2 — Δσ vs e-graph NN-size ratio

![Sigma shift vs neutral-network ratio](plot_sigma_shift_vs_nn_ratio.png)

The mechanism's second commitment: even if shifts vary, their magnitudes should track the independent e-graph measurement of how much more redundant the rich primitive set is on each target. They don't. The point cloud is essentially flat-to-negatively-sloped (Spearman ρ = −0.30 across n=5; bootstrap CI is enormous and brackets zero by a wide margin). The two targets where σ moved in the predicted direction (T2, T4) actually have **higher** NN ratios than two of the targets that moved the wrong way (T3, T5) — but the targets with the largest NN ratios are not the targets with the largest positive Δσ. There is no signal here.

### Test 3 — Pooled late-phase plateau-duration tails

![CCDF by target](plot_tail_ccdf_by_target.png)

CCDFs of pooled late-phase plateau durations, faceted by target, colored by condition. Visually: the curves are similar within each panel. Where the rich condition's tail is heavier (T2, T4) it matches the σ table; where it's lighter (T1, T5) it also matches. **The CSN log-normal-vs-exponential comparison favors log-normal significantly in 3/5 rich-condition targets** (T1, T4, T5; p < 0.05) — meaning plateau durations *are* heavy-tailed enough to be distinguishable from exponential in most cases, satisfying criterion (iii). So the tails are real; the **directional dependence on primitive-set richness** is what's absent.

## Discussion

The mechanism as stated predicted three things in sequence: (a) plateau durations heavy-tailed, (b) heavy-tailedness inflated by primitive-set enrichment, (c) per-target shift magnitude tracks the e-graph NN-size ratio. Part (a) survives (criterion iii passes in 3/5 cells). Parts (b) and (c) do not. The cleanest interpretation: at this scale, plateau-duration tails are a generic GP-dynamics phenomenon (failure mode M2 in the plan), not a reflection of the induced fitness landscape's neutral-network structure as measured by e-graph equivalence classes. The 10–14× inflation of mean class size in the rich condition (a large and consistent effect on the independent variable) produced no consistent σ-response — so it isn't that the manipulation was too weak. Two follow-ups worth considering before scaling up: (1) check whether the per-target final-fitness gap (rich often *worse* final best on T1, T3, T5 — e.g. T1_rich mean = 0.124 vs T1_min = 0.011) is contaminating the late-phase window definition, since "late phase" may be measuring different optimisation regimes across conditions; (2) test alternative structural axes (e.g. arity, depth budget) before declaring the whole mechanism dead — only the primitive-set axis was probed here.

## Caveats

- **n=5 targets** is the binding constraint on statistical power for both Wilcoxon and Spearman; the bootstrap CI on ρ is essentially uninformative. The preregistered design accepted this — power on these tests is target-count-limited, not replicate-count-limited.
- **Late-phase σ depends on the late-phase window definition** baked into the experiment code. The rich condition sometimes terminates with much worse final fitness (notably T1, where ~4/40 rich replicates ended near final_best ≈ 1.0 instead of ~0.01), so "late-phase" plateaus may sample a different regime across conditions.
- **CSN test in criterion (iii)** is computed from the `csn_late.R` / `csn_late.p` fields in the result rows (signed log-likelihood ratio, two-sided p) — this is `powerlaw.distribution_compare('lognormal', 'exponential')` as preregistered. Sign convention: R > 0 favours log-normal.
- **e-graph "neutral_proxy"** is a sampling estimate (`mean_class_size_weighted` over `n_valid=2000` random trees); it is consistent within condition but is a proxy for the true NN-size ratio.
- This is the **scaled-down** version of tk_9bdd20b5 (40 reps × pop_size 80 × 1500 gens). A scale-up could in principle reveal an effect this run missed, but the *direction* being mixed (not just the magnitude small) argues against simple under-powering as the explanation.
