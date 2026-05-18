# tk_4f5c6894 — Neutral-network mechanism for GP plateau heavy-tails: **NULL**

## TL;DR

The pre-registered prediction was that enriching the GP primitive set with identity-equivalent / mul-by-constant operators would lengthen plateau-duration tails (higher log-normal σ in the late phase), and that the per-target Δσ would correlate positively with an independent e-graph neutral-network-size ratio. **Neither pattern holds.** Δσ direction is split 2-positive / 3-negative across the 5 targets (paired one-sided Wilcoxon p = 0.500), the Spearman ρ between Δσ and NN-ratio is essentially zero with a CI spanning ±1, and the CSN log-normal-vs-exponential model comparison only favours log-normal at α=0.05 in 1 of 5 rich-condition targets (and actively favours **exponential** in 2). All three null-result criteria fire — the mechanism, as operationalised, is not demonstrated.

## Headline figure

![Sigma shift vs NN-ratio](plot_sigma_shift_vs_nn_ratio.png)

The hypothesis predicts a positive slope: bigger neutral networks (rich/minimal NN-ratio on x) → larger plateau-tail inflation (Δσ on y). The 5 targets show no such relationship — Δσ straddles zero independent of NN-ratio.

## Results table

| Pre-registered test | Observed | Null-result criterion | Pass? |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal), n=5 targets | V = 6, **p = 0.500**; mean Δσ = +0.003, direction 2+/3− | p > 0.05 **or** mixed direction → null | ✗ FAIL (both) |
| Spearman ρ(Δσ, NN-ratio), bootstrap 95% CI | **ρ = 0.300**, CI = [−1.000, 1.000] | CI crosses 0 → null | ✗ FAIL |
| CSN log-normal vs exponential (rich, per target) | 1/5 favours log-normal at α=0.05 (T4); 2/5 favour **exponential** (T1, T5); 2/5 indistinguishable | Not significant in **any** target → null | ✗ FAIL (≥3 targets) |

All three preregistered null criteria fire. **Mechanism not demonstrated.**

## Per-test figures

### Test 1 — Paired σ_late, rich vs minimal

![Paired sigma by target](plot_sigma_paired_by_target.png)

Within-target lines should slope upward if rich primitives inflate plateau-tail σ. Three of five targets slope **down** (T1, T3, T5); two slope up (T2, T4). The mean within-target shift is ≈ 0 (+0.003) and the one-sided Wilcoxon p = 0.500. Not consistent with the mechanism's directional prediction.

### Test 2 — Δσ vs e-graph NN-size ratio

![Sigma shift vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

NN-ratios are all large (9–14×) — the manipulation **did** inflate neutral-network sizes by an order of magnitude per the independent e-graph proxy, so the manipulation check works. But Δσ does not track NN-ratio: Spearman ρ = 0.300 with bootstrap CI [−1.000, 1.000] (n=5 is too small to constrain ρ regardless, but the point estimate is also small). The proposed coupling between NN size and plateau-tail heaviness is not detectable here.

### Test 3 — CSN log-normal vs exponential (rich condition, late phase)

![CSN model comparison](plot_csn_lognormal_vs_exp.png)

Per-target Vuong-style log-likelihood ratio R (positive ⇒ log-normal preferred, negative ⇒ exponential preferred), with the α=0.05 significance bands shaded. Only **T4_deep_mul** clears the bar for log-normal (R=+3.21, p=.001). **T1_poly_sep** (R=−3.92, p<.001) and **T5_transcend** (R=−2.95, p=.003) significantly prefer the exponential. T2 and T3 are indistinguishable. The CSN criterion required log-normal preference at α=0.05 in *every* target — it fails in 3 of 5, and in 2 of those the preferred distribution is the *thin-tailed* alternative.

### Diagnostic — pooled CCDFs

![CCDF by target](plot_tail_ccdf_by_target.png)

Visual check on the pooled late-phase plateau durations: the rich (orange) vs minimal (red) CCDFs largely overlap within each target panel. No systematic rightward shift of the rich condition's tail.

## Discussion

The data clearly disconfirm the operationalised hypothesis at this scale. The e-graph manipulation check passed — rich primitives inflated semantic-equivalence-class sizes by ~10× over minimal — but that inflation did **not** translate into the predicted plateau-tail lengthening, nor did the per-target magnitude track NN-ratio. Worse for the CSN sub-hypothesis: pooled late-phase durations in the rich condition are *not* uniformly heavy-tailed; in two targets (T1, T5) an exponential fits significantly better than a log-normal. This is more consistent with failure-mode **M2** (heavy-tail appearance in GP progress curves is not driven by landscape neutrality in the way the e-graph proxy measures it) than with **M3** (real shift but uncorrelated). The mechanism as written — *NN-size ratio → plateau-σ ratio* — is not supported. Plausible follow-ups: (a) examine whether the e-graph proxy captures the *relevant* neutral structure (it counts semantic equivalence at a fixed sample, not the operationally accessible neutral neighbourhood under the GP's actual variation operators); (b) re-examine the plateau-detection / phase-split choices, since `n_late_events` varies 72–147 across cells and several plateaus saturate at the 750-generation late-phase ceiling.

## Caveats

- **n_targets = 5** for the paired Wilcoxon and Spearman. Both tests are underpowered by design (the plan acknowledges this); a clean positive would still have been visible as a monotone pattern, which is absent. The Spearman CI [−1, 1] is uninformative on its own — interpretation rests on the point estimate (0.300) and the visibly non-monotone scatter.
- **Right-censoring at 750 generations** affects a non-trivial fraction of late-phase plateaus (visible as the spike at duration=750 in the CCDFs). σ estimates from a log-normal fit on censored data are biased; treat absolute σ values as comparative only.
- **NN-ratio proxy** is `mean_class_size_weighted` from random-tree sampling at fixed size_cap. It measures semantic redundancy of the primitive set on random genotypes, not the neutral neighbourhood actually traversed by the GP. A null on this proxy doesn't rule out a different operationalisation of "neutral-network size."
- **CSN R/p in the result JSON** is the powerlaw library's log-normal-vs-exponential comparison (Vuong-style); positive R favours log-normal. I report it verbatim — no re-fitting in R.
- Scaled-down run (pop=80, gen=1500, reps=40). Original tk_9bdd20b5 spec was 300/4000/120. Effect sizes here may not reflect the full-scale regime, though the CSN tail-event bar (≥200) is cleared by all cells via pooling.
