# GP Plateau Heavy-Tails vs. Primitive-Set Richness — Null Result

**Ticket tk_4f5c6894** · scaled-down first pass of tk_9bdd20b5

## TL;DR

The hypothesis — that enriching the GP primitive set with identity-equivalent / mul-by-constant operators lengthens plateau-duration tails (higher log-normal σ) on the same targets, and that the per-target shift tracks an e-graph neutral-network-size ratio — **fails on this dataset**. Direction of Δσ is mixed across the 5 targets (2 positive, 3 negative); the preregistered one-sided paired Wilcoxon does not reject (p ≈ 0.81); and Spearman ρ between Δσ and the NN-size ratio is **negative** (ρ = −0.30) with a wide bootstrap CI crossing zero. The CSN log-normal-vs-exponential comparison passes in 3 of 5 rich-condition cells but only weakly. **All three null criteria fire ⇒ mechanism not demonstrated at this scale.**

## Headline figure

![Headline](plot_sigma_paired.png)

Within-target paired σ_late, minimal → rich. The mechanism predicts every line slopes up. Three out of five slope **down**.

## Results table

| Test | Observed | Null criterion (any → null) | Pass mechanism? |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided) | V = 5, p ≈ 0.8125; 2/5 targets positive | p > 0.05 or mixed direction | **✗ FAIL** |
| Spearman ρ(Δσ, NN-ratio) | ρ = −0.30, 95% bootstrap CI ≈ [−1.0, +0.9] | CI crosses 0 | **✗ FAIL** |
| CSN log-normal vs exponential (rich, per target) | Significant (p<0.05) for T1, T2 (borderline), T4, T5; T3 not significant | Not significant in **any** rich target | **✓ PASS** (3–4/5 significant) |

Aggregate: ≥1 null criterion triggered ⇒ **mechanism not demonstrated**.

## Per-test figures

### 1. Paired σ_late by target

![Paired sigma](plot_sigma_paired.png)

The hypothesis predicts a monotone rich > minimal shift within every target. Observed: T2 and T4 move in the predicted direction (+0.27, +0.45); T1, T3, T5 move the **opposite** way (−0.32, −0.05, −0.33). A paired one-sided Wilcoxon (V = 5, n = 5, p ≈ 0.81) cannot reject the null of no rich-greater effect. With n = 5 targets, power is modest — but the direction isn't even consistent, which is the more damning observation.

### 2. Δσ vs. neutral-network-size ratio

![Delta sigma vs NN ratio](plot_shift_vs_nn.png)

If neutral-network inflation drives plateau heavy-tails, Δσ should rise with the e-graph NN-size ratio. Empirically the relationship is weakly *negative* (Spearman ρ = −0.30, bootstrap-resample 95% CI on n = 5 spans essentially the whole range). The NN-ratio is also a poor discriminator on its own — it ranges only from ~9 to ~14× across targets, while the sign of Δσ flips. M3 (heavy tails real but not from this neutrality mechanism) is consistent with what we see.

### 3. Pooled late-phase plateau-duration CCDFs

![CCDFs](plot_ccdf_facets.png)

Log-log CCDFs of pooled late-phase plateau durations, faceted by target. Both conditions show clear curvature (heavy-tailed, not exponential) and the rich/minimal curves substantially overlap within each panel. There is no visually obvious systematic separation in tail weight by condition — consistent with the paired-σ result. The tails are dominated by run-truncation events at gen = 750 (the late-phase budget), visible as the vertical spike near x = 750.

### 4. CSN log-normal vs. exponential, rich condition

![CSN](plot_csn_lr.png)

Per-target log-likelihood ratio (powerlaw-package R statistic) for log-normal vs. exponential on pooled late-phase durations in the rich condition. Positive R favours log-normal. T1, T4, T5 show significant log-normal preference; T2 is borderline (p ≈ 0.07); T3 is not significant. So plateau distributions *are* heavier-tailed than exponential in most rich-condition cells — heavy-tailedness itself isn't the problem, the *mechanism linking it to primitive-set richness* is.

## Discussion

The proposed mechanism (target → induced landscape neutral structure → heavy plateau tails, modulated by primitive-set richness) makes two predictions; the data confirms neither. Plateau durations are heavy-tailed (CSN passes most cells), so the phenomenon is real, but the manipulated variable does not move σ in the predicted direction, and the e-graph NN-ratio does not order the per-target shifts. Two things this experiment **cannot** distinguish: (i) the mechanism is wrong; (ii) the mechanism is right but the e-graph proxy used here (random-tree class-size histogram, n=2000 trees) is too coarse / not the relevant neutrality measure during search. A useful follow-up would be measuring neutrality *along realised search trajectories* rather than from random samples — the mismatch between the static NN proxy and the dynamic σ could itself be the finding.

Worth noting: the rich condition shows occasional catastrophic non-convergence (final-best ~1 for T1, ~22 for T4 in some replicates) absent in minimal. The richer primitive set may be *hurting* search on some targets via deceptive intermediate fitness, which would shorten effective plateaus rather than lengthen them — a plausible explanation for the negative Δσ signs.

## Caveats

- n = 5 targets gives the paired Wilcoxon very limited power; even a true moderate-effect mechanism could plausibly be missed. The null is **direction-of-effect mixed**, not "tight CI around zero" — read as "no signal at this scale," not "definitively no effect."
- Plateau durations are right-censored at the late-phase budget (~750 generations). A non-trivial fraction of late events are budget-truncated — biases σ downward and compresses tail discrimination.
- The neutral_proxy is a static random-tree sample, not a measurement of neutrality on realised search paths. The hypothesis as written referenced "e-graph-measured neutral-network sizes"; this is a proxy.
- Scale-down: pop_size = 80, max_gen = 1500, 40 reps. The original tk_9bdd20b5 plan called for 300/4000/120. If the true effect requires deeper search to manifest, this experiment cannot see it.
- CSN R-statistic significance is from the powerlaw package's `distribution_compare` (no bootstrap CI reported per-target).
