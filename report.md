## TL;DR

The mechanism is **not demonstrated**. Across 5 target equations, enriching the GP primitive set with identity-equivalent / multiplicative-constant operators did **not** systematically lengthen plateau-duration tails: only 2/5 targets showed positive Δσ, paired Wilcoxon p = 0.59375 (one-sided rich > minimal). The e-graph neutral-network ratio jumped ~10–14× in the rich condition for every target, yet the per-target Δσ ↔ NN-ratio Spearman ρ = 0.8 has bootstrap 95% CI [0.111, 1.000] — suggestive but rests on n=5. CSN log-normal-vs-exponential preference was only significant in 1/5 rich cells. **Verdict: any one of the three null criteria triggers a fail; here Test 1 and Test 3 both fail.**

## Headline figure

![Per-target σ shift vs neutral-network-size ratio](sigma_shift_vs_nn_ratio.png)

## Results table

| Test | Observed | Null criterion | Verdict |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal) | V = 7, p = 0.59375; 2/5 positive Δσ; mean Δσ = 0.00335 | p > 0.05 OR mixed direction → fail | ✗ |
| Spearman Δσ vs NN-ratio (bootstrap CI) | ρ = 0.80, 95% CI [0.111, 1.000] | CI crosses 0 → fail | ✓ (but n=5) |
| CSN log-normal vs exponential (rich) | 1/5 targets significant & log-normal-preferred (T4); 2/5 significant but exponential-preferred (T1, T5); 2/5 not distinguishable (T2, T3) | not significant in any target → fail | ✗ |

Pre-registered rule: **any one null criterion firing = mechanism not demonstrated.** Two fired.

## Test 1 — Paired Wilcoxon on σ_late

![Within-target paired sigma](sigma_paired_by_target.png)

The lines do not slope consistently upward. Only T2_coupled (Δσ = 0.266) and T4_deep_mul (Δσ = 0.454) move in the predicted direction; T1_poly_sep (Δσ = -0.321), T3_rational (Δσ = -0.050), and T5_transcend (Δσ = -0.332) move the wrong way. The one-sided paired Wilcoxon gives V = 7, p = 0.59375 — nowhere near α = 0.05. **The mechanism's primary prediction — that enriching the primitive set inflates plateau-duration heavy-tails on the same target — does not survive.**

## Test 2 — Δσ vs e-graph NN-size ratio

The Spearman correlation between per-target Δσ and the rich/minimal NN-size ratio is ρ = 0.8, 95% bootstrap CI [0.111, 1.000], not crossing zero. Taken alone, this is the *one* preregistered test that "passes." But interpretation is delicate: with n = 5 targets and a discrete rank statistic, the bootstrap distribution is coarse, and the rank ordering is driven entirely by the fact that the two targets with the largest NN-ratios (T2: 14.1, T4: 13.9) happen to also be the only two with positive Δσ. The correlation says "*if* there's a shift, it tracks NN-ratio order" — but Test 1 says the shift in level is not there on average.

## Test 3 — CSN log-normal vs exponential on rich-condition late-phase tails

![Tail CCDFs by target](tail_ccdf_by_target.png)

Per-target Vuong-style log-likelihood ratios (rich condition): T1_poly_sep R = -3.92, p = 0.000088 (significantly **exponential**-preferred); T2_coupled R = 1.82, p = 0.0688 (not significant); T3_rational R = 1.51, p = 0.130 (not significant); T4_deep_mul R = 3.21, p = 0.00131 (significantly log-normal-preferred); T5_transcend R = -2.95, p = 0.00322 (significantly exponential-preferred). Only T4 matches the mechanism's prediction. The plan required significance "in any target" reading literally — that single hit (T4) technically satisfies the wording, but the spirit of the criterion (heavy-tail-shape supports the log-normal mechanism) is contradicted by T1 and T5 going significantly the other way.

## Discussion

Mechanism M predicted (i) longer plateau tails under richer primitives, (ii) ordering of that shift by independently measured neutral-network ratio, and (iii) log-normal tail shape in the rich condition. We observed (i) no average shift, (ii) rank correlation in the predicted direction but with n=5 and one-sided dependence on two targets, and (iii) heterogeneous tail shapes — one log-normal-significant target, two exponential-significant. Read together: the e-graph measurement confirmed primitive-set enrichment did inflate semantic-class sizes ~10× across the board, but that inflation did not translate into a coherent change in GP plateau-duration statistics at this scale. The mechanism as stated — "neutral-network size controls heavy-tailedness" — is **not** the right summary at pop=80, gen=1500.

Plausible follow-ups: (a) the "rich" condition occasionally **hurts** convergence (T1 rich has mean final-best ~0.12 vs ~0.011 minimal — several replicates got stuck), so the comparison is partly contaminated by run-quality differences, not just landscape neutrality; (b) σ_late as a single scalar may be the wrong observable — Test 3 shows the tail family itself differs across targets; (c) the full-scale tk_9bdd20b5 design (pop=300, gen=4000, n=120) may surface effects the scaled-down version smooths over.

## Caveats

- **n = 5 targets** is the unit of the paired Wilcoxon and Spearman; the bootstrap CI on ρ = 0.80 is wide ([0.111, 1.000]) and the test is inherently underpowered for fine effect sizes.
- T1 (rich) and T5 (rich) had several replicates that failed to converge (final-best ≫ minimal), inflating early-tail mass and potentially deflating σ_late by changing which durations land in the "late" phase. This is a confound between primitive-set richness and run quality.
- The CCDF plot pools across replicates; right-censoring at gen 1500 produces the flat shelves at duration ≈ 750 visible in every panel — the apparent tail mass there is not informative about the underlying distribution.
- Neutral-network "proxy" is the e-graph weighted mean semantic-class size from 2000 sampled trees, not an exhaustive count — directionally correct but a point estimate without CI.
- Bootstrap CI for Spearman ρ uses 5000 resamples of 5 targets; 12 resamples were degenerate (constant column) and dropped — n_bootstrap = 4988.
