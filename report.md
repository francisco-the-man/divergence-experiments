# GA plateau heavy-tails vs. neutral-network structure — scaled-down first pass (tk_4f5c6894)

## TL;DR

The preregistered mechanism check **fails 2 of 3 tests**. The paired Wilcoxon on log-normal σ (rich vs minimal primitives) is null and direction-mixed (V = 7, p = 0.606, 3 of 5 targets move the *wrong* way; mean Δσ = 0.003). The CSN log-normal-vs-exponential comparison in the rich condition is significant in only 3 of 5 targets. The one positive signal is the per-target Spearman between Δσ and the e-graph NN-size ratio (ρ = 0.80, bootstrap 95% CI [0.11, 1.00], CI excludes zero) — but with n=5 targets and the primary contrast itself null, that correlation is not enough to rescue the mechanism. **Verdict: mechanism not demonstrated at this scale.**

## Headline figure

![Paired σ shift by target](sigma_paired_by_target.png)

## Results table

| Test | Observed | Null criterion (fail if…) | Verdict |
|---|---|---|---|
| Paired Wilcoxon, σ_late (rich > minimal, one-sided) | V = 7, p = 0.606; 2/5 positive, 3/5 negative | p > 0.05 or direction mixed | ✗ |
| Spearman ρ(Δσ, NN-ratio), 95% bootstrap CI | ρ = 0.80, CI [0.11, 1.00] | CI crosses 0 | ✓ |
| CSN log-normal vs exponential, rich condition, each target | 3/5 significant at α=0.05 | any target non-significant | ✗ |
| **Overall mechanism passes all three** | — | — | **✗** |

## Per-test figures

### σ shifts are direction-mixed across targets

![Paired σ shift by target](sigma_paired_by_target.png)

Within-target lines connect the minimal- and rich-primitive σ_late estimates. The mechanism predicts every line slopes upward (rich > minimal). Instead, only T2_coupled (Δσ = 0.266) and T4_deep_mul (Δσ = 0.454) move in the predicted direction; T1_poly_sep (Δσ = -0.321), T3_rational (Δσ = -0.050), and T5_transcend (Δσ = -0.332) move the wrong way. Mean Δσ is essentially zero (0.003). The one-sided Wilcoxon p = 0.606 reflects this — there is no detectable population-level shift toward heavier tails in the rich condition.

### Δσ correlates with NN-size ratio, but the underlying shifts are noisy

![Δσ vs NN ratio](sigma_shift_vs_nn_ratio.png)

Per-target Δσ correlates positively with the e-graph NN-size ratio (ρ = 0.80, bootstrap 95% CI [0.11, 1.00]). The CI excludes zero, so this test passes the preregistered criterion in isolation. But two caveats sharply limit the inference: (i) with n=5 the Spearman discretizes to only a few possible ρ values and the bootstrap CI is very wide; (ii) the y-axis variable (Δσ) is the same one that came up null in the Wilcoxon — a correlation between a null-mean shift and an independent ranking can arise from rank-coincidence rather than mechanism. NN ratios themselves span a narrow band (9.1× to 14.1×), so this is a weak ordinal contrast.

### Log-normal vs exponential: rich-condition tails are not uniformly heavy

![CSN log-normal vs exponential by target](csn_lognormal_vs_exponential.png)

The CSN log-likelihood-ratio test (R > 0 favors log-normal, R < 0 favors exponential) is significant at α=0.05 in only 3 of 5 targets in the rich condition, and crucially two of those significant cases (T1_poly_sep R = -3.92, T5_transcend R = -2.95) point the *wrong way* — exponential fits the late-phase plateau distribution better than log-normal. Only T4_deep_mul (R = 3.21, p = 0.001) shows the predicted "heavy-tailed log-normal beats exponential" pattern significantly; T2_coupled and T3_rational are not distinguishable. This means the rich condition does not produce uniformly heavy-tailed plateau distributions even before comparison with minimal.

## Discussion

The mechanism — *enriching the primitive set inflates neutral-network sizes → fatter plateau-duration tails* — predicts a coherent within-target shift that we don't see. The intervention *does* substantially inflate the e-graph NN-size proxy (rich/minimal ratio ~9–14×), so the manipulation worked at the landscape level. What didn't follow is the dynamical consequence: σ_late moves in different directions on different targets, and on the two targets where it moves *down* with enrichment (T1, T5), the late-phase rich-condition distribution is in fact closer to exponential than log-normal. This is consistent with failure mode M3 in the plan (heavy tails exist but are not driven by the proposed neutrality mechanism) and partially with M2 (some heavy-tail behavior may be generic GP dynamics). The surviving Spearman is suggestive but, given n=5 and the null primary contrast, not load-bearing. A natural follow-up is to test whether the rank ordering of Δσ vs NN ratio replicates at the full scale of tk_9bdd20b5 — if it does with the same direction, the mechanism might be real but require finer measurement; if it doesn't, the Spearman here was a coincidence over five points.

## Caveats

- **n = 5 targets** is the binding constraint on the paired Wilcoxon and Spearman; scaling up replicates (40 → 120) wouldn't change either. To strengthen, add targets, not reps.
- The "neutral-network size" measure is a coarse e-graph proxy (mean weighted equivalence-class size over a 2000-tree sample), not a true count of the induced fitness landscape's neutral networks.
- CSN R-values are reported as continuous; the preregistered criterion treats *any* per-target non-significant LR as a fail, which is strict — a less strict rule (e.g. "majority of targets significant in the right direction") would still fail here (only 1/5 targets, T4, is significantly log-normal>exponential).
- One-sided Wilcoxon (the preregistered alternative) — the two-sided p would be larger still given the mixed direction.
- T3_rational rich runs show much better final fitness than minimal (median 0.006 vs 0.34), so the late-phase windows are not on equivalent fitness landscapes — late-phase pooling assumes the search has reached a comparable regime, which may be violated.
- Bootstrap CI for Spearman is over only 5 paired points; some resamples are degenerate (4469 of 5000 retained).
