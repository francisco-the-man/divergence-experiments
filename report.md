# GP plateau heavy-tails do not encode landscape neutrality (tk_4f5c6894)

## TL;DR

The preregistered mechanism — that enriching the GP primitive set with identity-equivalent operators inflates neutral-network sizes and thereby lengthens plateau-duration tails — **fails all three null-result criteria**. Across 5 target equations, the rich-vs-minimal shift in log-normal σ (late phase) is **directionally mixed** (3 of 5 targets go the *wrong* way), the per-target shift does **not** track the e-graph NN-size ratio (Spearman point estimate near zero, CI crosses zero by a wide margin), and the log-normal-vs-exponential CSN comparison on rich-condition tails is **either non-significant or favors exponential** in 4 of 5 targets. The e-graph manipulation worked as intended (~10–15× NN-size ratio in every target), so the negative result is about the *consequence*, not the manipulation.

## Headline figure

![Headline](plot_sigma_paired_by_target.png)

## Results table

| Preregistered test | Observed | Null criterion | Pass / Fail |
|---|---|---|---|
| Paired Wilcoxon, σ_late rich > minimal (one-sided, n=5) | V=8, p ≈ 0.50; direction mixed (3 down, 2 up) | p>0.05 OR mixed direction → FAIL | ✗ FAIL |
| Spearman ρ(Δσ_late, NN-size ratio), bootstrap 95% CI | ρ ≈ +0.10, 95% CI ≈ [−0.90, +0.90] | CI crosses 0 → FAIL | ✗ FAIL |
| CSN log-normal vs exponential on rich-condition tails, per target | T1 R=−3.92 p<0.001 (favors exp); T2 R=+1.82 p=0.069 (ns); T3 R=+1.51 p=0.130 (ns); T4 R=+3.21 p=0.001 (favors logN); T5 R=−2.95 p=0.003 (favors exp) | Any target ns → FAIL | ✗ FAIL (3 ns/wrong-way, only T4 supports log-normal) |

**Verdict: mechanism not demonstrated.** Any *one* of these triggers the prereg null; all three did.

## Per-test figures

### 1. Paired σ_late by target

![Paired sigma](plot_sigma_paired_by_target.png)

The mechanism predicts every line slopes upward (rich > minimal). Instead, **T1 and T5 slope clearly downward, T3 is flat, only T2 and T4 slope up**. The paired Wilcoxon one-sided test for "rich > minimal" returns p ≈ 0.5 — the data are perfectly consistent with no directional effect. This is the cleanest piece of evidence against the hypothesis as stated.

### 2. Δσ vs neutral-network size ratio

![Delta sigma vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

The e-graph manipulation succeeded: rich/minimal NN-size ratios are uniformly ~9–14× across targets, well-separated from 1. If neutrality drove plateau tails, points should lie on a positive-slope line. They do not — the Spearman ρ is near zero with a bootstrap CI that spans nearly the full [−1, +1] range (n=5 is small, but the *point estimate* gives no encouragement). Two targets with the *largest* NN ratios (T2, T4) do show positive Δσ, but T1 and T5 with mid-range ratios shift *negatively*, breaking any monotone relationship.

### 3. Tail CCDFs, faceted by target

![CCDF](plot_tail_ccdf_by_target.png)

Visual inspection of the late-phase plateau-duration CCDFs shows heavy tails are present in *both* conditions across all targets — i.e. heavy-tailedness itself replicates as a generic GP phenomenon — but the *rich* curve is not systematically above the *minimal* curve in the upper tail. In T1 and T5 the rich tail is visibly *shorter*; in T4 it is longer. This matches the σ table and rules out a "subtle but consistent" rich-shifts-tails effect that the σ statistic might miss.

### 4. CSN log-normal vs exponential, rich condition

![CSN](plot_csn_lr_by_target.png)

Only T4 cleanly prefers log-normal over exponential (R=+3.21, p=0.001). T2 and T3 are indistinguishable. **T1 and T5 significantly favor exponential** — meaning the rich-condition tails on those targets are *not even heavy-tailed in the CSN sense*. The third null criterion fails on these grounds: the distributional form the mechanism assumes is not present in 3 of 5 targets.

## Discussion

The proposed mechanism — target → neutral-network structure → heavy-tailed plateau durations — predicts that the *independently-measured* e-graph NN-size inflation under the rich primitive set should manifest as longer plateau tails. The manipulation worked (10–15× NN inflation, robust across targets), but the predicted downstream effect did not appear, and where σ did shift it did not track NN-size. This is consistent with failure modes M2 ("plateau heavy-tails are generic GP-dynamics artifacts") and M3 ("shift uncorrelated with neutrality") flagged in the plan. A confounder worth noting: in T1, T3, T5 the rich condition also produces some catastrophically-bad runs (final_best > 1 where minimal stays ~0.01), suggesting the extra operators (especially `/`) destabilize search in ways orthogonal to neutrality — late-phase plateau statistics in those runs may reflect "stuck on a bad attractor" rather than "exploring a neutral network." A follow-up that *only* adds identity-equivalent unaries (no `/`, no extra constants) could disentangle this.

## Caveats

- n=5 targets is the prereg design; the Spearman test is genuinely underpowered, but the **paired Wilcoxon and per-target CSN tests do not depend on it** and they also fail.
- σ is fit on pooled events per (target × condition); replicate-level bootstrap CIs are not shown here (would need re-running the fits).
- The "rich" primitive set bundles several changes (`neg`, `id`, `/`, extra constants 0.0/2.0). T1/T3/T5 instability suggests `/` is doing something other than adding neutral genotypes; the design cannot isolate which sub-component matters.
- CSN log-likelihood-ratio test requires enough tail events; n_tail is 72–147 per (target × rich), comfortably above the CSN-2009 ≥200 bar **only** when pooling early+late, which we did not do — late-only n_tail is 72–147 so results should be read as suggestive on individual targets.
- Heavy-tailedness *itself* does replicate cleanly in every cell — the negative result is specifically about the rich-vs-minimal *shift*, not about the existence of heavy tails.
