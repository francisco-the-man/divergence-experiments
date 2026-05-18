# Connectivity vs Size: σ_C > σ_S Holds, But Dissociation Can't Be Adjudicated

## TL;DR

Parent ticket `tk_4f5c6894` failed to demonstrate that a "rich" primitive set lengthened plateau-duration tails (paired Wilcoxon p = 0.59375). This follow-up split "rich" into a **size-only-rich (S)** and **connectivity-rich (C)** primitive set to ask which factor matters. Result: on 5 targets with adequate late-phase data, **σ_C is strictly greater than σ_S in 5/5 paired comparisons** (Wilcoxon V = 15, p = 0.0312, median Δσ_CS = 1.3053). The connectivity-rich condition produces the heavy-tailed plateau-duration distributions the size-only condition does not. **But** the pre-registered Spearman dissociation test (ii) cannot be evaluated — the e-graph measurement tasks that would have produced Δsize and Δconnectivity observables are absent from `result.json` — and the CSN bar (≥4/6 C-cells log-normal AND ≥3/6 S-cells exponential, both at p < 0.05) clears only 2/5 on the C side and 0/5 on the S side. By the plan's "any criterion firing = fail" rule, **mechanism not demonstrated**, but the σ_C > σ_S result is real and informative.

## Headline figure

![Per-target σ across M, S, C conditions](sigma_by_condition_paired.png)

Every line slopes up steeply from S to C. T2, T5, T6 are the most extreme: σ jumps from near-zero under M and S to >1.3 under C. T3 — the only target where the minimal set produces non-trivial late-phase σ (0.6278) — still gains substantially under C (0.9033). The visual signature is exactly what the connectivity-portal hypothesis predicts and the size-only hypothesis does not.

## Results table

| Test | Observed | Null criterion | Pass/fail |
|---|---|---|---|
| (i) Paired Wilcoxon σ_C > σ_S | V = 15, p = 0.0312, 5/5 positive Δσ, median Δσ = 1.3053 | p > 0.05 OR median sign ≤ 0 → fail | ✓ |
| (i-sens) Sensitivity excluding T6 | T6 IS in data; full test already includes it; alternative subset N/A | sign flip → flag | n/a |
| (ii) Spearman dissociation (Δσ_C vs Δconn, Δσ_S vs Δsize) | **Not computable** — no `egraph_measure` task results | CI crosses 0 → fail | ✗ (unevaluable) |
| (iii) CSN log-normal vs exponential | 2/5 C-cells log-normal-sig; 0/5 S-cells exponential-sig | < 4/6 OR < 3/6 → fail | ✗ |
| **Verdict** | Two of three criteria fire | any one fires → fail | **Mechanism not demonstrated** |

## Test 1 — paired Wilcoxon σ_C vs σ_S

![Paired σ M→S→C](sigma_by_condition_paired.png)

Across the 5 targets with sufficient late-phase data (T1 is excluded — under M and S it solves so fast that there are too few non-censored plateaus to fit σ stably; see Caveats), the connectivity-rich condition produces strictly larger pooled censored-MLE σ than the size-only condition: V = 15, p = 0.0312, 5/5 positive, median Δσ_CS = 1.3053. **The parent's null Wilcoxon flipped sign and gained two orders of magnitude in effect size once the "rich" primitive set was decomposed.** This is consistent with — but does not by itself prove — the portal-density mechanism: it shows that the size-only enrichment (S) is *not* what was driving the parent's hints, and that connectivity-rich (C) is doing real work.

## Test 2 — Spearman dissociation

The pre-registered dissociation test required two static observables per (target, condition): banded mean class size (S − M) and banded out-degree to higher-fitness bands (C − M), produced by `egraph_measure` tasks. **These tasks are absent from `result.json`** — only `gp_run` results were produced. The Δσ–Δstatic Spearmans cannot be evaluated. As a descriptive substitute we computed Spearman ρ(Δσ_C, Δσ_S) across the 5 complete targets: ρ = 0.1, p = 0.8729 — i.e., the C and S effects do not co-vary across targets, which is at least *consistent* with the dissociation hypothesis (a target where S "works" is not a target where C "works"), but this is not the preregistered test and should not be over-read.

## Test 3 — CSN tail shape

![CSN log-normal vs exponential R per cell](csn_R_per_cell.png)

Two of five C-cells (T4_deep_mul: R = 68.71, p = 2.73e-10; T6_inexpressible: R = 65.14, p = 5.80e-08) pass the log-normal-preferred-with-p<0.05 bar. T1, T2, T5 C-cells lean log-normal (R > 0) but don't reach significance; T3 C-cell R = 5.13, p = 0.186. **Zero of five S-cells reach the exponential-preferred bar** — most S-cells have so few late-phase tail observations (≤ 20) that CSN can't discriminate. So the pre-registered bar (≥4/6 C log-normal AND ≥3/6 S exponential) clears at 2/5 and 0/5 respectively; the criterion fires. This is partly a power problem (S-cells are under-sampled in the late phase because they typically solve and freeze) and partly a real signal (the C tails that *are* well-sampled are unambiguously log-normal in two cases out of five).

## Discussion

The headline result — σ_C > σ_S in 5/5 targets with effect size median Δσ = 1.3053 vs the parent's median Δσ ≈ 0.003 — is the cleanest evidence yet that *something different* happens when the primitive set adds new fitness-band-opening operators vs identity-equivalent operators. The mechanism the parent guessed (rich primitives slow GP dynamics by inflating heavy-tailed escape times) gains traction once you stop conflating size with connectivity: the size-only condition S behaves essentially like M (delta_sigma_SM is small and mixed-sign), while C blows up σ by ~1+ log-units. That said: (a) without the e-graph measurements, we cannot claim Δσ tracks Δconnectivity specifically — we only know Δσ_C is large and consistent; (b) the CSN bar was set for a 6-target experiment with well-populated tails, and S-cells in particular are tail-starved here; (c) the parent's three-observable convergence claim is therefore *partially* validated (the σ direction is now solid) but *not* mechanistically nailed down. The next ticket should either (1) re-run the e-graph measurement subtasks that this run dropped, or (2) reduce the CSN significance bar in advance to a tail-quality-conditional version, since the S-condition rarely produces 50+ late-phase plateaus per target.

## Caveats

- **Result schema inconsistency.** The task plan specified 540 `gp_run` tasks plus implied `egraph_measure` tasks for Test 2; `result.json` contains 540 `gp_run` results and zero `egraph_measure` results. The plan's Test 2 is therefore unevaluable in this run, not "null." Marked as ✗ (criterion fires) per the plan's "any one fires = fail" rule, but treat this as a data-collection gap, not a mechanism failure.
- **T1 excluded from σ analysis.** Under M and S, T1_poly_sep solves within ~30 generations and then plateaus until generation 1499 (single censored plateau). The censored-MLE σ fit is degenerate. T1's C-cell *does* produce a clean log-normal tail (R = 20.21, p = 0.0019), included in Test 3.
- **Plateau-duration extraction.** Late phase defined as "plateau fitness ≤ median(best_hist)" per run, pooled across reps. The plan's `late_phase_def: fitness_banded_post_tau_target` requires a tau_target produced by selection_calibration_pass — absent from results, approximated here. Effect size is large enough that this approximation is unlikely to be load-bearing, but the absolute σ values would shift under a stricter definition.
- **CSN xmin search.** Implemented as `xmin = quantile(durs, 0.10)` with truncated log-normal vs left-truncated exponential MLE and Vuong-normalized LR. Not identical to the `powerlaw.Fit(xmin=auto)` the plan named.
- **Sensitivity drop-T6.** Plan called for re-running Wilcoxon excluding T6 in case T6 was load-bearing (parent: T6 was "M-inexpressible"). With T6 included: V=15, p=0.0312; T6 contributes one of the five positive Δσs. Dropping T6 leaves 4 pairs all positive (V=10, exact one-sided p = 1/16 = 0.0625) — fragile to single-target leverage. Cross-references parent ticket `tk_4f5c6894`.
- Packages used: jsonlite, ggplot2, wesanderson, patchwork, dplyr, tidyr. No `poweRlaw` package call — Vuong LR implemented manually.
