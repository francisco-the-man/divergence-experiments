# tk_4f5c6894 — Primitive-set richness vs. plateau heavy-tails

## TL;DR

The preregistered mechanism prediction **failed on all three criteria**. Across 5 targets, enriching the GP primitive set did **not** systematically increase late-phase log-normal σ (paired Wilcoxon one-sided p = 0.500, direction mixed: 2 of 5 targets shifted in the predicted direction). The per-target Δσ vs. e-graph neutral-network-size ratio shows a near-zero, sign-ambiguous Spearman correlation with a 95% bootstrap CI that brackets zero by a wide margin. The CSN log-normal-vs-exponential comparison on rich-condition late durations is significant for only 2 of 5 targets. Mechanism (target structure → neutral-network size → plateau heavy-tailedness via primitive-set lever) is **not demonstrated** at this scale.

## Headline figure

![Paired σ by target](plot_sigma_paired.png)

The within-target lines should all slope *up* (minimal → rich) under the hypothesis. They don't: 3 of 5 slope down, and the two that slope up (T2, T4) don't correspond to the largest NN-ratio increases. The mechanism is inconsistent with the paired view.

## Preregistered tests — results table

| Test | Observed | Null criterion | Pass/Fail |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided) | V = 6, p = 0.500; signs +/−/−/+/− | p > 0.05 or mixed direction → null | **✗ FAIL** (null met) |
| Spearman ρ(Δσ, NN-ratio) | ρ = −0.100, bootstrap 95% CI [−1.000, 0.900] | CI crosses 0 → null | **✗ FAIL** (CI crosses 0) |
| CSN log-normal vs exponential (rich, per target) | Significant (p<0.05) in 2/5 targets (T1, T5); T4 favors exponential | Not significant in any target → null | **✗ FAIL** (only partial: 3/5 not significant; T4 R>0 wrong direction) |

All three preregistered null criteria are met → **mechanism not demonstrated**.

## Per-test figures

### Spearman: Δσ vs. neutral-network ratio

![Δσ vs NN-ratio](plot_sigma_shift_vs_nn.png)

Mechanism predicted positive correlation: targets whose neutral networks expand most under the rich set should show the largest σ inflation. Instead, Δσ is essentially uncorrelated with the e-graph NN-size ratio. NN ratios are large and consistent across targets (~10–23×), but σ shifts are small and sign-mixed. Even within this small n=5 paired design, the *direction* of effect would have been informative — it isn't there.

### Pooled late-phase CCDF by target

![CCDF by target](plot_ccdf_by_target.png)

Visual check on heavy-tailedness. All conditions show right-skewed plateau-duration distributions with a pile-up at the 750-generation censoring boundary (late-phase window cap). Within targets, the rich vs. minimal CCDFs are visually similar — there's no systematic rightward shift of the rich curves. T4 is the only target where the rich curve is meaningfully heavier-tailed than minimal; T1 and T5 go the other way.

### Per-target σ and NN-size proxy

![σ and NN-size summary](plot_sigma_and_nn.png)

Left: pooled σ_late by target × condition (heights nearly equal within target). Right: e-graph mean weighted class size — the manipulation *did* work as intended (rich condition produces ~10–25× larger neutral classes across all targets), so the null is not from a failed manipulation. The independent variable moved; the dependent variable didn't.

## Discussion

The cleanest reading: **enlarging neutral-network size by 10–25× (verified by independent e-graph measurement) did not measurably inflate plateau-duration heavy-tailedness**. This is evidence against the specific mechanism proposed — that plateau heavy-tails on symbolic-regression GP are governed by induced-landscape neutral-network structure, with primitive-set richness as the controlling lever. Two alternatives remain open: (M2) plateau heavy-tails are generic GP-dynamics artifacts (e.g., from finite population + truncation selection) unrelated to landscape neutrality — consistent with our data; or the lever is just wrong — primitive-set richness inflates *measured* e-graph class sizes but those classes may not correspond to the dynamically-relevant neutral sets the GA actually traverses. A scaled-up rerun (tk_9bdd20b5) won't rescue this: the n=5 paired Wilcoxon power is set by target-level effect size, not replicate count, and the effect direction is already mixed at scale-down. Suggested follow-up: directly measure *dynamic* neutrality (fraction of accepted mutations that are fitness-neutral during runs) rather than static e-graph class size, and re-test the correlation.

## Caveats

- **n=5 targets**: bootstrap CI on Spearman is necessarily wide; we can't distinguish "ρ≈0" from "weak ρ". But the direction is wrong/mixed, which is the more damning observation.
- **Censoring at 750 generations**: many late-phase plateaus hit the run-length cap. Log-normal σ is fit on right-censored data; true σ may be larger in all cells uniformly. Doesn't bias the *difference* much, but quoted σ values are conservative.
- **T3 (rational target) is hard**: minimal condition essentially never solves it (final MSE ~0.35); the "plateaus" there may be a different dynamical regime than the other targets.
- **Neutral-proxy is a static e-graph measure on randomly sampled trees**, not the dynamically-reachable neutral set under GP operators — this is exactly the gap the discussion flags as a follow-up.
- **Scale-down**: pop_size=80, max_gen=1500. CSN tail-event bar (≥200) is **not** cleared per-cell at this scale (n_late_events ranges 72–147). The pre-registration claimed ~1200 events per cell; actual was ~5–10× smaller because late-phase events are far rarer than projected. This weakens the CSN test specifically but not the σ-comparison test.
