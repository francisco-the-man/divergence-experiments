# tk_4f5c6894 — Primitive-set richness vs. plateau heavy-tailedness

## TL;DR

The mechanism **did not survive contact with the data**. Across 5 symbolic-regression targets, enriching the GP primitive set inflated the e-graph neutral-network proxy by ~10× as intended, but the late-phase plateau-duration σ shifted in **mixed directions** (rich > minimal on 2/5 targets, rich < minimal on 2/5, ≈ on 1/5). All three preregistered null criteria fail. The CSN comparison further shows log-normal is not even consistently preferred over exponential in the rich condition — only 1 of 5 targets supports it.

## Headline figure

![Paired sigma by target](plot_sigma_paired_by_target.png)

If the proposed mechanism (richer primitive set → larger neutral networks → heavier-tailed plateau durations) were operating, every line should slope up. Two slope up, two slope down, one is flat. This is the headline null.

## Results table

| Test | Observed | Null criterion | Pass/Fail |
|---|---|---|---|
| Paired Wilcoxon σ_late (rich > minimal, one-sided) | V=7, p ≈ 0.81 (direction mixed: 2↑/2↓/1≈) | p > 0.05 OR mixed direction ⇒ null | **✗ FAIL** |
| Spearman ρ(Δσ_late, NN-ratio) | ρ ≈ -0.10, 95% bootstrap CI spans 0 (≈ [-1, +1] with n=5) | CI crosses 0 ⇒ null | **✗ FAIL** |
| CSN log-normal vs exponential (rich, per target) | T1 R=-3.92 p<.001 (exp); T2 R=1.82 p=.07 (ns); T3 R=1.51 p=.13 (ns); T4 R=3.21 p=.001 (LN); T5 R=-2.95 p=.003 (exp) | Any target with LR not significant at α=0.05 favouring log-normal ⇒ null | **✗ FAIL** (3/5 ns or favour exp) |

Mechanism not demonstrated on any of the three preregistered axes.

## Per-test figures

### 1. Paired σ_late by target

![Paired sigma](plot_sigma_paired_by_target.png)

Lines connect the same target across conditions. Hypothesis predicts uniformly positive slopes. Observed: T2 and T4 go up, T1 and T5 go down, T3 ≈ flat. This is consistent with M2 in the plan (heavy tails are a generic GP-dynamics property, not modulated by the primitive-set-richness manipulation in a directional way) or with the manipulation interacting with target structure in a way the simple "more neutral = heavier tail" story doesn't capture.

### 2. Δσ vs e-graph NN-size ratio

![Shift vs NN ratio](plot_sigma_shift_vs_nn_ratio.png)

The manipulation hit the intended independent variable hard — every target's NN-proxy ratio is ~9-14× (rich/min). If neutrality drove tail heaviness, Δσ should scale with that ratio. It doesn't: the scatter is essentially flat with no monotone trend (Spearman ρ ≈ -0.1). This rules out M3-style "tails exist but track an unmeasured property correlated with primitives": the e-graph measurement is good and the relationship simply isn't there.

### 3. CCDF of late-phase plateau durations

![Tail CCDFs](plot_tail_ccdf_by_target.png)

Visual check on the σ summary. Tails are heavy in both conditions and the rich/minimal CCDFs are close to overlapping for T1, T3, T5 (where σ goes the "wrong" way or flat), with rich notably heavier only for T4. Note the prominent right-edge mass — many late-phase plateaus run to the generation cap (750 gens late phase). This censoring affects both conditions symmetrically but inflates the apparent tail.

## Discussion

The proposed mechanism — that primitive-set-induced neutral-network inflation causes heavier plateau-duration tails — is **not supported** at this scale. The manipulation cleanly inflated the e-graph NN proxy (the "M3 escape hatch" of an unmeasurable independent variable is closed: we measured it, it moved, σ didn't track it). The directional inconsistency across targets (T2/T4 up, T1/T5 down) suggests an **interaction with target structure** the original framing didn't anticipate: e.g. for separable polynomial targets (T1) the rich set's extra constants/identities may let runs *escape* plateaus faster via more exit-mutation paths, *shortening* tails rather than lengthening them. The log-normal-vs-exponential CSN comparison is also informative: for 2/5 targets exponential is the better fit in the rich condition, so "log-normal σ" is not even a well-grounded primary observable for those cells. A follow-up that scales up shouldn't just re-run this — it should first revisit whether σ is the right summary and whether the "more neutrality ⇒ heavier tails" link needs to be replaced by a signed model that accounts for plateau-exit dynamics.

## Caveats

- **n=5 targets is the binding power constraint.** The Wilcoxon and Spearman are weak by construction; we report them per pre-registration but the visual mixed-direction is the real evidence.
- **Right-censoring at 1500 generations** affects late-phase durations. ~30-40% of late plateaus hit the cap in both conditions; σ fits treat these as observed values.
- **CSN n_tail values (72-147) clear the 200-event bar only when pooled across phases**; the per-(target × condition) late-phase counts are below 200 for some cells, marginal for others. Scale-up would help here.
- **`final_best` distributions differ between conditions** (rich has more catastrophic failures on T1 — see `per_rep_final_best_mean`=0.124 vs 0.011). This means rich and minimal runs aren't exploring identical fitness regions, which complicates the "same target, just more neutrality" framing.
- **The e-graph "neutral_proxy" is a sampled mean class size**, not a true NN-size measurement. It's the same proxy for both conditions so the *ratio* is meaningful, but absolute values shouldn't be over-interpreted.
- T3 minimal runs barely make progress (`final_best` median 0.34 vs init ~0.34) — the minimal primitive set may be effectively unable to express `x0/(1+x1²)`, so the T3 comparison is closer to "stuck vs. searching" than "two valid searches with different landscapes."
