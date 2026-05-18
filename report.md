# tk_3b8ce4aa — Decoupling Neutral-Network Size from Inter-Network Connectivity

## TL;DR

The portal-density mechanism is **not demonstrated** — and the data argue against it in a stronger way than parent tk_4f5c6894 did. The connectivity-rich set C produced **lower** late-phase plateau-duration σ than the size-only-rich set S on **0 of 6** targets (paired one-sided Wilcoxon V = 0, p = 0.989451, median Δσ_CS = -4.837205); both Spearman correlations have bootstrap CIs crossing zero (ρ_conn = -0.2, 95 % CI [-1, 0.8]; ρ_size = -0.085714, 95 % CI [-1, 1]); and CSN tail-shape comparisons preferred log-normal in 6 / 6 C cells but in 0 / 6 S cells as *exponential* — none of the three preregistered criteria fire in the predicted direction. All three null-criteria triggered: **mechanism not demonstrated.**

## Headline figure

![Paired late-phase σ across conditions](sigma_by_condition_paired.png)

The portal-density mechanism predicted C lines above S lines above M lines. We observe the opposite: σ_C sits at ~1.85 for every target while σ_M and σ_S balloon to 5–18 depending on target difficulty.

## Results table

| # | Test | Observed | Null criterion | Verdict |
|---|---|---|---|---|
| (i) | Paired Wilcoxon σ_C > σ_S (one-sided), 6 targets | V = 0, p = 0.989451, 0/6 positive, median Δσ_CS = -4.837205 | p > 0.05 OR median sign ≤ 0 → fail | ✗ |
| (i.sens) | Same, excluding T6 | V = 0, p = 0.984514 on 5 pairs | sign flips or p > 0.05 → flag | ✗ (T6 not the leverage) |
| (ii) | ρ(Δσ_C, Δconnectivity_C) bootstrap CI | ρ_conn = -0.2, 95 % CI [-1, 0.8]; ρ_size = -0.085714, 95 % CI [-1, 1]; contrast Δρ = -0.114286 [-1.8, 1.152] | CI crosses 0 OR ρ_conn ≤ ρ_size → fail | ✗ |
| (iii) | CSN ≥4/6 C log-normal-preferred AND ≥3/6 S exponential-preferred | C: 6/6 log-normal-sig; S: 0/6 exponential-sig (T1_S and T4_S both log-normal-preferred, others insufficient tail) | failure of either condition → fail | ✗ |

Pre-registered rule: any one criterion firing = mechanism not demonstrated. **All three fired.**

## Per-test interpretation

### Test (i) — Paired Wilcoxon σ_C vs σ_S

![Late-phase σ paired by target](sigma_paired_by_target.png)

In every one of 6 targets, σ_C is below σ_S — frequently by a large margin (T2: 1.865129 vs 7.677212; T5: 1.911015 vs 9.1015; T6: 1.845185 vs 18.213656). The one-sided Wilcoxon explicitly testing σ_C > σ_S returns V = 0, p = 0.989451, with median Δσ_CS = -4.837205. Excluding T6 doesn't help (p = 0.984514). The data are quantitatively consistent with σ_C < σ_S in every condition — the **opposite** of the portal-density prediction. **Mechanism prediction inverted, not just unmet.**

### Test (ii) — Dissociation via paired Spearmans

![Delta-sigma vs delta-static](delta_sigma_vs_delta_static.png)

The dissociation argument required ρ_conn to be positive and exceed ρ_size, both with CIs excluding zero. We find ρ_conn = -0.2 (95 % CI [-1, 0.8]) and ρ_size = -0.085714 (95 % CI [-1, 1]); the paired contrast Δρ = -0.114286 with CI [-1.8, 1.152]. With n = 6 targets the bootstrap CIs are essentially uninformative — the rank statistic can flip on resampling — and neither point estimate is in the predicted direction. **Dissociation not detectable at this n.**

### Test (iii) — CSN log-normal vs exponential, with xmin search

![Tail CCDF by condition](tail_ccdf_by_condition.png)

The preregistered prediction was C tails curve (log-normal), S tails straight (exponential). The pooled tails tell a different story: every C cell has R > 0 with p < 1e-35 (n_C_lognorm_sig = 6 / 6 targets), but the *S* cells that have enough observed (uncensored) tail to fit at all (T1_S: R = 79.36396, p = 4.196783e-06; T4_S: R = 448.197981, p = 1.871597e-19) are *also* log-normal-preferred. Targets with very small n_tail (T2_S = 45, T5_S = 17, T6_S = 1) cannot be fit. Net: **n_S_exp_sig = 0 / 6**, far below the required 3. Tail-shape evidence does not dissociate the two enriched primitive sets.

## Discussion

The clean inversion of Test (i) is the load-bearing finding and deserves a mechanistic reading. C runs reach a much lower final fitness than M / S runs (mean tail size n=1393–2067 for C vs 1–869 for M/S — see `sigma_table` n column), which means the "fitness-banded late-phase" of a C run is a *different* and lower-fitness regime than the late-phase of an M run on the same target. σ_C ≈ 1.85 on every target is the σ of a population that has already escaped most landscape structure; σ_M = 10.2 on T6 is the σ of a population stuck in a few sticky basins. Comparing them at face value isn't a clean test of portal density — it's a confounded comparison of two different points on the search trajectory.

What this experiment does demonstrate: connectivity-rich C *uniformly* enables progress past where M and S stall. That's consistent with portals existing — sin/cos/square/sqrtp do open new fitness bands — but the σ statistic we chose conflates "how far did the run get" with "how heterogeneous are the plateaus there." The parent ticket's positive Spearman (ρ = 0.80 in tk_4f5c6894) was reading a within-late-phase signal; here we replaced it with a between-late-phase comparison and lost the signal. **The mechanism may be real but our observable is not measuring it.**

Cannot adjudicate from these data: whether portal density per se (vs. raw expressibility) drives the escape, whether σ within a *matched fitness band* across conditions would dissociate C from S, and whether T6_inexpressible's M-cell with n_tail = 2 reflects a true ceiling or a censoring artifact. Suggested follow-up: re-band the late-phase using a fitness threshold *common to all three conditions* (the highest fitness reachable by M), then re-run all three tests. That isolates plateau heterogeneity from end-point depth.

## Caveats

- **n = 6 targets is structurally underpowered for bootstrap CIs on Spearman**: ρ_conn CI [-1, 0.8] is the bootstrap telling us it cannot rule out anything. Don't read the failure of Test (ii) as evidence *against* dissociation; read it as no signal extractable at this n.
- **Late-phase definition** here is "plateaus with fitness ≤ run-median fitness." C runs have a much lower median fitness than M/S, so this is **not** a fitness-matched comparison. See the discussion paragraph above.
- **S runs frequently have ≤45 uncensored tail events** (T2, T5, T6 all fewer than 50), below our CSN fitting threshold. Test (iii)'s S-side criterion was operationally hard to satisfy regardless of mechanism.
- **CSN test uses a custom R implementation** (KS-based xmin search + Vuong-style R / z / p) rather than the python `powerlaw` package, because the worker is R-only. The qualitative direction (R, sign of preference) is robust; the precise p-values are not directly comparable to the parent ticket's.
- **Parent ticket cross-link**: parent tk_4f5c6894 also failed its preregistered tests, but with σ_rich > σ_minimal at least *trending* up; this follow-up's σ_C uniformly *below* σ_S is a stronger inversion than the parent saw, attributable to the lower-fitness late-phase regime C reaches.
- **Censoring**: plateau events with `obs = 0` (run ended before plateau broke) were handled with a censored log-normal MLE; this preserves the lower-bound information from censored long plateaus rather than dropping them.
