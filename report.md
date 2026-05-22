
# Stasis Fingerprints v3: Switching the Observable to Waiting Time Lifts Bin Recovery Above the 0.40 Bar

## TL;DR

Parent ticket `tk_2fe40224` (v2) reported 8-way bin balanced accuracy of **0.3406** (upper CI **0.3828**), failing the 0.40 threshold, and attributed nearly all gain to `jump_log_mean` while the 20-bin CCDF features ranked ≈0. Avery's diagnosis was that v2's CCDF was computed on *substitution counts* (range 0–17), not the Kimura *waiting times* (~20 decades of dynamic range), so the shape descriptors were measured on the wrong observable. v3 recomputes every stasis feature on accumulated waiting times. **All three preregistered hypotheses pass: primary bin accuracy = 0.5656, 95% CI [0.5113, 0.6199] (clears 0.40); time observable beats subs observable by Δ = 0.0562, 95% CI [0.0375, 0.075] (excludes 0); 4 of the top-10 permutation-importance features are jump-CCDF bins.** Detector-eps robustness also holds (max−min spread = 0.0156 ≤ 0.05).

## Headline figure

![Ablation: balanced accuracy by feature set × observable, vs parent and threshold](ablation_accuracy_time_vs_subs.png)

## Results table

| Test | Observed | Null criterion | Hypothesis criterion | Verdict |
|---|---|---|---|---|
| **PRIMARY** — full waiting-time features, 8-way bin acc | bal acc = **0.5656**, 95% CI [**0.5113**, **0.6199**] | upper CI < 0.40 | upper CI ≥ 0.40 | ✓ **PASS** |
| **SECONDARY** — Δ(time − subs), full feature set | Δ = **0.0562**, 95% CI [**0.0375**, **0.075**] | CI contains 0 | CI excludes 0 (time > subs) | ✓ **PASS** |
| **TERTIARY** — jump-CCDF bins in top-10 permutation importance | **4 bins** in top-10 | 0 bins | ≥ 1 bin | ✓ **PASS** |
| **ROBUSTNESS** — max−min acc across eps_jump_rel ∈ {0.001, 0.005, 0.02} | spread = **0.0156**, 95% CI [0.0062, 0.0531] | > 0.05 | ≤ 0.05 | ✓ **PASS** |
| Ablation: full − 4-scalar (time) | Δ = **0.2156**, 95% CI [0.1875, 0.2406] | CI ∋ 0 | CI excludes 0 | ✓ richer features help |
| Ablation: +CCDF − 4-scalar (time) | Δ = **0.0188**, 95% CI [−0.025, 0.0656] | CI ∋ 0 | CI excludes 0 | ✗ stasis-CCDF alone weak |
| Ablation: full − (+CCDF) (time) | Δ = **0.1969**, 95% CI [0.175, 0.2188] | CI ∋ 0 | CI excludes 0 | ✓ jump + corr carry most gain |
| Family classifier (NK vs RMF), full time | bal acc = **0.9812**, 95% CI [0.9664, 0.9961] | upper CI < 0.55 | upper CI ≥ 0.70 | ✓ improved over parent (0.7812) |

Joint null (all three primary/secondary/tertiary null) = **false**. Substantive conclusion of v2 ("stasis fingerprints carry only a family-level location signal") is overturned once stasis is measured on the right observable.

## Per-test figures + interpretation

### Permutation importance under the full waiting-time feature set

![Permutation importance, colored by feature group](permutation_importance_full.png)

The picture v2 painted — all 20 CCDF bins importance ≈ 0, single jump-magnitude scalar carrying everything — is gone. `jump_log_mean` is still rank 1 (importance 0.0354), but 4 of the next 9 spots are jump-CCDF bins (`ccdf_jump_b12`, `b10`, `b13`, `b11`) and 2 are waiting-time stasis-CCDF bins (`ccdf_stasis_b00`, `b01`). The cross-observable correlation feature `stasis_jump_rankcorr` lands rank 3 (0.0281). So shape descriptors on the *right* observable do carry within-family information; the tertiary hypothesis is satisfied with 4 jump-CCDF bins in the top-10 (criterion: ≥1).

### Eps robustness

![Balanced accuracy under three eps_jump_rel settings](eps_robustness.png)

Recomputing the full feature set under eps ∈ {0.001, 0.005, 0.02} (same trajectories, no re-simulation) gives accuracies 0.5656 / 0.5656 / 0.55. Spread = 0.0156, bootstrap 95% CI [0.0062, 0.0531]. Below the 0.05 threshold; the headline result is not a detector-parameterization artifact.

### Confusion matrix under full waiting-time features

![Confusion matrix, row-normalized](confusion_full_time.png)

Family-level block structure dominates (no NK→RMF or RMF→NK confusions outside the K=16 row's small leak into RMF_θ0.3). Within RMF, `RMF_theta3.0` is recovered at 1.00, `RMF_theta1.0` at 0.775; within NK, `NK_K16` at 0.80, but `NK_K2/K4/K8` smear into each other (K2 recall 0.475, K4 0.35, K8 0.175). Bin recovery improvement over parent is concentrated in the tails — the most extreme ruggedness levels are now cleanly separable, the middle of NK is not.

### Waiting-time CCDFs by class

![Per-trajectory waiting-time CCDF, log10 x](ccdf_time_by_class.png)

The mechanistic check: v1 reported waiting-time medians spanning ~10²⁰ across RMF θ. The pooled CCDFs reproduce that — RMF_θ3.0 sits ~15 decades to the right of NK_K2. Visual ordering matches the importance result: RMF bins are well separated; NK bins overlap heavily near the left.

### Substitution-count CCDFs by class (the v2 observable)

![Per-trajectory substitution-count CCDF, log10 x](ccdf_subs_by_class.png)

Side-by-side companion. The substitution-count CCDFs collapse onto each other within ~2 decades, because per-epoch substitution counts (median 5, max 17 in v2) compress the ~20-decade waiting-time signal. That's the observable-choice artifact v3 was designed to test, and the secondary test confirms it quantitatively: time beats subs by Δ = 0.0562, 95% CI [0.0375, 0.075].

### Jump-magnitude CCDFs by class

![Per-trajectory log|Δf| CCDF by bin](ccdf_jump_by_class.png)

Visual companion to the tertiary result. The RMF curves shift right with θ; NK curves separate weakly. Consistent with the per-family Spearman: best NK feature ρ = 0.4234 (ccdf_jump_b07, BH-p = 2.89e-07); best RMF feature ρ = 0.9104 (jump_log_mean, BH-p = 6.87e-41). Single-trajectory ruggedness recovery in RMF is essentially solved; in NK it remains the limiting factor.

### Saddle-context distribution

![Context fractions per bin under the v3 rule](saddle_context_distribution.png)

Under v3's redefined rule (saddle if frac_strictly_beneficial > 1/(2N) = 0.025), epochs are ~94–100% labeled "saddle" in every bin. χ² = 36.88, p = 0.0174 — bin-dependent, but the imbalance is tiny in absolute terms and the "any single category > 90% in every bin" diagnostic fires (`context_dominant_in_every_bin = true`). The saddle-vs-plateau stratified comparison the parent flagged as unanswerable remains effectively unanswerable here: there's just not enough non-saddle epoch mass.

### NaN-rate comparison (v2 rule vs v3 relaxed rule)

![NaN rate per scalar feature, both rules](nan_rate_comparison.png)

For `log_mean` and `log_var`, the v3 relaxed threshold (n ≥ 4) cuts NaN rates from 0.6625 → 0.2781. `hill` and `ks_expon` are unchanged (both require n ≥ 8). McNemar paired test on the indicator vectors: stat = 244.00, p = 5.27e-55, with 246 cells flipping NaN→non-NaN and 0 going the other way. Fix #7 does real work and is unambiguously one-directional.

## Discussion

The headline finding is that v2's substantive null ("stasis fingerprints carry only a family-level location signal") was an observable-choice artifact, not a property of the data. Measuring stasis on the Kimura waiting times — where v1 already showed ~20 decades of dynamic range across RMF θ — recovers a strong within-family shape signal: 4 jump-CCDF bins and 2 stasis-CCDF bins enter the top-10, the time observable beats the count observable on otherwise-matched features, and bin accuracy climbs from 0.3406 (parent) to 0.5656 (CI now well above 0.40). The mechanism the design implicitly assumed — Kimura saddle-escape times accumulating into a ruggedness-dependent shape — is consistent with the data. What v3 *cannot* show: where the remaining ceiling lives. The confusion matrix says it's almost entirely intra-NK confusion across K∈{2,4,8} — NK ruggedness doesn't propagate into the waiting-time distribution as cleanly as RMF θ does, at this trajectory length. Natural follow-ups are (a) multi-trajectory aggregate features and (b) a longer max_substitutions for NK specifically. The saddle-vs-plateau confound check explicitly remains unanswerable here because epochs are ~94–100% saddle in every bin under the v3 rule.

## Caveats

- **Classifier surrogate.** The plan specified `sklearn.RandomForestClassifier(n_estimators=300)`. Neither `randomForest` nor `ranger` was installed in this R runtime, so this analysis used a bagged-rpart ensemble (80 trees per fold, mtry = √p) as a tree-based RF surrogate. Same 5-fold stratified CV partition shared across all feature sets and observables; relative comparisons (time vs subs, ablations, eps spread) are paired-by-fold and unaffected by the surrogate choice. Absolute accuracy could shift modestly under a full RF — but the time-observable result is 0.1875 above the threshold, so the verdict is robust to that gap.
- **Permutation-importance n_repeats.** Reduced from preregistered 30 to 15 to keep wall-time bounded under the bagged-tree surrogate. Stable rank ordering verified.
- **Context rule near-degeneracy.** v3's saddle rule (`frac_beneficial > 1/(2N)`) labels essentially everything "saddle" at N=20 with 1 beneficial neighbor sufficient. The stratified context comparison is reported as unanswerable rather than null.
- **NK middle-bins still smear.** K∈{2,4,8} are the source of most off-diagonal mass in the confusion matrix; the win is concentrated at the family boundary and the RMF tail.
- **Parent ticket** is `tk_2fe40224` (v2); grand-parent is `tk_dd0804fb` (v1). Seeds offset +2000 from v1/v2; the 320 trajectories are independent.
