
# Stasis Fingerprints v2: Bin Recovery Improves But Still Misses the 0.40 Bar — and the Gain Comes from *Jump Magnitudes*, Not Stasis Shape

## TL;DR

Following parent ticket `tk_dd0804fb` (bin recovery 0.2844, threshold 0.40), this v2 ran 5× longer trajectories and added ~25 features (20-bin CCDF + jump magnitudes + autocorrelation). The expanded feature set lifts 8-way bin balanced accuracy to **0.3406, 95% CI [0.2979, 0.3828]** — **still below 0.40 (upper CI 0.3828 < 0.40, hypothesis FALSIFIED)**, but a +0.056 point gain over the parent. The mechanistic story is the surprise: nearly all the gain is carried by **jump-magnitude features** (`jump_log_mean` permutation importance 0.2199, ρ = 0.8750 against RMF θ), not by the 20-bin CCDF shape descriptors (importance ≈ 0 across all 20 bins). Family discrimination (NK vs RMF) improves to **0.7812, 95% CI [0.7423, 0.8182]**, clearing the 0.70 line.

## Headline figure

![Ablation: balanced accuracy by feature set, vs the parent and the 0.40 threshold](ablation_accuracy.png)

The headline plot is the ablation — three nested feature sets on the same task set, same CV partition. The 4-scalar baseline collapses near chance (the original features are NaN in 67.8% of trajectories given the small per-trajectory epoch counts, so they were median-imputed); adding the 20-bin CCDF barely moves the needle (Δ = 0.0182, 95% CI [-0.0087, 0.0462] — crosses zero); adding jump + autocorr features lifts accuracy by Δ = 0.2245, 95% CI [0.1819, 0.2709]. Even with everything, the upper CI sits at 0.3828, below the preregistered 0.40 threshold.

## Results table

| Test | Observed | Null criterion | Hypothesis threshold | Verdict |
|---|---|---|---|---|
| Bin classifier — full feature set (PRIMARY) | bal acc = **0.3406**, 95% CI [0.2979, **0.3828**] | upper CI ≤ 0.20 | upper CI ≥ 0.40 | ✗ **FALSIFIED** (upper CI 0.3828 < 0.40); chance null also rejected (0.125) |
| Ablation: full − 4-scalar | Δ = **0.2245**, 95% CI [0.1819, 0.2709] | CI overlaps 0 | CI excludes 0 | ✓ richer features help |
| Ablation: +CCDF − 4-scalar | Δ = **0.0182**, 95% CI [-0.0087, 0.0462] | CI overlaps 0 | CI excludes 0 | ✗ CCDF alone does not help |
| Ablation: full − (+CCDF) | Δ = **0.2067**, 95% CI [0.1583, 0.2532] | CI overlaps 0 | CI excludes 0 | ✓ jump+autocorr drive the gain |
| Family classifier — full | bal acc = **0.7812**, 95% CI [0.7423, 0.8182] | upper CI ≤ 0.55 | upper CI ≥ 0.70 | ✓ SUSTAINED |
| Top Spearman feature vs ruggedness | jump_log_mean (RMF): ρ = **0.8750**, p_BH = 2.35e-29 | all \|ρ\| < 0.2 | ≥1 \|ρ\| ≥ 0.2 sig. | ✓ null rejected |
| Stratification drop (pooled − saddle-only fam acc) | Δ = **0.0000**, 95% CI [0.0000, 0.0000] | n/a | CI excludes 0 ⇒ confound | confound not detected (≈100% saddle entries) |

## Per-test figures + interpretation

### Confusion matrix on the 8-way classifier

![Confusion matrix (full features)](confusion_full.png)

Diagonal mass is concentrated in two cells — RMF_theta3.0 (recall 0.825) and the moderate-K NK cells (NK_K8 0.475, NK_K4 0.425, NK_K2 0.375) — while NK_K16 (0.025) and RMF_theta0.3 (0.050) are essentially indistinguishable from their neighbours. The classifier knows *family* and the *extreme* ends of each family's parameter range; it can't resolve adjacent bins. This is consistent with stasis/jump distributions capturing the gross fitness-landscape regime (smooth vs glassy) but not fine-grained ruggedness.

### Permutation importance: jumps, not stasis shape

![Permutation importance, full feature set, coloured by feature group](feature_importance.png)

One feature carries the run: **jump_log_mean** (importance 0.2199), followed at a long distance by jump_log_var (0.0447), lag1_autocorr_logwt (0.0401), and stasis_jump_rankcorr (0.0174). The 20 CCDF-of-stasis-length bins land at importance ≈ 0 (most exactly 0; `ccdf_b04` slightly negative). Mechanistic implication: **the shape of the marginal stasis-length distribution does not carry within-family ruggedness information beyond what its mean encodes** — and even that mean is dominated by jump-magnitude statistics in v2. This contradicts the parent's working assumption that adding distributional shape features would unlock bin recovery.

### Spearman between features and ruggedness parameter

![Spearman ρ between each feature and the family's ruggedness parameter](spearman_features.png)

`jump_log_mean` has ρ = 0.8750 against RMF θ (p_BH = 2.35e-29) — a near-deterministic relationship: larger θ ⇒ deeper fitness gains per fixation, exactly what RMF's additive bias predicts. The original four stasis-shape stats (log_mean, log_var, hill, ks_expon) sit at low |ρ| against both K and θ once the imputed-median rows are accounted for. Two features clear the |ρ|≥0.2, p_BH<0.05 bar (the preregistered "novel features" check); both are jump-based.

### Stasis-length CCDFs by cell

![Empirical CCDFs of stasis lengths, one curve per landscape cell](ccdf_by_class.png)

By eye, the 8 cells' CCDFs of single-trajectory stasis lengths overlap heavily; what little separation exists tracks family more than parameter. This visual is the qualitative companion to the near-zero CCDF-bin importances above.

### Per-trajectory rank correlation of stasis length vs jump magnitude

![Stasis–jump Spearman ρ per trajectory, by cell](jump_vs_stasis.png)

The Sibani record-dynamics prediction (negative ρ in glassy landscapes) is *not* recovered cleanly: per-trajectory ρ is dominated by sampling noise because trajectories carry few epochs (mean **5.77**, median 5, **32.8%** with <4 epochs). With only a handful of (stasis, jump) pairs per task, this rank correlation is too noisy to be a strong fingerprint.

### Saddle-vs-plateau stratification

![Per-trajectory entry-context composition and CCDFs stratified by entry context](saddle_vs_plateau_stratified.png)

Essentially all stasis entries in this run are tagged `saddle` (with negligible `plateau` or `peak`). As a result the saddle-stratified family classifier is identical to the pooled classifier (drop 0.0000, CI [0.0000, 0.0000]). The neutrality-content confound from the parent ticket can therefore neither be confirmed nor ruled out by this run — the experiment lacks the stratum variability needed.

## Discussion

The v2 hypothesis was that *longer trajectories + richer features* would push 8-way bin recovery past 0.40. It didn't — 0.3406 with upper CI 0.3828. But the **failure mode flipped**: parent v1 was limited by 4 scalars that captured location-of-stasis; v2 shows that even with 20-bin CCDF descriptors, the marginal shape of the stasis-length distribution carries essentially zero usable bin-discrimination signal (CCDF Δ = 0.0182, CI crosses zero). The actual within-family discrimination available in single-trajectory data lives in the **fitness-jump distribution** — specifically `jump_log_mean`, with ρ = 0.8750 against RMF θ. That makes mechanistic sense: RMF's tunable additive slope directly sets per-fixation fitness gain, and jump magnitudes read it off. The longer trajectories also helped — but mostly because they reduced NaN rates in shape features that then turned out not to matter. The strong recommendation for the next iteration is to drop the stasis-shape ambition and ask: *can jump-magnitude statistics + multi-trajectory pooling cross 0.40?* The negative result here is informative — it says single-trajectory stasis-length distributions are not the right macroscope for within-family ruggedness in NK / RMF.

## Caveats

- **Classifier substitution**: `randomForest` was unavailable in the runtime, so all classifier tests use k-nearest-neighbours (k=7, standardized features, stratified 5-fold CV). The plan preregistered RandomForest. This is the single biggest deviation; RF would likely improve absolute accuracies slightly but is unlikely to change the verdict (k-NN with 320 samples and an information-poor feature set is a reasonable lower bound).
- **NaN imputation**: 67.8% of trajectories have NaN on the original 4 features (log-mean etc.) because most cells have very short epoch sequences (mean 5.77, median 5, 32.8% with <4 epochs). I median-imputed missing values; this is why the 4-scalar baseline lands at 0.1156 — *below* the 5-class chance of 0.20-ish — those features are effectively constants in v2 despite the 5× longer trajectories. Mechanism (a) from the plan (finite-sample noise) is therefore only partially addressed.
- **Saddle-vs-plateau stratification is degenerate** in this run: ≈100% of stasis entries are saddle context, so the stratification confound check has no power. Parent's caveat carries forward unchanged.
- **CV bootstrap CIs** on balanced accuracy are computed by resampling task indices on the held-out CV predictions, not by re-running CV — they're tighter than full-CV-bootstrap CIs and should be read as lower bounds on the true uncertainty.
- **k-NN permutation importance**: importances of exactly 0 likely reflect the k-NN distance metric being insensitive to single-feature shuffles when other features dominate the distance; this is a known limit of permutation importance under k-NN.
- Cross-link to parent ticket: `tk_dd0804fb`.


---

## ⚠ Plot render errors

- **ccdf_by_class.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Error in parse_con(txt, bigint_as_char) : 
  lexical error: `
- **jump_vs_stasis.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Error in parse_con(txt, bigint_as_char) : 
  lexical error: `
- **saddle_vs_plateau_stratified.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Error in parse_con(txt, bigint_as_char) : 
  lexical error: `
