# Stasis-Length Distributions as a Landscape Fingerprint — Partial Signal

## TL;DR

Across 320 origin-fixation trajectories on NK (K ∈ {2, 4, 8, 16}) and RMF (θ ∈ {0.1, 0.3, 1.0, 3.0}) landscapes, four distributional summary statistics of stasis-length distributions distinguish landscape **family** (NK vs RMF) at balanced accuracy 0.7062 [95% CI 0.6563, 0.7562] — right at the preregistered hypothesis threshold of 0.70 — but **fail** to recover the within-family ruggedness bin (balanced accuracy 0.2844 [0.2349, 0.3338], threshold 0.40, chance 0.125). Spearman rank correlations between summary stats and ruggedness parameters are significant and moderately large for RMF (max |ρ| = 0.6382) but weaker for NK. The classifier signal **survives** restricting to saddle-context stasis (drop 0.1031, bootstrap 95% CI [-0.0628, 0.1695] crosses zero), so the family signature isn't just a neutrality-content artifact.

## Headline figure

![Four-panel summary](headline.png)

The status is **PARTIAL / MIXED** (`stats$hypothesis_status`): family discrimination passes, parameter-bin recovery does not.

## Results table

| Test | Observed | Null criterion | Hypothesis threshold | Verdict |
|---|---|---|---|---|
| Family classifier (NK vs RMF) | balanced acc = 0.7062, 95% CI [0.6563, 0.7562] | upper CI ≤ 0.55 | upper CI ≥ 0.70 | ✓ hypothesis sustained (upper CI 0.7562 ≥ 0.70); null rejected |
| Cell classifier (8-way) | balanced acc = 0.2844, 95% CI [0.2349, 0.3338], chance = 0.125 | upper CI ≤ 0.20 | upper CI ≥ 0.40 | ✗ hypothesis falsified (upper CI 0.3338 < 0.40); above chance but null not met either |
| Spearman summary-stat vs ruggedness | max \|ρ\| = 0.6382 (RMF log-mean), any \|ρ\| ≥ 0.2 with p<0.05: TRUE | all \|ρ\| < 0.2 | ≥1 \|ρ\| ≥ 0.2 sig. | ✓ null rejected |
| Stratification drop (pooled − saddle-only) | 0.7062 − 0.6031 = 0.1031, bootstrap 95% CI [-0.0628, 0.1695] | n/a (diagnostic) | CI excludes 0 ⇒ neutrality confound | confound **not** confirmed: CI includes 0 |

## Per-test interpretation

### Family discrimination works — barely

![Family classifier ROC-style summary and CIs](family_classifier.png)

NK vs RMF classification reaches balanced accuracy 0.7062 with the lower CI bound (0.6563) safely above the 0.55 null threshold and the upper bound (0.7562) above the 0.70 hypothesis threshold. The point estimate sits essentially **on** the preregistered "success" line — a result the plan called "supported" only because of CI inclusion. Mechanistically this is consistent with M1 ("topology imprints on stasis distribution") for the family axis: per-cell aggregates differ markedly — NK log-mean ranges 13.8 → 36.5 across K, while RMF log-mean spans 35.5 → -0.08 as θ increases.

### Parameter-bin recovery fails

![8-class confusion matrix](confusion.png)

The 8-way classifier reaches 0.2844 balanced accuracy — roughly **2.3× chance** (0.125) and statistically distinguishable from chance, but well below the preregistered 0.40 hypothesis threshold (upper CI 0.3338 < 0.40). The confusion matrix shows family blocks (NK rows confused mostly with NK columns, RMF with RMF) but adjacent ruggedness bins are routinely mixed. Mechanistically this favors M2 ("family-level structure registers; ruggedness-axis fine structure does not") over M1.

### Summary stats correlate with ruggedness — asymmetrically

![Spearman correlations per family and stat](spearman.png)

For **RMF**, log-mean shows ρ = -0.6382 (p = 7.07e-08, n = 58) — a strong, highly significant decrease in mean log-stasis as θ increases (RMF becomes more single-peaked, escape rate rises, stasis shrinks). log-var also correlates (ρ = -0.3595, p = 0.0056). Hill and KS-to-exponential show essentially no signal in RMF. For **NK**, the picture is weaker and more mixed: log-mean rises with K (ρ = 0.416, p ≈ 0, n = 91), but the other three stats hover at |ρ| < 0.2 with p > 0.08. This explains the bin-classifier ceiling: log-mean carries most of the K/θ signal, but a single feature can't separate four bins within a family at 0.40 accuracy.

### Stratification: neutrality confound not confirmed

![Pooled vs saddle-only classifier accuracy](stratification.png)

Restricting to saddle-context stasis (the lit-reviewer's recommended stratum) drops family balanced accuracy from 0.7062 to 0.6031, a delta of 0.1031. The bootstrap 95% CI on the drop is [-0.0628, 0.1695] — it **includes zero**. We can't conclude the pooled-distribution signal is a neutrality-content artifact (M2's strong form): saddle-only stratum still discriminates families above chance, and the observed drop is within bootstrap noise. The conservative reading is that ~10 percentage points of family-discrimination signal *might* come from neutral-plateau content, but we lack power to confirm.

## Discussion

The data say: **landscape family leaves a faint but real fingerprint on origin-fixation stasis-length distributions, while ruggedness-parameter bins within a family do not — at least not via these four summary statistics on N=20, 400-substitution trajectories.** The mechanism in the plan (saddle-escape rates inheriting from k-distributions) predicts a monotone trend in mean log-stasis with ruggedness; that prediction is confirmed (NK ρ = 0.416, RMF ρ = -0.6382), but the higher-moment / tail features (Hill, KS) don't carry independent signal, capping multi-bin discrimination well below 0.40. For inverse inference from real punctuated trajectories (paradigm shifts, speciation), this is sobering: a 70% NK-vs-RMF readout from one trajectory of ~14 epochs is roughly the practical ceiling here. The natural follow-up is to either (a) use more samples per trajectory (longer runs) or (b) engineer richer features (full empirical CDF, jump-magnitude joint distribution).

## Caveats

- **Saddle fraction is high but variable** (mean across runs = 0.7469): most stasis epochs are saddle-context, so the "pooled" and "saddle-only" feature sets share most of their data. This biases the stratification drop toward zero by construction.
- **Time-scale extremes**: median stasis-length-in-time spans ~10⁰ (RMF θ=3.0: median 1.17) to ~10²⁰ (RMF θ=0.3: median 3.02e+20). Logarithm helps, but the dynamic range strains any classifier; some runs have inverse-fixation-probabilities that diverge numerically.
- **One trajectory per landscape draw, 400 substitutions, N=20**: small genome, modest epoch counts (mean 19.975 for NK_K16 down to 7.725 for RMF_theta0.1). Trajectories with very few epochs yield noisy summary stats — the per-cell ns vary from 309 (RMF_theta0.1) to 5296 (NK_K4) pooled stasis events.
- **Class imbalance was a coding concern earlier in the pipeline but resolved**: all 8 cells × 40 reps are present in the final result.json (320 rows).
- **Logistic regression** was used for the binary family classifier and **multinomial logistic** (via `nnet::multinom`) for the 8-way cell classifier; a tree-based learner might extract more from these features, but that's outside the preregistered plan.
- **No NaN/Infinity scrubbing in the original JSON**: `result.json` contained bare `NaN` literals that aren't valid JSON; analysis.R preprocessed them to `null` before parsing. A handful of trajectories with degenerate summary stats were imputed at the column median for classification.
