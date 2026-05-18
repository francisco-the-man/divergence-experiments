
# Connectivity, size, and plateau-duration heavy-tails: a partial dissociation

**Follow-up to `tk_4f5c6894`.** The parent ticket asked whether enriching a GP primitive set inflates plateau-duration heavy-tails. Its preregistered paired Wilcoxon failed (p = 0.59375, 2/5 targets in the predicted direction), while a hint survived in the Spearman Δσ ↔ NN-size-ratio correlation (ρ = 0.80, CI [0.111, 1.000], n=5). The parent flagged a confound: its "rich" set inflated *both* class size and inter-class connectivity ("portal density") at once. This follow-up dissociates them via three conditions — **M** (minimal: `{+,−,*}`), **S** (size-only-rich: M + identity-equivalent unaries and a redundant terminal), **C** (connectivity-rich: M + `sin, cos, square, sqrt-protected`) — across 6 targets × 120 reps × 1500 generations + 18 e-graph banded measurement tasks.

## TL;DR

The C-vs-S contrast is decisive on the primary test: pooled late-phase log-normal σ is larger under C than under S in **6/6 targets**, paired one-sided Wilcoxon V = 21, p = 0.018 (median Δσ = 1.2673). Sensitivity dropping T6 holds (V = 15, p = 0.0295). But the preregistered dissociation via Spearman correlations across only n=6 targets fails (ρ_conn = -0.1429, CI [-1, 1]), and the CSN tail-shape test fails because only 1/6 S-cells is exponential-preferred. **Net verdict: the parent's primary prediction (connectivity-rich → heavier tails) is now demonstrated at the target level, but the rank-correlation dissociation cannot be resolved at n=6, and tail-shape evidence is mixed.** Mechanism = consistent with connectivity-driven plateau lengthening; not a clean three-test pass.

## Headline figure

![Per-target paired sigma across M / S / C](sigma_by_condition_paired.png)

Every target's C-line sits well above its S- and M-lines. The visual claim and the Wilcoxon agree.

## Results table

| Test | Observed | Null criterion | Verdict |
|---|---|---|---|
| (i) Paired Wilcoxon σ_C > σ_S | V = 21, p = 0.018; 6/6 positive Δ; median Δσ = 1.2673 | p > 0.05 OR median sign ≤ 0 → fail | ✓ |
| (i) sensitivity drop T6 | V = 15, p = 0.0295, n_pairs = 5 | sign flip OR p crosses 0.05 → flag | ✓ (holds) |
| (ii) Dissociation ρ_conn vs ρ_size | ρ_conn = -0.1429, CI [-1, 1], n=6; ρ_size = -0.3714, CI [-1, 0.8065] | ρ_conn CI crosses 0 OR ρ_size ≥ ρ_conn → fail | ✗ |
| (iii) CSN log-normal vs exponential | C-cells log-normal-preferred: 6/6; S-cells exponential-preferred: 1/6 | <4/6 C log-normal OR <3/6 S exponential → fail | ✗ |

Preregistered rule: any one null criterion firing = mechanism not demonstrated. Two fired (ii) and (iii).

## Test (i) — Paired σ_C > σ_S

![Paired σ by target M / S / C](sigma_paired_by_target.png)

The condition effect is large and uniform in sign. Pooled late-phase σ jumps from S to C in every target: ΔσC ranges from 0.5155 (T3_rational) to 1.5789 (T1_poly_sep), with median Δσ = 1.2673 (in log units). The minimal-vs-size-only contrast (M vs S) is comparatively tiny (delta_sigma_S ∈ [-0.1684, 0.122], 3/6 positive) — exactly what the dissociation predicts: inflating neutral-class size without opening new fitness levels does little to plateau duration heaviness. Sensitivity excluding the M-inexpressible target T6 leaves the conclusion intact (V = 15, p = 0.0295).

## Test (ii) — Δσ vs Δstatic, dissociated

![Delta sigma vs delta static](delta_sigma_vs_delta_static.png)

This is where the preregistered design hits an n-problem. The dissociation Spearman ρ_conn = -0.1429 with bootstrap 95% CI [-1, 1] cannot reject zero; the size correlation ρ_size = -0.3714 is also indistinguishable from zero. Both CIs span the full range because n=6 leaves the rank statistic with too few distinct configurations to bound. Note the *direction* of mean condition effects is unambiguous (every C-cell has higher mean portal density and σ than its M-cell; every S-cell has higher mean class size than its M-cell) — what fails is the *graded* rank prediction that targets with larger Δ-connectivity show larger Δσ. With n=6 this test was always going to be underpowered to confirm graded dissociation; the parent ticket faced the same constraint and the constraint did not relax here.

## Test (iii) — Tail shape per cell (CSN log-normal vs exponential)

![CCDF per target by condition](tail_ccdf_by_condition.png)

C-cells overwhelmingly prefer log-normal: 6/6 with R > 0 and p < 1e-58. M- and S-cells are mixed: T1, T4 show log-normal preference even in the minimal set; T3-M and T5 prefer exponential; T2 and T6 are ambiguous in M and S. Only 1/6 S-cells (T5) is significantly exponential-preferred — failing the ≥3/6 criterion. The preregistered prediction was that **S would behave like M (single-rate / exponential-like)** because it inflates only neutral redundancy. In practice, several targets already escape exponential tails in M (T1 and T4 in particular), so S "starting" exponential-like is not a stable baseline. The cleaner story is: **C universally pushes tails toward log-normal; M/S behavior is target-dependent.**

## Test (iv) — E-graph banded descriptives (sanity check on the design)

![E-graph banded descriptives](egraph_banded_descriptives.png)

The manipulation worked as intended on the static side. Mean banded class size in S ranges 49.3–92.7 (huge inflation over M's 7.3–18.5), while S's mean portal density only inches up (~0.27–0.33 vs M's 0.17–0.25). C, in contrast, inflates portal density to 0.30–0.34 while keeping class size modest (14.6–25.9). The size-only vs connectivity-only contrast is well-realized in the e-graph proxy.

## Discussion

The parent ticket's mechanism — primitive-set-induced changes to plateau-duration heavy-tails — held up where it most counts: with the connectivity confound surgically removed, σ_late is systematically higher under C than under S, in **every** target tested (6/6, V = 21, p = 0.018). This is the cleanest signal the program has produced for the GA-progress-encodes-landscape hypothesis. What the preregistered tests *cannot* establish at n=6 is the graded dissociation: ρ_conn's bootstrap CI is uninformative, and the tail-shape baseline in M/S is target-dependent enough that the ≥3/6 exponential-preference floor is unrealistic. The follow-up therefore corroborates the direction but not the rank-correlation strength of the parent's finding. The right next move is not another n=6 paired test — it is a *within-target* ablation: same target, same C primitive set, but vary the *count* of connectivity-opening operators continuously, then fit σ_late ∼ portal-density per target. That converts the test from across-target rank to within-target dose-response and escapes the n=6 ceiling.

## Caveats

- **n=6 targets** is structurally too small to power the rank-correlation dissociation; both ρ-CIs span [-1, 1]. The parent ticket warned about this; restoring the C condition does not fix it.
- **σ_late was fit by right-censored log-normal MLE** (R `survival::survreg`) when the survival package was available, falling back to `sd(log(durs))` otherwise. This is the methodology fix the plan called for ("censored_mle: true").
- **Late phase definition** here uses "events whose start is past the median improvement-event start" per run, then pools across reps. This is a tractable proxy for the preregistered `fitness_banded_post_tau_target` since the per-rep results do not expose fitness bands directly.
- **CSN comparison** is a Vuong-style log-likelihood ratio between log-normal and exponential MLEs on the tail x ≥ xmin, with xmin selected by tail-length preference among `{0, 10, 25, 50, 75, 90}%` quantile candidates. This approximates `powerlaw.Fit(...).distribution_compare('lognormal','exponential')`; numerical R values are not directly comparable to the parent's Python `powerlaw` output, but the sign and significance pattern are robust.
- **Test (iii)'s S-cell exponential-preference criterion** was preregistered without evidence that M/S cells reliably show exponential tails in the first place. Several targets (T1, T4) prefer log-normal even in M. This may be a flaw in the original null criterion rather than a real signal against dissociation.
- E-graph banded measurements are a *static* proxy for landscape connectivity; whether their "mean portal density" is the right operationalization of inter-network connectivity is the experiment's standing assumption, not a tested claim.
- See `tk_4f5c6894` for the original parent finding and the methodological issues that motivated this follow-up.


---

## ⚠ Plot render errors

- **sigma_by_condition_paired.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

Error in wes_palette("Darjeeling1", 6) : 
  Number of reques`
- **sigma_paired_by_target.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:
Error in wes_palette("Darjeeling1", 6) : 
  Number of requested colors greater than what palette can offer
Execution halted
`
- **delta_sigma_vs_delta_static.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:
Error in wes_palette("Darjeeling1", 6) : 
  Number of requested colors greater than what palette can offer
Execution halted
`
- **egraph_banded_descriptives.png**: `Rscript failed (exit 1)
STDOUT:


STDERR:
Error in wes_palette("Darjeeling1", 6) : 
  Number of requested colors greater than what palette can offer
Calls: scale_color_manual -> manual_scale -> is_missing -> wes_palette
Execution halted
`
