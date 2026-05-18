# CLT Smoke Test — Pipeline Verification

## TL;DR

Pipeline works end-to-end. Sample mean and stddev of N(0,1) both converge as predicted by CLT: at n=10,000 the median sample mean is ~0.006 (criterion: |·| < 0.05) and median sample stddev is ~0.999 (criterion: |·−1| < 0.05). Both preregistered tests **pass**. Error visually scales as ~1/√n across the three orders of magnitude tested.

## Headline figure

![Convergence of sample mean and stddev](plot_convergence.png)

## Results table

| Test | Observed (n=10,000) | Null criterion | Pass/Fail |
|---|---|---|---|
| `mean_converges` | median(sample_mean) = 0.00631 → \|·\| = 0.00631 | \|median\| < 0.05 | ✓ PASS |
| `stddev_converges` | median(sample_stddev) = 0.99855 → \|·−1\| = 0.00145 | \|median−1\| < 0.05 | ✓ PASS |

## Per-test figures

### Mean convergence

![Sample mean vs n](plot_mean.png)

At n=100 the three seeds scatter across ~[−0.07, 0.08]; by n=10,000 they cluster within ±0.013 of zero. The spread shrinks consistent with the CLT-predicted 1/√n rate (each 10× in n should shrink SE by ~3.16×; observed shrinkage from n=100→10,000 is roughly 10×, matching √100). Consistent with the hypothesis.

### Stddev convergence

![Sample stddev vs n](plot_stddev.png)

Sample stddev approaches 1 from below at small n (seed 1 at n=100 is 0.856) and tightens to within 0.006 of unity by n=10,000. Same 1/√n tightening story as the mean. Consistent with the hypothesis.

## Discussion

This was a smoke test, not a scientific question — the goal was to confirm the Divergence pipeline (PI → coder → worker → reporter → email/GitHub) executes cleanly on a trivial computation with a known answer. It does. The CLT result itself is unsurprising: `numpy.random.default_rng` produces N(0,1) samples whose moments converge at the textbook rate. No follow-up is warranted on the statistical content; any pipeline issues surfaced (or not) by this run are the actual signal.

## Caveats

- n=3 seeds per condition — "median" is the middle of three points, not a robust estimate. Fine for a smoke test, would be inadequate for a real claim.
- Only three values of n_samples; the 1/√n claim is asserted visually, not fit.
- `ddof=1` used for sample stddev; at n≥100 this is numerically indistinguishable from ddof=0.
- No exploratory analyses added beyond the preregistered plan.
