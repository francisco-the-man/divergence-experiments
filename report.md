# Stub experiment: Monte Carlo π estimate

**Ticket:** `tk_753c562a`
**Submitted:** 2026-05-17T23:38:04+00:00
**Original idea:** task 6 test — github + email pipeline

This is the task-5 stub. It runs a hardcoded Monte Carlo π estimation
across 10 parallel Modal containers to prove the compute
+ cost-capture path before real experiments and agents land.

## Setup

- Function: `pi_estimator_task` (0.25 vCPU / 128 MB per container)
- Containers in parallel: 10
- Samples per container: 200,000
- Total samples: 2,000,000

## Per-container results

| seed | samples | inside | estimate |
|------|--------:|-------:|---------:|
| 0 | 200,000 | 157,144 | 3.142880 |
| 1 | 200,000 | 156,860 | 3.137200 |
| 2 | 200,000 | 157,245 | 3.144900 |
| 3 | 200,000 | 156,891 | 3.137820 |
| 4 | 200,000 | 156,955 | 3.139100 |
| 5 | 200,000 | 157,025 | 3.140500 |
| 6 | 200,000 | 156,801 | 3.136020 |
| 7 | 200,000 | 157,102 | 3.142040 |
| 8 | 200,000 | 157,471 | 3.149420 |
| 9 | 200,000 | 157,064 | 3.141280 |

## Aggregate

- π estimate: **3.141116**  (true value: 3.141593)
- Absolute error: 0.000477

## Cost & timing

- Wall time: 2.87 s
- Estimated Modal cost: **$0.000290** USD
