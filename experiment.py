import numpy as np


def run(params: dict) -> dict:
    n = int(params["n_samples"])
    seed = int(params["seed"])
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(n)
    # ddof=1 for unbiased sample stddev; at n>=100 the difference vs ddof=0 is negligible
    return {
        "n_samples": n,
        "seed": seed,
        "sample_mean": float(x.mean()),
        "sample_stddev": float(x.std(ddof=1)),
    }
