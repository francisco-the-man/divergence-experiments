
"""
One task = one (target x primitive_set) cell of the GP plateau-statistics experiment.

Returns per-cell observables that the aggregator will pair across conditions:
  - log-normal sigma fits (early/late plateaus, replicate-level + pooled)
  - CSN power-law vs log-normal/exponential LR tests (pooled, late phase)
  - neutral-network proxy: count of distinct semantic equivalence classes among
    random size-bounded trees whose training MSE falls below each of a fixed
    set of TARGET-DETERMINED fitness checkpoints (so rich and minimal cells
    for the same target are probed at matched MSE thresholds).
"""

import math
import numpy as np
from scipy import stats
import powerlaw
import sympy as sp


# Target-only seeds for training data: identical (X,y) across rich/minimal cells.
TARGET_DATA_SEED = {
    "T1_poly_sep": 101,
    "T2_coupled":  102,
    "T3_rational": 103,
    "T4_deep_mul": 104,
    "T5_transcend": 105,
}

# Fraction-of-var(y) checkpoints; identical across conditions per target.
CHECKPOINT_FRACTIONS = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001]

N_VARS = 3
N_POINTS = 60


# ---------- training data ----------

def make_data(target_id, expr_str):
    rng = np.random.default_rng(TARGET_DATA_SEED[target_id])
    X = rng.uniform(-2.0, 2.0, size=(N_POINTS, N_VARS))
    xs = sp.symbols("x0 x1 x2")
    f = sp.lambdify(xs, sp.sympify(expr_str), modules="numpy")
    y = f(X[:, 0], X[:, 1], X[:, 2])
    y = np.asarray(y, dtype=float)
    return X, y


# ---------- tree representation ----------
# Nodes: ("op", name, [children]) or ("var", idx) or ("const", value)

BINARY_OPS = {
    "add": np.add,
    "sub": np.subtract,
    "mul": np.multiply,
}
UNARY_OPS = {
    "neg": np.negative,
    "sin": np.sin,
    "cos": np.cos,
    "square": np.square,
    # 'rich' identity-equivalent / multiplicative-constant operators:
    "id":   lambda a: a,
    "mul1": lambda a: a * 1.0,
    "add0": lambda a: a + 0.0,
}

PROTECTED_DIV = lambda a, b: np.where(np.abs(b) > 1e-6, a / b, 1.0)
BINARY_OPS["pdiv"] = PROTECTED_DIV


def primitive_sets(kind):
    # Minimal: add, sub, mul, pdiv, sin, square, neg
    # Rich:    minimal + id, mul1, add0  (identity-equivalent ops --> inflate neutral nets)
    minimal_bin = ["add", "sub", "mul", "pdiv"]
    minimal_un  = ["sin", "square", "neg"]
    if kind == "minimal":
        return minimal_bin, minimal_un
    return minimal_bin, minimal_un + ["id", "mul1", "add0"]


def random_tree(rng, depth, bin_ops, un_ops, max_depth=4):
    # Grow method with early termination at terminals.
    if depth >= max_depth or rng.random() < 0.3:
        if rng.random() < 0.7:
            return ("var", int(rng.integers(0, N_VARS)))
        return ("const", float(rng.uniform(-2.0, 2.0)))
    choices = []
    for op in bin_ops:
        choices.append(("bin", op))
    for op in un_ops:
        choices.append(("un", op))
    arity, name = choices[int(rng.integers(0, len(choices)))]
    if arity == "bin":
        return ("op2", name,
                random_tree(rng, depth + 1, bin_ops, un_ops, max_depth),
                random_tree(rng, depth + 1, bin_ops, un_ops, max_depth))
    return ("op1", name, random_tree(rng, depth + 1, bin_ops, un_ops, max_depth))


def tree_size(t):
    if t[0] in ("var", "const"):
        return 1
    if t[0] == "op1":
        return 1 + tree_size(t[2])
    return 1 + tree_size(t[2]) + tree_size(t[3])


def tree_eval(t, X):
    k = t[0]
    if k == "var":
        return X[:, t[1]]
    if k == "const":
        return np.full(X.shape[0], t[1])
    if k == "op1":
        return UNARY_OPS[t[1]](tree_eval(t[2], X))
    return BINARY_OPS[t[1]](tree_eval(t[2], X), tree_eval(t[3], X))


def all_nodes(t, path=()):
    yield path, t
    if t[0] == "op1":
        yield from all_nodes(t[2], path + (0,))
    elif t[0] == "op2":
        yield from all_nodes(t[2], path + (0,))
        yield from all_nodes(t[3], path + (1,))


def replace_at(t, path, new_sub):
    if not path:
        return new_sub
    i, rest = path[0], path[1:]
    if t[0] == "op1":
        return ("op1", t[1], replace_at(t[2], rest, new_sub))
    if i == 0:
        return ("op2", t[1], replace_at(t[2], rest, new_sub), t[3])
    return ("op2", t[1], t[2], replace_at(t[3], rest, new_sub))


def mutate(t, rng, bin_ops, un_ops, max_depth=4):
    nodes = list(all_nodes(t))
    path, _ = nodes[int(rng.integers(0, len(nodes)))]
    new_sub = random_tree(rng, 0, bin_ops, un_ops, max_depth=max_depth)
    return replace_at(t, path, new_sub)


def crossover(a, b, rng):
    nodes_a = list(all_nodes(a))
    nodes_b = list(all_nodes(b))
    pa, _ = nodes_a[int(rng.integers(0, len(nodes_a)))]
    pb, sub_b = nodes_b[int(rng.integers(0, len(nodes_b)))]
    return replace_at(a, pa, sub_b)


def fitness(t, X, y, size_cap):
    if tree_size(t) > size_cap:
        return np.inf
    try:
        with np.errstate(all="ignore"):
            yhat = tree_eval(t, X)
        yhat = np.asarray(yhat, dtype=float)
        if not np.all(np.isfinite(yhat)):
            return np.inf
        return float(np.mean((yhat - y) ** 2))
    except (ValueError, FloatingPointError, ZeroDivisionError):
        return np.inf


# ---------- one GP run ----------

def run_gp(seed, X, y, pop_size, max_gens, size_cap, kind):
    rng = np.random.default_rng(seed)
    bin_ops, un_ops = primitive_sets(kind)
    pop = [random_tree(rng, 0, bin_ops, un_ops) for _ in range(pop_size)]
    fits = np.array([fitness(t, X, y, size_cap) for t in pop])
    best_history = []
    for gen in range(max_gens):
        best_history.append(float(np.min(fits)))
        # Tournament select + crossover + mutation, elitism of 1.
        new_pop = []
        elite_idx = int(np.argmin(fits))
        new_pop.append(pop[elite_idx])
        while len(new_pop) < pop_size:
            # tournament of 3
            idx = rng.integers(0, pop_size, size=3)
            p1 = pop[idx[np.argmin(fits[idx])]]
            idx = rng.integers(0, pop_size, size=3)
            p2 = pop[idx[np.argmin(fits[idx])]]
            child = crossover(p1, p2, rng) if rng.random() < 0.7 else p1
            if rng.random() < 0.3:
                child = mutate(child, rng, bin_ops, un_ops)
            new_pop.append(child)
        pop = new_pop
        fits = np.array([fitness(t, X, y, size_cap) for t in pop])
    best_history.append(float(np.min(fits)))
    return np.array(best_history)


def plateaus_from_history(hist, atol=1e-12):
    # Strictly-improving events (allow tiny numerical slack).
    improvements = []   # list of (gen, fitness_after)
    cur = hist[0]
    improvements.append((0, cur))
    for g in range(1, len(hist)):
        if hist[g] < cur - atol:
            cur = hist[g]
            improvements.append((g, cur))
    # Durations = gaps between consecutive improvement events (in generations).
    gens = [g for g, _ in improvements]
    durations = np.diff(gens).astype(float)
    fits_at = np.array([f for _, f in improvements[:-1]])  # fitness DURING that plateau
    return durations, fits_at, improvements


# ---------- log-normal sigma fit ----------

def lognormal_sigma(durations):
    d = np.asarray(durations, dtype=float)
    d = d[d > 0]
    if len(d) < 5:
        return float("nan"), len(d)
    logs = np.log(d)
    return float(np.std(logs, ddof=1)), int(len(d))


# ---------- CSN diagnostics (pooled, late phase) ----------

def csn_diagnostics(durations):
    d = np.asarray(durations, dtype=float)
    d = d[d > 0]
    if len(d) < 50:
        return {"n_tail": int(len(d)), "alpha": None, "xmin": None,
                "R_exp": None, "p_exp": None, "R_lognorm": None, "p_lognorm": None}
    # suppress powerlaw's verbose stdout
    fit = powerlaw.Fit(d, discrete=True, verbose=False)
    alpha = float(fit.power_law.alpha)
    xmin = float(fit.power_law.xmin)
    n_tail = int(np.sum(d >= xmin))
    R_exp, p_exp = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)
    R_ln, p_ln = fit.distribution_compare("power_law", "lognormal", normalized_ratio=True)
    return {"n_tail": n_tail, "alpha": alpha, "xmin": xmin,
            "R_exp": float(R_exp), "p_exp": float(p_exp),
            "R_lognorm": float(R_ln), "p_lognorm": float(p_ln)}


# ---------- semantic-equivalence-class proxy (e-graph stand-in) ----------
# We approximate "neutral network size at fitness level F" by:
#   sample random trees of size <= size_cap from the primitive set,
#   keep those with MSE <= F on training data,
#   count distinct semantic equivalence classes via a behavioral signature
#   (rounded output vector on a fixed probe input set).
# More classes at a matched F => richer neutral-network structure.
# Checkpoints are TARGET-determined (var_y * fraction), so rich vs minimal
# are probed at IDENTICAL fitness thresholds for the same target.

def neutral_proxy(X, y, checkpoints, kind, seed, n_samples=4000, size_cap=40):
    rng = np.random.default_rng(seed)
    bin_ops, un_ops = primitive_sets(kind)
    sigs_below = {c: set() for c in checkpoints}
    counts_below = {c: 0 for c in checkpoints}
    for _ in range(n_samples):
        t = random_tree(rng, 0, bin_ops, un_ops, max_depth=4)
        if tree_size(t) > size_cap:
            continue
        f = fitness(t, X, y, size_cap)
        if not np.isfinite(f):
            continue
        # signature: rounded output vector
        with np.errstate(all="ignore"):
            yhat = tree_eval(t, X)
        yhat = np.asarray(yhat, dtype=float)
        if not np.all(np.isfinite(yhat)):
            continue
        sig = tuple(np.round(yhat, 4).tolist())
        for c in checkpoints:
            if f <= c:
                sigs_below[c].add(sig)
                counts_below[c] += 1
    return {
        "n_classes": {f"{c:.6g}": len(sigs_below[c]) for c in checkpoints},
        "n_members": {f"{c:.6g}": counts_below[c] for c in checkpoints},
    }


# ---------- main task ----------

def run(params: dict) -> dict:
    target_id = params["target_id"]
    expr = params["expr"]
    kind = params["primitives"]
    n_reps = int(params["n_replicates"])
    max_gens = int(params["max_generations"])
    pop_size = int(params["pop_size"])
    size_cap = int(params["size_cap"])
    seed_base = int(params["seed_base"])

    X, y = make_data(target_id, expr)
    var_y = float(np.var(y))
    checkpoints = [var_y * f for f in CHECKPOINT_FRACTIONS]

    # Per-replicate stats
    rep_sigma_early = []
    rep_sigma_late = []
    rep_n_early = []
    rep_n_late = []
    rep_final_fit = []
    rep_n_improvements = []

    # Pooled durations across replicates (for CSN late-phase test).
    pooled_late = []
    pooled_early = []

    for r in range(n_reps):
        seed = seed_base + r
        hist = run_gp(seed, X, y, pop_size, max_gens, size_cap, kind)
        durations, fits_at, improvements = plateaus_from_history(hist)
        rep_final_fit.append(float(hist[-1]))
        rep_n_improvements.append(int(len(improvements)))
        if len(durations) < 2:
            rep_sigma_early.append(float("nan"))
            rep_sigma_late.append(float("nan"))
            rep_n_early.append(0); rep_n_late.append(0)
            continue
        # Phase split by median final fitness of THIS run's improvement trail.
        # Late = plateaus at fitness BELOW the run's median improvement fitness
        # (closer-to-solved regime; theoretically cleanest for neutral-network claim).
        median_fit = float(np.median(fits_at))
        late_mask = fits_at <= median_fit
        early_mask = ~late_mask
        s_e, n_e = lognormal_sigma(durations[early_mask])
        s_l, n_l = lognormal_sigma(durations[late_mask])
        rep_sigma_early.append(s_e); rep_n_early.append(n_e)
        rep_sigma_late.append(s_l); rep_n_late.append(n_l)
        pooled_early.extend(durations[early_mask].tolist())
        pooled_late.extend(durations[late_mask].tolist())

    pooled_early_arr = np.array(pooled_early, dtype=float)
    pooled_late_arr = np.array(pooled_late, dtype=float)
    pooled_sigma_early, n_pooled_early = lognormal_sigma(pooled_early_arr)
    pooled_sigma_late, n_pooled_late = lognormal_sigma(pooled_late_arr)

    # CSN diagnostics on pooled late-phase tail (primary heavy-tail test).
    csn_late = csn_diagnostics(pooled_late_arr)

    # Replicate-level bootstrap CI for sigma_late.
    rep_sigma_late_arr = np.array(rep_sigma_late, dtype=float)
    finite = rep_sigma_late_arr[np.isfinite(rep_sigma_late_arr)]
    boot_seed_rng = np.random.default_rng(seed_base + 7777)
    if len(finite) >= 5:
        B = 2000
        boots = np.array([
            np.mean(boot_seed_rng.choice(finite, size=len(finite), replace=True))
            for _ in range(B)
        ])
        ci_lo, ci_hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))
        mean_sigma_late = float(np.mean(finite))
    else:
        ci_lo = ci_hi = mean_sigma_late = float("nan")

    # E-graph-style neutral-network proxy, at TARGET-determined checkpoints.
    np_proxy = neutral_proxy(
        X, y, checkpoints, kind,
        seed=seed_base + 99999,  # condition-distinct sampling seed is fine; checkpoints/data identical
        n_samples=4000, size_cap=size_cap,
    )

    return {
        "target_id": target_id,
        "primitives": kind,
        "n_replicates_run": n_reps,
        "var_y": var_y,
        "checkpoint_fractions": CHECKPOINT_FRACTIONS,
        "checkpoints": checkpoints,
        # Replicate-level sigmas (for paired Wilcoxon at aggregator).
        "rep_sigma_early": rep_sigma_early,
        "rep_sigma_late": rep_sigma_late,
        "rep_n_early_events": rep_n_early,
        "rep_n_late_events": rep_n_late,
        "rep_final_fit": rep_final_fit,
        "rep_n_improvements": rep_n_improvements,
        # Pooled (primary statistic per plan).
        "pooled_sigma_early": pooled_sigma_early,
        "pooled_sigma_late":  pooled_sigma_late,
        "n_pooled_early": n_pooled_early,
        "n_pooled_late":  n_pooled_late,
        # Bootstrap CI on the replicate-mean of sigma_late.
        "sigma_late_mean": mean_sigma_late,
        "sigma_late_ci95": [ci_lo, ci_hi],
        # CSN diagnostics (secondary).
        "csn_late": csn_late,
        # Neutral-network proxy (per-condition; aggregator divides rich/minimal).
        "neutral_proxy": np_proxy,
    }
