
"""GP symbolic regression: plateau-duration heavy-tail analysis under
minimal vs rich primitive sets, with independent e-graph-style neutral-
network size estimate per target.

Returns observables needed by plan.null_result_criteria:
  - sigma_lognormal_late          (primary: late = start_f <= median(final_fit))
  - sigma_lognormal_late_genmid   (secondary: late = gen >= median gen)
  - sigma_lognormal_early
  - per_rep_sigma_late_list       (replicate-bootstrap input)
  - csn_loglik_ratio_exp_vs_lognormal, csn_p_exp_vs_lognormal
  - csn_loglik_ratio_exp_vs_powerlaw, csn_p_exp_vs_powerlaw
  - csn_xmin, csn_alpha, n_tail
  - neutral_network_size_estimate (target-level; same for both primitive
    conditions of a given target, but reported per task so the pairing
    in the meta-analysis is unambiguous)
"""

import math
import random
from typing import List, Tuple

import numpy as np
import sympy as sp


# ---------- target / dataset ----------

def detect_n_vars(expr_str: str) -> int:
    expr = sp.sympify(expr_str)
    idxs = [int(str(s)[1:]) for s in expr.free_symbols if str(s).startswith("x")]
    return max(idxs) + 1 if idxs else 1


def make_dataset(expr_str: str, n_vars: int, n_points: int, rng):
    X = rng.uniform(-2.0, 2.0, size=(n_points, n_vars))
    expr = sp.sympify(expr_str)
    syms = [sp.Symbol(f"x{i}") for i in range(n_vars)]
    f = sp.lambdify(syms, expr, modules="numpy")
    y = f(*[X[:, i] for i in range(n_vars)])
    y = np.asarray(y, dtype=float)
    return X, y


# ---------- GP tree ----------

# A node is either:
#   ("var", i)            terminal variable x_i
#   ("const", v)          terminal constant
#   (op_name, *children)  internal

MINIMAL_OPS = {
    "add": (2, lambda a, b: a + b),
    "sub": (2, lambda a, b: a - b),
    "mul": (2, lambda a, b: a * b),
}

# Rich set: identity-like and mul-by-constant operators that inflate
# semantic-equivalence-class sizes without changing expressive power
# enough to alter target attainability.
RICH_OPS = {
    "add": (2, lambda a, b: a + b),
    "sub": (2, lambda a, b: a - b),
    "mul": (2, lambda a, b: a * b),
    "neg": (1, lambda a: -a),
    "id":  (1, lambda a: a),                       # pure identity
    "add0": (1, lambda a: a + 0.0),                # identity-equivalent
    "mul1": (1, lambda a: a * 1.0),                # identity-equivalent
    "scale2": (1, lambda a: a * 2.0),
    "scale_half": (1, lambda a: a * 0.5),
}


def tree_size(t):
    if t[0] in ("var", "const"):
        return 1
    return 1 + sum(tree_size(c) for c in t[1:])


def tree_depth(t):
    if t[0] in ("var", "const"):
        return 1
    return 1 + max(tree_depth(c) for c in t[1:])


def random_tree(rng, ops, n_vars, max_depth, depth=0, terminal_p=0.3):
    if depth >= max_depth or (depth > 0 and rng.random() < terminal_p):
        if rng.random() < 0.7:
            return ("var", int(rng.integers(0, n_vars)))
        return ("const", float(rng.uniform(-2.0, 2.0)))
    op = list(ops.keys())[int(rng.integers(0, len(ops)))]
    arity = ops[op][0]
    return (op,) + tuple(random_tree(rng, ops, n_vars, max_depth, depth + 1, terminal_p)
                         for _ in range(arity))


def eval_tree(t, X, ops):
    tag = t[0]
    if tag == "var":
        return X[:, t[1]]
    if tag == "const":
        return np.full(X.shape[0], t[1])
    fn = ops[tag][1]
    args = [eval_tree(c, X, ops) for c in t[1:]]
    return fn(*args)


def fitness(t, X, y, ops):
    with np.errstate(all="ignore"):
        yhat = eval_tree(t, X, ops)
    yhat = np.asarray(yhat, dtype=float)
    if not np.all(np.isfinite(yhat)):
        return 1e12
    err = yhat - y
    mse = float(np.mean(err * err))
    if not math.isfinite(mse):
        return 1e12
    return mse


def all_nodes(t, path=()):
    yield path, t
    if t[0] not in ("var", "const"):
        for i, c in enumerate(t[1:]):
            yield from all_nodes(c, path + (i,))


def get_subtree(t, path):
    for i in path:
        t = t[1 + i]
    return t


def set_subtree(t, path, new):
    if not path:
        return new
    i = path[0]
    children = list(t[1:])
    children[i] = set_subtree(children[i], path[1:], new)
    return (t[0],) + tuple(children)


def crossover(a, b, rng):
    nodes_a = list(all_nodes(a))
    nodes_b = list(all_nodes(b))
    pa = nodes_a[int(rng.integers(0, len(nodes_a)))][0]
    pb = nodes_b[int(rng.integers(0, len(nodes_b)))][0]
    sub_b = get_subtree(b, pb)
    return set_subtree(a, pa, sub_b)


def mutate(t, rng, ops, n_vars, max_depth):
    nodes = list(all_nodes(t))
    p = nodes[int(rng.integers(0, len(nodes)))][0]
    new_sub = random_tree(rng, ops, n_vars, max_depth=3)
    return set_subtree(t, p, new_sub)


def tournament(pop_fits, rng, k=3):
    n = len(pop_fits)
    idxs = rng.integers(0, n, size=k)
    best = idxs[0]
    for i in idxs[1:]:
        if pop_fits[i] < pop_fits[best]:
            best = i
    return int(best)


# ---------- one GP run ----------

def gp_run(seed, expr_str, primitives, max_generations, pop_size, size_cap,
           wall_budget_s, n_vars, X, y):
    import time
    rng = np.random.default_rng(seed)
    random.seed(seed)

    ops = MINIMAL_OPS if primitives == "minimal" else RICH_OPS
    pop = [random_tree(rng, ops, n_vars, max_depth=4) for _ in range(pop_size)]
    fits = [fitness(t, X, y, ops) for t in pop]

    best_curve = []
    t0 = time.time()
    for gen in range(max_generations):
        b = float(min(fits))
        best_curve.append(b)
        if time.time() - t0 > wall_budget_s:
            break
        new_pop = []
        # elitism
        elite_idx = int(np.argmin(fits))
        new_pop.append(pop[elite_idx])
        while len(new_pop) < pop_size:
            i = tournament(fits, rng)
            j = tournament(fits, rng)
            child = crossover(pop[i], pop[j], rng)
            if rng.random() < 0.3:
                child = mutate(child, rng, ops, n_vars, max_depth=4)
            if tree_size(child) > size_cap or tree_depth(child) > 8:
                child = pop[i]
            new_pop.append(child)
        pop = new_pop
        fits = [fitness(t, X, y, ops) for t in pop]
    return np.array(best_curve, dtype=float)


# ---------- plateau extraction ----------

def extract_plateaus(curve, tol=1e-9):
    """Plateau = consecutive generations where best-fitness is unchanged
    (within tol). Returns list of (duration, start_gen, start_fit)."""
    plateaus = []
    if len(curve) < 2:
        return plateaus
    start = 0
    for i in range(1, len(curve)):
        if curve[i] < curve[start] - tol:
            plateaus.append((i - start, start, float(curve[start])))
            start = i
    # exclude the trailing plateau (right-censored)
    return plateaus


# ---------- log-normal sigma on a set of durations ----------

def fit_lognormal_sigma(durations):
    arr = np.asarray([d for d in durations if d >= 1], dtype=float)
    if len(arr) < 5:
        return float("nan"), len(arr)
    logs = np.log(arr)
    return float(np.std(logs, ddof=1)), len(arr)


# ---------- e-graph-style neutral-network estimate ----------
# We estimate, for the target's fitness landscape, the average semantic-
# equivalence-class size at fitness levels achievable by trees of size
# <= size_cap, separately for the two primitive sets. The ratio rich:min
# is the per-target observable used in the meta-analysis correlation.

def semantic_signature(t, X, ops, n_probe=32):
    """Behavioral fingerprint of a tree: rounded output on a probe set."""
    with np.errstate(all="ignore"):
        y = eval_tree(t, X[:n_probe], ops)
    y = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(y)):
        return None
    return tuple(np.round(y, 4).tolist())


def neutral_network_estimate(expr_str, n_vars, primitives, size_cap, n_samples, seed):
    """Sample random trees up to size_cap, group by semantic signature on
    a fixed probe set; report mean class size (proxy for NN size)."""
    rng = np.random.default_rng(seed)
    ops = MINIMAL_OPS if primitives == "minimal" else RICH_OPS
    probe_rng = np.random.default_rng(12345)  # fixed probe across conditions
    Xp = probe_rng.uniform(-2.0, 2.0, size=(64, n_vars))

    from collections import Counter
    sigs = Counter()
    drawn = 0
    while drawn < n_samples:
        depth = int(rng.integers(2, 6))
        t = random_tree(rng, ops, n_vars, max_depth=depth)
        if tree_size(t) > size_cap:
            continue
        s = semantic_signature(t, Xp, ops, n_probe=32)
        if s is None:
            continue
        sigs[s] += 1
        drawn += 1
    if not sigs:
        return float("nan")
    counts = np.array(list(sigs.values()), dtype=float)
    # size-biased mean: probability a random sample lands in a class of
    # size c is proportional to c, so E[class size | sample] = sum c^2 / sum c
    return float((counts * counts).sum() / counts.sum())


# ---------- CSN comparison ----------

def csn_analysis(durations):
    import powerlaw
    arr = np.asarray([d for d in durations if d >= 1], dtype=float)
    out = {"csn_xmin": float("nan"), "csn_alpha": float("nan"),
           "n_tail": 0,
           "csn_loglik_ratio_exp_vs_lognormal": float("nan"),
           "csn_p_exp_vs_lognormal": float("nan"),
           "csn_loglik_ratio_exp_vs_powerlaw": float("nan"),
           "csn_p_exp_vs_powerlaw": float("nan")}
    if len(arr) < 30:
        return out
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fit = powerlaw.Fit(arr, discrete=True, verbose=False)
        except Exception:
            return out
        out["csn_xmin"] = float(fit.xmin) if fit.xmin else float("nan")
        out["csn_alpha"] = float(fit.alpha) if fit.alpha else float("nan")
        out["n_tail"] = int((arr >= (fit.xmin or 0)).sum())
        try:
            R, p = fit.distribution_compare("exponential", "lognormal",
                                            normalized_ratio=True)
            out["csn_loglik_ratio_exp_vs_lognormal"] = float(R)
            out["csn_p_exp_vs_lognormal"] = float(p)
        except Exception:
            pass
        try:
            R, p = fit.distribution_compare("exponential", "power_law",
                                            normalized_ratio=True)
            out["csn_loglik_ratio_exp_vs_powerlaw"] = float(R)
            out["csn_p_exp_vs_powerlaw"] = float(p)
        except Exception:
            pass
    return out


# ---------- main task entry point ----------

def run(params: dict) -> dict:
    target_id = params["target_id"]
    expr_str = params["expr"]
    primitives = params["primitives"]
    n_replicates = int(params["n_replicates"])
    max_generations = int(params["max_generations"])
    pop_size = int(params["pop_size"])
    size_cap = int(params["size_cap"])
    seed_base = int(params["seed_base"])

    # Budget management: 1 vCPU, 30 min. Reserve a few min for analysis +
    # NN estimate. Distribute the rest across replicates. With 120 reps
    # this gives ~12-13 s per rep; reps will be generation-capped or
    # wall-capped depending on which bites first.
    total_budget_s = 26 * 60
    nn_budget_s = 90
    per_rep_budget_s = max(8.0, (total_budget_s - nn_budget_s) / n_replicates)

    n_vars = detect_n_vars(expr_str)

    # Dataset: fixed across reps within a task (same target, same fit
    # landscape); seeded from seed_base so reps share an X,y but each
    # rep's GP search uses its own rng stream.
    data_rng = np.random.default_rng(seed_base)
    X, y = make_dataset(expr_str, n_vars, n_points=64, rng=data_rng)

    import time
    t_start = time.time()

    all_durations: List[Tuple[int, int, float, int]] = []  # (dur, start_gen, start_fit, rep)
    final_fits = []
    curve_lens = []
    completed = 0

    for r in range(n_replicates):
        if time.time() - t_start > total_budget_s - nn_budget_s:
            break
        seed = seed_base + r + 1
        curve = gp_run(seed, expr_str, primitives, max_generations, pop_size,
                       size_cap, per_rep_budget_s, n_vars, X, y)
        if len(curve) < 50:
            continue
        curve_lens.append(int(len(curve)))
        final_fits.append(float(curve[-1]))
        plats = extract_plateaus(curve)
        for d, sg, sf in plats:
            all_durations.append((d, sg, sf, r))
        completed += 1

    # ---- aggregate ----
    if not all_durations:
        return {
            "target_id": target_id, "primitives": primitives,
            "n_replicates_completed": completed,
            "sigma_lognormal_late": float("nan"),
            "sigma_lognormal_late_genmid": float("nan"),
            "sigma_lognormal_early": float("nan"),
            "per_rep_sigma_late_list": [],
            "neutral_network_size_estimate": float("nan"),
            "n_plateaus_total": 0,
            "error": "no_plateaus",
        }

    # Phase split — PRIMARY: by start-fitness vs median final fit.
    median_final = float(np.median(final_fits)) if final_fits else float("nan")
    late_durs_primary = [d for (d, _sg, sf, _r) in all_durations if sf <= median_final]
    early_durs = [d for (d, _sg, sf, _r) in all_durations if sf > median_final]

    sigma_late, n_late = fit_lognormal_sigma(late_durs_primary)
    sigma_early, n_early = fit_lognormal_sigma(early_durs)

    # SECONDARY: by gen-midpoint, per-rep.
    median_gen_per_rep = {}
    for d, sg, _sf, r in all_durations:
        median_gen_per_rep.setdefault(r, []).append(sg)
    median_gen_per_rep = {r: float(np.median(g)) for r, g in median_gen_per_rep.items()}
    late_durs_genmid = [d for (d, sg, _sf, r) in all_durations
                        if sg >= median_gen_per_rep.get(r, 0)]
    sigma_late_genmid, n_late_genmid = fit_lognormal_sigma(late_durs_genmid)

    # Per-rep sigma (late, primary) — for replicate-level bootstrap.
    per_rep_sigma_late = []
    reps_present = sorted({rr for _, _, _, rr in all_durations})
    for r in reps_present:
        rep_late = [d for (d, _sg, sf, rr) in all_durations
                    if rr == r and sf <= median_final]
        s, n = fit_lognormal_sigma(rep_late)
        if math.isfinite(s) and n >= 5:
            per_rep_sigma_late.append(s)

    # CSN on late (primary) durations.
    csn = csn_analysis(late_durs_primary)

    # ---- neutral-network estimate ----
    nn_size = neutral_network_estimate(
        expr_str, n_vars, primitives, size_cap=size_cap,
        n_samples=3000, seed=seed_base + 99991)

    out = {
        "target_id": target_id,
        "primitives": primitives,
        "expr": expr_str,
        "n_replicates_completed": completed,
        "n_replicates_requested": n_replicates,
        "n_plateaus_total": len(all_durations),
        "n_plateaus_late_primary": n_late,
        "n_plateaus_late_genmid": n_late_genmid,
        "n_plateaus_early": n_early,
        "median_final_fitness": median_final,
        "mean_curve_len": float(np.mean(curve_lens)) if curve_lens else 0.0,
        "sigma_lognormal_late": sigma_late,                     # primary
        "sigma_lognormal_late_genmid": sigma_late_genmid,       # secondary
        "sigma_lognormal_early": sigma_early,
        "per_rep_sigma_late_list": per_rep_sigma_late,
        "neutral_network_size_estimate": nn_size,
        "wall_time_s": float(time.time() - t_start),
    }
    out.update(csn)
    return out
