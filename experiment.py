
import numpy as np
import powerlaw
import warnings
warnings.filterwarnings("ignore")


def _safe_div(a, b):
    return np.where(np.abs(b) < 1e-6, 1.0, a / np.where(np.abs(b) < 1e-6, 1.0, b))

BINOPS = {"+": (np.add, 2), "-": (np.subtract, 2),
          "*": (np.multiply, 2), "/": (_safe_div, 2)}
UNOPS = {"neg": (np.negative, 1), "id": (lambda x: x, 1), "sin": (np.sin, 1)}

def primitive_set(kind, n_vars):
    # Minimal: small closed set under +,-,* with one constant.
    # Rich: adds identity-equivalent unaries (neg, id), division (often a no-op
    # when divisor=1), and extra constants {0, 2} — these create many distinct
    # genotypes mapping to the same phenotype, inflating neutral-network sizes
    # of the induced fitness landscape. This is the manipulated variable.
    if kind == "minimal":
        funcs = [("+", 2), ("-", 2), ("*", 2)]
        terms = [f"x{i}" for i in range(n_vars)] + ["1.0"]
    else:
        funcs = [("+", 2), ("-", 2), ("*", 2), ("/", 2),
                 ("neg", 1), ("id", 1)]
        terms = [f"x{i}" for i in range(n_vars)] + ["1.0", "0.0", "2.0"]
    return funcs, terms


def gen_tree(rng, funcs, terms, max_depth, depth=0):
    if depth >= max_depth or (depth > 0 and rng.random() < 0.3):
        t = terms[rng.integers(len(terms))]
        if t.startswith("x"):
            return ("var", t)
        return ("const", float(t))
    op, arity = funcs[rng.integers(len(funcs))]
    children = [gen_tree(rng, funcs, terms, max_depth, depth+1) for _ in range(arity)]
    return (op, *children)


def tree_size(t):
    if t[0] in ("var", "const"):
        return 1
    return 1 + sum(tree_size(c) for c in t[1:])


def tree_nodes(t, acc=None):
    if acc is None:
        acc = []
    acc.append(t)
    if t[0] not in ("var", "const"):
        for c in t[1:]:
            tree_nodes(c, acc)
    return acc


def replace_node(t, target, replacement):
    if t is target:
        return replacement
    if t[0] in ("var", "const"):
        return t
    return (t[0], *[replace_node(c, target, replacement) for c in t[1:]])


def evaluate(t, X):
    op = t[0]
    if op == "var":
        return X[:, int(t[1][1:])]
    if op == "const":
        return np.full(X.shape[0], t[1])
    if op in BINOPS:
        f, _ = BINOPS[op]
        return f(evaluate(t[1], X), evaluate(t[2], X))
    if op in UNOPS:
        f, _ = UNOPS[op]
        return f(evaluate(t[1], X))
    raise ValueError(op)


def fitness(t, X, y):
    with np.errstate(all="ignore"):
        yp = evaluate(t, X)
    if not np.all(np.isfinite(yp)):
        return 1e10
    err = yp - y
    mse = float(np.mean(err * err))
    if not np.isfinite(mse):
        return 1e10
    return mse


def crossover(rng, p1, p2, size_cap):
    nodes1 = tree_nodes(p1)
    nodes2 = tree_nodes(p2)
    a = nodes1[rng.integers(len(nodes1))]
    b = nodes2[rng.integers(len(nodes2))]
    child = replace_node(p1, a, b)
    if tree_size(child) > size_cap:
        return p1
    return child


def mutate(rng, p, funcs, terms, size_cap, max_depth=3):
    nodes = tree_nodes(p)
    a = nodes[rng.integers(len(nodes))]
    sub = gen_tree(rng, funcs, terms, max_depth)
    child = replace_node(p, a, sub)
    if tree_size(child) > size_cap:
        return p
    return child


def run_gp(rng, target_fn, n_vars, kind, pop_size, max_gen, size_cap, noise_sigma):
    funcs, terms = primitive_set(kind, n_vars)
    X = rng.uniform(-2, 2, size=(60, n_vars))
    y_clean = target_fn(X)
    # Additive y-noise creates a non-zero MSE floor. Without it, easy targets
    # (T1, T2) hit exact-zero error and the late-phase plateau collapses to a
    # single block, killing the plateau-duration statistic. With noise, the GP
    # keeps making small strict-improvement jumps even after convergence.
    y = y_clean + rng.normal(0, noise_sigma * (np.std(y_clean) + 1e-9),
                              size=y_clean.shape)

    pop = [gen_tree(rng, funcs, terms, max_depth=4) for _ in range(pop_size)]
    fits = np.array([fitness(t, X, y) for t in pop])
    best_hist = np.empty(max_gen)
    tour_k = 4
    px = 0.85

    for g in range(max_gen):
        new_pop = [pop[int(np.argmin(fits))]]   # elitism: best survives
        new_fits = [fits.min()]
        while len(new_pop) < pop_size:
            idx = rng.integers(0, pop_size, size=tour_k)
            p1 = pop[idx[np.argmin(fits[idx])]]
            if rng.random() < px:
                idx2 = rng.integers(0, pop_size, size=tour_k)
                p2 = pop[idx2[np.argmin(fits[idx2])]]
                child = crossover(rng, p1, p2, size_cap)
            else:
                child = mutate(rng, p1, funcs, terms, size_cap)
            if rng.random() < 0.15:
                child = mutate(rng, child, funcs, terms, size_cap)
            new_pop.append(child)
            new_fits.append(fitness(child, X, y))
        pop = new_pop
        fits = np.array(new_fits)
        best_hist[g] = fits.min()
    return best_hist


def plateau_durations(best_hist):
    """Inter-improvement times. Elitism makes best_hist non-increasing, so a
    strict decrease ends a plateau. Returns one duration per plateau, including
    the final unterminated one."""
    durs = []
    cur = 1
    for i in range(1, len(best_hist)):
        if best_hist[i] < best_hist[i-1]:
            durs.append(cur)
            cur = 1
        else:
            cur += 1
    durs.append(cur)
    return np.array(durs)


def neutral_network_proxy(kind, n_vars, n_samples=2000, max_depth=4, seed=0):
    """Independent measurement of phenotypic redundancy of the primitive set.
    Sample random expressions, fingerprint each by its output vector on fixed
    probe points (semantic equivalence). Return the membership-weighted mean
    class size = E[|class containing a random sample|] = sum(s^2)/sum(s).
    This is the e-graph stand-in: probe-point semantic equivalence is the same
    canonicalization principle e-graphs use (Kronberger/de França 2024), just
    via numerical signature instead of explicit rewriting."""
    rng = np.random.default_rng(seed)
    funcs, terms = primitive_set(kind, n_vars)
    probe = rng.uniform(-1.5, 1.5, size=(12, n_vars))
    fps = {}
    for _ in range(n_samples):
        t = gen_tree(rng, funcs, terms, max_depth)
        with np.errstate(all="ignore"):
            try:
                v = evaluate(t, probe)
            except Exception:
                continue
        if not np.all(np.isfinite(v)):
            continue
        key = tuple(np.round(v, 4).tolist())
        fps[key] = fps.get(key, 0) + 1
    sizes = np.array(list(fps.values()))
    if len(sizes) == 0:
        return {"mean_class_size_weighted": 0.0, "n_classes": 0, "n_valid": 0,
                "max_class_size": 0}
    return {
        "mean_class_size_weighted": float((sizes * sizes).sum() / sizes.sum()),
        "n_classes": int(len(sizes)),
        "n_valid": int(sizes.sum()),
        "max_class_size": int(sizes.max()),
    }


def fit_lognormal(durs):
    d = durs[durs > 0].astype(float)
    if len(d) < 10:
        return None
    ld = np.log(d)
    return {"mu": float(ld.mean()),
            "sigma": float(ld.std(ddof=1)),
            "n": int(len(d))}


def csn_compare(durs):
    """CSN log-likelihood-ratio test, lognormal vs exponential. normalized_ratio
    returns standardized R (Vuong-style); p tests H0: R == 0."""
    d = durs[durs >= 1].astype(float)
    if len(d) < 30:
        return {"R": None, "p": None, "n_tail": int(len(d)),
                "err": "too few events for CSN comparison"}
    try:
        fit = powerlaw.Fit(d, discrete=True, verbose=False, xmin=1)
        R, p = fit.distribution_compare("lognormal", "exponential",
                                         normalized_ratio=True)
        return {"R": float(R), "p": float(p), "n_tail": int(len(d))}
    except Exception as e:
        return {"R": None, "p": None, "n_tail": int(len(d)), "err": str(e)}


def run(params: dict) -> dict:
    target_id = params["target_id"]
    expr = params["expr"]
    kind = params["primitives"]
    n_reps = int(params["n_replicates"])
    max_gen = int(params["max_generations"])
    pop_size = int(params["pop_size"])
    size_cap = int(params["size_cap"])
    seed_base = int(params["seed_base"])
    noise_sigma = float(params.get("noise_sigma", 0.05))

    import sympy as sp
    syms = sorted({s for s in ["x0", "x1", "x2"] if s in expr})
    sym_objs = sp.symbols(syms)
    sym_map = dict(zip(syms, sym_objs))
    sym_expr = sp.sympify(expr, locals=sym_map)
    f_lambd = sp.lambdify(sym_objs, sym_expr, "numpy")
    n_vars = len(syms)

    def target_fn(X):
        args = [X[:, i] for i in range(n_vars)]
        out = f_lambd(*args)
        return np.broadcast_to(out, (X.shape[0],)).astype(float)

    early_durs_all = []
    late_durs_all = []
    per_rep = []
    for r in range(n_reps):
        rng = np.random.default_rng(seed_base + r)
        hist = run_gp(rng, target_fn, n_vars, kind, pop_size, max_gen,
                      size_cap, noise_sigma)
        half = max_gen // 2
        early_durs_all.append(plateau_durations(hist[:half]))
        late_durs_all.append(plateau_durations(hist[half:]))
        per_rep.append({
            "final_best": float(hist[-1]),
            "init_best": float(hist[0]),
            "n_plateaus_early": int(len(early_durs_all[-1])),
            "n_plateaus_late": int(len(late_durs_all[-1])),
        })

    pooled_early = np.concatenate(early_durs_all) if early_durs_all else np.array([])
    pooled_late = np.concatenate(late_durs_all) if late_durs_all else np.array([])

    ln_early = fit_lognormal(pooled_early)
    ln_late = fit_lognormal(pooled_late)
    csn_late = csn_compare(pooled_late)

    nn = neutral_network_proxy(kind, n_vars, n_samples=2000,
                                max_depth=4, seed=seed_base + 99999)

    def _sample(arr, k=2000):
        if len(arr) <= k:
            return arr.tolist()
        return np.random.default_rng(0).choice(arr, k, replace=False).tolist()

    return {
        "target_id": target_id,
        "primitives": kind,
        "expr": expr,
        "n_replicates": n_reps,
        # PRIMARY: log-normal σ on pooled late-phase plateau durations.
        # Hypothesis predicts rich > minimal (paired across the 5 targets).
        "pooled_sigma_late": ln_late["sigma"] if ln_late else None,
        "pooled_mu_late": ln_late["mu"] if ln_late else None,
        "n_late_events": ln_late["n"] if ln_late else 0,
        # Secondary descriptive: early-phase σ/μ.
        "pooled_sigma_early": ln_early["sigma"] if ln_early else None,
        "pooled_mu_early": ln_early["mu"] if ln_early else None,
        "n_early_events": ln_early["n"] if ln_early else 0,
        # CSN lognormal-vs-exponential on pooled late durations.
        "csn_late": csn_late,
        # Independent NN-size measurement (e-graph proxy).
        "neutral_proxy": nn,
        # Per-rep aggregates for replicate-level bootstrap.
        "per_rep_final_best_mean": float(np.mean([p["final_best"] for p in per_rep])),
        "per_rep_final_best_median": float(np.median([p["final_best"] for p in per_rep])),
        "per_rep_n_plateaus_late_mean": float(np.mean([p["n_plateaus_late"] for p in per_rep])),
        "per_rep_n_plateaus_early_mean": float(np.mean([p["n_plateaus_early"] for p in per_rep])),
        "per_rep_records": per_rep,
        # Raw pooled durations (subsampled to ≤2000 each) for replots / refits.
        "pooled_late_durations_sample": _sample(pooled_late),
        "pooled_early_durations_sample": _sample(pooled_early),
        "max_generations": max_gen,
        "pop_size": pop_size,
        "noise_sigma": noise_sigma,
    }
