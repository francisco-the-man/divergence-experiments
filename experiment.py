
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# -- primitives -------------------------------------------------------------
# Changes vs parent:
#  * Dropped _safe_div entirely (PI: "drop _safe_div"); division removed.
#  * Three primitive sets M, S, C (was {minimal, rich}):
#      M  baseline:           {+, -, *} on x_i ∪ {1.0}
#      S  size-only-rich:     M + identity-equivalent unaries {neg, id}
#                             + a redundant identity terminal (extra "1.0").
#                             Inflates neutral-class size but adds no new
#                             expressible functional form -> no new fitness
#                             bands ("portals").
#      C  connectivity-rich:  M + {sin, cos, square, sqrt-protected}.
#                             Opens new fitness bands (genuine portals to
#                             phenotypes M cannot reach).

def _sqrtp(x):
    return np.sqrt(np.abs(x))

def _square(x):
    return x * x

BINOPS = {"+": np.add, "-": np.subtract, "*": np.multiply}
UNOPS_S = {"neg": np.negative, "id": (lambda x: x)}
UNOPS_C = {"sin": np.sin, "cos": np.cos, "sq": _square, "sqrtp": _sqrtp}
ALL_UN = {**UNOPS_S, **UNOPS_C}


def primitive_set(cond, n_vars):
    base_funcs = [("+", 2), ("-", 2), ("*", 2)]
    base_terms = [f"x{i}" for i in range(n_vars)] + ["1.0"]
    if cond == "M":
        return base_funcs, base_terms
    if cond == "S":
        return base_funcs + [("neg", 1), ("id", 1)], base_terms + ["1.0"]
    if cond == "C":
        return base_funcs + [("sin", 1), ("cos", 1), ("sq", 1), ("sqrtp", 1)], base_terms
    raise ValueError(cond)


# -- tree ops (kept from parent) -------------------------------------------
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
        return BINOPS[op](evaluate(t[1], X), evaluate(t[2], X))
    if op in ALL_UN:
        return ALL_UN[op](evaluate(t[1], X))
    raise ValueError(op)


def fitness(t, X, y):
    with np.errstate(all="ignore"):
        yp = evaluate(t, X)
    if not np.all(np.isfinite(yp)):
        return 1e10
    err = yp - y
    mse = float(np.mean(err * err))
    if not np.isfinite(mse) or mse > 1e10:
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


def tree_to_str(t):
    if t[0] == "var":
        return t[1]
    if t[0] == "const":
        return repr(t[1])
    return "(" + t[0] + " " + " ".join(tree_to_str(c) for c in t[1:]) + ")"


# -- GP loop ---------------------------------------------------------------
# Changes vs parent:
#   * tournament k = 2 (was 4): weaker selection allows true neutral drift.
#   * pop_size and max_gen are params-driven (defaults 200 / 1500).
#   * collect trajectory-visited best trees (sub-sampled at log-spaced gens)
#     for downstream banded e-graph measurement; serialized as strings so the
#     payload stays JSON-clean.
def run_gp(rng, target_fn, n_vars, cond, pop_size, max_gen, size_cap, noise_sigma):
    funcs, terms = primitive_set(cond, n_vars)
    X = rng.uniform(-2, 2, size=(60, n_vars))
    y_clean = target_fn(X)
    # additive y-noise to keep a non-zero MSE floor (parent's rationale)
    y = y_clean + rng.normal(0, noise_sigma * (np.std(y_clean) + 1e-9), size=y_clean.shape)

    pop = [gen_tree(rng, funcs, terms, max_depth=4) for _ in range(pop_size)]
    fits = np.array([fitness(t, X, y) for t in pop])
    best_hist = np.empty(max_gen)
    tour_k = 2
    px = 0.85
    visited_samples = []

    for g in range(max_gen):
        new_pop = [pop[int(np.argmin(fits))]]
        new_fits = [float(fits.min())]
        while len(new_pop) < pop_size:
            idx = rng.integers(0, pop_size, size=tour_k)
            p1 = pop[idx[int(np.argmin(fits[idx]))]]
            if rng.random() < px:
                idx2 = rng.integers(0, pop_size, size=tour_k)
                p2 = pop[idx2[int(np.argmin(fits[idx2]))]]
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

        if g < 64 or (g & (g - 1)) == 0:
            best_idx = int(np.argmin(fits))
            visited_samples.append((g, float(fits[best_idx]), tree_to_str(pop[best_idx])))

    return best_hist, visited_samples


# -- plateau extraction with censoring -------------------------------------
# A plateau is a maximal run where best-of-population is constant (elitism
# preserves it). A strict decrease ends a plateau; the trailing plateau is
# right-censored (we never observed it end). The fitness *level* of each
# plateau is its constant best-fit value; banding to a coarse log-fitness
# bin is done in post-processing where per-target ranges are pooled across
# conditions (so bands are comparable across M, S, C).
def plateau_events(best_hist):
    events = []
    n = len(best_hist)
    start = 0
    for i in range(1, n):
        if best_hist[i] < best_hist[i-1] - 1e-15:
            events.append((int(start), int(i - 1), int(i - start),
                           float(best_hist[start]), 1))  # observed
            start = i
    events.append((int(start), int(n - 1), int(n - start),
                   float(best_hist[start]), 0))          # censored
    return events


# -- banded e-graph proxy --------------------------------------------------
# New task kind 'egraph_banded'. Sample random expressions from the primitive
# set, evaluate them on a noiseless probe set, bin by fitness band, and per
# band compute:
#   * mean class size = membership-weighted E[|class|] via semantic
#     fingerprint (rounded output vector) -> "within-network size" proxy
#   * out-degree portal density = P(one-node random mutation lands in a
#     strictly better band) -> "between-network connectivity" proxy
# The two measurements are independent of the GP dynamics (static landscape
# probe). Hypothesis predicts S inflates class size; C inflates portal
# density; M is the floor on both.
def egraph_banded_measure(rng, target_fn, n_vars, cond, n_samples, max_depth, n_bands):
    from collections import Counter
    funcs, terms = primitive_set(cond, n_vars)
    X = rng.uniform(-2, 2, size=(40, n_vars))
    y = target_fn(X)
    trees = [gen_tree(rng, funcs, terms, max_depth) for _ in range(n_samples)]
    fits, sigs = [], []
    for t in trees:
        with np.errstate(all="ignore"):
            yp = evaluate(t, X)
        if not np.all(np.isfinite(yp)):
            fits.append(np.inf)
            sigs.append(None)
            continue
        fits.append(float(np.mean((yp - y) ** 2)))
        sigs.append(tuple(np.round(yp, 4).tolist()))
    fits = np.array(fits)
    finite = np.isfinite(fits)
    if finite.sum() < 10:
        return {"bands": [], "size": [], "portal": [], "edges": [],
                "n_finite": int(finite.sum())}
    log_f = np.log10(fits[finite] + 1e-12)
    edges = np.quantile(log_f, np.linspace(0, 1, n_bands + 1))
    edges[-1] += 1e-6
    band_of = np.full(len(trees), -1)
    band_of[finite] = np.clip(np.digitize(log_f, edges) - 1, 0, n_bands - 1)

    out_bands, out_size, out_portal, out_n = [], [], [], []
    for b in range(n_bands):
        idx = np.where(band_of == b)[0]
        if len(idx) < 5:
            continue
        sig_counts = Counter(sigs[i] for i in idx if sigs[i] is not None)
        sizes = np.array(list(sig_counts.values()))
        mean_class = float((sizes ** 2).sum() / sizes.sum()) if sizes.sum() else 0.0
        n_probe = min(len(idx), 100)
        probe_idx = rng.choice(idx, size=n_probe, replace=False)
        portals = 0
        for i in probe_idx:
            mut = mutate(rng, trees[i], funcs, terms, 30, max_depth=3)
            mf = fitness(mut, X, y)
            if not np.isfinite(mf):
                continue
            mb_log = np.log10(mf + 1e-12)
            mb = int(np.clip(np.digitize([mb_log], edges)[0] - 1, 0, n_bands - 1))
            if mb < b:
                portals += 1
        out_bands.append(int(b))
        out_size.append(mean_class)
        out_portal.append(portals / n_probe)
        out_n.append(int(len(idx)))
    return {"bands": out_bands, "size": out_size, "portal": out_portal,
            "n_in_band": out_n, "edges": edges.tolist(),
            "n_finite": int(finite.sum())}


# -- target dispatch -------------------------------------------------------
def make_target(expr, n_vars):
    table = {
        "x0**2 + x1**2 + x2": lambda X: X[:, 0]**2 + X[:, 1]**2 + X[:, 2],
        "x0*x1 + x2*x0":      lambda X: X[:, 0]*X[:, 1] + X[:, 2]*X[:, 0],
        "x0 / (1 + x1**2)":   lambda X: X[:, 0] / (1 + X[:, 1]**2),
        "(x0 + x1) * (x2 + x0) * (x1 + x2)":
            lambda X: (X[:, 0]+X[:, 1])*(X[:, 2]+X[:, 0])*(X[:, 1]+X[:, 2]),
        "sin(x0) + x1**2":    lambda X: np.sin(X[:, 0]) + X[:, 1]**2,
        "sin(x0)*x1 + x2":    lambda X: np.sin(X[:, 0])*X[:, 1] + X[:, 2],
    }
    return table[expr]


def run(params: dict) -> dict:
    kind = params.get("kind", "gp_run")
    target = params["target"]
    target_id = target["id"]
    cond = params["condition"]
    seed = int(params["seed"])
    rng = np.random.default_rng(seed)
    target_fn = make_target(target["expr"], target["n_vars"])

    if kind == "egraph_banded":
        out = egraph_banded_measure(
            rng, target_fn, target["n_vars"], cond,
            n_samples=params.get("n_samples", 4000),
            max_depth=params.get("max_depth", 4),
            n_bands=params.get("n_bands", 8),
        )
        return {"kind": "egraph_banded", "target_id": target_id, "condition": cond,
                "seed": seed, **out}

    # gp_run: one rep of one (target, condition)
    pop_size = params.get("pop_size", 200)
    max_gen = params.get("max_generations", 1500)
    size_cap = params.get("size_cap", 30)
    noise_sigma = params.get("noise_sigma", 0.05)

    best_hist, visited = run_gp(rng, target_fn, target["n_vars"], cond,
                                 pop_size, max_gen, size_cap, noise_sigma)
    events = plateau_events(best_hist)
    # subsample best_hist for payload size
    n = len(best_hist)
    if n <= 600:
        bh_idx = list(range(n))
    else:
        bh_idx = sorted(set(list(range(0, 200)) +
                            list(range(200, n, max(1, (n - 200) // 400)))))
    bh_sample = [float(best_hist[i]) for i in bh_idx]

    return {
        "kind": "gp_run",
        "target_id": target_id,
        "condition": cond,
        "rep_idx": int(params.get("rep_idx", 0)),
        "seed": seed,
        "n_gens": int(max_gen),
        "best_hist_idx": bh_idx,
        "best_hist": bh_sample,
        "final_best": float(best_hist[-1]),
        "plateau_events": events,         # (start, end, dur, fit_level, observed)
        "visited_samples": visited[:200], # (gen, fit, tree_str)
    }
