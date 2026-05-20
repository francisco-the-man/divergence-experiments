
import numpy as np
from scipy import stats as sstats


def _nk_fitness_table(N, K, rng):
    # Each locus depends on itself + K other loci. Random epistatic table per locus.
    neighbors = np.zeros((N, K), dtype=np.int64)
    for i in range(N):
        choices = [j for j in range(N) if j != i]
        neighbors[i] = rng.choice(choices, size=K, replace=False)
    contribs = rng.random((N, 1 << (K + 1)))
    return neighbors, contribs


def _nk_fitness(genome, neighbors, contribs, K):
    N = genome.shape[0]
    idx = genome.astype(np.int64).copy()
    for k in range(K):
        idx = idx + (genome[neighbors[:, k]].astype(np.int64) << (k + 1))
    return float(contribs[np.arange(N), idx].mean())


def _make_rmf(N, theta, rng):
    # Rough Mt. Fuji: F(g) = -theta * d_H(g, ref) + eta(g), eta iid N(0,1) per genotype.
    # We lazily fill eta on demand and cache by genome bytes (small relative to 2^20).
    ref = rng.integers(0, 2, size=N, dtype=np.int8)
    cache = {}
    field_rng = np.random.default_rng(int(rng.integers(0, 2**31 - 1)))
    return ref, cache, theta, field_rng


def _rmf_fitness(genome, ref, cache, theta, field_rng):
    d = int(np.sum(genome != ref))
    key = genome.tobytes()
    eta = cache.get(key)
    if eta is None:
        eta = float(field_rng.standard_normal())
        cache[key] = eta
    return -theta * d + eta


def _origin_fixation_step(genome, f_cur, fitness_fn, N, Ne, rng):
    # Enumerate N single-mutant neighbors, compute Kimura haploid fixation prob per neighbor,
    # sample waiting time ~ Exp(sum pi) in units of 1/(mu_per_locus * Ne) (since mu*Ne cancels
    # in relative rates; we record dimensionless waiting times), then pick which to fix.
    fits = np.empty(N)
    for i in range(N):
        genome[i] ^= 1
        fits[i] = fitness_fn(genome)
        genome[i] ^= 1
    denom = max(abs(f_cur), 1e-6)
    s = (fits - f_cur) / denom
    pi = np.empty(N)
    for i in range(N):
        si = s[i]
        if abs(si) < 1e-9:
            pi[i] = 1.0 / Ne
        else:
            x = 2.0 * si
            y = 2.0 * Ne * si
            # Safe forms of pi(s) = (1-e^{-2s}) / (1-e^{-2 Ne s}) in three regimes.
            if y > 50:
                pi[i] = 1.0 - np.exp(-x) if x < 50 else 1.0
            elif y < -50:
                pi[i] = np.exp(y) * (np.exp(x) - 1.0) if abs(x) < 50 else 0.0
            else:
                num = 1.0 - np.exp(-x)
                den = 1.0 - np.exp(-y)
                pi[i] = num / den if den != 0 else 1.0 / Ne
        if pi[i] < 0:
            pi[i] = 0.0
    total = pi.sum()
    if total <= 0:
        return None
    wt = float(rng.exponential(1.0 / total))
    idx = int(rng.choice(N, p=pi / total))
    genome[idx] ^= 1
    return float(fits[idx]), wt, idx


def _hill_estimator(x, k_frac=0.25):
    # Hill tail index estimator on the top k_frac of values (one-sided heavy-tail).
    # Returns the inverse-shape estimate; NaN if too few positive samples.
    x = np.sort(np.asarray(x, dtype=float))
    x = x[x > 0]
    n = len(x)
    if n < 8:
        return float('nan')
    k = max(3, int(k_frac * n))
    if k >= n:
        k = n - 1
    top = x[-(k + 1):]
    xk = top[0]
    if xk <= 0:
        return float('nan')
    return float(np.mean(np.log(top[1:] / xk)))


def _ks_to_exponential(x):
    # KS statistic of x against Exp(mean=x.mean()). Small => more exponential-like.
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if len(x) < 8:
        return float('nan')
    scale = x.mean()
    if scale <= 0:
        return float('nan')
    stat, _ = sstats.kstest(x, 'expon', args=(0, scale))
    return float(stat)


def _summary_stats(lengths):
    x = np.asarray(lengths, dtype=float)
    x = x[x > 0]
    if len(x) < 8:
        return {'log_mean': float('nan'), 'log_var': float('nan'),
                'hill': float('nan'), 'ks_expon': float('nan'), 'n': int(len(x))}
    lx = np.log(x)
    return {
        'log_mean': float(lx.mean()),
        'log_var': float(lx.var()),
        'hill': _hill_estimator(x),
        'ks_expon': _ks_to_exponential(x),
        'n': int(len(x)),
    }


def _parse_cell(cell):
    if cell.startswith('NK_K'):
        return 'NK', int(cell[4:])
    if cell.startswith('RMF_theta'):
        return 'RMF', float(cell[len('RMF_theta'):])
    raise ValueError(f"unknown cell: {cell}")


def _neighbor_context(genome, f_cur, fitness_fn, N, eps):
    nf = np.empty(N)
    for i in range(N):
        genome[i] ^= 1
        nf[i] = fitness_fn(genome)
        genome[i] ^= 1
    nb = int(np.sum(nf > f_cur + eps))
    nn = int(np.sum(np.abs(nf - f_cur) <= eps))
    if nb > 0:
        return 'saddle'
    if nn > 0:
        return 'plateau'
    return 'peak'


def run(params: dict) -> dict:
    """One origin-fixation trajectory on one landscape; returns stasis-length distribution
    and four distributional summary stats. Aggregator across tasks trains the classifier."""
    seed = int(params['seed'])
    cell = params['landscape_family_and_param']
    N = int(params.get('N', 20))
    Ne = int(params.get('pop_size', 1000))
    max_subs = int(params.get('max_substitutions', 400))
    # A fixation counts as a beneficial "jump" (ends a stasis) when Δf > eps_jump_rel * |f_cur|.
    # Scale-relative so the threshold transfers across NK (f in [0,1]) and RMF (f unbounded).
    eps_jump_rel = float(params.get('eps_jump_rel', 0.005))

    rng = np.random.default_rng(seed)
    family, param_value = _parse_cell(cell)

    if family == 'NK':
        K = int(param_value)
        neighbors, contribs = _nk_fitness_table(N, K, rng)
        def fitness_fn(g, _n=neighbors, _c=contribs, _K=K):
            return _nk_fitness(g, _n, _c, _K)
    else:
        theta = float(param_value)
        ref, cache, _, field_rng = _make_rmf(N, theta, rng)
        def fitness_fn(g, _r=ref, _c=cache, _t=theta, _fr=field_rng):
            return _rmf_fitness(g, _r, _c, _t, _fr)

    genome = rng.integers(0, 2, size=N, dtype=np.int8)
    f_cur = float(fitness_fn(genome))
    f_init = f_cur

    # Saddle-aware stasis detector: stasis epoch = consecutive fixations with Δf ≤ eps
    # terminating when a beneficial (Δf > eps) fixation occurs. Entry context is recorded
    # at the start of each epoch by examining single-mutant neighbors (saddle if any
    # beneficial neighbor exists, plateau if only neutral neighbors exist, peak otherwise).
    stasis_lengths_subs = []
    stasis_lengths_time = []
    stasis_contexts = []
    jump_magnitudes = []

    cur_subs = 0
    cur_time = 0.0
    pending_context = _neighbor_context(
        genome, f_cur, fitness_fn, N, eps_jump_rel * max(abs(f_cur), 1e-6)
    )

    for _step in range(max_subs):
        out = _origin_fixation_step(genome, f_cur, fitness_fn, N, Ne, rng)
        if out is None:
            break
        f_new, wt, _idx = out
        eps = eps_jump_rel * max(abs(f_cur), 1e-6)
        df = f_new - f_cur
        cur_time += wt
        cur_subs += 1
        if df > eps:
            stasis_lengths_subs.append(cur_subs)
            stasis_lengths_time.append(cur_time)
            stasis_contexts.append(pending_context)
            jump_magnitudes.append(float(df))
            cur_subs = 0
            cur_time = 0.0
            pending_context = _neighbor_context(
                genome, f_new, fitness_fn, N, eps_jump_rel * max(abs(f_new), 1e-6)
            )
        f_cur = f_new

    # Drop trailing partial epoch (right-censored) — already done by not appending after loop.

    stats_subs = _summary_stats(stasis_lengths_subs)
    stats_time = _summary_stats(stasis_lengths_time)
    saddle_lengths_subs = [L for L, c in zip(stasis_lengths_subs, stasis_contexts) if c == 'saddle']
    saddle_lengths_time = [L for L, c in zip(stasis_lengths_time, stasis_contexts) if c == 'saddle']
    stats_saddle_subs = _summary_stats(saddle_lengths_subs)
    stats_saddle_time = _summary_stats(saddle_lengths_time)

    return {
        'seed': seed,
        'family': family,
        'cell': cell,
        'param_value': float(param_value),
        'n_epochs': int(len(stasis_lengths_subs)),
        'n_substitutions': int(max_subs),
        'initial_fitness': f_init,
        'final_fitness': float(f_cur),
        'stasis_lengths_subs': [int(x) for x in stasis_lengths_subs],
        'stasis_lengths_time': [float(x) for x in stasis_lengths_time],
        'stasis_contexts': stasis_contexts,
        'jump_magnitudes': jump_magnitudes,
        'summary_stats': stats_subs,
        'summary_stats_time': stats_time,
        'summary_stats_saddle': stats_saddle_subs,
        'summary_stats_saddle_time': stats_saddle_time,
        'n_saddle_epochs': int(len(saddle_lengths_subs)),
        'eps_jump_rel': eps_jump_rel,
    }
