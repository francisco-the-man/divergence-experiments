
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
    # Rough Mt. Fuji: F(g) = -theta * d_H(g, ref) + eta(g), eta iid N(0,1), cached by genome.
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
    # sample waiting time ~ Exp(sum pi) (dimensionless), pick which neighbor fixes.
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
    # Hill tail index on the top k_frac of values. NaN if too few positive samples.
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
    # KS statistic against Exp(mean=x.mean()); small => more exponential-like.
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if len(x) < 8:
        return float('nan')
    scale = x.mean()
    if scale <= 0:
        return float('nan')
    stat, _ = sstats.kstest(x, 'expon', args=(0, scale))
    return float(stat)


# Fixed log10(stasis-length) edges for CCDF features so the per-task feature dimension is
# constant across families and seeds. Range -2..6 covers single-step epochs through very
# long stases without depending on per-task data.
_CCDF_LOG10_EDGES = np.linspace(-2.0, 6.0, 21)


def _ccdf_features(lengths):
    # 20 values: empirical P(log10(stasis) > edge) at each interior edge.
    x = np.asarray(lengths, dtype=float)
    x = x[x > 0]
    if len(x) < 1:
        return [float('nan')] * 20
    lx = np.log10(x)
    n = float(len(lx))
    return [float(np.sum(lx > e) / n) for e in _CCDF_LOG10_EDGES[1:]]


def _jump_log_stats(j):
    j = np.asarray(j, dtype=float)
    j = j[j > 0]
    if len(j) < 4:
        return float('nan'), float('nan')
    lj = np.log(j)
    return float(lj.mean()), float(lj.var())


def _stasis_jump_rankcorr(s, j):
    # Sibani record-dynamics signature: negative ρ in glassy regimes.
    s = np.asarray(s, dtype=float)
    j = np.asarray(j, dtype=float)
    m = (s > 0) & (j > 0)
    if m.sum() < 6:
        return float('nan')
    rho, _ = sstats.spearmanr(np.log(s[m]), np.log(j[m]))
    return float(rho) if np.isfinite(rho) else float('nan')


def _lag1_autocorr_log(x):
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if len(x) < 6:
        return float('nan')
    lx = np.log(x) - np.log(x).mean()
    d = float((lx * lx).sum())
    if d <= 0:
        return float('nan')
    return float((lx[:-1] * lx[1:]).sum() / d)


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
    """One origin-fixation trajectory on one landscape; returns stasis-length distribution,
    jump magnitudes, waiting times, entry contexts, and the expanded per-trajectory feature
    vector (4 originals + 20-bin CCDF + jump log-mean/log-var + stasis-jump rank-corr +
    lag-1 autocorr of log waiting times). Aggregator across tasks trains the classifier."""
    seed = int(params['seed'])
    cell = params['landscape_family_and_param']
    N = int(params.get('N', 20))
    Ne = int(params.get('pop_size', 1000))
    max_subs = int(params.get('max_substitutions', 2000))
    # Scale-relative jump threshold so it transfers across NK (f in [0,1]) and RMF (f unbounded).
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

    # Saddle-aware stasis detector: an epoch is a run of fixations whose relative fitness
    # change vs epoch-start stays within eps_jump_rel; it ends when a fixation produces a
    # relative gain > eps_jump_rel (the "jump"). Stasis length = #subs in epoch (>=1).
    # We record jump magnitude (|Δf|), accumulated waiting time, and the entry context
    # (saddle/plateau/peak) at epoch start for the confound-check plot.
    stasis_lengths = []
    jump_mags = []
    waiting_times = []
    entry_contexts = []

    epoch_start_f = f_cur
    epoch_subs = 0
    epoch_wait = 0.0
    entry_ctx = _neighbor_context(genome, f_cur, fitness_fn, N,
                                  max(abs(f_cur), 1e-6) * eps_jump_rel)

    for _ in range(max_subs):
        step = _origin_fixation_step(genome, f_cur, fitness_fn, N, Ne, rng)
        if step is None:
            break
        f_new, wt, _idx = step
        epoch_subs += 1
        epoch_wait += wt
        denom = max(abs(epoch_start_f), 1e-6)
        rel = (f_new - epoch_start_f) / denom
        if rel > eps_jump_rel:
            stasis_lengths.append(epoch_subs)
            jump_mags.append(abs(f_new - epoch_start_f))
            waiting_times.append(epoch_wait)
            entry_contexts.append(entry_ctx)
            f_cur = f_new
            epoch_start_f = f_new
            epoch_subs = 0
            epoch_wait = 0.0
            entry_ctx = _neighbor_context(genome, f_cur, fitness_fn, N,
                                          max(abs(f_cur), 1e-6) * eps_jump_rel)
        else:
            f_cur = f_new
    # Trailing unclosed epoch dropped to avoid right-censoring bias (consistent with parent).

    base = _summary_stats(stasis_lengths)
    ccdf = _ccdf_features(stasis_lengths)
    jm_lmean, jm_lvar = _jump_log_stats(jump_mags)
    sj_corr = _stasis_jump_rankcorr(stasis_lengths, jump_mags)
    # Plan asks for lag-1 autocorr of log waiting times; fall back to log stasis lengths
    # if waiting-times sample is too small (rare-trajectory regimes).
    autocorr = _lag1_autocorr_log(waiting_times if len(waiting_times) >= 6 else stasis_lengths)

    sl = np.asarray(stasis_lengths)
    ec = np.asarray(entry_contexts)
    stratified = {ctx: (sl[ec == ctx].tolist() if (ec == ctx).any() else [])
                  for ctx in ('saddle', 'plateau', 'peak')}

    features = {
        'log_mean': base['log_mean'],
        'log_var': base['log_var'],
        'hill': base['hill'],
        'ks_expon': base['ks_expon'],
    }
    for i, v in enumerate(ccdf):
        features[f'ccdf_b{i:02d}'] = v
    features['jump_log_mean'] = jm_lmean
    features['jump_log_var'] = jm_lvar
    features['stasis_jump_rankcorr'] = sj_corr
    features['lag1_autocorr_logwt'] = autocorr

    # Feature-group tags for the aggregator's ablation and permutation-importance plots.
    feature_groups = {'log_mean': 'original', 'log_var': 'original',
                      'hill': 'original', 'ks_expon': 'original',
                      'jump_log_mean': 'jump', 'jump_log_var': 'jump',
                      'stasis_jump_rankcorr': 'jump', 'lag1_autocorr_logwt': 'autocorr'}
    for i in range(20):
        feature_groups[f'ccdf_b{i:02d}'] = 'ccdf'

    return {
        'seed': seed,
        'cell': cell,
        'family': family,
        'param_value': float(param_value) if family == 'RMF' else int(param_value),
        'n_epochs': int(len(stasis_lengths)),
        'n_substitutions': int(sum(stasis_lengths)),
        'stasis_lengths': [int(x) for x in stasis_lengths],
        'jump_magnitudes': [float(x) for x in jump_mags],
        'waiting_times': [float(x) for x in waiting_times],
        'entry_contexts': entry_contexts,
        'stasis_lengths_stratified': {k: [int(x) for x in v] for k, v in stratified.items()},
        'features': features,
        'feature_groups': feature_groups,
        # Echo the bin label so the aggregator can build the 8-way classification target
        # without re-parsing the cell string.
        'param_bin_label': cell,
        'family_label': family,
    }
