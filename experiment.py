
import numpy as np
from scipy import stats as sstats


# ---------- landscapes (unchanged from parent) ----------

def _nk_fitness_table(N, K, rng):
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


# ---------- Kimura origin-fixation step (parent's body; extra return value: pre-state
# neighbor-fitness vector, which the v3 context rule and eps-sensitivity re-derivation need) ----------

def _origin_fixation_step(genome, f_cur, fitness_fn, N, Ne, rng):
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
        return None, None, None, fits
    wt = float(rng.exponential(1.0 / total))
    idx = int(rng.choice(N, p=pi / total))
    genome[idx] ^= 1
    return float(fits[idx]), wt, idx, fits


# ---------- v3 saddle-aware context rule (NEW; replaces parent's _neighbor_context) ----------
# "saddle if frac_strictly_beneficial > 1/(2N); peak if all neighbors strictly worse;
#  plateau if majority within ±eps; else mixed"

def _neighbor_context_v3(nf, f_cur, N, eps_abs):
    frac_ben = float(np.sum(nf > f_cur)) / N
    if frac_ben > 1.0 / (2.0 * N):
        return 'saddle'
    if np.all(nf < f_cur):
        return 'peak'
    if float(np.sum(np.abs(nf - f_cur) <= eps_abs)) / N > 0.5:
        return 'plateau'
    return 'mixed'


# ---------- feature primitives ----------

def _hill_estimator(x, k_frac=0.25):
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
    x = np.asarray(x, dtype=float)
    x = x[x > 0]
    if len(x) < 8:
        return float('nan')
    scale = x.mean()
    if scale <= 0:
        return float('nan')
    stat, _ = sstats.kstest(x, 'expon', args=(0, scale))
    return float(stat)


def _summary_stats_v3(x_pos, n_log_min=4, n_hill_min=8):
    """v3 relaxed thresholds: log_mean / log_var require n>=4; hill and KS require n>=8."""
    x = np.asarray(x_pos, dtype=float)
    x = x[x > 0]
    n = len(x)
    out = {'log_mean': float('nan'), 'log_var': float('nan'),
           'hill': float('nan'), 'ks_expon': float('nan'), 'n': int(n)}
    if n >= n_log_min:
        lx = np.log(x)
        out['log_mean'] = float(lx.mean())
        out['log_var'] = float(lx.var())
    if n >= n_hill_min:
        out['hill'] = _hill_estimator(x)
        out['ks_expon'] = _ks_to_exponential(x)
    return out


def _summary_stats_v2_rule(x_pos):
    """v2's rule (all four require n>=8). Used only for the NaN-rate validation test."""
    x = np.asarray(x_pos, dtype=float)
    x = x[x > 0]
    n = len(x)
    if n < 8:
        return {'log_mean': float('nan'), 'log_var': float('nan'),
                'hill': float('nan'), 'ks_expon': float('nan'), 'n': int(n)}
    lx = np.log(x)
    return {'log_mean': float(lx.mean()), 'log_var': float(lx.var()),
            'hill': _hill_estimator(x), 'ks_expon': _ks_to_exponential(x), 'n': int(n)}


def _ccdf_features(values_pos, log10_lo, log10_hi, nbins):
    """Fixed-edge CCDF on log10(values), returning nbins values P(log10(v) > edge_i) at
    nbins interior edges spanning [log10_lo, log10_hi]. NaN vector if no positive samples."""
    x = np.asarray(values_pos, dtype=float)
    x = x[x > 0]
    if len(x) == 0:
        return [float('nan')] * nbins
    lx = np.log10(x)
    edges = np.linspace(log10_lo, log10_hi, nbins + 1)[1:]
    n = float(len(lx))
    return [float(np.sum(lx > e) / n) for e in edges]


def _ccdf_features_natural_log(values_pos, log_lo, log_hi, nbins):
    """Same shape as _ccdf_features but on natural-log axis — jump-magnitude range is
    given in ln per the plan's ccdf_jump_log_range = [-8, 4]."""
    x = np.asarray(values_pos, dtype=float)
    x = x[x > 0]
    if len(x) == 0:
        return [float('nan')] * nbins
    lx = np.log(x)
    edges = np.linspace(log_lo, log_hi, nbins + 1)[1:]
    n = float(len(lx))
    return [float(np.sum(lx > e) / n) for e in edges]


def _jump_log_stats(j):
    j = np.asarray(j, dtype=float)
    j = j[j > 0]
    if len(j) < 4:
        return float('nan'), float('nan')
    lj = np.log(j)
    return float(lj.mean()), float(lj.var())


def _stasis_jump_rankcorr(s, j):
    s = np.asarray(s, dtype=float)
    j = np.asarray(j, dtype=float)
    m = (s > 0) & (j > 0)
    if m.sum() < 6:
        return float('nan')
    # Guard against constant inputs (e.g. all stasis_lengths == 1 in glassy regimes);
    # spearmanr returns nan in that case anyway, but we suppress the noisy warning.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
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


def _parse_cell(cell):
    if cell.startswith('NK_K'):
        return 'NK', int(cell[4:])
    if cell.startswith('RMF_theta'):
        return 'RMF', float(cell[len('RMF_theta'):])
    raise ValueError(f"unknown cell: {cell}")


# ---------- v3 feature builder (NEW; combines stasis-shape + jump-shape + corr/autocorr) ----------

def _build_features(stasis_obs, jump_mags, waiting_times,
                    ccdf_time_log10_range, ccdf_time_nbins,
                    ccdf_jump_log_range, ccdf_jump_nbins,
                    nan_thresholds):
    """One feature dict on ONE choice of 'stasis observable' (waiting times for the v3
    primary view, substitution counts for the matched secondary view). Jump-magnitude
    features are identical across both views — they don't depend on the observable."""
    n_lm = int(nan_thresholds.get('log_mean', 4))
    n_hill = int(nan_thresholds.get('hill', 8))
    base = _summary_stats_v3(stasis_obs, n_log_min=n_lm, n_hill_min=n_hill)
    ccdf_stasis = _ccdf_features(stasis_obs,
                                 ccdf_time_log10_range[0], ccdf_time_log10_range[1],
                                 ccdf_time_nbins)
    jm_lmean, jm_lvar = _jump_log_stats(jump_mags)
    jump_ccdf = _ccdf_features_natural_log(jump_mags,
                                           ccdf_jump_log_range[0], ccdf_jump_log_range[1],
                                           ccdf_jump_nbins)
    jump_hill = _hill_estimator(np.asarray(jump_mags, dtype=float))
    jump_ks = _ks_to_exponential(np.asarray(jump_mags, dtype=float))
    sj_corr = _stasis_jump_rankcorr(stasis_obs, jump_mags)
    autocorr = _lag1_autocorr_log(waiting_times if len(waiting_times) >= 6 else stasis_obs)

    feats = {
        'log_mean': base['log_mean'],
        'log_var': base['log_var'],
        'hill': base['hill'],
        'ks_expon': base['ks_expon'],
    }
    groups = {'log_mean': 'scalar', 'log_var': 'scalar',
              'hill': 'scalar', 'ks_expon': 'scalar'}
    for i, v in enumerate(ccdf_stasis):
        k = f'ccdf_stasis_b{i:02d}'
        feats[k] = v
        groups[k] = 'ccdf_stasis'
    feats['jump_log_mean'] = jm_lmean
    feats['jump_log_var'] = jm_lvar
    feats['jump_hill'] = jump_hill
    feats['jump_ks_expon'] = jump_ks
    groups.update({'jump_log_mean': 'jump_scalar', 'jump_log_var': 'jump_scalar',
                   'jump_hill': 'jump_scalar', 'jump_ks_expon': 'jump_scalar'})
    for i, v in enumerate(jump_ccdf):
        k = f'ccdf_jump_b{i:02d}'
        feats[k] = v
        groups[k] = 'ccdf_jump'
    feats['stasis_jump_rankcorr'] = sj_corr
    feats['lag1_autocorr_logwt'] = autocorr
    groups['stasis_jump_rankcorr'] = 'corr'
    groups['lag1_autocorr_logwt'] = 'autocorr'
    return feats, groups, base


# ---------- main task ----------

def run(params: dict) -> dict:
    """One origin-fixation trajectory on one landscape. v3 changes vs parent:
      - Build features twice at the primary eps: once on waiting times (PRIMARY observable)
        and once on substitution counts (matched SECONDARY observable). Same epoch
        decomposition; only the 'stasis quantity' differs.
      - 25-bin CCDF on log10(time) over [-2, 22] (fixed log10 edges); 20-bin CCDF on
        natural log |Δf| over [-8, 4] for the jump observable (NEW).
      - Add jump-shape scalars (jump_hill, jump_ks_expon) alongside jump_log_mean/var.
      - Saddle context follows v3's 4-category rule (saddle/plateau/peak/mixed).
      - Relaxed NaN thresholds: log_mean/log_var require n>=4; hill/KS require n>=8.
      - eps sensitivity ({0.001, 0.005, 0.02}) re-derived from the SAVED fixation
        trajectory — no re-simulation. The single simulated trajectory drives all three.
      - Report per-feature NaN indicators under v2's and v3's rules for fix-#7 validation.
    Aggregator across tasks trains the 8-way RandomForest classifier and runs the bootstrap /
    permutation-importance / chi-square / Spearman analyses."""
    seed = int(params['seed'])
    cell = params['landscape_family_and_param']
    N = int(params.get('N', 20))
    Ne = int(params.get('pop_size', 1000))
    max_subs = int(params.get('max_substitutions', 2000))
    eps_primary = float(params.get('eps_jump_rel_primary', 0.005))
    eps_sensitivity = list(params.get('eps_jump_rel_sensitivity', [0.001, 0.005, 0.02]))

    ccdf_time_range = list(params.get('ccdf_time_log10_range', [-2.0, 22.0]))
    ccdf_time_nbins = int(params.get('ccdf_time_nbins', 25))
    ccdf_jump_range = list(params.get('ccdf_jump_log_range', [-8.0, 4.0]))
    ccdf_jump_nbins = int(params.get('ccdf_jump_nbins', 20))
    nan_thr = dict(params.get('nan_thresholds',
                              {'log_mean': 4, 'log_var': 4, 'hill': 8, 'ks_expon_min_n': 8}))

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

    # Simulate the FULL fixation trajectory once; cache per-step neighbor-fit vector
    # so eps sensitivity (and v3 context evaluation at each epoch start) needs no re-sim.
    # Memory: 2000 steps * (8 floats N=20 + a few scalars) ≈ 350 KB; well under 1 GB.
    traj = []
    for _ in range(max_subs):
        prev_f = f_cur
        f_new, wt, _idx, nf_pre = _origin_fixation_step(genome, f_cur, fitness_fn, N, Ne, rng)
        if f_new is None:
            break  # trajectory dead-ended (Σ π ≤ 0; consistent with parent)
        traj.append({'prev_f': float(prev_f), 'f_new': float(f_new), 'wt': float(wt),
                     'neighbor_fits_at_entry': nf_pre.copy()})
        f_cur = f_new

    def derive_at_eps(eps):
        """Re-derive epochs from the saved trajectory at the given eps_jump_rel.
        Same per-step rule as the parent (relative gain vs epoch-start fitness > eps closes
        the epoch); trailing unclosed epoch dropped to avoid right-censoring bias."""
        sl, jm, wt_acc, ec = [], [], [], []
        if not traj:
            return sl, jm, wt_acc, ec
        epoch_start_idx = 0
        epoch_start_f = traj[0]['prev_f']
        for i, step in enumerate(traj):
            denom = max(abs(epoch_start_f), 1e-6)
            rel = (step['f_new'] - epoch_start_f) / denom
            if rel > eps:
                subs = i - epoch_start_idx + 1
                wt_sum = sum(traj[j]['wt'] for j in range(epoch_start_idx, i + 1))
                sl.append(subs)
                jm.append(abs(step['f_new'] - epoch_start_f))
                wt_acc.append(wt_sum)
                # Context evaluated AT the epoch-start state (neighbor fits cached there).
                nf = traj[epoch_start_idx]['neighbor_fits_at_entry']
                eps_abs = max(abs(epoch_start_f), 1e-6) * eps
                ec.append(_neighbor_context_v3(nf, epoch_start_f, N, eps_abs))
                if i + 1 < len(traj):
                    epoch_start_idx = i + 1
                    epoch_start_f = traj[i + 1]['prev_f']
                else:
                    epoch_start_idx = len(traj)
                    epoch_start_f = step['f_new']
        return sl, jm, wt_acc, ec

    # ----- primary eps: build BOTH the waiting-time feature view AND the matched
    # substitution-count view from the same epoch decomposition (key v3 adjudicator). -----
    sl_p, jm_p, wt_p, ec_p = derive_at_eps(eps_primary)

    feats_time_primary, groups_time, base_time = _build_features(
        wt_p, jm_p, wt_p,
        ccdf_time_range, ccdf_time_nbins,
        ccdf_jump_range, ccdf_jump_nbins, nan_thr,
    )
    feats_subs_primary, _, base_subs = _build_features(
        sl_p, jm_p, wt_p,
        ccdf_time_range, ccdf_time_nbins,
        ccdf_jump_range, ccdf_jump_nbins, nan_thr,
    )

    # ----- eps sensitivity: full waiting-time feature set at each eps (ROBUSTNESS test). -----
    eps_features = {}
    for eps in eps_sensitivity:
        sl_e, jm_e, wt_e, ec_e = derive_at_eps(eps)
        feats_e, _, _ = _build_features(
            wt_e, jm_e, wt_e,
            ccdf_time_range, ccdf_time_nbins,
            ccdf_jump_range, ccdf_jump_nbins, nan_thr,
        )
        eps_features[f'{eps}'] = {
            'features': feats_e,
            'n_epochs': int(len(wt_e)),
            'context_counts': {ctx: int(ec_e.count(ctx))
                               for ctx in ('saddle', 'plateau', 'peak', 'mixed')},
        }

    # ----- NaN-rate validation (fix #7): per-scalar indicators under v2's vs v3's rules. -----
    base_v2_time = _summary_stats_v2_rule(wt_p)
    base_v3_time = _summary_stats_v3(
        wt_p, n_log_min=int(nan_thr.get('log_mean', 4)),
        n_hill_min=int(nan_thr.get('hill', 8)),
    )
    nan_indicator_time = {
        k: {'v2_rule_nan': not np.isfinite(base_v2_time[k]),
            'v3_rule_nan': not np.isfinite(base_v3_time[k])}
        for k in ('log_mean', 'log_var', 'hill', 'ks_expon')
    }
    base_v2_subs = _summary_stats_v2_rule(sl_p)
    base_v3_subs = _summary_stats_v3(
        sl_p, n_log_min=int(nan_thr.get('log_mean', 4)),
        n_hill_min=int(nan_thr.get('hill', 8)),
    )
    nan_indicator_subs = {
        k: {'v2_rule_nan': not np.isfinite(base_v2_subs[k]),
            'v3_rule_nan': not np.isfinite(base_v3_subs[k])}
        for k in ('log_mean', 'log_var', 'hill', 'ks_expon')
    }

    context_counts_primary = {ctx: int(ec_p.count(ctx))
                              for ctx in ('saddle', 'plateau', 'peak', 'mixed')}

    return {
        'seed': seed,
        'cell': cell,
        'family': family,
        'param_value': float(param_value) if family == 'RMF' else int(param_value),
        'param_bin_label': cell,
        'family_label': family,
        'n_fixations': int(len(traj)),
        'n_epochs_primary': int(len(wt_p)),
        # Raw observables — for CCDF plots and per-family Spearman tests.
        'stasis_lengths': [int(x) for x in sl_p],
        'jump_magnitudes': [float(x) for x in jm_p],
        'waiting_times': [float(x) for x in wt_p],
        'entry_contexts': ec_p,
        # PRIMARY-test feature view (waiting times).
        'features_time_primary': feats_time_primary,
        # Matched SECONDARY-test feature view (substitution counts, same epochs).
        'features_subs_primary': feats_subs_primary,
        'feature_groups': groups_time,
        # ROBUSTNESS-test inputs.
        'eps_sensitivity': eps_features,
        'eps_primary': float(eps_primary),
        # NaN-rate validation inputs (fix #7).
        'nan_indicator_time': nan_indicator_time,
        'nan_indicator_subs': nan_indicator_subs,
        # Context-distribution chi-square inputs (v3 4-category rule).
        'context_counts_primary': context_counts_primary,
        # Scalar summaries for direct reporting.
        'base_time_v3': base_v3_time,
        'base_subs_v3': base_v3_subs,
    }
