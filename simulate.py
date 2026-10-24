"""Two independent continuous outcomes; fixed reproducible scenario grid.
Run: python simulate.py --workers 4. Requires numpy scipy pandas and g++.
No patient data. All effects are mean differences in population RMS-SD units.
"""
from pathlib import Path
import argparse
import concurrent.futures
import ctypes
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd
import scipy
from scipy import stats
ROOT = Path(__file__).resolve().parent
PUBLISHED = ROOT / 'results' / 'published'
OUTPUT = ROOT / 'results' / 'reproduced'
LIBRARY = ROOT / '.build' / 'resample.so'
SEED = 20261024
NS = [10, 20, 30, 50, 100]
DISTS = ['normal', 'lognormal', 't3']
DESIGNS = {'equal': (1.0, 1.0, 1), 'hetero_balanced': (2.0, 1.0, 1), 'hetero_unbalanced': (2.0, 1.0, 2)}
DELTAS = [0.0, 0.5, 0.8]
R = 10000
RR = 10000
B = 1999

def generate(rng, shape, dist):
    """Generate centered observations with population variance one."""
    if dist == 'normal':
        return rng.normal(size=shape)
    if dist == 't3':
        return rng.standard_t(3, size=shape) / np.sqrt(3)
    if dist == 'lognormal':
        return (rng.lognormal(0, 1, size=shape) - np.exp(0.5)) / np.sqrt((np.e - 1) * np.e)
    raise ValueError(f'Unknown distribution: {dist}')

def analytic(x, y):
    """Mean difference, standard errors, Welch df, and two-sided t p-values."""
    nx = x.shape[-1]
    ny = y.shape[-1]
    d = x.mean(-1) - y.mean(-1)
    v1 = x.var(-1, ddof=1)
    v2 = y.var(-1, ddof=1)
    se = np.sqrt(v1 / nx + v2 / ny)
    df = (v1 / nx + v2 / ny) ** 2 / ((v1 / nx) ** 2 / (nx - 1) + (v2 / ny) ** 2 / (ny - 1))
    sp = np.sqrt(((nx - 1) * v1 + (ny - 1) * v2) / (nx + ny - 2) * (1 / nx + 1 / ny))
    return (d, se, df, sp, 2 * stats.t.sf(np.abs(d / sp), nx + ny - 2), 2 * stats.t.sf(np.abs(d / se), df))

def rate_metrics(reject):
    p = float(np.mean(reject))
    return {'rejection': p, 'rejection_mcse': float(np.sqrt(p * (1 - p) / len(reject)))}

def interval_metrics(lo, hi, d, target):
    cover = (lo <= target) & (hi >= target)
    p = float(cover.mean())
    return {'coverage': p, 'coverage_mcse': float(np.sqrt(p * (1 - p) / len(d))), 'width': float(np.mean(hi - lo)), 'bias': float(d.mean() - target), 'bias_mcse': float(d.std(ddof=1) / np.sqrt(len(d)))}

def scenario(item):
    """Apply seven procedures to the same independently seeded datasets."""
    idx, (n, dist, design, delta) = item
    ss = np.random.SeedSequence([SEED, idx])
    data_ss, resample_ss = ss.spawn(2)
    rng = np.random.default_rng(data_ss)
    s1, s2, mult = DESIGNS[design]
    ny = n * mult
    target = delta * np.sqrt((s1 * s1 + s2 * s2) / 2)
    x = np.ascontiguousarray(s1 * generate(rng, (R, n), dist) + target)
    y = np.ascontiguousarray(s2 * generate(rng, (R, ny), dist))
    d, se, df, sp, ps, pw = analytic(x, y)
    pm = stats.mannwhitneyu(x, y, axis=1, alternative='two-sided', method='asymptotic', use_continuity=True).pvalue
    base = {'scenario': idx, 'n1': n, 'n2': ny, 'distribution': dist, 'design': design, 'variance_ratio': s1 * s1 / (s2 * s2), 'delta': delta, 'target': target}
    rows = []
    for name, p in [('Student', ps), ('Welch', pw), ('Mann_Whitney', pm)]:
        row = base | {'method': name, 'repetitions': R} | rate_metrics(p < 0.05)
        if name != 'Mann_Whitney':
            margin = stats.t.ppf(0.975, n + ny - 2) * sp if name == 'Student' else stats.t.ppf(0.975, df) * se
            row |= interval_metrics(d - margin, d + margin, d, target)
        rows.append(row)
    lib = ctypes.CDLL(str(LIBRARY))
    f = lib.resample
    arr = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    f.argtypes = [arr, arr, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_uint64, arr, arr, arr, arr]
    f.restype = None
    pdiff = np.empty(RR)
    pt = np.empty(RR)
    bd = np.empty((RR, B))
    bt = np.empty((RR, B))
    seed = int(resample_ss.generate_state(1, dtype=np.uint64)[0])
    f(x[:RR], y[:RR], RR, n, ny, B, seed, pdiff, pt, bd, bt)
    if not np.isfinite(bt).all():
        raise ValueError('Degenerate bootstrap sample')
    for name, p in [('Permutation_raw', pdiff), ('Permutation_studentized', pt)]:
        rows.append(base | {'method': name, 'repetitions': RR} | rate_metrics(p < 0.05))
    bq = np.quantile(bd, [0.025, 0.975], axis=1, method='linear')
    tq = np.quantile(bt, [0.025, 0.975], axis=1, method='linear')
    for name, lo, hi in [('Bootstrap_percentile', bq[0], bq[1]), ('Bootstrap_t', d[:RR] - tq[1] * se[:RR], d[:RR] - tq[0] * se[:RR])]:
        rows.append(base | {'method': name, 'repetitions': RR} | rate_metrics((lo > 0) | (hi < 0)) | interval_metrics(lo, hi, d[:RR], target))
    return rows

def corrections(p):
    """Unadjusted, Bonferroni, Holm, and BH decisions; hypotheses along axis 1."""
    m = p.shape[1]
    order = np.argsort(p, axis=1)
    s = np.take_along_axis(p, order, axis=1)
    bon = s < 0.05 / m
    holm = np.logical_and.accumulate(s < 0.05 / np.arange(m, 0, -1), axis=1)
    ranks = np.arange(1, m + 1)
    k = np.max(np.where(s < 0.05 * ranks / m, ranks, 0), axis=1)
    bh = ranks[None, :] <= k[:, None]
    inv = np.argsort(order, axis=1)
    return {'None': p < 0.05, 'Bonferroni': np.take_along_axis(bon, inv, axis=1), 'Holm': np.take_along_axis(holm, inv, axis=1), 'BH': np.take_along_axis(bh, inv, axis=1)}

def multiple():
    """Simulate independent endpoint families under global and partial nulls."""
    rows = []
    m = 20
    reps = 10000
    for j, (dist, n, design, k) in enumerate(((di, nn, de, kk) for di in ['normal', 'lognormal'] for nn in [10, 30] for de in ['equal', 'hetero_unbalanced'] for kk in [0, 4])):
        rng = np.random.default_rng(np.random.SeedSequence([SEED, 999, j]))
        a, b, mult = DESIGNS[design]
        ny = n * mult
        x = a * generate(rng, (reps, m, n), dist)
        y = b * generate(rng, (reps, m, ny), dist)
        x[:, :k, :] += 0.8 * np.sqrt((a * a + b * b) / 2)
        *_, ps, pw = analytic(x, y)
        for method, p in [('Student', ps), ('Welch', pw)]:
            for correction, reject in corrections(p).items():
                v = reject[:, k:].sum(1)
                tot = reject.sum(1)
                fdr = v / np.maximum(tot, 1)
                fw = v > 0
                power = reject[:, :k].mean(1) if k else np.full(reps, np.nan)
                rows.append({'distribution': dist, 'n1': n, 'n2': ny, 'design': design, 'nonnull': k, 'm': m, 'method': method, 'correction': correction, 'repetitions': reps, 'FWER': float(fw.mean()), 'FWER_mcse': float(fw.std(ddof=1) / np.sqrt(reps)), 'FDR': float(fdr.mean()), 'FDR_mcse': float(fdr.std(ddof=1) / np.sqrt(reps)), 'power': float(power.mean())})
    return pd.DataFrame(rows)

def output_directory(value=None, smoke=False):
    """Keep shipped reference results immutable during reproduction."""
    output = Path(value).resolve() if value else OUTPUT
    if smoke:
        output = output / 'smoke'
    if output == PUBLISHED or PUBLISHED in output.parents:
        raise ValueError('Choose an output directory outside results/published')
    output.mkdir(parents=True, exist_ok=True)
    return output

def build_resampler():
    """Build the native C++17 resampling routine; shell execution is never used."""
    LIBRARY.parent.mkdir(parents=True, exist_ok=True)
    configured = os.environ.get('CXX')
    compiler = shlex.split(configured) if configured else [shutil.which('g++') or shutil.which('clang++') or '']
    if not compiler[0]:
        raise RuntimeError('A C++17 compiler is required: install g++ or clang++')
    linking = '-dynamiclib' if sys.platform == 'darwin' else '-shared'
    command = compiler + ['-O3', '-std=c++17', linking, '-fPIC', str(ROOT / 'resample.cpp'), '-o', str(LIBRARY)]
    subprocess.run(command, check=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--output-dir', type=Path, help='Default: results/reproduced; smoke runs use a smoke subdirectory')
    args = ap.parse_args()
    if args.workers < 1:
        ap.error('--workers must be positive')
    output = output_directory(args.output_dir, args.smoke)
    build_resampler()
    grid = [(n, dist, de, delta) for n in NS for dist in DISTS for de in DESIGNS for delta in DELTAS]
    if args.smoke:
        grid = grid[:1]
    start = time.monotonic()
    allrows = []
    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        jobs = {pool.submit(scenario, item): item[0] for item in enumerate(grid)}
        for fut in concurrent.futures.as_completed(jobs):
            allrows.extend(fut.result())
            checkpoint = output / 'single_tests.partial.csv'
            pd.DataFrame(allrows).sort_values(['scenario', 'method']).to_csv(checkpoint, index=False)
            print(f'{len(allrows) // 7}/{len(grid)} scenarios; {time.monotonic() - start:.1f}s', flush=True)
    checkpoint.replace(output / 'single_tests.csv')
    if not args.smoke:
        multiple().to_csv(output / 'multiple_tests.csv', index=False)
    meta = {'seed': SEED, 'R_analytic': R, 'R_resampling': RR, 'B': B, 'scenarios': len(grid), 'smoke': args.smoke, 'numpy': np.__version__, 'scipy': scipy.__version__, 'pandas': pd.__version__, 'python': platform.python_version(), 'elapsed_seconds': time.monotonic() - start}
    (output / 'run_metadata.json').write_text(json.dumps(meta, indent=2))
    print(meta, flush=True)
if __name__ == '__main__':
    main()
