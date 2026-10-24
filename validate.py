"""Independent checks of numerical formulae and resampling implementation."""
import argparse
import ctypes
import itertools
import json
from pathlib import Path
import numpy as np
from scipy import stats
from simulate import analytic, corrections, LIBRARY, build_resampler, output_directory
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    output = output_directory(args.output_dir)
    build_resampler()
    checks = []
    rng = np.random.default_rng(345983)
    for nx, ny in [(10, 10), (10, 20), (30, 50)]:
        x = rng.normal(size=(100, nx)) * 2
        y = rng.normal(size=(100, ny))
        *_, ps, pw = analytic(x, y)
        assert np.allclose(ps, stats.ttest_ind(x, y, axis=1, equal_var=True).pvalue, rtol=1e-12, atol=1e-14)
        assert np.allclose(pw, stats.ttest_ind(x, y, axis=1, equal_var=False).pvalue, rtol=1e-12, atol=1e-14)
    checks.append('Student and Welch p-values match scipy for 300 independent datasets')
    p = rng.uniform(size=(100, 6))
    p[:, :2] **= 6
    ours = corrections(p)
    assert np.array_equal(stats.false_discovery_control(p, axis=1, method='bh') < 0.05, ours['BH'])
    closed = np.ones_like(p, dtype=bool)
    for size in range(1, 7):
        for subset in itertools.combinations(range(6), size):
            local = p[:, subset].min(axis=1) * size < 0.05
            closed[:, subset] &= local[:, None]
    assert np.array_equal(closed, ours['Holm'])
    test = corrections(np.array([[0.001, 0.015, 0.04, 0.2]]))
    assert np.array_equal(test['Bonferroni'], [[True, False, False, False]])
    assert np.array_equal(test['Holm'], [[True, True, False, False]])
    checks.append('BH matches scipy; Holm matches exhaustive closed Bonferroni testing for 100 six-hypothesis families; fixed Bonferroni/Holm reference example matches')
    lib = ctypes.CDLL(str(LIBRARY))
    f = lib.resample
    arr = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    f.argtypes = [arr, arr, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_uint64, arr, arr, arr, arr]
    f.restype = None
    
    def call(x, y, B=100000, seed=76543):
        x = np.array([x], dtype=np.float64)
        y = np.array([y], dtype=np.float64)
        pd = np.empty(1)
        pt = np.empty(1)
        bd = np.empty((1, B))
        bt = np.empty((1, B))
        f(x, y, 1, x.shape[1], y.shape[1], B, seed, pd, pt, bd, bt)
        return (pd[0], pt[0], bd[0], bt[0])
    x = np.array([-0.4, 0.1, 0.8, 1.2, 2.0])
    y = np.array([-0.8, -0.3, 0.2, 0.6, 1.1, 1.5])
    z = np.r_[x, y]
    nx = len(x)
    ny = len(y)
    obs = x.mean() - y.mean()
    obs_t = stats.ttest_ind(x, y, equal_var=False).statistic
    cd = ct = 0
    total = 0
    for comb in itertools.combinations(range(len(z)), nx):
        a = z[list(comb)]
        b = np.delete(z, comb)
        dp = a.mean() - b.mean()
        tp = stats.ttest_ind(a, b, equal_var=False).statistic
        cd += abs(dp) >= abs(obs) - 1e-12
        ct += abs(tp) >= abs(obs_t) - 1e-12
        total += 1
    pd, pt, bd, bt = call(x, y)
    assert abs(pd - cd / total) < 0.008 and abs(pt - ct / total) < 0.008
    assert abs(bd.mean() - obs) < 0.01
    assert np.isfinite(bt).all()
    re = call(x, y)
    assert pd == re[0] and pt == re[1] and np.array_equal(bd, re[2]) and np.array_equal(bt, re[3])
    checks.append('C++ permutation probabilities agree within 0.008 with exhaustive enumeration of 462 partitions')
    checks.append('Bootstrap mean agrees with observed difference within 0.01; studentized draws finite')
    checks.append('Fixed seed reproduces all resampling outputs bit for bit')
    result = {'status': 'PASS', 'checks': checks, 'exact_raw_p': cd / total, 'sampled_raw_p': pd, 'exact_studentized_p': ct / total, 'sampled_studentized_p': pt}
    (output / 'validation.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
