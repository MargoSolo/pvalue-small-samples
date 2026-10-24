"""Selection bias of the same arithmetic mean-difference estimator.
No extra simulation outcomes are sampled; regenerate the identical streams.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.integrate import quad
from simulate import output_directory, SEED, NS, DISTS, DESIGNS, DELTAS, R, generate, analytic

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    output = output_directory(args.output_dir)
    rows = []
    grid = [(n, di, de, delta) for n in NS for di in DISTS for de in DESIGNS for delta in DELTAS]
    for idx, (n, dist, design, delta) in enumerate(grid):
        data_ss, _ = np.random.SeedSequence([SEED, idx]).spawn(2)
        rng = np.random.default_rng(data_ss)
        a, b, mult = DESIGNS[design]
        ny = n * mult
        target = delta * np.sqrt((a * a + b * b) / 2)
        x = a * generate(rng, (R, n), dist) + target
        y = b * generate(rng, (R, ny), dist)
        d, *_, ps, pw = analytic(x, y)
        for method, p in [('Student', ps), ('Welch', pw)]:
            chosen = d[p < 0.05]
            rows.append({'scenario': idx, 'n1': n, 'n2': ny, 'distribution': dist, 'design': design, 'delta': delta, 'target': target, 'method': method, 'unselected_mean': d.mean(), 'unselected_bias': d.mean() - target, 'selected_n': len(chosen), 'selected_mean': chosen.mean(), 'selected_bias': chosen.mean() - target, 'selected_mean_mcse': chosen.std(ddof=1) / np.sqrt(len(chosen))})
    pd.DataFrame(rows).to_csv(output / 'selection_bias.csv', index=False)
    mean = np.exp(0.5)

    def integrand(z):
        threshold = 2 * np.exp(z) - mean
        return stats.norm.pdf(z) * (1 if threshold <= 0 else stats.norm.sf(np.log(threshold)))
    theta, error = quad(integrand, -12, 12, epsabs=1e-11, points=[np.log(mean / 2)])
    (output / 'estimands.json').write_text(json.dumps({'lognormal_heteroscedastic_null_mean': 0, 'P_X_less_Y': theta, 'quadrature_absolute_error_estimate': error, 'explanation': 'X=2*(L1-E[L])/SD[L]; Y=(L2-E[L])/SD[L]; L1,L2 iid Lognormal(0,1). Equal means do not imply P(X<Y)=0.5.'}, indent=2))
if __name__ == '__main__':
    main()
