"""Regression-check the public code against shipped reference results."""
import concurrent.futures
import json

import pandas as pd

import simulate as sim


def normalized(frame):
    return frame.sort_values(['scenario', 'method']).reset_index(drop=True).sort_index(axis=1)


def main():
    sim.build_resampler()
    single = pd.read_csv(sim.PUBLISHED / 'single_tests.csv')
    grid = [(n, dist, design, delta) for n in sim.NS for dist in sim.DISTS
            for design in sim.DESIGNS for delta in sim.DELTAS]
    checks = []
    indices = [0, 15, 78, 123]
    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        for index, rows in zip(indices, pool.map(sim.scenario, [(i, grid[i]) for i in indices])):
            pd.testing.assert_frame_equal(
                normalized(pd.DataFrame(rows)), normalized(single[single.scenario == index]),
                check_exact=False, rtol=1e-12, atol=1e-13,
            )
            checks.append(f'All seven procedures match reference scenario {index}')
            print(checks[-1], flush=True)
    families = pd.read_csv(sim.PUBLISHED / 'multiple_tests.csv', keep_default_na=False)
    families['power'] = pd.to_numeric(families['power'], errors='coerce')
    pd.testing.assert_frame_equal(
        sim.multiple().sort_index(axis=1), families.sort_index(axis=1),
        check_exact=False, rtol=1e-12, atol=1e-13,
    )
    checks.append('All 128 multiple-testing rows match the reference')
    output = sim.output_directory()
    if not (output / 'selection_bias.csv').exists():
        raise FileNotFoundError('Run python supplementary.py before this regression check')
    pd.testing.assert_frame_equal(
        pd.read_csv(output / 'selection_bias.csv'),
        pd.read_csv(sim.PUBLISHED / 'selection_bias.csv'),
        check_exact=False, rtol=1e-12, atol=1e-13,
    )
    checks.append('All 270 selection-bias rows match the reference')
    try:
        sim.output_directory(sim.PUBLISHED)
    except ValueError:
        checks.append('Published-output overwrite protection works')
    else:
        raise AssertionError('Published directory was not protected')
    record = {'status': 'PASS', 'relative_tolerance': 1e-12,
              'absolute_tolerance': 1e-13, 'checks': checks}
    (output / 'repository_checks.json').write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
