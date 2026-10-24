"""Generate a summary and an English figure from complete CSV results."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from simulate import PUBLISHED, output_directory


def row(frame, **conditions):
    selected = frame
    for key, value in conditions.items():
        selected = selected[selected[key] == value]
    if len(selected) != 1:
        raise ValueError(f'Expected one result, got {len(selected)} for {conditions}')
    return selected.iloc[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=PUBLISHED)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    output = output_directory(args.output_dir)
    single = pd.read_csv(args.input_dir / 'single_tests.csv')
    families = pd.read_csv(args.input_dir / 'multiple_tests.csv', keep_default_na=False)
    selected = pd.read_csv(args.input_dir / 'selection_bias.csv')
    if len(single) != 945 or len(families) != 128:
        raise ValueError('A full run is required; smoke results are incomplete')
    base = dict(n1=10, design='hetero_unbalanced', delta=0)
    selection = row(selected, method='Welch', n1=10, distribution='normal', design='equal', delta=0.5)
    family = dict(method='Welch', n1=30, distribution='normal', design='equal', nonnull=0)
    summary = ['# Results summary', '',
               'Every value below is read from the CSV results.', '',
               '| Scenario | Procedure | Quantity | Result |',
               '| --- | --- | --- | --- |']
    for dist in ['normal', 'lognormal']:
        for method in ['Student', 'Welch', 'Permutation_studentized', 'Bootstrap_t']:
            value = row(single, method=method, distribution=dist, **base).rejection
            summary.append(f'| {dist}; n₁=10, n₂=20; variances 4:1 | {method} | Type I error | {value:.1%} |')
    for correction in ['None', 'Holm']:
        value = row(families, correction=correction, **family).FWER
        summary.append(f'| Normal; 20 independent null endpoints; n₁=n₂=30 | Welch + {correction} | FWER | {value:.1%} |')
    summary += ['',
        f'With normal equal-variance groups of 10 and true mean difference 0.5, the mean estimate was {selection.unselected_mean:.3f} before selection and {selection.selected_mean:.3f} among the {int(selection.selected_n)} Welch-significant results (p < 0.05).', '',
        'For unequal designs, the smaller group has the larger variance. These results do not establish a universal ranking. Mann–Whitney rejection under differing distributions is not automatically type I error for equal means.']
    (output / 'summary.md').write_text('\n'.join(summary) + '\n', encoding='utf-8')

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False})
    figure, axes = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True)
    methods = [('Student', '-', 'o'), ('Welch', '--', 's'),
               ('Permutation_studentized', ':', '^'), ('Bootstrap_t', '-.', 'D')]
    for axis, dist, title in zip(axes, ['normal', 'lognormal'], ['Normal', 'Lognormal']):
        for method, style, marker in methods:
            data = single[(single.method == method) & (single.distribution == dist) &
                          (single.design == 'hetero_unbalanced') & (single.delta == 0)].sort_values('n1')
            axis.plot(data.n1, 100 * data.rejection, color='black', linestyle=style,
                      marker=marker, markersize=4, markerfacecolor='white', linewidth=1,
                      label=method.replace('_', ' '))
        axis.axhline(5, color='gray', linewidth=0.8)
        axis.set_title(title)
        axis.set_xlabel('Smaller group size n (other group: 2n)')
        axis.set_xticks([10, 30, 50, 100])
        axis.set_ylim(0, 16)
        axis.grid(axis='y', alpha=0.2)
    axes[0].set_ylabel('Type I error (%)')
    figure.suptitle('Equal means; variance ratio 4:1; smaller group has larger variance', fontsize=10)
    figure.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=2,
                  frameon=False, fontsize=8)
    figure.subplots_adjust(left=0.08, right=0.99, top=0.82, bottom=0.25, wspace=0.15)
    figure.savefig(output / 'type1_error_en.png', dpi=300, bbox_inches='tight')
    plt.close(figure)
    # A focused figure for the selection example; the mean-difference scale
    # and zero baseline are identical for the selected and unselected bars.
    figure, axis = plt.subplots(figsize=(8, 7))
    figure.subplots_adjust(left=0.14, right=0.96, bottom=0.25, top=0.76)
    values = [selection.unselected_mean, selection.selected_mean]
    colors = ['#8C9BA8', '#176B87']
    axis.bar([0, 1], values, width=0.55, color=colors, zorder=3)
    axis.axhline(0.5, color='#303A43', linestyle='--', linewidth=1.3, zorder=4)
    axis.text(1.41, 0.51, 'True effect = 0.5', ha='right', va='bottom', fontsize=11)
    axis.set_xticks([0, 1], ['All simulated studies\nn = 10,000 studies',
                           f'Only Welch p < 0.05\nn = {int(selection.selected_n):,} studies'])
    axis.tick_params(axis='x', labelsize=12, pad=10)
    axis.set_ylabel('Average estimated mean difference', fontsize=12)
    axis.set_ylim(0, 1.32)
    axis.set_xlim(-0.55, 1.55)
    axis.grid(axis='y', alpha=0.15, zorder=0)
    for position, value in enumerate(values):
        axis.text(position, value + 0.035, f'{value:.3f}', ha='center', va='bottom',
                  fontsize=23, fontweight='bold', color=colors[position])
    figure.text(0.08, 0.92, 'Keep only “significant” results.\nWatch the estimated effect grow.',
                fontsize=23, fontweight='bold', color='#162B3A', va='top')
    figure.text(0.08, 0.805, 'Normal data · 10 observations per group · Equal variances',
                fontsize=12, color='#47545F')
    figure.text(0.08, 0.075,
                'The true mean difference is 0.5 in every simulated study.\n'
                'The second bar averages only results selected by p < 0.05.\n'
                'Synthetic data · Reproducible Python/C++ analysis',
                fontsize=10.5, color='#47545F', linespacing=1.55)
    figure.savefig(output / 'selection_bias.png', dpi=200, facecolor='white')
    plt.close(figure)
    print(output / 'summary.md')
    print(output / 'type1_error_en.png')
    print(output / 'selection_bias.png')


if __name__ == '__main__':
    main()
