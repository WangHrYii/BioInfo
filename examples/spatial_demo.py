"""Synthetic spatial autocorrelation exercise, not a real tissue analysis."""
import argparse
from pathlib import Path
import csv
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def morans_i(values, weights):
    centered = values - values.mean()
    denominator = centered @ centered
    if denominator == 0 or weights.sum() == 0:
        raise ValueError('Moran I needs variation and nonzero spatial weights')
    return len(values) / weights.sum() * (centered @ weights @ centered) / denominator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('results/spatial'))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    x, y = np.meshgrid(np.arange(15), np.arange(15))
    coordinates = np.column_stack([x.ravel(), y.ravel()])
    signal = rng.poisson(2 + 15 * np.exp(-((coordinates[:, 0] - 5) ** 2 +
                                          (coordinates[:, 1] - 7) ** 2) / 15))
    shuffled = rng.permutation(signal)
    distance = np.abs(coordinates[:, None, :] - coordinates[None, :, :]).sum(axis=2)
    weights = (distance == 1).astype(float)  # Symmetric, binary four-neighbor graph.
    observed = morans_i(signal, weights)
    null = np.array([morans_i(rng.permutation(signal), weights) for _ in range(999)])
    pvalue = (1 + np.count_nonzero(null >= observed)) / (len(null) + 1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, values, title in zip(axes, [signal, shuffled], ['Spatial pattern', 'Same values, shuffled']):
        points = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=values,
                            s=75, marker='s', cmap='viridis', vmin=0, vmax=signal.max())
        ax.set(title=title, xlabel='x (grid units)', ylabel='y (grid units)', aspect='equal')
        ax.invert_yaxis()
    fig.colorbar(points, ax=axes, label='Synthetic counts')
    fig.suptitle('Synthetic illustration; not experimental evidence')
    for extension in ['png', 'pdf']:
        fig.savefig(args.out / ('spatial_pattern.' + extension), dpi=180, bbox_inches='tight')
    plt.close(fig)
    with (args.out / 'spots.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['x', 'y', 'synthetic_counts', 'shuffled_counts'])
        writer.writerows(np.column_stack([coordinates, signal, shuffled]).tolist())
    result = {'morans_i': float(observed), 'shuffled_i': float(morans_i(shuffled, weights)),
              'one_sided_permutation_p': float(pvalue), 'permutations': len(null),
              'graph': 'symmetric binary four-neighbor grid', 'synthetic': True}
    (args.out / 'statistics.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
