"""Draw reproducible synthetic chart examples; values are not research results."""
import argparse
from pathlib import Path
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def bh_adjust(pvalues):
    order = np.argsort(pvalues)
    ranked = pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def save(fig, out, name):
    fig.suptitle('Synthetic teaching examples - not experimental evidence', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for extension in ['png', 'pdf']:
        fig.savefig(out / (name + '.' + extension), dpi=160, bbox_inches='tight')
    plt.close(fig)


def km_curve(time, event):
    xs, ys = [0.0], [1.0]
    survival = 1.0
    for t in np.unique(time):
        at_risk = np.sum(time >= t)
        deaths = np.sum((time == t) & event)
        if deaths:
            survival *= 1 - deaths / at_risk
        xs.append(t)
        ys.append(survival)
    return xs, ys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('results/figures'))
    out = parser.parse_args().out
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    colors = ['#0072B2', '#D55E00']
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    ax = axes.flat[0]
    groups = [rng.normal(4, 0.7, 12), rng.normal(5, 0.8, 12)]
    ax.boxplot(groups, showfliers=False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['Control', 'Treated'])
    for i, values in enumerate(groups):
        ax.scatter(rng.normal(i + 1, 0.05, len(values)), values, c=colors[i], s=22)
    ax.set(title='A  Box plot + individual samples', ylabel='Synthetic expression')
    n = 300
    lfc = rng.normal(0, 1.2, n)
    # Artificial p-values illustrate axes and BH correction, not an expression test.
    p = rng.uniform(size=n)
    p[np.abs(lfc) > 1.7] *= 0.0001
    padj = bh_adjust(p)
    significant = (padj < 0.05) & (np.abs(lfc) >= 1)
    point_colors = np.where(significant, colors[1], '#999999')
    ax = axes.flat[1]
    ax.scatter(lfc, -np.log10(np.maximum(padj, 1e-300)), c=point_colors, s=12)
    for x in [-1, 1]:
        ax.axvline(x, ls='--', color='#777777', lw=0.8)
    ax.axhline(-np.log10(0.05), ls='--', color='#777777', lw=0.8)
    ax.set(title='B  Volcano (artificial statistics)', xlabel='log2 fold change', ylabel='-log10 adjusted p')
    mean = np.exp(rng.uniform(0, 9, n))
    ax = axes.flat[2]
    ax.scatter(mean, lfc, c=point_colors, s=12)
    ax.axhline(0, color='#777777', lw=0.8)
    ax.set(xscale='log', title='C  MA-style plot', xlabel='Mean normalized count', ylabel='log2 fold change')
    expression = rng.normal(5, 1, (30, 8))
    expression[:10, 4:] += 2
    z = (expression - expression.mean(axis=1, keepdims=True)) / expression.std(axis=1, keepdims=True)
    ax = axes.flat[3]
    im = ax.imshow(z, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax.set(title='D  Row Z-score heatmap', xlabel='Sample', ylabel='Gene')
    fig.colorbar(im, ax=ax, label='Within-gene Z-score')
    centered = expression.T - expression.T.mean(axis=0)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    score = u * singular
    explained = singular ** 2 / np.sum(singular ** 2)
    ax = axes.flat[4]
    for group in range(2):
        idx = slice(group * 4, group * 4 + 4)
        ax.scatter(score[idx, 0], score[idx, 1], color=colors[group], label=['Control', 'Treated'][group])
    ax.set(title='E  PCA', xlabel='PC1 ({:.1%})'.format(explained[0]), ylabel='PC2 ({:.1%})'.format(explained[1]))
    ax.legend(frameon=False)
    ax = axes.flat[5]
    counts = np.array([4, 6, 8, 12])
    ratios = counts / 50
    dots = ax.scatter(ratios, np.arange(4), s=counts * 18, c=[0.04, 0.025, 0.01, 0.002], cmap='viridis_r')
    ax.set_yticks(np.arange(4))
    ax.set_yticklabels(['Set A', 'Set B', 'Set C', 'Set D'])
    ax.set(title='F  Enrichment-style dot plot', xlabel='GeneRatio; size = gene count')
    fig.colorbar(dots, ax=ax, label='Adjusted p (artificial)')
    save(fig, out, 'chart_core')
    np.savetxt(out / 'artificial_de.csv', np.column_stack([mean, lfc, p, padj]), delimiter=',',
               header='baseMean,log2FoldChange,pvalue,padj', comments='')

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    chromosome = np.repeat(np.arange(1, 5), 100)
    position = np.arange(400)
    gwas_p = rng.uniform(size=400)
    gwas_p[170:175] = np.linspace(1e-10, 2e-8, 5)
    ax = axes.flat[0]
    ax.scatter(position, -np.log10(gwas_p), c=np.where(chromosome % 2, colors[0], '#777777'), s=9)
    ax.axhline(-np.log10(5e-8), color=colors[1], ls='--')
    ax.set_xticks([49.5, 149.5, 249.5, 349.5])
    ax.set_xticklabels(['1', '2', '3', '4'])
    ax.set(title='A  Manhattan (artificial p-values)', xlabel='Chromosome', ylabel='-log10 p')
    ax = axes.flat[1]
    expected = -np.log10((np.arange(1, 401) - 0.5) / 400)
    observed = -np.log10(np.sort(gwas_p))
    ax.scatter(expected, observed, s=10)
    ax.plot([0, max(expected)], [0, max(expected)], '--', color='#777777')
    ax.set(title='B  QQ plot', xlabel='Expected -log10 p', ylabel='Observed -log10 p')
    ax = axes.flat[2]
    distance = np.linspace(-2000, 2000, 120)
    for i, amplitude in enumerate([2, 4]):
        signal = amplitude * np.exp(-(distance / 450) ** 2) + 0.2
        ax.plot(distance, signal, color=colors[i], label=['Control', 'Treated'][i])
    ax.set(title='C  Mean signal around TSS', xlabel='Distance to TSS (bp)', ylabel='Normalized signal')
    ax.legend(frameon=False)
    ax = axes.flat[3]
    marker_mean = np.array([[3, 0.3, 0.1], [0.2, 4, 0.1], [0.1, 0.3, 3]])
    fraction = marker_mean / 5
    xx, yy = np.meshgrid(np.arange(3), np.arange(3))
    dots = ax.scatter(xx.ravel(), yy.ravel(), s=fraction.ravel() * 450, c=marker_mean.ravel(), cmap='viridis')
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Marker A', 'Marker B', 'Marker C'])
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['Type A', 'Type B', 'Type C'])
    ax.set(title='D  Marker dot plot', xlabel='Dot area = expressing fraction')
    fig.colorbar(dots, ax=ax, label='Mean expression')
    ax = axes.flat[4]
    x, y = np.meshgrid(np.arange(12), np.arange(12))
    values = np.exp(-((x - 4) ** 2 + (y - 7) ** 2) / 12)
    dots = ax.scatter(x, y, c=values, s=45, cmap='viridis')
    ax.invert_yaxis()
    ax.set(title='E  Spatial expression', xlabel='x (grid units)', ylabel='y (grid units)', aspect='equal')
    fig.colorbar(dots, ax=ax, label='Synthetic intensity')
    ax = axes.flat[5]
    nodes = {'TF': (0, 0), 'Gene A': (-1, 1), 'Gene B': (1, 1), 'Gene C': (1, -1)}
    for name, point in nodes.items():
        ax.scatter(*point, s=700, c=colors[0] if name == 'TF' else '#CCCCCC')
        ax.text(*point, name, ha='center', va='center', fontsize=8)
        if name != 'TF':
            ax.annotate('', xy=point, xytext=(0, 0), arrowprops={'arrowstyle': '->', 'shrinkA': 17, 'shrinkB': 17})
    ax.set(title='F  Candidate regulatory network', xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
    ax.axis('off')
    save(fig, out, 'chart_omics')

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    ax = axes.flat[0]
    for group in range(2):
        death = rng.exponential(18 + group * 8, 40)
        censor = rng.uniform(12, 36, 40)
        duration = np.minimum(death, censor)
        event = death <= censor
        xs, ys = km_curve(duration, event)
        ax.step(xs, ys, where='post', color=colors[group], label='Group {}'.format(group + 1))
        for t in duration[~event]:
            survival = ys[np.searchsorted(xs, t, side='right') - 1]
            ax.plot(t, survival, '+', color=colors[group], ms=5)
    ax.set(title='A  Kaplan-Meier (+ = censored)', xlabel='Months', ylabel='Survival estimate', ylim=(0, 1.05))
    ax.legend(frameon=False)
    ax = axes.flat[1]
    hr, lower, upper = np.array([0.7, 1.2, 1.8]), np.array([0.4, 0.8, 1.1]), np.array([1.1, 1.8, 2.9])
    ax.errorbar(hr, range(3), xerr=np.vstack([hr - lower, upper - hr]), fmt='o', capsize=3)
    ax.axvline(1, ls='--', color='#777777')
    ax.set_yticks(range(3))
    ax.set_yticklabels(['Factor A', 'Factor B', 'Factor C'])
    ax.set(xscale='log', title='B  Forest plot (artificial estimates)', xlabel='Hazard ratio and 95% CI')
    truth = rng.integers(0, 2, 300)
    scores = 1 / (1 + np.exp(-(rng.normal(0, 1, 300) + truth * 1.5 - 0.75)))
    order = np.argsort(scores)[::-1]
    tpr = np.r_[0, np.cumsum(truth[order]) / np.sum(truth)]
    fpr = np.r_[0, np.cumsum(1 - truth[order]) / np.sum(1 - truth)]
    ax = axes.flat[2]
    auc = np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2)
    ax.plot(fpr, tpr, label='AUC={:.2f}'.format(auc))
    ax.plot([0, 1], [0, 1], '--', color='#777777')
    ax.set(title='C  Binary ROC', xlabel='False positive rate', ylabel='True positive rate')
    ax.legend(frameon=False)
    ax = axes.flat[3]
    for i in range(5):
        mask = (scores >= i / 5) & (scores < (i + 1) / 5)
        if mask.any():
            ax.scatter(scores[mask].mean(), truth[mask].mean(), s=mask.sum(), color=colors[0])
    ax.plot([0, 1], [0, 1], '--', color='#777777')
    ax.set(title='D  Calibration; area = bin size', xlabel='Mean predicted probability', ylabel='Observed event fraction', xlim=(0, 1), ylim=(0, 1))
    ax = axes.flat[4]
    sets = [set(range(0, 7)), set(range(3, 10)), set(range(5, 12))]
    patterns = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    universe = set.union(*sets)
    counts = [sum(tuple(int(item in group) for group in sets) == pattern for item in universe) for pattern in patterns]
    ax.bar(range(7), counts, color=colors[0])
    ax.set_xticks(range(7))
    ax.set_xticklabels(['A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC'])
    ax.set(title='E  Exclusive intersections (UpSet idea)', xlabel='Exact set membership', ylabel='Number of elements')
    ax = axes.flat[5]
    x, y = np.meshgrid(np.linspace(-1, 1, 9), np.linspace(-1, 1, 9))
    ax.quiver(x, y, -0.3 * y, 0.3 * x, angles='xy', scale_units='xy', scale=1)
    ax.set(title='F  Vector field schematic', xlabel='Embedding 1', ylabel='Embedding 2', aspect='equal')
    save(fig, out, 'chart_clinical')
    (out / 'README.json').write_text(json.dumps({
        'synthetic': True, 'seed': 42, 'numpy': np.__version__, 'matplotlib': matplotlib.__version__,
        'note': 'Illustrations only. Artificial p-values and estimates are not analysis of experimental data.'
    }, indent=2), encoding='utf-8')
    print('Saved three chart sheets as PNG/PDF in {}'.format(out.resolve()))


if __name__ == '__main__':
    main()
