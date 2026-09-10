"""Public PBMC3k Scanpy teaching workflow; first run downloads the dataset."""
from pathlib import Path
import importlib.metadata
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc


def main():
    out = Path('results/scrna')
    out.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 2
    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()
    adata.obs['sample'] = 'PBMC3k'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata.obs.to_csv(out / 'qc_before_filter.csv')
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.2, multi_panel=True, show=False)
    plt.savefig(out / 'qc.png', dpi=180, bbox_inches='tight')
    plt.close('all')
    # Historical PBMC3k teaching thresholds, not universal QC thresholds.
    adata = adata[(adata.obs.n_genes_by_counts >= 200) &
                  (adata.obs.n_genes_by_counts < 2500) &
                  (adata.obs.pct_counts_mt < 5), :].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()  # Log-normalized expression, NOT raw counts.
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000)
    representation = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(representation, max_value=10)
    sc.tl.pca(representation, n_comps=40, svd_solver='arpack', random_state=42)
    adata.obsm['X_pca'] = representation.obsm['X_pca'].copy()
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30, use_rep='X_pca', random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, resolution=0.5, random_state=42, flavor='leidenalg',
                 n_iterations=2, directed=False, key_added='cluster')
    sc.tl.rank_genes_groups(adata, 'cluster', method='wilcoxon', use_raw=True, pts=True)
    sc.get.rank_genes_groups_df(adata, group=None).to_csv(out / 'cluster_markers.csv', index=False)
    sc.pl.umap(adata, color=['cluster', 'sample'], show=False)
    plt.savefig(out / 'umap.png', dpi=180, bbox_inches='tight')
    plt.close('all')
    markers = ['CD3D', 'MS4A1', 'LYZ', 'NKG7', 'PPBP']
    markers = [gene for gene in markers if gene in adata.raw.var_names]
    sc.pl.dotplot(adata, markers, groupby='cluster', use_raw=True, show=False)
    plt.savefig(out / 'markers.png', dpi=180, bbox_inches='tight')
    plt.close('all')
    adata.write_h5ad(out / 'pbmc3k_processed.h5ad')
    versions = {name: importlib.metadata.version(name)
                for name in ['scanpy', 'anndata', 'numpy', 'scipy', 'igraph', 'leidenalg']}
    (out / 'versions.json').write_text(json.dumps(versions, indent=2), encoding='utf-8')
    print('Retained {} cells and {} genes. Outputs: {}'.format(*adata.shape, out.resolve()))


if __name__ == '__main__':
    main()
