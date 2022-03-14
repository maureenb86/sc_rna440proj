import numpy as np
import pandas as pd
import scanpy as sc

# set up scanpy to function as intended
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = 'write/440.h5ad'  # the file that will store the analysis results

# read in the scRNA-seq data to adata
adata = sc.read(
    "medullo_immune_cells/matrix.csv",
    cache=True)

# create top twenty genes figure and save
twenty_top_genes = sc.pl.highest_expr_genes(adata,
                                            n_top=20,
                                            save="twenty_top_genes.png")

# filter out bad outlier cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,
                           qc_vars=['mt'],
                           percent_top=None,
                           log1p=False,
                           inplace=True)

# create qc violin plots and save
sc.pl.violin(adata,
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4,
             multi_panel=True,
             save="violin_qc.png")

# plot mitochondrial DNA plots (no save)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save="mit_percent_qc.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
