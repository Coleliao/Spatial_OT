#script.py matirx.gem ssDNA.tif outdir
import warnings
warnings.filterwarnings("ignore")
import spateo as st
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.pyplot as plt
import os
import sys
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
#import dynamo as dyn

fpath=sys.argv[1] # matirx.gem
spath=sys.argv[2] # ssDNA.tif
outdir=sys.argv[3]


st.config.n_threads = 1
# data input
adata = st.io.read_bgi_agg(
    fpath, spath,
)

adata

os.chdir(outdir)
if not os.path.exists("figures"):
    os.mkdir("figures")
os.chdir("./figures")

st.pl.imshow(adata, 'stain')
plt.savefig("fig1_stain.png",dpi=600)
plt.close()

# alignment refine
'''
before = adata.layers['stain'].copy()
st.cs.refine_alignment(adata, mode='rigid', transform_layers=['stain'])
fig, axes = plt.subplots(ncols=2, figsize=(8, 4), tight_layout=True)
axes[0].imshow(before)
st.pl.imshow(adata, 'unspliced', ax=axes[0], alpha=0.6, cmap='Reds', vmax=2, use_scale=False, save_show_or_return='return')
axes[0].set_title('before alignment')
st.pl.imshow(adata, 'stain', ax=axes[1], use_scale=False, save_show_or_return='return')
st.pl.imshow(adata, 'unspliced', ax=axes[1], alpha=0.6, cmap='Reds', vmax=2, use_scale=False, save_show_or_return='return')
axes[1].set_title('after alignment')
plt.savefig("fig2_refine_alignment.png",dpi=600)
plt.close()
'''

# segmentation
st.cs.mask_nuclei_from_stain(adata)
st.pl.imshow(adata, 'stain_mask')
plt.savefig("fig3_stain_mask.png",dpi=600)
plt.close()
# seg1 by watershed
st.cs.find_peaks_from_mask(adata, 'stain', 7)
st.cs.watershed(adata, 'stain', 5, out_layer='watershed_labels')
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'watershed_labels', labels=True, alpha=0.5, ax=ax)
plt.savefig("fig4_waterfall.png",dpi=600)
plt.close()
# seg2 by deep learning with stardist
st.cs.stardist(adata, tilesize=-1, equalize=2.0, out_layer='stardist_labels')
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'stardist_labels', labels=True, alpha=0.5, ax=ax)
plt.savefig("fig5_stardist.png",dpi=600)
plt.close()
# Augment labels
st.cs.augment_labels(adata, 'watershed_labels', 'stardist_labels', out_layer='augmented_labels')
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'augmented_labels', labels=True, alpha=0.5, ax=ax)
plt.savefig("fig6_argmented_labels.png",dpi=600)
plt.close()

# Expand from nuclei to cytoplasm
st.cs.mask_cells_from_stain(adata, out_layer='stain_cell_mask')
st.cs.watershed(
    adata, 'stain',
    mask_layer='stain_cell_mask',
    markers_layer='augmented_labels',
    out_layer='cell_labels',
) # first is to use watershed method with a lenient threshold to expand nuclei to cell level
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'cell_labels', labels=True, alpha=0.5, ax=ax)
plt.savefig("fig7_expand_cell.png",dpi=600)
plt.close()
'''
st.cs.expand_labels(adata, 'augmented_labels', distance=5, max_area=400)
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'augmented_labels_expanded', labels=True, alpha=0.5, ax=ax)
plt.savefig("fig8_augmented_labels_expanded.png",dpi=600)
plt.close()
'''

# Obtain a cell x gene matrix
st.cs.expand_labels(
    adata, 'cell_labels', distance=10, out_layer='cell_labels_expanded',max_area=1200
)
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'cell_labels_expanded', labels=True, alpha=0.5, ax=ax)
plt.savefig("fig9_final_cell_labels_expanded.png",dpi=700)
plt.close()


cell_adata = st.io.read_bgi(
    fpath,   # matrix path
    segmentation_adata=adata,
    labels_layer='cell_labels_expanded',
)
cell_adata    # final result of cell segmentation by spateo

cell_adata.write('../cell_adata.h5ad')

###############part2 start##############
print("part2 is going to start!\n")

data=cell_adata

#filter
plt.figure(dpi=600,figsize=(12,16)) #adjust dpi and figsize
plt.xticks(fontsize=10) #adjust fontsize
sc.pl.highest_expr_genes(data, n_top=20)
plt.savefig("Highest_expr_genes.pdf")
plt.close()
sc.pp.filter_cells(data, min_genes=200)
sc.pp.filter_genes(data, min_cells=3)

sc.pp.calculate_qc_metrics(data, percent_top=None, log1p=False, inplace=True)
sc.pl.violin(data, ['n_genes_by_counts', 'total_counts'],jitter=0.4, multi_panel=True)
plt.savefig("QC_violin.pdf")
data = data[data.obs.n_genes < 2500, :]
print(data.obs['n_genes_by_counts'].median(),'\n')
print(data.obs['total_counts'].median(),"\n")


#pre-processiing and pick up HVGs
sc.pp.normalize_total(data, target_sum=1e4) 
sc.pp.log1p(data)
data.raw = data
sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(data)
plt.savefig("highly_variable_genes.pdf")
data = data[:, data.var.highly_variable]
sc.pp.regress_out(data, ['total_counts'])
sc.pp.scale(data, max_value=10)

#PCA and dimension reduction
sc.tl.pca(data, svd_solver='arpack')
sc.pl.pca_variance_ratio(data, log=True)
plt.savefig("PCA.pdf")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.leiden(data) #leiden for clutering
sc.tl.umap(data)
sc.pl.umap(data, color=['leiden'])
plt.savefig("Umap.pdf")
data.write('results_file.h5ad')
data.obsm['X_umap']=data.obsm['spatial']
sc.pl.umap(data, color=['leiden'])
plt.savefig("Umap2.pdf")

# different expression genes
sc.tl.rank_genes_groups(data, 'leiden', method='t-test')
sc.pl.rank_genes_groups(data, n_genes=25, sharey=False)
plt.savefig('diff_feature.png',dpi=600)



