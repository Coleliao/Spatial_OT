import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
import anndata
import argparse
import spateo as st
def parse_args():
    description='try to run spateo with your gem matrix!'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-I',"--input",help='input the path of your gem matrix')
    parser.add_argument('-B',"--bins",help='bins you want',default=50)
    parser.add_argument('-O','--output',help='output path',default=".")
    return parser.parse_args()

args=parse_args()
mtx_path=args.input
bs=int(args.bins)
out=args.output

os.chdir(out)
plt.figure(dpi=600,figsize=(24,18))

#adata=anndata.read('./cell_adata.h5ad')
adata = st.io.read_bgi(
    mtx_path,
    binsize=bs
)
print(adata)
'''
sc.pp.filter_cells(adata, min_genes=200) # filter low-quality bins
print(adata)
'''
if not os.path.exists("figures"):
    os.mkdir("figures")
os.chdir("./figures")

#Set coordinates
x_array=adata.obsm['spatial'][:,0].tolist()
y_array=adata.obsm['spatial'][:,1].tolist()
x_pixel=adata.obsm['spatial'][:,0].tolist()
y_pixel=adata.obsm['spatial'][:,1].tolist()
adata.obs["x_array"]=x_array
adata.obs["y_array"]=y_array
adata.obs["x_pixel"]=x_array
adata.obs["y_pixel"]=y_array

#Calculate adjacent matrix
s=1
b=49
adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
#preprocessing
adata.var_names_make_unique()
#sc.pp.filter_cells(adata, min_genes=200) # filter low-quality bins
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
print(adata)
#Normalize and take log for UMI
adata.layers['counts']=adata.X # copy the raw count matrix
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
#Set hyper-parameters
p=0.5
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)  # Find the l value given p
n_clusters=28 # prospective cluster number you want
r_seed=t_seed=n_seed=100 #Set seed
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
# Run SpaGCN
clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
#Save results
meta=adata.obs
meta.to_csv("metadata.xls",sep="\t",float_format = '%.0f')
#adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
try:
    adata.write('../adata_sp_clustering.h5ad')
    print(adata)
except:
    print(adata)
    print('adata is not able to be stored.\n')

#visualization
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785","#FFFF99", "#B15928","#7FC97F", "#BEAED4", "#FDC086","#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462","#D9D9D9"]
#Plot spatial domains
domains="pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig("spGCN_pred.png", dpi=600)
plt.close()
#Plot refined spatial domains
domains="refined_pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig("spGCN_refined_pred.png", dpi=600)
plt.close()

#Plot split refine spatial domains
adata.obsm['umap']=adata.obsm['spatial']
celltype=adata.obs['refined_pred'].unique().tolist()
celltype=map(int,celltype) #tranfer into numeric
celltype=list(celltype)
celltype.sort()
celltype= [str(x) for x in celltype] #transfer back into str

adata.obs.refined_pred=adata.obs.refined_pred.astype('str')
print(celltype)
if len(celltype)>=36:
    celltype=celltype[0:36] # just abandon the rest of it, 'cause they're not important
print(celltype)
if len(celltype)>=25:
    n=6
else:
    n=5
fig, axes = plt.subplots(n, n, figsize=(24, 18), tight_layout=True,dpi=600)
x=0;y=0
for i in celltype:
    fig=sc.pl.umap(adata, color='refined_pred', groups=[i], ax=axes[y][x],show=False)
    x = x +1 if x <(n-1) else 0
    y = y +1 if y <(n-1) and x ==0 else y
plt.show()
plt.savefig('./spGCN_split.png')

