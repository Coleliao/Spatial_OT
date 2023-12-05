#/*************************************************************************
# File Name: h5ad_to_h5.py
# Author:liaokuo
# Created Time: Fri 01 Jul 2022 04:21:23 PM CST
# ************************************************************************/
# turn h5ad to Seurat object in RDS format
# script.py -I input_data_path -O outpath
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import h5py
import argparse
import os

def parse_args():
    description='try to turn h5ad to Seurat object in RDS format!'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-I',"--input",help='input the path of your h5ad')
    parser.add_argument('-O','--output',help='output path',default=".")
    return parser.parse_args()

args=parse_args()
inpath=args.input
out=args.output
os.chdir(out)

data=anndata.read(inpath)
is_numpy_array = isinstance(data.layers['counts'], np.ndarray) 
if is_numpy_array: 
    mat=pd.DataFrame(data=data.layers['counts'],index=data.obs_names,columns=data.var_names)
else:
    mat=pd.DataFrame(data=data.layers['counts'].todense(),index=data.obs_names,columns=data.var_names)
mat.to_hdf("mat.h5","mat")
data.obs['x']=data.obsm['spatial'][:,0]
data.obs['y']=data.obsm['spatial'][:,1]
#data.obs['bin50_x']=data.obsm['spatial_50'][:,0]
#data.obs['bin50_y']=data.obsm['spatial_50'][:,1]
meta=pd.DataFrame(data=data.obs)
meta.to_csv('metadata.tsv',sep="\t")
