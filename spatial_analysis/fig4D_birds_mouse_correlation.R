library(Seurat)
# mouse data
allen=readRDS('./allen_mouse_ctx_hip_singlecell/vis_ssp_mop/allen_VIS_SSp_MOp_downsample_30k_0519.rds')
head(allen@meta.data)
type=table(allen$label)
type=type[type>=100] # less 100 omit
type=names(type)
Idents(allen)='label'
allen=subset(allen,idents=type)
Idents(allen)='class_label'
glu=subset(allen,idents='Glutamatergic')
saveRDS(glu,'./allen_vis_ssp_mop_glu_filter_0522.rds')
gaba=subset(allen,idents='GABAergic')
saveRDS(gaba,'./allen_vis_ssp_mop_gaba_filter_0522.rds')

# mine
mine=readRDS('./data_subclutering_EXN_2023-05-17.rds')

glu$celltype1=paste0(glu$region_label,glu$label)
m=rownames(glu)
rownames(glu@assays$RNA@data)=toupper(m)
rownames(glu@assays$RNA@counts)=toupper(m)

source('../functions/tools.R')
mtr=corSeurat(mine,glu,idents1='seurat_clusters',idents2='region_label')


# visualization
pdf('cor_hvg_tf_heatmap.pdf')
corheatmap(mtr$r,mtr$p)
corheatmap(mtr2$r,mtr2$p)
dev.off()



