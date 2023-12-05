## correlation between two species
source('../functions/process.R')
source('../functions/tools.R')
source('../functions/Visualization.R')
library(Seurat)
library(ComplexHeatmap)
library(psych)
Idents(data)='species'
pg=subset(data,idents='pigeon')
zf=subset(data,idents='zebrafinch')
deg=read.csv('./deg_OT_domain2_0729.CSV',row.names=1)
deg=subset(deg,deg$p_val_adj<0.05&deg$avg_log2FC>0.5)
genes=unique(deg$gene)
zf$domain2=factor(zf$domain2,levels=c('SOp','SGF','SGC','SAC','SGP','SFP','BCS','Ipc','Imc'))
pg$domain2=factor(pg$domain2,levels=c('SOp','SGF','SGC','SAC','SGP','SFP','BCS','Ipc','Imc'))
mtr=corSeurat(zf,pg,idents1='domain2',idents2='domain2',plot=T,filename='heatmap_zf_pg_cor_deg_0729.pdf',genes=genes,cluster_row=F,cluster_col=F)
top50 <- deg %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 50)
top50=top50$gene
mtr=corSeurat(zf,pg,idents1='domain2',idents2='domain2',plot=T,filename='heatmap_zf_pg_cor_deg_top50_0729.pdf',genes=top50,cluster_row=F,cluster_col=F)


