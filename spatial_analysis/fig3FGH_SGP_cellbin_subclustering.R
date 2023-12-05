library(Seurat)
library(dplyr)
## 0806 subcluter SGP+SFP
Idents(data)='domain'
sub=subset(data,idents=c('SGP','SFP'))
source('../functions/process.R')
meta=sub@meta.data
write.csv(meta,'./cellbin_sgp_sfp_metadata_0808.CSV')
# DBSCAN filtered outliers
meta=read.csv('./cellbin_sgp_sfp_metadata_adj_0808.CSV',row.names=1)
sub=sub[,rownames(meta)]
sub=quickSeurat(sub,method='sct')
pdf('dimplot_SGP_SFP_subcluster_0808.pdf') 
DimPlot(sub)
DimPlot(sub,group.by='celltype2')
dev.off()
pdf('featureplot_SGP_SFP_subcluster_0806.pdf')
genes=c('SLC17A6','GAD2','PLP1','PDGFRA','FLT1','RGS5','GJA1','NRN1','PDGFRA','DCN','OLIG2','AQP4')
genes=intersect()
for(i in genes){
        print(FeaturePlot(sub,features=i,order=T))
}
dev.off()
DefaultAssay(sub)='Spatial'
sub=NormalizeData(sub)
deg=FindAllMarkers(sub,only.pos=T)
top5 <- deg %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 5)
##
saveRDS(sub,'cellbin_sgp_sfp_0807.rds')
write.csv(deg,'cellbin_sgp_sfp_deg_0807.CSV')
pdf('featureplot_SGP_SFP_subcluster_2_0806.pdf')
FeaturePlot(sub,features=c('PENK','CCK','NPTX2','NTS','NMB','CRHBP'))
dev.off()
DoDimPlot(sub,'spatial_dimplot_cellbin_SGP_SFP_0807.pdf')
cells=CellsByIdentities(sub)
DoSpotHighlight(sub,cell=cells[[7]],filename='spatial_dimplot_cellbin_SGP_SFP_clu6_0807.pdf',pt.size=2,sizes.highlight=2)
DoSpotHighlight(sub,cell=cells[[8]],filename='spatial_dimplot_cellbin_SGP_SFP_clu7_0807.pdf',pt.size=2,sizes.highlight=3)
mtr=scheatmap(sub,features=top5$gene,plot=T,filename='heatmap_SGP_SFP_deg_0807.pdf')
# deg filter
deg=deg[-grep('A306|transcript',deg$gene),]
top5 <- deg %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 5)
mtr=scheatmap(sub,features=top5$gene,plot=T,filename='heatmap_SGP_SFP_deg_2_0807.pdf')
deg=subset(deg,deg$p_val_adj<0.05)
write.csv(deg,'cellbin_sgp_sfp_deg_filter_0807.CSV')
meta=sub@meta.data
write.csv(meta,'./cellbin_sgp_sfp_metadata_0808.CSV')

