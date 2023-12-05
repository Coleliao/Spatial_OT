# 2023/03/22 analyse the cellbin data of pigeon D6 by Seurat
library(Seurat)
data=readRDS('../data.rds')
data <- data %>%
  subset(subset = nFeature_Spatial > 150 & nFeature_Spatial < 2500) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures=3500) %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>% FindClusters()
png('Dimplot_seurat_clustering_0322.png')
DimPlot(data,label=T)
dev.off()
genes=c('GAD1','GAD2','SLC17A7','SLC17A6','SLC1A2','FLT1','APBB1IP','PTPRC','MOG','PLP1','PDGFRA','VIP','RELN','PDGFRB','SST','PVALB','LAMP5')
pdf('Featureplot_markers_1_0322.pdf',18,12)
FeaturePlot(data,features=genes[1:11],ncol=3,pt.size=0.05)
dev.off()
pdf('Featureplot_markers_2_0322.pdf',14,10)
FeaturePlot(data,features=genes[12:17],ncol=3,pt.size=0.05)
dev.off()
saveRDS(data,'D6_data_clustering_0322.RDS')

# subset OT
ot=subset(data,subset = y>=9000)
ot$bin100.y=ot$x
ot$bin100.x=ot$y
DoDimPlot(ot,filename='ot_dimplot_0322.pdf')
ot <- ot %>%
  subset(subset = nFeature_Spatial > 200 & nFeature_Spatial < 2500) %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures=3500) %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>%
  RunUMAP(dims = 1:50) %>% FindClusters(resolution=1.5)
saveRDS(ot,'ot_d6_clutering_0322.rds')

pdf('ot_dimplot_clustering_0323.pdf')
DimPlot(ot,label=T)
dev.off()
DefaultAssay(ot)='Spatial'
# DCN--VLMC GJA1/AQP4--AST NRN1/SLC17A6--EXN GAD2--INN PLP1/MBP--OLI OLIG2--OPC
genes=c('SNAP25','SLC17A6','NRN1','GAD2','PLP1','MBP','GJA1','AQP4','GFAP','OLIG2','DCN')
pdf('ot_featureplot_markers_final_0327.pdf',10,10)
for(i in genes){
    print(FeaturePlot(ot,features=i,pt.size=0.01,order=T))
}
dev.off()
pdf('ot_featureplot_markers_2_0322.pdf',14,10)
FeaturePlot(ot,features=genes[12:17],ncol=3,pt.size=0.05)
dev.off()
DoDimPlot(ot,filename='ot_dimplot_spatial_0324.pdf',pt.size=0.25)

# Annotation
ano <- data.frame(seurat_clusters=c(0:22,25),# 23/24 with too less cells, omit
                      celltype=c('INN','EXN','INN','AST','INN','OLIG','OLIG','OLIG',
                                 'INN','OLIG','EXN','OLIG','OLIG','OLIG','EXN','OLIG',
                                 'VLMC','EXN','AST','OPC','INN','OLIG','INN','EXN'))
Idents(ot)='seurat_clusters'
ot=subset(ot,idents=c(23,24),invert=T)
ot$celltype1 <- ano[match(ot$seurat_clusters,ano$seurat_clusters),'celltype']
Idents(ot)='celltype1'
pdf('ot_dimplot_clustering_anno_0328.pdf',10,8)
library(RColorBrewer)
col <- c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"))
DimPlot(ot,label=T,cols=col)
dev.off()
saveRDS(ot,'./ot_d6_clutering_anno_0328.rds')
ot$celltype1=factor(ot$celltype1,levels=c('AST','OLIG','VLMC','EXN','INN','OPC'))
DoDimPlot(ot,filename='ot_dimplot_anno_spatial_0328.pdf',pt.size=0.25,group='celltype1')

# fig3D
cellbin=readRDS('./D6_OT_cellbin_with_domain_info_0402.rds')
meta=cellbin@meta.data
meta$celltype1=as.vector(meta$celltype1)
meta[which(meta$seurat_clusters==14),'celltype1']='CHN'
cellbin$celltype2=unname(meta$celltype1)
cellbin$celltype2=factor(cellbin$celltype2,levels=c('AST','OLIG','VLMC','EXN','INN','OPC','CHN'))
DoDimPlot(cellbin,filename='ot_dimplot_anno_spatial_0615.pdf',pt.size=0.25,group='celltype2',col=col37[10:16])
# fig3E cell proportion
meta=cellbin@meta.data
meta <- subset(meta,meta$domain%in%c('SOp','SGF','SGC','SAC','SGP','SFP','BCS','Imc','Ipc'))
meta <- mutate(meta,domain=factor(domain,levels = rev(c('SOp','SGF','SGC','SAC','SGP','SFP','BCS','Imc','Ipc'))))
meta$celltype1 <- factor(meta$celltype2,levels=c('AST','OLIG','VLMC','EXN','INN','OPC','CHN'))
p <- ggplot(data=meta,aes(x=domain,y=..count../sum(..count..),fill=celltype2))
source('../functions/Visualization.R')
pdf('./D6_OT_domain_cell_proportion_0615.pdf',10,6)
p + geom_bar(stat = 'count',position = 'fill')+ 
  theme(panel.background=element_rect(fill='transparent',
                                      color ="gray"),
        axis.text.x = element_text(angle = 90, hjust = 0.5,
                                   vjust = 0.5,color = "black",size=12),
        axis.text.y = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size=14))+
  labs(title = 'Cell proportion in OT', x = 'Domain',y="Cell proportion")+
  scale_fill_manual(values = col37[10:16]) +
  ylim(0,1)+
  coord_flip()
dev.off()
saveRDS(cellbin,'./D6_OT_cellbin_with_domain_info_0615.rds')

