# set section Z1 from zebrafinch as example, all 5 sections followed.

library(Seurat)
library(dplyr)
source('../functions/Visualization.R')

data=readRDS('./Z1_cellbin_data_0720.rds')
data <- data %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures=3500) %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>% FindClusters(resolution=1.0)
library(RColorBrewer)
col <- c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"))
data$bin100.x=data$x
data$bin100.y=data$y
DoDimPlot(data,pt.size=3.15,filename='Dimplot_Z1_zf_spGCN_0720.pdf',group='refined_pred')
sam=data
Idents(sam)='refined_pred'
sam@reductions$umap@cell.embeddings[,"UMAP_1"] <- as.numeric(sam$bin100.x)
sam@reductions$umap@cell.embeddings[,"UMAP_2"] <- as.numeric(sam$bin100.y)
cell=CellsByIdentities(object=sam)
library(ggplot2)
p <- list()
for(i in 1:length(levels(sam))){
  p[[i]] <- DimPlot(sam, cells.highlight=cell[[i]],cols.highlight = c("#DE2D26", "grey90"),pt.size=0.05,sizes.highlight=0.05)
  p[[i]] <- p[[i]] + NoAxes()+NoLegend()+scale_x_reverse()
}
pdf("Dimplot_zf_Z1_spatial_split_spGCN_0720.pdf",16,16)
print(cowplot::plot_grid(plotlist=p,labels = names(cell),label_x = 0.5,label_y=1))
dev.off()

genes=c("ANLN","NTS","ST18","NMB","CALB2","HAPLN2","SCN4B","CABP7","CBLN2","SYT13","CBLN1","CPLX3","ENKUR","OTX2","CRHBP","SLC13A3","APOD","MLPH",'CTGF')
DoFeaturePlot(data,features=genes,pt.size=3)
genes2=c('TCF7L2','LHX9','CRHBP','PENK')
DoFeaturePlot(data,features=genes2,pt.size=3,filename='featureplot2.pdf')

# annotation, omit domain1 here
anno <- data.frame(clusters=c(9,7,0,10,3,8,1,4,11,12,14,6,15,13,5,2),
                   #domain1=c('L1','L2','L3','L4','L4','L5','L6','L7','L8','L9','Imc','MB1','MB2','MB3',
                   #          'MB4','MB5','MB6','BCS','CB1','CB2','CB3','CB4','CB5','CB6','CB7','CB8','Other1'),
                   domain2=c('SOp','SOp','SGF','SGF','SGF','SGF','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc','Mld','MB1'))
data$domain2 <- anno[match(data$refined_pred,anno$clusters),'domain2']
DoDimPlot(data,pt.size=3.15,filename='Dimplot_Z1_spGCN_ann_domain2_0720.pdf',group='domain2')


# UMI spatial distribution for figS1B
library(ggplot2)
source('/hwfssz1/ST_SUPERCELLS/P22Z10200N0664/liaokuo/shanghai_monkeybrain/Visualization.R')
Z1=readRDS('./zebrafinch_Z1_anno_count_0727.rds')
col=colorRampPalette(c('#041C54','#5ecc74', "#e59e36",'#DC143C'))(100)
DoFeaturePlot(Z1,features='nCount_Spatial',col='stereo',pt.size = 3,filename='spatial_UMI_distribution_Z1_1130.pdf')


