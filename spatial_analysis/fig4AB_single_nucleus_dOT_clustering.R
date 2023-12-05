library(Seurat)

samlist=c('./02.zebrafinch/ZF2_2206_SGP/','./02.zebrafinch/ZF7_2206_SGP/')
datalist <- list()
for(i in 1:length(samlist)){
  file=dir(samlist[i])
  for(j in 1:length(file)){
    path <- paste0(samlist[i],file[j])
    data <- Read10X(data.dir = path,gene.column = 1)
    data <- CreateSeuratObject(data,project = 'turtle',min.cells = 3, min.features = 200)
    data$library <- file[j] 
    if(j==1){object=data}else{object=merge(object,data)}
  }
  info <- strsplit(samlist[i],split = '/')
  object$sample_areaID <- info[[1]][10]
  tmp=strsplit(info[[1]][10],split='_')
  object$sample=tmp[[1]][[1]] 
  object$areID=tmp[[1]][[3]] 
  object <- subset(object,subset = nFeature_RNA > 300 & nFeature_RNA < 3000 ) # QC
  if(ncol(object)>=200){ # cell number is too less to be integrated
    object <- SCTransform(object,verbose = FALSE,assay='RNA') # SCT
    object=list(object)
    names(object)=info[[1]][10]
    datalist=c(datalist,object)
  }
}
features <- SelectIntegrationFeatures(object.list = datalist)
datalist <- PrepSCTIntegration(object.list = datalist, anchor.features = features) # 这步是SCT做整合需要的 就和SCT做差异分析一样
datalist <- FindIntegrationAnchors(object.list = datalist, normalization.method = "SCT",anchor.features = features)
sam.combined <- IntegrateData(anchorset = datalist, normalization.method = "SCT")
sam.combined <- RunPCA(sam.combined, verbose = FALSE)
sam.combined <- RunUMAP(sam.combined, reduction = "pca", dims = 1:30)
sam.combined <- FindNeighbors(sam.combined, reduction = "pca", dims = 1:30)
sam.combined <- FindClusters(sam.combined, resolution = 1.0)
saveRDS(sam.combined,"./pigeon_sgp_integrated_data_raw_0322.rds")
print(sam.combined)


# visualization
library(ggplot2)
library(RColorBrewer)
pdf('dimplot_inte_sgp_sc_0322.pdf',10,7)
color=c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"))
DimPlot(sam.combined,label = T,cols=color,repel=T)
DimPlot(sam.combined,label = F,cols=color)+NoLegend()
DimPlot(sam.combined,group.by = 'sample',label=F)
DimPlot(sam.combined,group.by = 'areID',label = F)
DimPlot(sam.combined,group.by = 'library',label = F)
dev.off()
DefaultAssay(sam.combined)='RNA'
# ! do NormalizedData on RNA assay first, or else, Featureplot will use counts slot instead.
sam.combined=NormalizeData(sam.combined)
sam.combined=FindVariableFeatures(sam.combined)
genes=c('GAD1','GAD2','SLC17A7','SLC17A6','SLC1A2','FLT1','APBB1IP','PTPRC','MOG','PLP1','PDGFRA','C1QB','VIP','RELN','PDGFRB','SST','PVALB','LAMP5')
pdf('Featureplot_markers_1_0314.pdf',18,12)
FeaturePlot(sam.combined,features=genes[1:11],ncol=3,pt.size=0.05)
dev.off()
pdf('Featureplot_markers_2_0314.pdf',14,10)
FeaturePlot(sam.combined,features=genes[12:17],ncol=3,pt.size=0.05)
dev.off()

# 0328 main cell type annotation 
genes=c('HBAD','HBAA','C1QC','C1QB',"AGT","CX3CR1","DRAXIN","GFAP",'ZCCHC24')
sam.combined=readRDS('./pigeon_sgp_integrated_data_raw_0322.rds')
DefaultAssay(sam.combined)='RNA'
sam.combined=NormalizeData(sam.combined)
sam.combined=FindVariableFeatures(sam.combined)
pdf('Featureplot_markers_3_0328.pdf',18,12)
FeaturePlot(sam.combined,features=genes,ncol=3,pt.size=0.05)
dev.off()
ano <- data.frame(seurat_clusters=c(0:29),# 23/25/29 omit
                      celltype=c('EXN','EXN','EXN','INN','AST','INN','INN','INN',
                                 'OLIG','OLIG','INN','EXN','INN','INN','INN','INN',
                                 'EXN','EXN','INN','EXN','OPC','INN','Blood','INN',
                                 'EXN','EXN','AST','MIC','OLIG','INN')) 
sam.combined$celltype1 <- ano[match(sam.combined$seurat_clusters,ano$seurat_clusters),'celltype']
pdf('dimplot_inte_sgp_sc_anno_0328.pdf',10,7)
DimPlot(sam.combined,label = T,cols=color,repel=T)
DimPlot(sam.combined,label = T,cols=color,group.by='celltype1',repel=T)
dev.off()
saveRDS(sam.combined,'./pigeon_sgp_integrated_data_anno_0328.rds')

# 0528 re-plot and markers dotplot
data=readRDS('./pigeon_sgp_integrated_data_anno_0328.rds')
genes=c('SLC17A6','GAD2','SLC1A2','PLP1','PDGFRA','C1QB','HBAD')
data$celltype1=factor(data$celltype1,levels=c('EXN','INN','AST','OLIG','OPC','MIC','Blood'))
Idents(data)='celltype1'
pdf('dotplot_main_celltype_markers_0528.pdf',8,4)
DotPlot(data,features=genes)+theme_bw()+scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))
dev.off()
# main celltype highlight
cells=CellsByIdentities(data)
pdf('cells_highlight_0528.pdf')
DimPlot(data, cells.highlight=cells$EXN,cols.highlight = c("#DE2D26", "grey90"),pt.size=0.05,sizes.highlight=0.05)
DimPlot(data, cells.highlight=cells$INN,cols.highlight = c("#DE2D26", "grey90"),pt.size=0.05,sizes.highlight=0.05)
dev.off()



# 0822 save data
data=readRDS('pigeon_sgp_integrated_data_anno_0328.rds')
count=as.matrix(data@assays$RNA@counts)

meta=data@meta.data
meta$`orig.ident`='zebrafinch'
meta$areID='deep_OT'
meta$sample_name='zebrafinch'
write.csv(meta,'zf_deep_OT_sc_qc_metadata.CSV')

pos=data@reductions$umap@cell.embeddings
pos=as.data.frame(pos)
pos$cluster=meta$seurat_clusters
pos$celltype=meta$celltype1
write.csv(pos,'zf_deep_OT_sc_qc_umap.CSV')



