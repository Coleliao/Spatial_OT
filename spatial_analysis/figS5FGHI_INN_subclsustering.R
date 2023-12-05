#sub INN
library(Seurat)
data=readRDS("./data_pigeon_mt_2023-05-17.rds")
DefaultAssay(data)=ifelse('RNA'%in%Assays(data),'RNA','Spatial')
sub=c('INN')
idents='celltype1'
split='sample'
setwd("./INN")
source('./functions/process.R')
source('./functions/Visualization.R')
sub=Subclutering(data,subset=sub,idents=idents,split=split,method='sct',dims=1:10,resolution=1)
pdf(paste0('Dimplot_subclutering_INN_r1_',Sys.Date(),'.pdf'))
DimPlot(sub,cols=col37,label=T)
DimPlot(sub,group.by=split)
dev.off()
saveRDS(sub,paste0('data_subclutering_INN_r1_',Sys.Date(),'.rds'))

#featureplot
pdf(paste0('Featureplot_subclutering_INN_r1_',Sys.Date(),'.pdf'))
FeaturePlot(sub,features = c('GAD1','GAD2'))
FeaturePlot(sub,features = c('SNCG','NRN1'))
dev.off()

#deg
Idents(sub)="seurat_clusters"
sub=NormalizeData(sub)
deg=FindMarkers(sub,ident.1=c(4,13,14),min.pct = 0.25)
write.csv(deg,paste0('deg_FindMarkers_INN_41314_r1_',Sys.Date(),'.CSV'))
deg2=deg %>% subset(deg$p_val_adj<0.05)
write.csv(deg,paste0('deg_FindMarkers_INN_41314_fitrate_r1_',Sys.Date(),'.CSV'))

#HAc
DefaultAssay(sub)="RNA"
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)
sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
gene=sub@assays$RNA@var.features
Idents(sub)="seurat_clusters"
mtr <- AverageExpression(sub,features=gene,assays='RNA')
mtr <- mtr[[1]] %>% t()
mtr <- scale(mtr)
d<-dist(mtr) 
fit1<-hclust(d,method = "average")
pdf('HAC_INN_r1_2000.pdf',10,10)
plot(fit1,hang =-1,cex=.8,main = "HAC of EXN in pigeon")
dev.off()


#dotplot
sub$seurat_clusters=factor(sub$seurat_clusters,levels=c(16,15,10,8,9,13,4,14,11,12,6,7,5,3,2,0,1))
Idents(sub)="seurat_clusters" 
gene <- c("FOXP1","FOXP2","ZFHX4","MEIS2","PBX3","TSHZ1","PVALB","ETV1","TRPS1","SST","SOX6","LHX6","MAF","ARX","SATB1","NPY","EFNA5","ELFN1","CALB1","LAMP5","NR2F2","ADARB2","CNR1","ZBTB16","PROX1","RELN","PENK","LHX8")
pdf(paste0('Dotplot_doublet_INN_gene_',Sys.Date(),'.pdf'))
DotPlot(sub,features=gene)+coord_flip()+theme_bw()+theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+ scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66'))
dev.off()


