library(Seurat)
data=readRDS("./pigeon_sgp_integrated_data_anno_0328.rds")
DefaultAssay(data)=ifelse('RNA'%in%Assays(data),'RNA','Spatial')
sub=c('EXN')
idents='celltype1'
split='sample'
source('../functions/process.R')
source('../functions/Visualization.R')
sub=Subclutering(data,subset=sub,idents=idents,split=split,method='sct',dims=1:10,resolution=1)
pdf(paste0('Dimplot_subclutering_EXN_r1',Sys.Date(),'.pdf'))
DimPlot(sub,cols=col37,label=T)
DimPlot(sub,group.by=split)
dev.off()
saveRDS(sub,paste0('data_subclutering_EXN_r1',Sys.Date(),'.rds'))


