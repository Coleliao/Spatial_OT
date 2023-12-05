# 2023/07/26
# integrate the spatial data of pigeon and zebrafinch
library(Seurat)
library(dplyr)
library(ggplot2)
source('../functions/shanghai_monkeybrain/process.R')
source('../functions/Visualization.R')


# pigeon
d6=readRDS('./pigeon_d6_anno_count_0727.rds')
d6$slide='D6'
d6$species='pigeon'
Idents(d6)='domain2'
d6=subset(d6,idents=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'))
c6=readRDS('./pigeon_c6_anno_count_0727.rds')
c6$slide='C6'
c6$species='pigeon'
Idents(c6)='domain2'
c6=subset(c6,idents=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'))
e3=readRDS('./pigeon_e3_anno_count_0727.rds')
e3$slide='E3'
e3$species='pigeon'
Idents(e3)='domain2'
e3=subset(e3,idents=c('SOp','SGF','SGC','SAC','BCS','SGP','Imc'))

# zebrafinch
f2=readRDS('./zebrafinch_f2_anno_count_0727.rds')
f2$slide='F2'
f2$species='zebrafinch'
Idents(f2)='domain2'
f2=subset(f2,idents=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'))
a6=readRDS('./zebrafinch_a6_anno_count_0727.rds')
a6$slide='A6'
a6$species='zebrafinch'
Idents(a6)='domain2'
a6=subset(a6,idents=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc'))

# QC
mlist=list(d6,c6,e3,f2,a6)
mlist=lapply(mlist,function(x){
x=subset(x,nFeature_Spatial>200)
return(x)
})
#min(colSums(d6@assays$Spatial@count))

data <- DoIntegration(mlist[[1]],mlist[[2]],mlist[[3]],mlist[[4]],mlist[[5]],nfeatures=2000,method='sct')
saveRDS(data,'spatial_bin100_zf_pigeon_integrated_0728.rds')
print(data)
pdf('dimplot_bin100_zf_pigeon_integrated_0728.pdf')
DimPlot(data,label=T,cols=col37)
DimPlot(data,label=F,cols=col37)
DimPlot(data,group.by='species')
DimPlot(data,group.by='slide')
DimPlot(data,group.by='domain2',label=T)
dev.off()
pdf('dimplot_bin100_zf_pigeon_integrated_split_0728.pdf',15,8)
DimPlot(data,group.by='domain2',split.by='species',label=T)
DimPlot(data,group.by='domain2',split.by='slide',label=T)
dev.off()
pdf('dimplot_bin100_zf_pigeon_integrated_split_2_0728.pdf',15,6)
#DimPlot(data,group.by='domain2',split.by='species',label=T)
DimPlot(data,group.by='domain2',split.by='slide',label=T,cols=col37)
dev.off()
# adjust color
data$domain2=factor(data$domain2,levels=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'))
col=c('#FDBF6F','#FB9A99',RColorBrewer::brewer.pal(10, "Paired"))%>%unique()
pdf('dimplot_bin100_zf_pigeon_integrated_0730.pdf')
DimPlot(data,group.by='species')
DimPlot(data,group.by='slide')
DimPlot(data,group.by='domain2',label=T,cols=col)
DimPlot(data,group.by='domain2',label=F,cols=col)
dev.off()










