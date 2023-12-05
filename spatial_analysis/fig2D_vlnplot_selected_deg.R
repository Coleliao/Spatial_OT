library(ggplot2)
library(Seurat)
data=readRDS('spatial_bin100_zf_pigeon_integrated_0728.rds')
DefaultAssay(data)='Spatial'
data=NormalizeData(data)
genes=c('TAC1','NRXN1','PENK','PCDH17','SNCG','CCK','NPTX2','SFRP2','SYNGR3')
library(reshape2)
library(ggplot2)
vln.df=as.data.frame(data[["Spatial"]]@data[genes,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene") # reshape2 dcast-melt
colnames(vln.df)[c(2,3)]=c("CB","exp") # CB is cell barcodes
data$CB=colnames(data)
anno=data@meta.data[,c("CB","domain2")]
vln.df=dplyr::inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = genes) # gene order
vln.df$domain2=factor(vln.df$domain2,levels = c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc')) # type order
p=ggplot(vln.df,aes(domain2,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
pdf('vlnplot_SGP_selected_degs_0731.pdf')
print(p)
dev.off()
#ggsave("vln2.pdf",width = 11,height = 22,units = "cm")

## QC Vlnplot
data=readRDS('spatial_bin100_zf_pigeon_integrated_0728.rds')
data$slide=factor(data$slide,levels=c('D6','C6','E3','F2','A6'))
pdf('vlnplot_QC_bin100_0921.pdf',12,10)
VlnPlot(data,group.by='slide',features='nCount_Spatial',pt.size=0,cols=col37)
VlnPlot(data,group.by='slide',features='nFeature_Spatial',pt.size=0,cols=col37)
dev.off()

