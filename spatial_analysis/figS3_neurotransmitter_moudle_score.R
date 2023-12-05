# neurotransmitter geneset score

# d6
d6=readRDS('./pigeon_d6_anno_0328.rds')
dopa <- c('DRD1','DRD5','DRD2','DRD3','DRD4')
gaba <- c('GABRA1', 'GABRA2', 'GABRA3', 'GABRA4', 'GABRA5', 'GABRA6','GABRB1', 
          'GABRB2', 'GABRB3','GABRG1', 'GABRG2', 'GABRG3','GABRD','GABRE','GABRP',
          'GABRQ','GABRR1','GABRR2', 'GABRR3','GABBR1','GABBR2')
glu <- c('GRIA1','GRIA2','GRIA3','GRIA4','GRIK1','GRIK2','GRIK3','GRIK4','GRIK5',
         'GRIN1','GRIN2A','GRIN2B','GRIN2C','GRIN2D','GRIN3A','GRIN3B')
ht5 <- c('HTR1A','HTR1B','HTR1D','HTR1E','HTR1F','HTR2A','HTR2B','HTR2C','HTR3A', 
         'HTR3B', 'HTR3C', 'HTR3D', 'HTR3E','HTR4','HTR5A', 'HTR5BP','HTR6','HTR7')
ach <- c('CHRM1','CHRM2','CHRM3','CHRM4','CHRM5','CHRNA1','CHRNA10','CHRNA2','CHRNA3',
         'CHRNA4','CHRNA5','CHRNA6','CHRNA7','CHRNA9','CHRNB1','CHRNB2','CHRNB3',
         'CHRNB4','CHRND','CHRNE','CHRNG')
genelist <- list(Dopamine=dopa,GABA=gaba,Glutamate=glu,Serotonin=ht5,Acetylcholine=ach)
d6 <- AddModuleScore(d6,features = genelist)
col <- colorRampPalette(viridis::inferno(n=20,direction = -1))(100)
pdf('score_neurotransmitter_receptors_0516.pdf')
for(i in paste0('Cluster',1:5)){
  print(FeaturePlot(d6,features = i,pt.size=0.4))
}
dev.off()

# c6
c6=readRDS('./pigeon_C6_anno_0330.rds')
c6 <- AddModuleScore(c6,features = genelist)
c6@reductions$umap@cell.embeddings[,1]=c6$bin100.x
c6@reductions$umap@cell.embeddings[,2]=c6$bin100.y
pdf('score_neurotransmitter_receptors_c6_0516.pdf')
for(i in paste0('Cluster',1:5)){
    #print(FeaturePlot(c6,features = i,pt.size=0.4)+scale_colour_gradientn(colours=col,name=paste0(i,'\nexpression')))
    print(FeaturePlot(c6,features = i,pt.size=0.4))
}
dev.off()

# e3
e3=readRDS('./pigeon_E3_anno_0330.rds')
e3 <- AddModuleScore(e3,features = genelist)
e3@reductions$umap@cell.embeddings[,1]=e3$bin100.x
e3@reductions$umap@cell.embeddings[,2]=e3$bin100.y
pdf('score_neurotransmitter_receptors_e3_0516.pdf')
for(i in paste0('Cluster',1:5)){
    print(FeaturePlot(e3,features = i,pt.size=0.4))
}
dev.off()



#contains pigeons and zebrafinches
data=readRDS('./spatail_zf_pg_integ/spatial_bin100_zf_pigeon_integrated_0728.rds')
DefaultAssay(data)='Spatial'
data=NormalizeData(data)
data=AddModuleScore(data,features = genelist)
data$domain2=factor(data$domain2,levels=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'))
cols <- c('#FDBF6F','#FB9A99',RColorBrewer::brewer.pal(10, "Paired"))
pdf('vlnplot_neuotransmitter_score_5slides_0801.pdf',6,12)
VlnPlot(data,features=paste0('Cluster',1:5),group.by='domain2',stack=T,flip=T,cols=cols)
dev.off()


glut <- c('SLC1A3','SLC1A2','SLC1A1','SLC1A6','SLC1A7','SLC17A7','SLC17A6','SLC17A8')
gabat <- c('SLC6A1','SLC6A13','SLC6A11','SLC6A12','SLC32A1')
monot <- c('SLC6A3','SLC6A2','SLC6A4','SLC18A1','SLC18A2') # 单胺类（儿茶酚胺包括多巴胺、去甲肾上腺素和肾上腺素。吲哚胺主要是5-羟色胺）
acht <- c('SLC18A3')
trans <- list(GABA=gabat,Glutamate=glut,Monoamine=monot,Acetylcholine=acht)
d6 <- AddModuleScore(d6,features = trans)
col <- colorRampPalette(viridis::inferno(n=20,direction = -1))(100)
pdf('score_neurotransmitter_transporters_0525.pdf')
for(i in paste0('Cluster',1:4)){
  print(FeaturePlot(d6,features = i,pt.size=0.4))
}
dev.off()
data=AddModuleScore(data,features = trans)
cols <- c('#FDBF6F','#FB9A99',RColorBrewer::brewer.pal(10, "Paired"))
pdf('vlnplot_neuotransmitter_transporters_score_5lides_0801.pdf',6,12)
VlnPlot(data,features=paste0('Cluster',1:4),group.by='domain2',stack=T,flip=T,cols=cols[2:5])
dev.off()


# Kruskal-Wallis Rank Sum Test
data=readRDS('./spatail_zf_pg_integ/spatial_bin100_zf_pigeon_integrated_0728.rds')
DefaultAssay(data)='Spatial'
data=NormalizeData(data)
data=AddModuleScore(data,features = genelist)
data$domain2=factor(data$domain2,levels=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'))
meta=data@meta.data
colnames(data@meta.data)=c(colnames(data@meta.data)[1:(ncol(meta)-5)],names(genelist))
write.csv(data@meta.data,'bin100_integrated_zf_pg_receptor_score_1201.CSV') # cluster1~5
# 
trans <- list(GABA=gabat,Glutamate=glut,Monoamine=monot,Acetylcholine=acht)
data=AddModuleScore(data,features = trans)
meta=data@meta.data
colnames(data@meta.data)=c(colnames(data@meta.data)[1:(ncol(meta)-4)],names(trans))
write.csv(data@meta.data,'bin100_integrated_zf_pg_transportor_score_1201.CSV') # cluster1~4

compare_means(GABA~domain2,data = trans,method = "kruskal.test")
compare_means(Glutamate~domain2,data = trans,method = "kruskal.test")
compare_means(Monoamine~domain2,data = trans,method = "kruskal.test")
compare_means(Acetylcholine~domain2,data = trans,method = "kruskal.test")




