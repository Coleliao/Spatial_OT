deg <- read.csv('./Exn/deg_FindMarkers_EXN_all_2023-05-17.CSV',row.names = 1)
deg <- subset(deg,deg$cluster=='3'&deg$p_val_adj<0.05)
deg <- deg[order(deg$avg_log2FC,decreasing = T),]

data <- readRDS('./spatail_zf_pg_integ/spatial_bin100_zf_pigeon_integrated_0728.rds')
DefaultAssay(data)='Spatial'
data=NormalizeData(data)

deg2 <- subset(deg,deg$gene%in%rownames(data))
gene1 <- deg2$gene[1:5]
gene2 <- deg2$gene[1:10]
gene3 <- deg2$gene[1:20]
gene4 <- deg2$gene[1:50]
genelist <- list(top5=gene1,top10=gene2,top20=gene3,top50=gene4)
data <- AddModuleScore(data,features = genelist)
pdf('vlnplot_EXN_clu3_addmodulescore_1203.pdf')
VlnPlot(data,features = c('Cluster1','Cluster2','Cluster3','Cluster4'),stack = T,group.by = 'domain2')
dev.off()

col <- data.frame(domain2=c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc'),
                  col= c('#FDBF6F','#FB9A99',RColorBrewer::brewer.pal(9, "Paired"))%>%unique(),
                  row.names = 1)
c('#A6CEE3','#B2DF8A')

library(dplyr)
library(ggpubr)
meta <- data@meta.data
meta <- read.csv('D:/下载/bin100_EXN_clu3_deg5-10-20_addmodulescore_1203.CSV',row.names = 1)
tmp <- meta %>% group_by(domain2) %>% 
  summarise(top5=median(Cluster1),top10=median(Cluster2),top20=median(Cluster3),top50=median(Cluster4)) %>% 
  as.data.frame()
ord1 <- tmp[order(tmp$top5,decreasing = T),]$domain2
ord2 <- tmp[order(tmp$top10,decreasing = T),]$domain2


# my_comparisons <- list(c("Imc","SGC","SGF","SAC","SOp","BCS",'Ipc','SFP','SGP'))
# my_comparisons <- list(c('Imc','SFP'))
# meta$domain2 <- factor(meta$domain2,levels = ord1)
# ggviolin(meta, x = "domain2", y = "Cluster1",
#          fill = "domain2",
#          ylab = "top10_deg_score", xlab = "group",add = "boxplot", add.params = list(fill="white"))+
#   theme_bw()+ scale_fill_manual(values = col[ord1,])+ 
#   stat_compare_means(label.y = 2)+     # Global p-value
#   stat_compare_means(comparisons = my_comparisons)
# dev.off()
meta$domain2 <- factor(meta$domain2,levels = ord1)
p1 <- ggplot(data = meta,aes( x = domain2, y = Cluster1,fill = domain2))+
  geom_violin(width=1,size=0)+
  theme_bw()+ scale_fill_manual(values = col[ord1,])+
  geom_boxplot(fill='white',outlier.colour = 'black',outlier.size = 0.5,width=0.1,size=0)+
  theme(axis.line = element_line(color = 'black'),
        axis.text = element_text(colour = 'black'))+
  labs(x = 'Spatial layers',y="top5_deg_score")+
  ylim(-0.3,1.25)

meta$domain2 <- factor(meta$domain2,levels = ord2)
p2 <- ggplot(data = meta,aes( x = domain2, y = Cluster2,fill = domain2))+
  geom_violin(width=1,size=0)+
  theme_bw()+ scale_fill_manual(values = col[ord2,])+
  geom_boxplot(fill='white',outlier.colour = 'black',outlier.size = 0.5,width=0.1,size=0)+
  theme(axis.line = element_line(color = 'black'),
        axis.text = element_text(colour = 'black'))+
  labs(x = 'Spatial layers',y="top10_deg_score")+
  ylim(-0.3,1.25)

pdf('vlnplot_EXN_clu3_deg5-10_addmodulescore_inBin100spatial_1203.pdf')
cowplot::plot_grid(p1,p2,ncol = 1)
dev.off()

# wilcox.test
meta$group <- ifelse(meta$domain2=='Imc','top','low')
compare_means(Cluster2~group,data = meta,method = "wilcox.test")
compare_means(Cluster2~domain2,data = meta,method = "wilcox.test")
ggplot(meta,aes(x=domain2,y=Cluster2))+geom_violin()
#  ***






