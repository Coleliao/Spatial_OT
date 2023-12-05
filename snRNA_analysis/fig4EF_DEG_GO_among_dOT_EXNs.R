# deg
library(Seurat)
data=readRDS("./data_subclutering_EXN_2023-05-17.rds")
DefaultAssay(data)="RNA"
Idents(data)='seurat_clusters'
deg=FindMarkers(data,ident.1=c(4,7,11,12),ident.2=c(3,0,1,8),min.pct = 0.25)
write.csv(deg,paste0('deg_FindMarkers_EXN_ideng12_',Sys.Date(),'.CSV'))
source('../functions/process.R')
deg2=deg %>% subset(deg$p_val_adj<0.05)
write.csv(deg,paste0('deg_FindMarkers_EXN_ideng12_filtrate_',Sys.Date(),'.CSV'))

# volcano plot
library(ggplot2)
library(ggrepel)
#ggplot(deg,aes(x=avg_log2FC,y=-log10(p_val)))+geom_point()
deg$threshold = factor(ifelse(deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) >= 0.5, ifelse(deg$avg_log2FC>= 0 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
deg$gene=rownames(deg)
ran<-runif(nrow(deg),min = 50,max = 100)
ran<-100^-ran
deg$p_val_adj<-deg$p_val_adj+ran 
setwd("E:/pigeon_EXN&INN/EXN_enrichGO_adjust")
pdf('volcanoplot_deg_pigeon_ident.1_ident.2_adjust_0525.pdf',8,6)
ggplot(deg,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#E41A1C","#03396B","#CDDBDB"))+
  geom_text_repel(
    data = deg[deg$p_val<0.05&abs(deg$avg_log2FC)>0.5,],
    aes(label = gene),
    size = 4,
    segment.color = "black", show.legend = FALSE )+
  theme_bw()+
  theme(
    legend.title = element_blank()
  )+
  ylab('-log10 (p_val)')+
  xlab('avg_log2FC')+
  geom_vline(xintercept=c(-0.5,0.5),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)
dev.off()

# GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
data<-read.csv("./deg_FindMarkers_EXN_ideng12_adjust2023-05-25.csv")
data$gene <- data$X
ids=bitr(data$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ids <- na.omit(ids)
data=merge(data,ids,by.x='gene',by.y='SYMBOL')
#by.x，by.y：指定依据哪些行合并数据框，默认值为相同列名的列
#data$degree=ifelse(data$avg_log2FC>0,'up','down')
data2 <- subset(data,(data$avg_log2FC)>0.25)#<-0.5 down >0.5up
enrich.go.BP = enrichGO(gene = data2$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",#这里可以选BP,MF,CC，ALL
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)
library(ggplot2)
dev.new()
pdf(file="GO_EXN_enrichGO_adjust_0.25.pdf")
dotplot(enrich.go.BP,font.size=9,title="GO of EXN 0.25 ")+theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5))
dev.off()
res <- enrich.go.BP@result   
item <- res$geneID
item_name <- c()
for(i in item){
  tmp <- strsplit(i,split = '/')
  tmp <- tmp[[1]]
  ids=bitr(tmp, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  gname <- ids$SYMBOL
  tmp <- paste(gname,collapse = '/')
  item_name <- append(item_name,tmp)
}
res$gene_name <- item_name
write.csv(res,"GO_EXN_enrich_adjust_0.25.csv")
# select key GO items
df <- res
go <- c('GO:0050804','GO:0050808','GO:0061564','GO:0097485','GO:0006836','GO:0050890','GO:0007611','GO:0001764','GO:0060560')
df2 <- subset(df,df$ID%in%unlist(go))
df2$Cluster <- 'Up'
df2$GeneRatio <- DOSE::parse_ratio(df2$GeneRatio)
df2$Description<- factor(df2$Description,levels = rev(unique(df2$Description)))
pdf('D:/Documents/BGI/B_cross_species/OT视顶盖/fig/fig4/EXN/红蓝比较/red_blue_up_GO_selected_0525.pdf',5,5)
ggplot(df2,aes(x=Cluster,y = Description,size = GeneRatio, color = p.adjust))+
  geom_point() +
  scale_y_discrete(position = "right",labels=function(x) stringr::str_wrap(x, width=80)) + 
  scale_color_gradientn(colours = rev(viridis::viridis(20))) + 
  cowplot::theme_cowplot() +
  ylab("") + xlab("") + theme_bw() + 
  theme(
    axis.text.x = element_text(size=10, angle=90, hjust=0.5, color="black"),
    axis.text.y = element_text(size=10, color="black",margin = c(5,1)),
    axis.title = element_text(size=14)
  )
dev.off()

