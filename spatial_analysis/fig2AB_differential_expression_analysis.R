library(Seurat)
library(ggplot2)
library(dplyr)
data=readRDS('spatial_bin100_zf_pigeon_integrated_0728.rds')
Idents(data)='domain2'
DefaultAssay(data)='Spatial'
data=NormalizeData(data)%>% FindVariableFeatures() %>% ScaleData()
deg=FindAllMarkers(data,only.pos=T,logfc.threshold = 0.1)
write.csv(deg,'./deg_OT_domain2_0729.CSV')
# visulization
deg$cluster <- factor(deg$cluster,levels = rev(c('SOp','SGF','SGC','SAC','BCS','SGP','SFP','Imc','Ipc')))
deg$label <- ifelse(deg$p_val_adj<0.05,"adjust P-val<0.05","adjust P-val>=0.05")
top10 <- deg %>% subset(deg$p_val_adj<0.05) %>% group_by(cluster) %>% slice_max(order_by = abs(avg_log2FC),n = 10)
deg$size <- case_when(!(deg$gene %in% top10$gene)~ 1,deg$gene %in% top10$gene ~ 2)
dt <- filter(deg,size==1)
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
ymax <- aggregate(deg$avg_log2FC,by=list(deg$cluster),max)
ymax <- ymax$x %>% signif(digits = 3)+0.2
ymin <- aggregate(deg$avg_log2FC,by=list(deg$cluster),min)
ymin <- ymin$x %>% signif(digits = 3)-0.2
dfbar1<-data.frame(x=unique(top10$cluster),y=ymax) # 统一坐标轴位置
dfbar2<-data.frame(x=unique(top10$cluster), y=ymin)
p <- ggplot()+
  geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar2,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
dfcol<-data.frame(x=unique(top10$cluster),y=0,label=unique(top10$cluster))
mycol <- RColorBrewer::brewer.pal(nrow(dfcol),'Set1')
p2 <- p+geom_tile(data = dfcol,
                  aes(x=x,y=y),
                  height=0.6,
                  color = "black",
                  fill = mycol,
                  alpha = 1,
                  show.legend = F)
library(ggrepel)
p3 <- p2+
  geom_text_repel(
    data=top10,
    aes(x=cluster,y=avg_log2FC,label=gene),
    force = 1.2,
    max.overlaps = 40
    #arrow = arrow(length = unit(0.008, "npc"), # arrow
    #              type = "open", ends = "last"),
  )+
  scale_color_manual(name=NULL, # color of dots
                     values = c("red","black"))
# axis label
p4 <- p3+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="black")+theme()
# theme adjustment
p5 <- p4+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               linewidth = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
pdf('./deg_birds_OT_layers_volcanoplot_0729.pdf',12,6)
p5
dev.off()

