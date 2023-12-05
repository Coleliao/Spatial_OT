library(ggplot2)
library(ggrepel)
library(clusterProfiler)
deg <- read.csv("./deg_OT_domain2_0729.CSV",row.names = 1)
deg <- subset(deg,deg$p_val_adj<0.05)
top50 <- deg %>% group_by(cluster) %>% slice_max(n=50,order_by = avg_log2FC)
table(deg$cluster)
genelist <- list()
deg <- as.data.frame(top50)
for(i in unique(deg$cluster)){
  gene <- deg[which(deg$cluster==i),'gene']
  ids=bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ids <- na.omit(ids)
  gene <- ids$ENTREZID
  genelist[[i]] <- gene
}
sam.go.BP <- compareCluster(genelist,
                            fun = "enrichGO",
                            OrgDb = "org.Hs.eg.db",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05
)
p <- dotplot(sam.go.BP,font.size=5,title="GO of stromal cells",showCategory = 10) 
df <- p$data
pdf('./fig2/OT_bin100_layers_GO_0921.pdf',9,14)
ggplot(df,aes(x=Cluster,y = Description,size = GeneRatio, color = p.adjust))+
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


# focus on SGP and plot
sgp <- subset(deg,deg$cluster=='SGP')
ids=bitr(sgp$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ids <- na.omit(ids)
gene <- ids$ENTREZID
sgp.go.BP <- enrichGO(gene,
                      OrgDb = "org.Hs.eg.db",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05
)
barplot(sgp.go.BP)
df <- sgp.go.BP@result
item <- df$geneID
item_name <- c()
for(i in item){
  tmp <- strsplit(i,split = '/')
  tmp <- tmp[[1]]
  ids=bitr(tmp, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  gname <- ids$SYMBOL
  tmp <- paste(gname,collapse = '/')
  item_name <- append(item_name,tmp)
}
df$gene_name <- item_name
write.csv(df,'./fig2/OT_SGP_GO_0731.csv')
selected_cate <- c('adult behavior','adult behavior','dopamine transport','sensory perception of pain',
                   'learning','response to epinephrine','catecholamine secretion','regulation of behavior',
                   'cognition','associative learning','midbrain development','feeding behavior','behavioral fear response')
sgp.go.BP@result$geneID=item_name
pdf('./fig2/OT_SGP_GO_selected_731.pdf',8,6)
cnetplot(sgp.go.BP,showCategory = selected_cate,layout = 'gem',color_category='#A50026',color_gene='#1A9850')
dev.off()
