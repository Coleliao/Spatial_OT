# Visualization
library(ggplot2)
col37=c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(12, "Paired"),RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1"))
neu_markers <- c('SLC1A2','FLT1','APBB1IP','PTPRC','VIP','RELN','PDGFRB','PDGFRA',
                 'MOG','GAD1','GAD2','SLC17A7','SST','PVALB','LAMP5')
#' Call a suitable size for DimPlot and FeaturePlot
#'
#' @param obj A spatial RNA seurat object, whcih must contain spatial position info in its embeddings.
#'
#' @return a ratio which is the height devided bu the width
#' @export
#'
#' @examples
call_size <- function(obj){
  wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
  hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
  rat <- wid/hei
  return(rat)}
chumap <- function(object,x,y){
  if(!'umap'%in%Reductions(object = object)){stop('Please do UMAP redution first!\n',call. = F)}
  if(!all(length(x)==length(y),length(x)==ncol(object))){stop('The length of x or y information does not macth the object!')}
  object@reductions$umap@cell.embeddings[,1]=unname(x)
  object@reductions$umap@cell.embeddings[,2]=unname(y)
  return(object)
}


#' Single factor highlight in Spatial
#' @example 
#' Idents(data)='seurat_clusters'
#' cell=CellsByIdentities(data, idents=9)
#' DoSpotHighlight(data,cell=cell[[1]])
DoSpotHighlight <- function(object,cell,filename='dimplot.pdf',pt.size=1,sizes.highlight=1,p2=TRUE){
  object@reductions$umap@cell.embeddings[,1] <- as.numeric(as.character(object$bin100.x))
  object@reductions$umap@cell.embeddings[,2] <- as.numeric(as.character(object$bin100.y))
  call_size <- function(obj){
    wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
    hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
    rat <- wid/hei
    return(rat)}
  rat <- call_size(object);wei=rat*10;hei=10
  outlier <- setdiff(cell,colnames(object))
  if(length(outlier)>=1){warning('There are some cells not contained in the object.\n',call. = F)}
  pdf(filename,wei,hei)
  p1 <- DimPlot(object, cells.highlight=cell,cols.highlight = c("#DE2D26", "#F2F2F2"),pt.size=pt.size,sizes.highlight=sizes.highlight)
  print(p1)
  if(p2){
    p2 <- p1+NoAxes()+NoLegend()
    print(p2)
  }
  dev.off()
}


# Multiple factors highlight in Spatial


DoFeaturePlot <- function(object,features,filename='featureplot.pdf',col='blue',p2=TRUE,pt.size=1,order=F){
  object@reductions$umap@cell.embeddings[,1] <- as.numeric(as.character(object$bin100.x))
  object@reductions$umap@cell.embeddings[,2] <- as.numeric(as.character(object$bin100.y))
  call_size <- function(obj){
    wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
    hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
    rat <- wid/hei
    return(rat)
  }
  col <- switch(col,
                blue=colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100),
                viridis=colorRampPalette(viridis::viridis(20))(100),
                red=colorRampPalette(viridis::inferno(n=20,direction = -1))(100),
                redblue=colorRampPalette(c("#3288BD","white","#D53E4F"))(50),
                stereo=colorRampPalette(c('#08306B','#5ecc74', "#FFFF33",'#DC143C'))(100)) # color
  rat <- call_size(object);wei=rat*10;hei=10
  features <- unique(features)
  gene=intersect(features,c(rownames(object),colnames(object@meta.data)))
  if(length(gene)!=length(features)){
    outlier <- setdiff(features,c(rownames(object),colnames(object@meta.data)))
    if(length(gene)==0){stop('No features were found in the object!\n',call. = F)}
    cat('Warning: ',outlier,'did not exist in the object.\n')
  }
  pdf(filename,wei,hei)
  for(i in gene){
  p1 <- FeaturePlot(object,features = i,pt.size = pt.size,order=order)+scale_colour_gradientn(colours=col,name=paste0(i,'\nexpression'))
  print(p1)
  if(p2){
    print(p1+NoLegend()+NoAxes())
    }
  }
  dev.off()
}

# spatial dimplot adjust for paper
DoDimPlot <- function(object,filename='dimplot.pdf',col=NULL,p2=TRUE,pt.size=1,group='seurat_clusters'){
  object@reductions$umap@cell.embeddings[,1] <- as.numeric(as.character(object$bin100.x))
  object@reductions$umap@cell.embeddings[,2] <- as.numeric(as.character(object$bin100.y))
  call_size <- function(obj){
    wid <- max(obj@reductions$umap@cell.embeddings[,1])-min(obj@reductions$umap@cell.embeddings[,1])
    hei <- max(obj@reductions$umap@cell.embeddings[,2])-min(obj@reductions$umap@cell.embeddings[,2])
    rat <- wid/hei
    return(rat)
  }
  rat <- call_size(object);wei=rat*10;hei=10
  pdf(filename,wei,hei)
  if(is.null(col)){
    col=c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(12, "Paired"),RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1"))
  }
  meta <- colnames(object@meta.data)
  group <- intersect(group,meta)
  if(length(group)==0){stop('Please input right group types!\n',call. = F)}
  for(i in group){
    p1 <- DimPlot(object = object,cols = col,pt.size = pt.size,group.by = i)
    print(p1)
    if(p2){
      print(p1+NoLegend()+NoAxes())
    }
  }
  dev.off()
}



# cell proportion heatmap
# cell proportion barplot


scheatmap <- function(object,features,assay=NULL,idents=NULL,slot='data',plot=F,filename='average_heatmap.pdf',w=5,h=10){
  if(is.null(assay)){assay <- ifelse('RNA'%in%Assays(object),'RNA','Spatial')}
  DefaultAssay(object) <- assay
  if(assay!='SCT'){object <- NormalizeData(object)}
  if(is.null(idents)){idents <- 'seurat_clusters'}
  Idents(object)=idents
  genes <- intersect(features,rownames(object))
  if(length(genes)==0){stop('No features is found in this object!\n',call. = F)}
  mtr <- AverageExpression(object = object,assays = assay,features = genes,slot=slot)
  mtr <- mtr[[1]]
  if(plot){
    require(ComplexHeatmap)
    require(dplyr)
    mtr=t(mtr)%>%scale()%>%t()
    pdf(filename,w,h)
    p <- Heatmap(mtr,cluster_columns = F,cluster_rows = F,show_column_names = T)
    print(p)
    dev.off()
    return(mtr)
  }else(return(mtr))
}

