# functions for some basic tools to accelerate the analysis process
# different from process.R which prefers to be complete data pipeline
# tools are likely to be some function
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)


#' @example 
#' data <- chumap(data,x=data$bin100.x,y=data$bin100.y)
chumap <- function(object,x,y){
  if(!'umap'%in%Reductions(object = object)){stop('Please do UMAP redution first!\n',call. = F)}
  if(!all(length(x)==length(y),length(x)==ncol(object))){stop('The length of x or y information does not macth the object!')}
  object@reductions$umap@cell.embeddings[,1]=unname(x)
  object@reductions$umap@cell.embeddings[,2]=unname(y)
  return(object)
}




#' @example 
#' res <- consist(meta,'seurat_clusters','sample')
consist <- function(df,level,consist,plot=F,col=NULL){
  if(!is.data.frame(df)){df <- as.data.frame(df)}
  if(!all(c(level,consist)%in%colnames(df))){stop('The level or consist is not contained in the obejct!\n',call. = F)}
  df <- df[,c(level,consist)]
  colnames(df) <- c('Level','Consist')
  df$count <- 1
  df=data.table::as.data.table(df)
  mtr <- df[, sum(count), by = .(Level,Consist)]
  mtr <- dcast(mtr,Level~mtr$Consist,value.var = 'V1') # long to wide
  mtr[is.na(mtr)]=0
  item <- colnames(mtr)[-1]
  for(i in item){ #
    mtr[,paste0(i,'_percent')] <- mtr[,i]/sum(mtr[,i])
  }
  for(i in paste0(item,'_percent')){ # 
    mtr[,paste0(i,'_related')] <- mtr[,i]/rowSums(mtr[,grep('percent$',colnames(mtr),perl = T)])
  }
  if(plot){
    tmp=melt(mtr,id.vars='Level',paste0(item,'_percent_related')) # wide to long
    if(is.null(col)){col <- c(RColorBrewer::brewer.pal(8, "Dark2"),RColorBrewer::brewer.pal(12, "Paired"),RColorBrewer::brewer.pal(8, "Set2"),brewer.pal(9, "Set1"))}
    pdf(paste0('proportion_',level,'_',consist,Sys.Date(),'.pdf'))
    print(ggplot(tmp,aes(x=Level,y=value))+geom_bar(aes(fill=variable),stat= 'identity',position='stack')+scale_fill_manual(values=col))
    dev.off()
  }
  colnames(mtr)[1] <- level
  return(mtr)
}




IsSeuratObject <- function(object){
  type <- class(object)[1]
  if(type=='Seurat'){return(TRUE)}else{return(FALSE)}
}
# Hierarchical clustering for Seurat object
hseurat <- function(object,genes=NULL,plot=T,n=2500){
  if(!IsSeuratObject(object)){stop('Please start with a Seurat Object!\n',call. = F)}
  assay <- ifelse('RNA'%in%Assays(object),'RNA','Spatial')
  DefaultAssay(object) <- assay
  object <- NormalizeData(object)
  if(is.null(genes)){
    object <- FindVariableFeatures(object,nfeatures=n)
    if(assay=='RNA'){hvg <- object@assays$RNA@var.features}else{
      hvg <- object@assays$Spatial@var.features
    }
    genes <- hvg
  }
  genes <- intersect(genes,row.names(object))
  if(length(genes)==0){stop('No genes is contained in the object!/n',call. = F)}
  mtr <- AverageExpression(object,features=gene,assays=assay)
  mtr <- mtr[[1]] %>% t()
  mtr <- scale(mtr)
  d<-dist(mtr) 
  fit1<-hclust(d,method = "average")
  pdf(paste0('HAC_plot',Sys.Date(),'.pdf'))
  plot(fit1,hang = -1,cex=.8,main = "HAC result")
  dev.off()
  return(fit1)
}




if(!require(psych)){
  cat('Notes:\n Function corSeurat is supportable because of the absence of package: \n')
}

corheatmap <- function(r,p){
  if(!require(pheatmap)){stop('No package called pheatmap!\n')}
  if (!is.null(p)){
    if(!all(dim(r)==dim(p))){stop('The matrix of R and p donot match!\n')}
    p[is.na(p)] <- 0
    sig1 <- p< 0.001
    p[sig1] <-'***'
    sig2 <- p> 0.001&p<0.01
    p[sig2] <- '**'
    sig3 <- p >0.01& p <0.05
    p[sig3] <- '*'
    p[!sig1&!sig2&!sig3]<- ''
  } else {
    p <- F
  }
  pheatmap(r,scale = "none",cluster_row = T, cluster_col = T, border=NA,
           display_numbers = p,fontsize_number = 12, number_color = "white",
           cellwidth = 20, cellheight =20) # color=
}
corSeurat <- function(object1,object2,idents1,idents2,genes=NULL,n=2500,method='spearman',plot=F,filename='cor_heatmap.pdf'){
  mlist <- list(object1,object2)
  if(!idents1%in%colnames(object1@meta.data)|!idents2%in%colnames(object2@meta.data)){
    stop('Wrong idents given!\n',call. = F)
  }
  mlist <- lapply(mlist, function(x){
    if(!IsSeuratObject(x)){stop('Only Seurat Object is supportable!\n',call. = T)}
    assay <- ifelse('RNA'%in%Assays(x),'RNA','Spatial')
    DefaultAssay(x) <- assay
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x,nfeatures=n) # default HVG is 2500
    return(x)
  })
  if(is.null(genes)){
    assay <- DefaultAssay(mlist[[1]])
    hvg1 <- mlist[[1]][[assay]]@var.features
    assay <- DefaultAssay(mlist[[2]])
    hvg2 <- mlist[[2]][[assay]]@var.features
    genes <- intersect(hvg1,hvg2)
  }else{
    genes <- intersect(genes,rownames(mlist[[1]]))
    genes <- intersect(genes,rownames(mlist[[2]]))
  }
  # genes <- Reduce(f = intersect,x = list(genes,rownames(mlist[[1]]),rownames(mlist[[2]])))
  if(length(genes)==0){stop('No shared genes in these two objects!\n')}
  cat(paste0('genes used in corrlation are ',length(genes),'\n'))
  mtr1 <- AverageExpression(mlist[[1]],features=genes,group.by = idents1)
  mtr1 <- mtr1[[1]]
  mtr2 <- AverageExpression(mlist[[2]],features=genes,group.by = idents2)
  mtr2 <- mtr2[[1]]
  res <- corr.test(mtr1,mtr2,use="pairwise",method=method, 
                   adjust="holm",alpha=0.05)
  if(!require(pheatmap)){plot=F}
  if(plot){
    pdf(filename)
    corheatmap(res$r,res$p)
    dev.off()
  }
  return(res)
}
