# functions for data processes or new analysis methods
# complete and well-rounded workflows are not contained here, which will be packaged into some independent scripts. (如果整理了的化)
library(Seurat)
library(ggplot2)
library(dplyr)


## 基础快速分析流程（标准化-降维聚类，不包含QC和regress out）
quickSeurat <- function(object,assay=NULL,method='log',dims=1:30,resolution=0.5){
  if(is.null(assay)){
    assay <- ifelse('RNA'%in%Assays(object),'RNA','Spatial')
    DefaultAssay(object) <- assay
  }else{
    DefaultAssay(object)=assay
  }
  if(!method%in%c('log','sct')){stop('Please input right methods for normalization!\n',call. = F)}
  if(method=='log'){
    object <- object %>% 
      NormalizeData() %>% FindVariableFeatures(nfeatures=3000) %>% ScaleData() %>% RunPCA() %>% 
      FindNeighbors(dims = dims) %>% RunUMAP(dims = dims) %>% FindClusters(resolution=resolution)
  }else{
    object <- object %>% SCTransform(assay = assay) %>% 
      RunPCA() %>% FindNeighbors(dims = dims) %>% 
      RunUMAP(dims = dims) %>% FindClusters(resolution=resolution)
  }
  return(object)
}


### 整合 ###
#' @example 
#' data <- DoIntegration(data1,data2,method='sct')
DoIntegration <- function(object1,object2,...,method='sct',nfeatures=3000,dims=1:30,resolution=1.0,gene_filter=T){
  object.list <- c(list(object1,object2),list(...))
  for(i in 1:length(object.list)){
    DefaultAssay(object.list[[i]]) <- ifelse('RNA'%in%Assays(object.list[[i]]),'RNA','Spatial') # 特别情况手动修改
  }
  # 
  if(gene_filter){
    gene_list <- lapply(object.list, rownames)
    genes <- Reduce(intersect,gene_list)
    if(length(genes)==0){stop('No genes shared in two object!\n',call. = F)}
    print(paste0('The number of genes shared in objects for integration is ',length(genes),'\n'))
    for(i in 1:length(object.list)){
      object.list[[i]] <- object.list[[i]][genes,]
    }
  }
  if(!method%in%c('sct','log')){stop('Please input right method for integaration!')}
  if(method=='sct'){
    object.list <- lapply(X = object.list, FUN = function(x){
      x <- SCTransform(x,assay = Assays(x)[1])
      return(x)
    })
    features <- SelectIntegrationFeatures(object.list = object.list,nfeatures=nfeatures)
    object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = features) # 这步是SCT做整合需要的 就和SCT做差异分析一样
    object.list <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT",anchor.features = features)
    sam.combined <- IntegrateData(anchorset = object.list, normalization.method = "SCT")
    sam.combined <- RunPCA(sam.combined, verbose = FALSE)
    sam.combined <- RunUMAP(sam.combined, reduction = "pca", dims = dims)
    sam.combined <- FindNeighbors(sam.combined, reduction = "pca", dims = dims)
    sam.combined <- FindClusters(sam.combined, resolution = resolution)
  }else{
    object.list <- lapply(X = object.list, FUN = function(x){
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, selection.method = "vst",  nfeatures = 3000, verbose = FALSE)
      return(x)
    })
    sam.anchors <- FindIntegrationAnchors(object.list=object.list, dims = 1:30)
    sam.combined <- IntegrateData(anchorset = sam.anchors, dims = 1:30)
    sam.combined <- sam.combined %>%
      ScaleData() %>% RunPCA() %>% FindNeighbors(dims = dims) %>%
      RunUMAP(dims = dims) %>% FindClusters(resolution=resolution)
  }
  return(sam.combined)
}

### 亚群重聚类 ###
# subset-(重整合)-重聚类  #目前只支持log-normalize
#' @examples
#' data <- Subclutering(data,subset=c(1,2,3),idents='seurat_clusters')
Subclutering <- function(data,subset,idents=NULL,split=NULL,assay='RNA',method='log',dims=1:30,resolution=1.0){
  DefaultAssay(data)=assay
  if(is.null(idents)){idents <- Idents(data)} # use default idents
  Idents(data)=idents
  if(!all(subset%in%levels(data))){stop('Please input right idents in the obeject',call. = F)}
  if(!method%in%c('sct','log')){stop('Please input right normalization method!',call. = F)}
  data=subset(data,idents=subset)
  print(data)
  if(is.null(split)){
    if(method=='log'){
      data=NormalizeData(data)%>% FindVariableFeatures(nfeatures=3500) %>% ScaleData()
    }else{
      data=SCTransform(data,assay=assay)
    }
      data=RunPCA(data) %>% FindNeighbors(dims = dims) %>%
       RunUMAP(dims = dims) %>% FindClusters(resolution=resolution)
  }else{
    sam.list=SplitObject(data,split.by=split)
    if(method=='log'){
      sam.list=lapply(sam.list, function(x){
        x=NormalizeData(x)%>% FindVariableFeatures(nfeatures=3500) %>% ScaleData()
        return(x)
      })
      sam.anchors <- FindIntegrationAnchors(object.list=sam.list, dims = 1:30)
      sam.combined <- IntegrateData(anchorset = sam.anchors, dims = 1:30)
      sam.combined <- sam.combined %>%
        ScaleData() %>% RunPCA() %>% FindNeighbors(dims = dims) %>%
        RunUMAP(dims = dims) %>% FindClusters(resolution=resolution)
      data <- sam.combined
    }else{
      sam.list=lapply(sam.list, function(x){
        x <- SCTransform(x,assay = Assays(x)[1])
        return(x)
      })
      features <- SelectIntegrationFeatures(object.list = sam.list,nfeatures=2500)
      sam.list <- PrepSCTIntegration(object.list = sam.list, anchor.features = features) # 这步是SCT做整合需要的 就和SCT做差异分析一样
      sam.list <- FindIntegrationAnchors(object.list = sam.list, normalization.method = "SCT",anchor.features = features)
      sam.combined <- IntegrateData(anchorset = sam.list, normalization.method = "SCT")
      sam.combined <- RunPCA(sam.combined, verbose = FALSE)
      sam.combined <- RunUMAP(sam.combined, reduction = "pca", dims = dims)
      sam.combined <- FindNeighbors(sam.combined, reduction = "pca", dims = dims)
      sam.combined <- FindClusters(sam.combined, resolution = resolution)
      data <- sam.combined
    }
  }
  return(data)
}



# 两个物种之间的同源基因转化

## doubletfinder去双胞（建议对文库或者同一样本使用，整合后的数据效果不佳的）
require(DoubletFinder)
rundf <- function(data){
  data=quickSeurat(data) # 这步是必须的
  sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  DoubletRate = ncol(data)*8*1e-6 # 每一千+0.8%计算
  homotypic.prop <- modelHomotypic(data$seurat_clusters)
  nExp_poi <- round(DoubletRate*ncol(data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, 
                           nExp = nExp_poi.adj, reuse.pANN = F, sct = T) # sct=F
  colnames(data@meta.data)[ncol(data@meta.data)]="DoubletFinder"
  Idents(data) <- 'DoubletFinder'
  data <- subset(data,idents='Singlet')
  return(data)
}



