.libPaths("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/huangbaoqian/software/R4.0/lib")
# turn h5ad to Seurat RDS
library(rhdf5)
library(Matrix)
library(Seurat)

parser =argparse:: ArgumentParser(description="to turn the h5 to RDS")
parser$add_argument('-I','--inputFiles', help='input h5 file')
parser$add_argument('-N','--name', help='Your sample name',default='data')
parser$add_argument('-O','--out', help='out directory')
args = parser$parse_args()


mydata <- h5read(args$inputFiles,"mat")
mat <- mydata$block0_values
rownames(mat) <- mydata$axis0
colnames(mat) <- mydata$axis1
mat <- Matrix(mat, sparse = TRUE)
setwd(args$out)
meta <- read.table('metadata.tsv',sep="\t",header=T,row.names=1)
data <- CreateSeuratObject(mat,assay='Spatial',meta.data=meta)
saveRDS(data,paste0(args$name,'.rds'))


library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
tissue_lowres_image <- matrix(1, max(data$y), max(data$x))
tissue_dataitions_list <- data.frame(row.names = colnames(data),
                                    tissue = 1,
                                    row = data$y, col = data$x,
                                    imagerow = data$y, imagecol = data$x)
scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))
mat <- data@assays$Spatial@counts

seurat_spatialObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=5, min.features=5)
generate_spatialObj <- function(image, scale.factors, tissue.dataitions, filter.matrix = TRUE)
{
  if (filter.matrix) {
    tissue.dataitions <- tissue.dataitions[which(tissue.dataitions$tissue == 1), , drop = FALSE]
  }

  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(image))
  return(new(Class = 'VisiumV1',
             image = image,
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                          fiducial = scale.factors$fiducial_diameter_fullres,
                                          hires = scale.factors$tissue_hires_scalef,
                                          lowres = scale.factors$tissue_lowres_scalef),
             coordinates = tissue.dataitions,
             spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image,
                                  scale.factors = fromJSON(scalefactors_json),
                                  tissue.dataitions = tissue_dataitions_list)

spatialObj <- spatialObj[Cells(seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'
seurat_spatialObj[['slice1']] <- spatialObj
seurat_spatialObj=AddMetaData(object=seurat_spatialObj,metadata=meta)
saveRDS(seurat_spatialObj,"./data_sp.rds")

