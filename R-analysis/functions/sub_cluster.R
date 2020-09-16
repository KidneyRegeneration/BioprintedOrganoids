
#Function to do basic re-transformation and clustering on a subset of cells

sub_cluster <- function(seurat){

#re-calc SCtransform on subset of cells
seurat <- SCTransform(seurat, verbose = F)

#vars.to.regress = c("S.Score","G2M.Score")

seurat <- RunPCA(seurat, assay = "SCT")
seurat <- RunUMAP(seurat, dims=1:30)

seurat <- FindNeighbors(seurat, dims = 1:30, verbose = F)
seurat <- FindClusters(seurat, resolution = seq(0.1, 1, 0.1), verbose = F)

return(seurat)

}


