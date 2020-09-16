plot_by_set <- function(seurat, clust, genes, ncol, pt.size){
  sub <- seurat %>% subset(idents = clust) %>% SetIdent(value = "orig.ident")
  plt <- VlnPlot(object = sub, features = genes, ncol = ncol, pt.size = pt.size)
  return(plt)
}