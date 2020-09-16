
get_props_table <- function(seurat){
  
  global.props <- as.data.frame(prop.table(table(seurat@meta.data$rep, seurat@active.ident), margin = 1))
  
  colnames(global.props) <- c("replicate","cluster","prop")
  
  global.props["condition"] <- sapply(as.character(global.props$replicate), 
                                      FUN = function(x){ return(strsplit(x, split = "_")[[1]][1]) } )
  
  global.props <- global.props %>% group_by(cluster) %>% dplyr::mutate('norm.prop' = prop / mean(prop))
  
  return(global.props)
  
}
