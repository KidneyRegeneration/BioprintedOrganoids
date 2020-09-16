
# export markers as an excel file with each sheet as a cluster.

export_marker_xls <- function(markers, xls.path){
  
  #get number of clusters - note first cluster is 0
  no.of.clusters <- length(unique(markers$cluster))
  
  #generate list split by cluster
  list.by.cluster <- split(markers, markers$cluster)
  
  for (c in unique(markers$cluster)){
    xlsx::write.xlsx(list.by.cluster[[c]], 
               xls.path, 
               sheetName = paste0("Cluster_",c), 
               append=TRUE
               )}
  } 