
#function to do paired DE analysis

do_pairwise_edgeR <- function(seurat){
  
  print("object contains active ident :")
  print(unique(seurat@active.ident))
  clust.sce <- as.SingleCellExperiment(seurat)
  clust.agg <- scater::sumCountsAcrossCells(clust.sce,
                                            ids=clust.sce@colData$rep)
  colnames(clust.agg)
  group = sapply(colnames(clust.agg), 
                 FUN = function(x){ return(strsplit(x, split = "_")[[1]][1]) } )
  y <- DGEList(counts = clust.agg, samples = colnames(clust.agg), group = group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  print(y$samples)
  design <- model.matrix(~group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf.r0_man <- glmQLFTest(fit, coef=2)
  qlf.r40_man <- glmQLFTest(fit, coef=3)
  qlf.r40_r0 <- glmQLFTest(fit, contrast = c(0,-1,1))
  qlf.any <- glmQLFTest(fit, coef = 2:3)
  return(list("r0_man" = topTags(qlf.r0_man, n = 500, p.value = 0.05),
              "r40_man" = topTags(qlf.r40_man, n = 500, p.value = 0.05),
              "r40_r0" = topTags(qlf.r40_r0, n = 500, p.value = 0.05),
              "any" = topTags(qlf.any, n = 500, p.value = 0.05)))
}

