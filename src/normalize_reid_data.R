library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
#setwd("~/Work/Writing/Papers/scPlasmodium/eLife/Supplementary Data/")


# Read in data
normalize.reid.data=function(counts.df,meta.df){
  samples=intersect(colnames(counts.df),rownames(meta.df))
  counts.mat=as.matrix(counts.df[,samples])
  meta.mat=as.matrix(meta.df[samples,])
  keep_feature <- rowSums(counts.mat> 0) > 0
  counts.mat <- counts.mat[keep_feature, ]
  sce <- SingleCellExperiment(
    assays = list(counts = counts.mat), 
    colData = meta.mat,metadata=meta.df[samples,])
  sce <- scater::calculateQCMetrics(object = sce)
  # Filter cells with low counts
  filter_by_total_counts <- (sce$total_counts > 25000)
  table(filter_by_total_counts)
  # Filter cells with low numbers of features detected
  filter_by_expr_features <- (sce$total_features > 1000)
  table(filter_by_expr_features)
  # Filter data
  sce$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
      # sufficient molecules counted
      filter_by_total_counts &
      # controls shouldn't be used in downstream analysis
      sce$is_control!='FALSE'
  )
  return(sce)
}
