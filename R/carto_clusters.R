#' Clustering of the consumers on their hedonic scores
#'
#' Hierarchical CLV clustering of the consumers (ClustVarLV, method "local",
#' without scaling), consolidated by CLV k-means from the partition into
#' ngroupe clusters (strategy "kplusone", rho = 0).
#'
#' @param hed the hedonic scores dataset (products x consumers)
#' @param ngroupe the number of clusters
#'
#' @returns the output of ClustVarLV::CLV_kmeans; its row 2 of clusters gives
#'   the cluster of each consumer
#' @noRd
carto_clusters <- function(hed, ngroupe) {
  cah <- ClustVarLV::CLV(hed, method = "local", sX = FALSE, graph = FALSE)
  ClustVarLV::CLV_kmeans(hed, method = "local", sX = FALSE,
                         clust = cah[[ngroupe]]$clusters[2, ],
                         iter.max = 1000, strategy = "kplusone", rho = sqrt(0))
}
