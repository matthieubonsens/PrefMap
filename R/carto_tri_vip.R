#' Sorting the sensory profils by importance in projection
#'
#' @param x the sensory profil matrix
#' @param vip the Variables Importance in Projection
#' @param noms the names of the descriptors
#'
#' @returns a matrix that contains the profils sorting by importance (VIP)
#'
#' @noRd
carto_tri_vip <- function(x, vip, noms) {
  x_sort <- c()
  noms_sort <- c()
  vip2 <- vip
  for (j in 1:ncol(x)) {
    max <- which.max(vip2[1, ])
    vip2[1, max] <- 0
    x_sort <- cbind(x_sort, x[, max])
    noms_sort <- cbind(noms_sort, noms[max])
  }
  x_sort <- as.matrix(x_sort)
  colnames(x_sort) <- noms_sort
  x_sort
}
