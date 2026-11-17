#' Back to the scale of the descriptors
#'
#' @param x sensory profil of a descriptor
#' @param carac the sensory caracterisation dataset
#' @param borne0 a logical value. If TRUE, negative values are set to 0
#'
#' @returns the sensory profil rescaling by his mean
#'
#' @noRd
carto_recentre <- function(x, carac, borne0 = TRUE) {
  for (j in 1:ncol(carac)) {
    x[j] <- x[j] + mean(carac[, j])
    if (borne0 && x[j] < 0) {
      x[j] <- 0
    }
  }
  x
}
