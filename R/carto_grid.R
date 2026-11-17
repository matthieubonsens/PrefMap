#' Grid building function for the map
#'
#' @param Dim1 the first dimension of the sensory analysis
#' @param Dim2 the second dimension of the sensory analysis
#' @param bornes the limits of the grid
#'
#' @returns a list corresponding to the grid
#'
#' @noRd
carto_grid <- function(Dim1, Dim2, bornes = NULL) {
  if (is.null(bornes)) {
    G1 <- seq(min(Dim1), max(Dim1), length.out = 100)
    G2 <- seq(min(Dim2), max(Dim2), length.out = 100)
  } else {
    G1 <- seq(min(min(Dim1), bornes$xmin), max(max(Dim1), bornes$xmax), length.out = 100)
    G2 <- seq(min(min(Dim2), bornes$ymin), max(max(Dim2), bornes$ymax), length.out = 100)
  }
  Gtot <- expand.grid(G1, G2)
  colnames(Gtot) <- c("Dim1", "Dim2")
  list(G1 = G1, G2 = G2, Gtot = Gtot)
}
