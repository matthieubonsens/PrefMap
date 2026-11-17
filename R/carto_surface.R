#' Preference surface on the grid
#'
#' Fits the regression model to the scores of each consumer of the cluster and
#' counts, at each grid point, the consumers whose predicted score exceeds
#' their mean score.
#'
#' @param hed.gr the hedonic scores of the consumers of the cluster
#' @param Dim1,Dim2 the coordinates of the products on the two dimensions
#' @param Gtot the grid, from carto_grid
#' @param regmod the regression model: 1 vectorial, 2 circular, 3 elliptic,
#'   4 full quadratic
#'
#' @returns Z, the number of consumers at each grid point
#' @noRd
carto_surface <- function(hed.gr, Dim1, Dim2, Gtot, regmod) {
  formules <- list(
    Liking ~ Dim1 + Dim2,                                         # vectoriel
    Liking ~ Dim1 + Dim2 + I(Dim1^2 + Dim2^2),                    # circulaire
    Liking ~ Dim1 + Dim2 + I(Dim1^2) + I(Dim2^2),                 # elliptique
    Liking ~ Dim1 + Dim2 + I(Dim1^2) + I(Dim2^2) + Dim1 * Dim2    # quadratique complet
  )

  Z <- rep(0, nrow(Gtot))
  if (regmod %in% 1:4) {
    for (j in seq_len(ncol(hed.gr))) {
      dta.mod <- cbind.data.frame(hed.gr[, j], Dim1, Dim2)
      colnames(dta.mod)[1] <- "Liking"
      mod.c <- stats::lm(formules[[regmod]], dta.mod)
      Yc <- stats::predict(mod.c, newdata = Gtot)
      m <- mean(hed.gr[, j])
      notes <- ifelse(Yc > m, 1, 0)
      Z <- Z + notes
    }
  }
  Z
}
