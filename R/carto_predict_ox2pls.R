#' Predicted hedonic scores of OX2PLS at given points of the map
#'
#' Each point is projected back into the descriptor space (regression of the
#' descriptors on the two PLS components), the intensities are bounded between
#' 0 and 10, the OX2PLS predictors are recomputed with the parameters estimated
#' on the products, and the scores are predicted with the coefficients of the
#' PLS2 model with two components.
#'
#' @param carac the sensory characterisation dataset (products x descriptors)
#' @param regpls the PLS2 model fitted on the OX2PLS predictors
#' @param Gtot the points of the map, a matrix or data frame with two columns
#'   (Dim1, Dim2): the grid of carto_grid or the ideal product
#' @param tens,centre,echelle the elements of the same names returned by
#'   carto_extension
#' @param cols the consumers to predict (columns of the Y of regpls): a logical
#'   vector, or TRUE for all
#'
#' @returns a points x consumers matrix of predicted scores, centred on the
#'   mean score of each consumer
#' @noRd
carto_predict_ox2pls <- function(carac, regpls, Gtot, tens, centre, echelle, cols = TRUE) {
  reslm <- stats::lm(as.matrix(carac) ~ regpls$scores[, 1:2])
  Gtot <- cbind(1, Gtot)
  p <- as.matrix(Gtot) %*% stats::coef(reslm)
  p[p < 0] <- 0
  p[p > 10] <- 10
  pr <- NULL
  nc <- ncol(carac)
  for (j in 1:nc) {
    zj <- cbind((p[, j]   - centre[1, j]) / echelle[1, j],
                (p[, j]^2 - centre[2, j]) / echelle[2, j])
    pr <- cbind(pr, zj %*% tens[, , j])
  }
  pr <- sweep(pr, 2, regpls$scale, "/")
  pr <- sweep(pr, 2, regpls$Xmeans, "-")
  pr %*% regpls$coefficients[, , 2][, cols, drop = FALSE]
}
