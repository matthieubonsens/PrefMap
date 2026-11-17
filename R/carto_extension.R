#' Sensory caracterization dataset augmented for X2PLS and OX2PLS
#'
#' @param carac the sensory caracterization dataset
#' @param method a character which define if we use X2PLS or OX2PLS
#'
#' @returns A list:
#' \itemize{
#'   \item \code{carac2}: the augmented matrix (products x 2p);
#'   \item \code{tens}: OX2PLS only, the 2 x 2 x p array of the inverses of the
#'     R factors, used to transform new profiles (NULL for X2PLS);
#'   \item \code{centre}, \code{echelle}: OX2PLS only, the 2 x p means and
#'     standard deviations of x and x^2 (NULL for X2PLS).
#' }
#'
#' @noRd
carto_extension <- function(carac, method) {
  tens <- NULL
  centre <- NULL
  echelle <- NULL
  if (method == "x2pls") {
    carac2 <- carac
    for (i in 1:ncol(carac)) {
      carac2 <- cbind(carac2, carac[, i]^2)
    }
  }

  if (method == "ox2pls") {
    tens <- array(0, dim = c(2, 2, ncol(carac)))
    centre <- matrix(0, nrow = 2, ncol = ncol(carac))
    echelle <- matrix(0, nrow = 2, ncol = ncol(carac))
    carac2 <- NULL
    for (i in 1:ncol(carac)) {
      X1 <- carac[, i]
      X2 <- carac[, i]^2

      X <- cbind(X1, X2)
      X <- scale(X)
      centre[, i] <- attr(X, "scaled:center")
      echelle[, i] <- attr(X, "scaled:scale")

      R <- qr(X)
      R <- qr.R(R)
      R <- R * sign(diag(R))
      R <- solve(R)
      tens[, , i] <- R
      carac2 <- cbind(carac2, X %*% R)
    }
  }
  list(carac2 = carac2, tens = tens, centre = centre, echelle = echelle)
}
