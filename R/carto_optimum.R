#' Function to find the optimum product
#'
#' @param Z the liking grid
#' @param Gtot the dimension of the map
#'
#' @returns the postion of the otpimal product
#'
#' @noRd
carto_optimum <- function(Z, Gtot) {
  maxi <- max(Z)
  ou.maxi <- which(Z == maxi)
  df.opt <- Gtot[ou.maxi, ]
  colMeans(df.opt)
}
