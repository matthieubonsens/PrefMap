#' Localization of the nearest optimal product
#'
#' @param coord_prod localization of the study product
#' @param coord_opt localization of the optimal product
#' @param n the number of study product
#'
#' @returns the nearest product of the optimal product
#'
#' @noRd
carto_plus_proche <- function(coord_prod, coord_opt, n) {
  d1 <- coord_prod[1:n, 1] - coord_opt[1]
  d2 <- coord_prod[1:n, 2] - coord_opt[2]
  which.min(sqrt(d1^2 + d2^2))
}
