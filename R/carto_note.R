#' Function to calculate the hedonic score of the optimal product
#'
#' @param hed.gr the hedonic scores dataset by cluster
#' @param Dim1 the first dimension of the map
#' @param Dim2 the second dimension of the map
#' @param coord_opt the position in the map of the optimal product
#'
#' @returns the hedonic score for the optimal product
#'
#' @noRd
carto_note <- function(hed.gr, Dim1, Dim2, coord_opt, regmod) {
  formules <- list(
    Liking ~ Dim1 + Dim2,
    Liking ~ Dim1 + Dim2 + I(Dim1^2 + Dim2^2),
    Liking ~ Dim1 + Dim2 + I(Dim1^2) + I(Dim2^2),
    Liking ~ Dim1 + Dim2 + I(Dim1^2) + I(Dim2^2) + Dim1 * Dim2)
  ideal <- data.frame(Dim1 = coord_opt[1], Dim2 = coord_opt[2])
  notes <- c()
  for (j in seq_len(ncol(hed.gr))) {
    dta.mod <- cbind.data.frame(hed.gr[, j], Dim1, Dim2)
    colnames(dta.mod)[1] <- "Liking"
    mod.c <- stats::lm(formules[[regmod]], dta.mod)
    notes <- c(notes, stats::predict(mod.c, newdata = ideal))
  }
  mean(notes)
}
