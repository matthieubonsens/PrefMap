#' Conversion of the datasets into numeric data frames
#'
#' @param X the sensory characterisation dataset (products x descriptors)
#' @param Y the hedonic scores dataset (products x consumers)
#'
#' @returns a list: carac, X as a data frame of numeric columns, and hed, Y as
#'   a data frame of numeric columns
#' @noRd
prepare_data <- function(X, Y) {
  carac <- X
  hed <- Y

  carac <- data.frame(carac)
  hed <- data.frame(hed)

  for (i in seq_len(ncol(carac))) {
    carac[, i] <- as.numeric(carac[, i])
  }
  for (i in seq_len(ncol(hed))) {
    hed[, i] <- as.numeric(hed[, i])
  }

  list(carac = carac, hed = hed)
}
