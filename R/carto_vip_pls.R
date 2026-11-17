#' Squared VIP of a PLS2 model for a cluster of consumers
#'
#' VIP2_k = p * sum_a w_ka^2 R2_a / sum_a R2_a over the two components, where
#' w are the loading weights and R2_a is the share of the variance of the
#' scores of the cluster explained by component a (each consumer weighted by
#' the variance of their scores).
#'
#' @param regpls the PLS2 model (pls::plsr) fitted on all the consumers
#' @param hed.gr the hedonic scores of the consumers of the cluster
#' @param p the value the squared VIP sum to: the number of predictors for
#'   X2PLS and OX2PLS, the number of descriptors for PLSS
#' @param cols logical vector of the consumers of the cluster among the columns
#'   of the Y of regpls (TRUE for all)
#'
#' @returns the vector of the squared VIP of the predictors (columns of the X
#'   of regpls)
#' @noRd
carto_vip_pls <- function(regpls, hed.gr, p, cols = TRUE) {
  rsq <- pls::R2(regpls, estimate = "train")$val[1, , ]
  rsq <- rsq[cols, 2:3, drop = FALSE]
  rsq[, 2] <- rsq[, 2] - rsq[, 1]
  w.rsq <- apply(hed.gr, 2, stats::var)
  varexp <- colSums(rsq * w.rsq) / sum(w.rsq)

  vipsq <- c()
  for (k in 1:nrow(regpls$loading.weights)) {
    loading.weights <- as.matrix(regpls$loading.weights[k, 1:2])^2
    rsqxlw <- c()
    for (j in 1:length(varexp)) {
      rsqxlw <- c(rsqxlw, loading.weights[j] * varexp[j])
    }
    rsqxlw <- sum(rsqxlw)
    vipsq <- c(vipsq, p * rsqxlw / sum(varexp))
  }
  vipsq
}
