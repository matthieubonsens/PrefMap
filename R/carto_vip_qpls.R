#' Squared VIP of the QPLS model
#'
#' Same formula as carto_vip_pls, with the X weights of QPLSniPond and the
#' linear part of the variance of the scores explained by each component.
#'
#' @param qpls2 the output of QPLSniPond
#' @param p the number of descriptors (the squared VIP sum to p)
#' @param hed.gr the hedonic scores of the consumers of the cluster; if NULL,
#'   the variance of Y explained by QPLSniPond is used
#'
#' @returns the vector of the p squared VIP
#' @noRd
carto_vip_qpls <- function(qpls2, p, hed.gr = NULL) {
  if (is.null(hed.gr)) {
    rsq <- qpls2$VarExpl[, 2]
  } else {
    Yc <- scale(as.matrix(hed.gr), scale = FALSE)
    rsq <- sapply(1:2, function(k)
      1 - sum(qr.resid(qr(cbind(1, qpls2$X_score[, 1:k])), Yc)^2) / sum(Yc^2))
  }
  rsq[2] <- rsq[2] - rsq[1]
  varexp <- rsq

  vipsq <- c()
  for (k in 1:nrow(qpls2$X_weight)) {
    loading.weights <- as.matrix(qpls2$X_weight[k, ])^2
    rsqxlw <- c()
    for (j in 1:length(varexp)) {
      rsqxlw <- c(rsqxlw, loading.weights[j] * varexp[j])
    }
    rsqxlw <- sum(rsqxlw)
    vipsq <- c(vipsq, p * rsqxlw / sum(varexp))
  }
  vipsq
}
