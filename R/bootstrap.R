#' Bootstrap of the weighted adjusted R-squared and standard deviation
#'
#' Resamples the consumers with replacement \code{boot} times. For each sample,
#' rebuilds the PLS spaces (OX2PLS, X2PLS, PLSS and QPLS; the PCA space does not
#' depend on the consumers) and recomputes the weighted adjusted R-squared of
#' the 20 method x model pairs.
#'
#' @param X the sensory characterisation dataset (products x descriptors)
#' @param Y the hedonic scores dataset (products x consumers), with the
#'   products in the same order as in \code{X}
#' @param ngroupe the number of consumer clusters
#' @param boot the number of bootstrap samples (default = 1000)
#'
#' @returns A list that contains:
#' \itemize{
#'   \item \code{r2aj}: the weighted adjusted R-squared on the whole sample,
#'     one data.frame per cluster (\code{cluster1}, \code{cluster2}, ...), with
#'     one row per method and one column per regression model
#'   \item \code{sd}: the standard deviation of the \code{boot} bootstrap values
#'     of the weighted adjusted R-squared, same layout as \code{r2aj}
#'   \item \code{echantillons}: the \code{boot} x consumers matrix of the
#'     resampled consumers (one row per sample, column numbers of \code{Y})
#'   \item \code{alpha}: the alpha of each QPLS component, estimated on the
#'     whole sample and used for every bootstrap sample
#' }
#' @export bootstrap
#'
#' @examples
#'
#' \donttest{
#' data(beer)
#'
#' set.seed(123)
#'
#' res <- bootstrap(senso.beer, hedo.beer, boot = 20)
#'
#' res$r2aj$cluster1
#' res$sd$cluster1
#' }
#'
bootstrap <- function(X, Y, ngroupe = 1, boot = 1000) {
  dat <- prepare_data(X, Y)
  carac <- dat$carac
  hed <- as.matrix(dat$hed)
  groupes <- carto_clusters(dat$hed, ngroupe)$clusters[2, ]

  methodes <- c(ACP = "acp", OX2PLS = "ox2pls", X2PLS = "x2pls", PLSS = "plss", QPLS = "qpls")
  espaces <- lapply(methodes, function(m) carto_espace(carac, hed, m))
  alpha <- espaces$QPLS$alpha
  r2aj <- carto_r2_groupes(lapply(espaces, `[[`, "scores"), hed, groupes, ngroupe)

  segments <- lapply(sort(unique(groupes)), function(g) which(groupes == g))
  echantillons <- matrix(0L, nrow = boot, ncol = length(groupes))
  for (b in seq_len(boot)) {
    echantillons[b, ] <- unlist(lapply(segments, function(idx)
      idx[sample.int(length(idx), length(idx), replace = TRUE)]))
  }

  valeurs <- array(NA_real_, dim = c(length(methodes), 4, ngroupe, boot))
  pb <- utils::txtProgressBar(min = 0, max = boot, style = 3)
  on.exit(close(pb))
  for (b in seq_len(boot)) {
    idx <- echantillons[b, ]
    hed.b <- hed[, idx, drop = FALSE]
    espaces.b <- lapply(methodes, function(m) {
      if (m == "acp") return(espaces$ACP$scores)
      carto_espace(carac, hed.b, m, alpha = if (m == "qpls") alpha)$scores
    })
    r2.b <- carto_r2_groupes(espaces.b, hed.b, groupes[idx], ngroupe)
    for (i in seq_len(ngroupe)) valeurs[, , i, b] <- as.matrix(r2.b[[i]])
    utils::setTxtProgressBar(pb, b)
  }

  ecarts <- lapply(seq_len(ngroupe), function(i) {
    s <- apply(valeurs[, , i, , drop = FALSE], c(1, 2), stats::sd)
    dimnames(s) <- dimnames(as.matrix(r2aj[[i]]))
    as.data.frame(s)
  })
  names(ecarts) <- names(r2aj)

  list(r2aj = r2aj, sd = ecarts, echantillons = echantillons, alpha = alpha)
}
