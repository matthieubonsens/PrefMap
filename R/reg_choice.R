#' Choice of the regression model for each perceptual space
#'
#' Returns the weighted adjusted R-squared of the
#' 20 method (PCA, QPLS, PLSS, X2PLS and OX2PLS) x model (Vectorial, Circular, Elliptic and Quadratic) pairs.
#' Each consumer's adjusted R-squared is weighted by the variance of his scores.
#'
#' The spaces are built on all the consumers, as in \code{PrefMap}
#'
#' @param X the sensory characterisation dataset (products x descriptors)
#' @param Y the hedonic scores dataset (products x consumers), with the
#'   products in the same order as in \code{X}
#' @param ngroupe the number of consumer clusters
#'
#' @returns A list with one element per cluster (\code{cluster1},
#'   \code{cluster2}, ...). Each element is a data.frame of weighted adjusted
#'   R-squared, with one row per method (\code{ACP}, \code{OX2PLS},
#'   \code{X2PLS}, \code{PLSS}, \code{QPLS}) and one column per regression model
#'   (\code{Vectorial}, \code{Circular}, \code{Elliptic},
#'   \code{Quadratic}). A model with as many parameters as products has no
#'   adjusted R-squared: its column is \code{NA}.
#'
#' @export reg_choice
#'
#' @examples
#'
#'\donttest{
#' data(beer)
#'
#' res <- reg_choice(senso.beer, hedo.beer)
#'
#' res$cluster1
#'}
#'
reg_choice <- function(X, Y, ngroupe = 1) {
  dat <- prepare_data(X, Y)
  carac <- dat$carac
  hed <- dat$hed
  groupes <- carto_clusters(hed, ngroupe)$clusters[2, ]

  methodes <- c(ACP = "acp", OX2PLS = "ox2pls", X2PLS = "x2pls", PLSS = "plss", QPLS = "qpls")
  espaces <- lapply(methodes, function(m) carto_espace(carac, hed, m)$scores)

  carto_r2_groupes(espaces, hed, groupes, ngroupe)
}


#' Weighted adjusted R-squared of each method, cluster by cluster
#' @noRd
carto_r2_groupes <- function(espaces, hed, groupes, ngroupe) {
  messages <- character(0)
  res <- stats::setNames(vector("list", ngroupe), paste0("cluster", seq_len(ngroupe)))
  for (i in seq_len(ngroupe)) {
    hed.gr <- hed[, groupes == i, drop = FALSE]
    r2 <- withCallingHandlers(
      t(sapply(espaces, function(s) carto_r2aj(hed.gr, s[, 1], s[, 2]))),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    res[[i]] <- as.data.frame(r2)
  }
  for (m in unique(messages)) warning(m, call. = FALSE)
  res
}


#' Two-component perceptual space of a PrefMap method
#' @noRd
carto_espace <- function(carac, hed, method, alpha = NULL) {
  if (method == "acp") {
    res.acp <- FactoMineR::PCA(carac, scale.unit = FALSE, ncp = 2, graph = FALSE)
    return(list(scores = res.acp$ind$coord[, 1:2], alpha = NULL))
  }
  if (method == "qpls") {
    qpls <- if (is.null(alpha)) {
      QPLSniPond(carac, hed, ncomp = 2)
    } else {
      QPLSniPond(carac, hed, ncomp = 2, alpha = alpha)
    }
    return(list(scores = qpls$X_score[, 1:2], alpha = qpls$aopt))
  }
  if (method == "plss") {
    B <- carto_base_plss(carac)
    regpls <- pls::plsr(as.matrix(hed) ~ B, ncomp = 2, scale = FALSE, method = "kernelpls")
  } else {
    carac2 <- as.matrix(carto_extension(carac, method)$carac2)
    regpls <- pls::plsr(as.matrix(hed) ~ carac2, ncomp = 2, scale = TRUE, method = "kernelpls")
  }
  list(scores = unclass(regpls$scores)[, 1:2], alpha = NULL)
}


#' B-spline basis of PLSS: degree 1, one knot at the median of each descriptor
#' @noRd
carto_base_plss <- function(carac) {
  B <- matrix(0, nrow = nrow(carac), ncol = 0)
  for (j in seq_len(ncol(carac))) {
    xj <- carac[, j]
    k  <- stats::median(xj)
    M <- splines::bs(xj, degree = 1, knots = k, intercept = TRUE)
    B <- cbind(B, M)
  }
  B
}


#' Weighted adjusted R-squared of the four regression models
#' @noRd
carto_r2aj <- function(hed.gr, Dim1, Dim2) {
  modeles <- list(
    Vectorial   = cbind(Dim1, Dim2),
    Circular  = cbind(Dim1, Dim2, Dim1^2 + Dim2^2),
    Elliptic  = cbind(Dim1, Dim2, Dim1^2, Dim2^2),
    Quadratic = cbind(Dim1, Dim2, Dim1^2, Dim2^2, Dim1 * Dim2)
  )
  Y <- as.matrix(hed.gr)
  n <- nrow(Y)
  w <- apply(Y, 2, stats::var)
  Y <- Y[, w > 0, drop = FALSE]
  w <- w[w > 0]
  sct <- colSums(scale(Y, center = TRUE, scale = FALSE)^2)

  sapply(modeles, function(M) {
    k <- ncol(M)
    if (n - k - 1 <= 0) {
      warning("Mod\u00e8le \u00e0 ", k + 1, " param\u00e8tres pour ", n,
              " produits : R\u00b2 ajust\u00e9 non d\u00e9fini (NA).")
      return(NA_real_)
    }
    scr <- colSums(qr.resid(qr(cbind(1, M)), Y)^2)
    r2aj <- 1 - (scr / (n - k - 1)) / (sct / (n - 1))
    sum(w * r2aj) / sum(w)
  })
}
