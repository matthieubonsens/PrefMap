#' External Preference Mapping
#'
#' Builds the two-dimensional perceptual space with the chosen method,
#' regresses the hedonic scores of every consumer on it with the chosen model
#' and draws, for each cluster of consumers, the preference map. The ideal product is the centre of the area
#' where this number is maximal. For the PLS methods, the function also gives
#' the VIP of the descriptors, the sensory profile of the ideal product and its
#' gap to the closest tested product.
#'
#' The perceptual space is built once, on all the consumers; with several
#' clusters, only the regressions and the outputs are computed cluster by
#' cluster. The plots are printed as they are built and returned in the list.
#'
#' @param X the sensory characterisation dataset: a matrix or data frame with
#'   the products in rows and the sensory descriptors in columns.
#' @param Y the hedonic scores dataset: a matrix or data frame with the same
#'   products, in the same order, in rows and the consumers in columns
#' @param method a character string defining the method to be used ("acp", "x2pls", "ox2pls", "qpls" and "plss")
#' @param regmod a value between 1 and 4 defining the regression model to be used : 1 for vectorial, 2 for circular, 3 for elliptic and 4 for quadratic
#' @param ngroupe the number of panelist cluster in the hedonic dataset with CLV_kmeans
#' @param var ogical. If \code{TRUE} (default), the descriptors are drawn as
#'   arrows on the map.
#' @param graph_carto logical. If \code{TRUE} (default), the preference map is
#'   drawn.
#' @param graph_profil_senso logical. If \code{TRUE} (default), the VIP of the
#'   descriptors, the sensory profile of the ideal product and its gap to the
#'   closest product are drawn. Not used with \code{method = "acp"}
#' @param ellipses logical. If \code{TRUE}, confidence ellipses of the
#'   products, computed from 1000 bootstrap samples of the consumers aligned by
#'   Procrustes rotation, are drawn on the map. Not used with
#'   \code{method = "acp"}. Default \code{FALSE}.
#'
#' @returns For \code{method = "acp"}, a list with:
#' \itemize{
#'   \item \code{acp}: the PCA, output of \code{FactoMineR::PCA};
#'   \item \code{coord_PI}: the coordinates of the ideal product, (Dim1, Dim2)
#'     for each cluster in turn;
#'   \item \code{coord_prod}: the coordinates of the products followed by the
#'     ends of the descriptor arrows, two columns per cluster.
#' }
#' For the other methods, a list with:
#' \itemize{
#'   \item \code{regpls}: the PLS2 model, or the output of
#'     \code{\link{QPLSniPond}} for \code{"qpls"};
#'   \item \code{alpha} (\code{"qpls"} only): the alpha selected for each of
#'     the two components (weight of the linear term in the QPLS criterion);
#'   \item \code{Bspline} (\code{"plss"} only): the B-spline basis;
#'   \item \code{profil_senso_PI}: the predicted sensory profile of the ideal
#'     product, one row per cluster, descriptors sorted by decreasing VIP;
#'   \item \code{coord_PI}: the coordinates of the ideal product, (Dim1, Dim2)
#'     for each cluster in turn;
#'   \item \code{delta}: the gap between the profile of the ideal product and
#'     the projected profile of the closest tested product, one row per
#'     cluster (positive: intensity to increase);
#'   \item \code{vip}: the squared VIP of the descriptors, one row per
#'     cluster; for \code{"x2pls"} and \code{"ox2pls"}, three rows per
#'     cluster: total, linear part and quadratic part;
#'   \item \code{Dim}: the coordinates of the products on the two components;
#'   \item \code{coord_prod}: the coordinates of the products (first rows),
#'     two columns per cluster; the following rows, when present, are the ends
#'     of the descriptor arrows;
#'   \item \code{graph}, \code{graph_vip}, \code{graph_profil_senso},
#'     \code{graph_delta}: lists with one element per cluster
#'     (\code{cluster1}, \code{cluster2}, ...): the preference map, the VIP
#'     barplot, the profile of the ideal product and its gap to the closest
#'     product (ggplot objects). An element is \code{NULL} when its plot is not
#'     drawn.
#' }
#'
#' @export PrefMap
#' @importFrom pls mvrValstats
#'
#' @examples
#'
#' data(beer)
#'
#' res_carto <- PrefMap(senso.beer, hedo.beer, method = "ox2pls", regmod = 1)
#'
PrefMap <- function(X, Y, method = c("ox2pls", "acp", "x2pls", "qpls", "plss"),
                  regmod = 2, ngroupe = 1, var = TRUE, graph_carto = TRUE,
                  graph_profil_senso = TRUE, ellipses = FALSE) {
  method <- match.arg(method)
  if (!regmod %in% 1:4) {
    stop("`regmod` doit valoir 1 (vectoriel), 2 (circulaire), 3 (elliptique) ou 4 (quadratique complet).",
         call. = FALSE)
  }

  # Structure du jeu de données
  dat <- prepare_data(X, Y)
  carac <- dat$carac
  hed <- dat$hed

  # Clustering
  km <- carto_clusters(hed, ngroupe)

  switch(method,
         acp = carto_acp(carac, hed, km, regmod = regmod, ngroupe = ngroupe, var = var,
                         graph_carto = graph_carto),
         x2pls = ,
         ox2pls = carto_x2pls(carac, hed, km, method = method, regmod = regmod, ngroupe = ngroupe,
                              var = var, graph_carto = graph_carto,
                              graph_profil_senso = graph_profil_senso, ellipses = ellipses),
         plss = carto_plss(carac, hed, km, regmod = regmod, ngroupe = ngroupe, var = var,
                           graph_carto = graph_carto, graph_profil_senso = graph_profil_senso,
                           ellipses = ellipses),
         qpls = carto_qpls(carac, hed, km, regmod = regmod, ngroupe = ngroupe, var = var,
                           graph_carto = graph_carto, graph_profil_senso = graph_profil_senso,
                           ellipses = ellipses)
  )
}
