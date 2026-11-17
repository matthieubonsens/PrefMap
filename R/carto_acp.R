#' PrefMap on PCA space
#'
#' @param carac the sensory characterisation dataset (products x descriptors)
#' @param hed the hedonic scores dataset (products x consumers)
#' @param km the clustering of the consumers, output of carto_clusters
#' @param regmod the regression model: 1 vectorial, 2 circular, 3 elliptic,
#'   4 full quadratic
#' @param ngroupe the number of clusters
#' @param var logical: draw the descriptors as arrows on the map
#' @param graph_carto logical: draw the map
#'
#' @returns A list that contain :
#'
#'\itemize{
#'   \item \code{acp}: the results of the PCA
#'   \item \code{coord_PI} : the optimal product coordinates
#'   \item \code{coord_product} : the product coordinates
#'}
#'
#'
#' @noRd
carto_acp <- function(carac, hed, km, regmod, ngroupe, var, graph_carto) {
  res.acp <- FactoMineR::PCA(carac, scale.unit = FALSE, ncp = Inf, graph = FALSE)

  coord_PI <- c()
  coord_product <- c()

  for (i in 1:ngroupe) {
    hed.gr <- hed[, km$clusters[2, ] == i]
    Dim1 <- res.acp$ind$coord[, 1]
    Dim2 <- res.acp$ind$coord[, 2]
    grille <- carto_grid(Dim1, Dim2)
    G1 <- grille$G1
    G2 <- grille$G2
    Gtot <- grille$Gtot

    Z <- carto_surface(hed.gr, Dim1, Dim2, Gtot, regmod)
    coord_opt <- carto_optimum(Z, Gtot)

    df.fleche <- carto_fleches(as.data.frame(res.acp$var$coord[, 1:2]), G1, G2, colnames(carac))
    resmod <- data.frame(G1 = Gtot[, 1], G2 = Gtot[, 2], Z = Z)
    coord_prod <- data.frame(cbind(res.acp$ind$coord[, 1], res.acp$ind$coord[, 2]))
    colnames(df.fleche) <- colnames(coord_prod)
    coord_prod <- data.frame(rbind(coord_prod, df.fleche))

    # Cartographie
    if (graph_carto) {
      g <- carto_plot_carte(resmod, coord_prod, df.fleche, coord_opt, n = nrow(carac), i = i,
                            var = var, style = "acp",
                            var_x = res.acp$eig[1:2, 2],
                            r2_hed = carto_r2_hed(hed.gr, Dim1, Dim2, regmod))
      print(g)
    }
    coord_PI <- c(coord_PI, coord_opt)
    coord_product <- data.frame(c(coord_product, coord_prod))
    row.names(coord_product) <- row.names(coord_prod)
  }

  list(acp = res.acp, coord_PI = coord_PI, coord_prod = coord_product)
}
