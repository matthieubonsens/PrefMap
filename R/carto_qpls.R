#' PrefMap on QPLS space
#'
#' @param carac the sensory characterisation dataset (products x descriptors)
#' @param hed the hedonic scores dataset (products x consumers)
#' @param km the clustering of the consumers, output of carto_clusters
#' @param method "x2pls" or "ox2pls"
#' @param regmod the regression model: 1 vectorial, 2 circular, 3 elliptic,
#'   4 full quadratic
#' @param ngroupe the number of clusters
#' @param var logical: draw the descriptors as arrows on the map
#' @param graph_carto logical: draw the map
#' @param graph_profil_senso logical: draw the VIP, the profile of the ideal
#'   product and its gap to the closest product
#' @param ellipses logical: draw the bootstrap confidence ellipses of the
#'   products
#'
#' @returns the list returned by PrefMap() for method = "qpls" (regpls, alpha,
#'   vip, profil_senso_PI, coord_PI, delta, Dim, coord_prod and the four lists
#'   of plots)
#'
#' @noRd
carto_qpls <- function(carac, hed, km, regmod, ngroupe, var, graph_carto,
                       graph_profil_senso, ellipses) {
  # Resultats
  profil_senso_PI <- array(0, dim = c(ngroupe, ncol(carac)),
                           dimnames = list(paste0("cluster", seq_len(ngroupe)), colnames(carac)))
  coord_PI <- c()
  res_delta <- array(0, dim = c(ngroupe, ncol(carac)),
                     dimnames = list(paste0("cluster", seq_len(ngroupe)), colnames(carac)))
  res_vip <- NULL
  coord_product <- c()
  noms_clusters <- paste0("cluster", seq_len(ngroupe))
  g_carte <- stats::setNames(vector("list", ngroupe), noms_clusters)
  g_vip_l <- stats::setNames(vector("list", ngroupe), noms_clusters)
  g_ps_l  <- stats::setNames(vector("list", ngroupe), noms_clusters)
  g_pp_l  <- stats::setNames(vector("list", ngroupe), noms_clusters)

  # Recherche d'alpha (QPLSniPond doit faire partie du package)
  qpls <- QPLSniPond(carac, hed, ncomp = 2)
  alpha <- qpls$aopt

  # Regression PLS
  qpls2 <- QPLSniPond(carac, hed, ncomp = 2, alpha = alpha)

  for (i in 1:ngroupe) {
    hed.gr <- hed[, km$clusters[2, ] == i]

    # Ellipses de confiance
    coord_prod <- data.frame(cbind(qpls2$X_score[, 1], qpls2$X_score[, 2]))
    row.names(coord_prod) <- row.names(carac)
    bornes <- list(xmin = 0, xmax = 0, ymin = 0, ymax = 0)
    ellips <- NULL
    if (ellipses) {
      ajuste <- function(hed.boot) {
        reg.boot <- QPLSniPond(carac, hed.boot, ncomp = 2, alpha = alpha)
        data.frame(cbind(reg.boot$X_score[, 1], reg.boot$X_score[, 2]))
      }
      res.ellips <- carto_ellipses(coord_prod, hed.gr, ajuste)
      ellips <- res.ellips$ellips
      bornes <- res.ellips[c("xmin", "xmax", "ymin", "ymax")]
    }

    # Construction du plan
    Dim1 <- qpls2$X_score[, 1]
    Dim2 <- qpls2$X_score[, 2]
    grille <- carto_grid(Dim1, Dim2, bornes)
    G1 <- grille$G1
    G2 <- grille$G2
    Gtot <- grille$Gtot

    Z <- carto_surface(hed.gr, Dim1, Dim2, Gtot, regmod)
    coord_opt <- carto_optimum(Z, Gtot)

    mod.load <- stats::lm(as.matrix(carac) ~ Dim1 + Dim2)
    load <- stats::coef(mod.load)
    load <- t(load[2:3, ])

    df.fleche <- carto_fleches(as.data.frame(load), G1, G2, colnames(carac))
    resmod <- data.frame(G1 = Gtot[, 1], G2 = Gtot[, 2], Z = Z)
    colnames(df.fleche) <- colnames(coord_prod)

    # Cartographie
    if (graph_carto) {
      g <- carto_plot_carte(resmod, coord_prod, df.fleche, coord_opt, n = nrow(carac), i = i,
                            var = var, ellips = ellips, style = "pls",
                            var_x = 100 * diff(c(0, qpls2$VarExpl[1:2, "X"])),
                            r2_hed = carto_r2_hed(hed.gr, Dim1, Dim2, regmod))
      if (var) {
        coord_prod <- data.frame(rbind(coord_prod, df.fleche))
      }
      print(g)
      g_carte[[i]] <- g
    }

    # VIP
    vipsq <- carto_vip_qpls(qpls2, ncol(carac), hed.gr)
    vipsq <- data.frame(vipsq)
    row.names(vipsq) <- colnames(carac)
    vipsq <- t(vipsq)

    # Profil sensoriel du point idéal
    carac_opt <- coord_opt %*% t(load)
    carac_opt <- data.frame(carac_opt)
    colnames(carac_opt) <- colnames(carac)

    note <- carto_note(hed.gr, Dim1, Dim2, coord_opt, regmod)

    carac_opt <- carto_recentre(carac_opt, carac, borne0 = TRUE)
    carac_opt <- as.matrix(carac_opt)
    carac_opt_sort <- carto_tri_vip(carac_opt, vipsq, colnames(carac_opt))

    # Comparaison avec le produit le plus proche
    occ_pp <- carto_plus_proche(coord_prod, coord_opt, nrow(carac))
    coord_pp <- coord_prod[occ_pp, ]
    coord_pp <- as.matrix(coord_pp)
    carac_pp <- coord_pp %*% t(load)
    carac_pp <- carac_pp[, 1:ncol(carac)]
    carac_pp <- data.frame(carac_pp)
    carac_pp <- t(carac_pp)

    carac_pp <- carto_recentre(carac_pp, carac, borne0 = TRUE)
    carac_pp_sort <- carto_tri_vip(carac_pp, vipsq, colnames(carac_opt))
    row.names(carac_pp_sort) <- row.names(coord_pp)
    delta <- carac_opt_sort - carac_pp_sort

    # Graphe VIP et profil sensoriel
    if (graph_profil_senso) {
      g_vip <- carto_plot_vip(vipsq)
      print(g_vip)
      g_vip_l[[i]] <- g_vip

      g_ps <- carto_plot_ideal(carac_opt_sort, i, sous_titre = paste("Note =", note))
      g_pp <- carto_plot_delta(delta, i, sous_titre = paste("Produit :", row.names(carac_pp_sort)))
      gridExtra::grid.arrange(g_ps, g_pp)
      g_ps_l[[i]] <- g_ps
      g_pp_l[[i]] <- g_pp
    }

    profil_senso_PI[i, ] <- carac_opt

    coord_PI <- c(coord_PI, coord_opt)

    res_delta[i, ] <- carac_opt - carac_pp

    res_vip <- rbind(res_vip, vipsq)
    colnames(res_vip) <- colnames(vipsq)

    coord_product <- data.frame(c(coord_product, coord_prod))
    row.names(coord_product) <- row.names(coord_prod)
  }

  list(regpls = qpls2, alpha = alpha, vip = res_vip, profil_senso_PI = profil_senso_PI,
       coord_PI = coord_PI, delta = res_delta, Dim = data.frame(Dim1, Dim2),
       coord_prod = coord_product, graph = g_carte, graph_vip = g_vip_l,
       graph_profil_senso = g_ps_l, graph_delta = g_pp_l)
}
