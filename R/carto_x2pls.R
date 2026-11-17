#' PrefMap on the X2PLS or OX2PLS space
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
#' @returns the list returned by PrefMap() for these methods (regpls,
#'   profil_senso_PI, coord_PI, delta, vip, Dim, coord_prod and the four lists
#'   of plots)
#' @noRd
carto_x2pls <- function(carac, hed, km, method, regmod, ngroupe, var, graph_carto,
                        graph_profil_senso, ellipses) {
  ext <- carto_extension(carac, method)
  carac2 <- ext$carac2
  tens <- ext$tens
  if (method == "x2pls") {
    lin  <- seq_len(ncol(carac))
    quad <- ncol(carac) + seq_len(ncol(carac))
  } else {
    lin  <- seq(1, ncol(carac2), by = 2)
    quad <- seq(2, ncol(carac2), by = 2)
  }

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

  # Regression PLS
  regpls <- pls::plsr(as.matrix(hed) ~ as.matrix(carac2), method = "kernelpls", scale = TRUE,
                      validation = "none", maxit = 1e12)

  for (i in 1:ngroupe) {
    hed.gr <- hed[, km$clusters[2, ] == i]

    # Ellipses de confiance
    coord_prod <- data.frame(cbind(regpls$scores[, 1], regpls$scores[, 2]))
    bornes <- list(xmin = 0, xmax = 0, ymin = 0, ymax = 0)
    ellips <- NULL
    if (ellipses) {
      ajuste <- function(hed.boot) {
        reg.boot <- pls::plsr(as.matrix(hed.boot) ~ as.matrix(carac2), method = "kernelpls",
                              scale = TRUE, validation = "none", maxit = 1e12)
        data.frame(cbind(reg.boot$scores[, 1], reg.boot$scores[, 2]))
      }
      res.ellips <- carto_ellipses(coord_prod, hed.gr, ajuste)
      ellips <- res.ellips$ellips
      bornes <- res.ellips[c("xmin", "xmax", "ymin", "ymax")]
    }

    # Construction du plan
    Dim1 <- regpls$scores[, 1]
    Dim2 <- regpls$scores[, 2]
    grille <- carto_grid(Dim1, Dim2, bornes)
    G1 <- grille$G1
    G2 <- grille$G2
    Gtot <- grille$Gtot

    if (method == "ox2pls" && regmod == 1) {
      Z <- carto_surface_ox2pls(carac, regpls, Gtot, tens, ext$centre, ext$echelle, cols = km$clusters[2,] == i)
    } else {
      Z <- carto_surface(hed.gr, Dim1, Dim2, Gtot, regmod)
    }
    coord_opt <- carto_optimum(Z, Gtot)

    if (method == "x2pls") {
      df.fleche <- as.data.frame(regpls$loadings[1:ncol(carac), 1:2])
    }
    if (method == "ox2pls") {
      df.fleche <- as.data.frame(regpls$loadings[seq(1, ncol(carac2), by = 2), 1:2])
    }
    df.fleche <- carto_fleches(df.fleche, G1, G2, colnames(carac))
    resmod <- data.frame(G1 = Gtot[, 1], G2 = Gtot[, 2], Z = Z)
    colnames(df.fleche) <- colnames(coord_prod)

    # Cartographie
    if (graph_carto) {
      g <- carto_plot_carte(resmod, coord_prod, df.fleche, coord_opt, n = nrow(carac), i = i,
                            var = var, ellips = ellips, style = "pls",
                            var_x = pls::explvar(regpls)[1:2],
                            r2_hed = carto_r2_hed(hed.gr, Dim1, Dim2, regmod))
      if (var) {
        coord_prod <- data.frame(rbind(coord_prod, df.fleche))
      }
      print(g)
      g_carte[[i]] <- g
    }

    # VIP
    vipsq <- carto_vip_pls(regpls, hed.gr, ncol(carac2),cols = km$clusters[2, ] == i)
    vipsq.lin <- vipsq[lin]
    vipsq.quad <- vipsq[quad]
    vipsq <- vipsq.lin + vipsq.quad
    vipsq <- data.frame(vipsq)
    row.names(vipsq) <- colnames(carac)
    vipsq <- cbind(vipsq, vipsq.lin, vipsq.quad)
    vipsq <- t(vipsq)

    descripteurs <- c()
    type <- c()
    vip <- c()
    for (j in 1:ncol(vipsq)) {
      descripteurs <- cbind(descripteurs, colnames(vipsq)[j], colnames(vipsq)[j])
      type <- cbind(type, "Lin\u00e9aire", "Quadratique")
      vip <- cbind(vip, vipsq[2, j], vipsq[3, j])
    }
    df.vip <- rbind(descripteurs, type, vip)
    row.names(df.vip) <- c("Descripteurs", "Type", "VIP")
    df.vip <- t(df.vip)
    df.vip <- data.frame(df.vip)
    df.vip$VIP <- as.numeric(df.vip$VIP)
    df.vip$Descripteurs <- as.factor(df.vip$Descripteurs)
    df.vip$Type <- factor(df.vip$Type, levels = c("Quadratique", "Lin\u00e9aire"))

    # Profil sensoriel du point ideal
    mod.load <- stats::lm(as.matrix(carac) ~ Dim1 + Dim2)
    load <- t(stats::coef(mod.load)[2:3, ])
    carac_opt <- coord_opt %*% t(load)
    colnames(carac_opt) <- colnames(carac)

    if (method == "ox2pls" && regmod == 1) {
      ideal <- data.frame(Dim1 = coord_opt[1], Dim2 = coord_opt[2])
      pr <- carto_predict_ox2pls(carac, regpls, ideal, tens, ext$centre, ext$echelle,
                                 cols = km$clusters[2, ] == i)
      note <- round(mean(pr + regpls$Ymeans[km$clusters[2, ] == i]), digits = 3)
    } else {
      note <- round(carto_note(hed.gr, Dim1, Dim2, coord_opt, regmod), digits = 3)
    }

    carac_opt <- carto_recentre(carac_opt, carac, borne0 = TRUE)
    carac_opt_sort <- carto_tri_vip(carac_opt, vipsq, colnames(carac_opt))

    # Comparaison avec le produit le plus proche
    occ_pp <- carto_plus_proche(coord_prod, coord_opt, nrow(carac))
    coord_pp <- coord_prod[occ_pp, ]
    coord_pp <- as.matrix(coord_pp)
    carac_pp <- coord_pp %*% t(load)

    carac_pp <- carto_recentre(carac_pp, carac, borne0 = TRUE)
    carac_pp_sort <- carto_tri_vip(carac_pp, vipsq, colnames(carac_opt))
    row.names(carac_pp_sort) <- row.names(coord_pp)
    delta <- carac_opt_sort - carac_pp_sort

    # Graphes profil sensoriel et VIP
    if (graph_profil_senso) {
      g_vip <- carto_plot_vip_lq(df.vip)
      print(g_vip)
      g_vip_l[[i]] <- g_vip

      g_ps <- carto_plot_ideal(carac_opt_sort, i, sous_titre = paste("Note =", note), taille_titre = 10)
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

  list(regpls = regpls, profil_senso_PI = profil_senso_PI, coord_PI = coord_PI,
       delta = res_delta, vip = res_vip, Dim = data.frame(Dim1, Dim2), coord_prod = coord_product,
       graph = g_carte, graph_vip = g_vip_l,
       graph_profil_senso = g_ps_l, graph_delta = g_pp_l)
}
