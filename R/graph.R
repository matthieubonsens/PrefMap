#' Preference map of a cluster
#'
#' @param resmod a data frame with the grid points (G1, G2) and the preference
#'   surface (Z)
#' @param coord_prod the coordinates of the products (first n rows), followed
#'   for the PCA by the ends of the descriptor arrows
#' @param df.fleche the ends of the descriptor arrows
#' @param coord_opt the coordinates of the ideal product
#' @param n the number of products
#' @param i the cluster number
#' @param var logical. If TRUE, the descriptors are drawn as arrows
#' @param ellips NULL, or the nboot x 2 x n array returned by carto_ellipses,
#'   to draw one ellipse per product
#' @param style "acp" for the PCA, "pls" for the other methods
#' @param var_x the percentages of variance of X explained by Dim1 and Dim2
#'   (NULL: axes without percentages)
#' @param r2_hed the R-squared of the hedonic scores returned by carto_r2_hed
#'   (NULL: no subtitle)
#'
#' @returns a ggplot object
#' @importFrom ggplot2 .data
#' @noRd
carto_plot_carte <- function(resmod, coord_prod, df.fleche, coord_opt, n, i, var,
                             ellips = NULL, style = c("pls", "acp"), var_x = NULL,
                             r2_hed = NULL) {
  style <- match.arg(style)

  g <- ggplot2::ggplot(data = resmod, ggplot2::aes(x = .data$G1, y = .data$G2, z = .data$Z))
  g <- g + ggplot2::geom_contour_filled(show.legend = FALSE, bins = 8)
  g <- g + ggplot2::scale_fill_brewer(palette = "RdBu", direction = -1)
  g <- g + ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed", col = "black")
  g <- g + ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed", col = "black")
  g <- g + ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(min(resmod$G1), max(resmod$G1)))
  g <- g + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(resmod$G2), max(resmod$G2)))
  g <- g + ggplot2::geom_point(ggplot2::aes(x = coord_prod[1:n, 1], y = coord_prod[1:n, 2], z = NULL),
                               data = data.frame(coord_prod[1:n, ]), color = "black")
  g <- g + ggplot2::annotate("point", x = coord_opt[1], y = coord_opt[2],
                             shape = 8, colour = "yellow")

  pct <- function(v) paste0(sub(".", ",", sprintf("%.1f", v), fixed = TRUE), " %")
  axes <- c("Dim1", "Dim2")
  if (!is.null(var_x)) axes <- paste0(axes, " (", pct(var_x[1:2]), ")")
  sous_titre <- NULL
  if (!is.null(r2_hed)) {
    sous_titre <- paste0("R\u00b2 of the hedonic scores (", names(r2_hed), " model) : ", pct(r2_hed))
  }
  titre <- ggplot2::labs(title = paste("Mapping cluster ", i), subtitle = sous_titre,
                         x = axes[1], y = axes[2])
  theme_carte <- ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
                                plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
                                axis.title = ggplot2::element_text(size = 15))

  if (style == "acp") {
    g <- g + carto_layer_labels(coord_prod)
    g <- g + titre + theme_carte
    g <- g + ggplot2::coord_fixed()
    if (var) {
      g <- g + carto_layer_fleches(df.fleche)
    }
    return(g)
  }

  g <- g + titre + theme_carte
  if (var) {
    coord_lab <- data.frame(rbind(coord_prod, df.fleche))
    g <- g + carto_layer_fleches(df.fleche)
    g <- g + carto_layer_labels(coord_lab)
  } else {
    g <- g + carto_layer_labels(coord_prod)
  }
  if (!is.null(ellips)) {
    for (prod in seq_len(n)) {
      g <- g + ggplot2::stat_ellipse(data = data.frame(ellips[, , prod]),
                                     mapping = ggplot2::aes(x = .data$X1, y = .data$X2, z = NULL),
                                     type = "norm", linetype = 5)
    }
  }
  g + ggplot2::coord_fixed()
}

#' Flèches des descripteurs
#' @noRd
carto_layer_fleches <- function(df.fleche) {
  ggplot2::geom_segment(data = df.fleche,
                        ggplot2::aes(x = 0, y = 0, z = NULL, xend = df.fleche[, 1], yend = df.fleche[, 2]),
                        arrow = grid::arrow(length = grid::unit(0.1, "cm"), type = "closed"),
                        colour = "steelblue4", linewidth = 0.5)
}

#' Labels : Produits, descripteurs, ...
#' @noRd
carto_layer_labels <- function(coord) {
  ggrepel::geom_label_repel(data.frame(coord),
                            mapping = ggplot2::aes(x = coord[, 1], y = coord[, 2], z = NULL,
                                                   label = rownames(coord)),
                            label.size = NA, size = 3.5, segment.size = 0.3,
                            segment.colour = "grey30", label.padding = 0,
                            min.segment.length = 1, max.overlaps = Inf)
}

#' VIP linéaire / quadratique
#' @noRd
carto_plot_vip_lq <- function(df.vip) {
  ggplot2::ggplot(df.vip, ggplot2::aes(fill = .data$Type, x = .data$Descripteurs, y = .data$VIP)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::labs(title = "VIP") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(size = 20),
                   axis.title = ggplot2::element_text(size = 20),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 15),
                   legend.text = ggplot2::element_text(size = 20))
}

#' VIP PLS
#' @noRd
carto_plot_vip <- function(vip) {
  ggplot2::ggplot(data = data.frame(vip[1, ]), ggplot2::aes(x = colnames(vip), y = vip[1, ])) +
    ggplot2::geom_col(color = "black", fill = "deepskyblue") +
    ggplot2::labs(title = "VIP by descriptors", x = "Descriptors", y = "VIP") +
    ggplot2::geom_text(label = colnames(vip), vjust = 0, hjust = "inward", color = "black",
                       size = 4, angle = 90) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(size = 15),
                   axis.title = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_blank())
}

#' Profil sensoriel du point idéal
#' @noRd
carto_plot_ideal <- function(carac_opt_sort, i, sous_titre = ggplot2::waiver(), taille_titre = 15) {
  ggplot2::ggplot(data.frame(carac_opt_sort[1, ]),
                  ggplot2::aes(x = factor(t(colnames(carac_opt_sort)), levels = t(colnames(carac_opt_sort))),
                               y = as.matrix(carac_opt_sort[1, ]))) +
    ggplot2::geom_col(color = "black", fill = "darkorange") +
    ggplot2::labs(title = paste("Score by descriptor for the ideal product for cluster", i),
                  subtitle = sous_titre, x = "Descriptors", y = "Scores") +
    ggplot2::geom_text(label = t(colnames(carac_opt_sort)), vjust = 0, hjust = "inward",
                       color = "black", size = 4, angle = 90) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = taille_titre, face = "bold", hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_blank())
}

#' Écart entre le point idéal et le produit le plus proche
#' @noRd
carto_plot_delta <- function(delta, i, sous_titre = ggplot2::waiver()) {
  ggplot2::ggplot(data.frame(delta[1, ]),
                  ggplot2::aes(x = factor(t(colnames(delta)), levels = t(colnames(delta))),
                               y = as.matrix(delta[1, ]))) +
    ggplot2::geom_col(color = "black", fill = "darkorange") +
    ggplot2::labs(title = paste("Scores differences between the ideal product and the nearest product for cluster", i),
                  subtitle = sous_titre, x = "Descriptors", y = "Scores differences") +
    ggplot2::geom_text(label = t(colnames(delta)), vjust = 0, hjust = "inward",
                       color = "black", size = 4, angle = 90) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_blank())
}

#' R2 de la regression du Liking
#' @noRd

carto_r2_hed <- function(hed.gr, Dim1, Dim2, regmod) {
  M <- switch(regmod,
              cbind(Dim1, Dim2),
              cbind(Dim1, Dim2, Dim1^2 + Dim2^2),
              cbind(Dim1, Dim2, Dim1^2, Dim2^2),
              cbind(Dim1, Dim2, Dim1^2, Dim2^2, Dim1 * Dim2))
  Y <- as.matrix(hed.gr)
  ssqr <- sum(qr.resid(qr(cbind(1, M)), Y)^2)
  ssqt <- sum(scale(Y, center = TRUE, scale = FALSE)^2)
  modeles <- c("vectorial", "circular", "elliptic", "quadratic")
  stats::setNames(100 * (1 - ssqr / ssqt), modeles[regmod])
}





