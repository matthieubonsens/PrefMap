#' Bootstrap confidence ellipses of the products
#'
#' Resamples the consumers of the cluster nboot times with replacement,
#' rebuilds the space on each sample with ajuste() and aligns it on the
#' reference configuration by a Procrustes rotation (vegan::procrustes,
#' without scaling).
#'
#' @param coord_prod the n x 2 coordinates of the products on the map
#' @param hed.gr the hedonic scores of the consumers of the cluster
#' @param ajuste a function that takes a resampled hedonic dataset and returns
#'   the n x 2 coordinates of the products in the space rebuilt on it
#' @param nboot the number of bootstrap samples
#'
#' @returns a list: ellips, the nboot x 2 x n array of the aligned coordinates
#'   (sample, axis, product), and xmin, xmax, ymin, ymax, the range of these
#'   coordinates (0 included), used to extend the grid
#' @noRd
carto_ellipses <- function(coord_prod, hed.gr, ajuste, nboot = 1000) {
  mat.ellips <- array(0, dim = c(nrow(coord_prod), ncol(coord_prod), nboot))

  sample.boot <- matrix(0, nrow = nboot, ncol = ncol(hed.gr))
  pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)
  on.exit(close(pb))

  for (iteration in 1:nboot) {
    sample.boot[iteration, ] <- sample(1:ncol(hed.gr), ncol(hed.gr), replace = TRUE)
  }
  for (iteration in 1:nboot) {
    hed.boot <- hed.gr[, sample.boot[iteration, ]]
    coord_prod_boot <- ajuste(hed.boot)
    coord_prod_boot <- vegan::procrustes(coord_prod, coord_prod_boot, scale = FALSE)$Yrot
    mat.ellips[, , iteration] <- as.matrix(coord_prod_boot)
    utils::setTxtProgressBar(pb, iteration)
  }

  # ellips[iteration, axe, produit]
  ellips <- aperm(mat.ellips, c(3, 2, 1))

  list(ellips = ellips,
       xmin = min(0, ellips[, 1, ]), xmax = max(0, ellips[, 1, ]),
       ymin = min(0, ellips[, 2, ]), ymax = max(0, ellips[, 2, ]))
}
