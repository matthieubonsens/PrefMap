#' Function to program the arrows for sensory descriptors on the map
#'
#' @param df.fleche the descriptors postions
#' @param G1 the first component of the map
#' @param G2 the second component of the map
#' @param noms_desc the names of the descriptors
#'
#' @returns the real position of the descriptors on the map
#'
#' @noRd
carto_fleches <- function(df.fleche, G1, G2, noms_desc) {
  rownames(df.fleche) <- noms_desc
  max.norm.prod <- min(max(G1), abs(min(G1)), max(G2), abs(min(G2)))
  max.norm.desc <- max(abs(df.fleche))
  expand <- max.norm.prod / max.norm.desc
  df.fleche * expand
}
