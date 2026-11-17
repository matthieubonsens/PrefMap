#' Preference surface of the vectorial OX2PLS model
#'
#' @param carac,regpls,Gtot,tens,centre,echelle,cols as in carto_predict_ox2pls
#'
#' @returns Z, the number of consumers whose predicted score (from
#'   carto_predict_ox2pls) exceeds their mean score, at each grid point
#' @noRd
carto_surface_ox2pls <- function(carac, regpls, Gtot, tens, centre, echelle, cols = TRUE) {
  pr <- carto_predict_ox2pls(carac, regpls, Gtot, tens, centre, echelle, cols)
  rowSums(pr > 0)
}
