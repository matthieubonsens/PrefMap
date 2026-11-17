#' PrefMap: External Preference Mapping for Sensory Data
#'
#' External preference mapping : relates the hedonic scores of consumers to the
#' sensory profiles of the products, in two steps: a two-dimensional perceptual
#' space is built from the sensory data, then the scores of each consumer are
#' regressed on it.
#'
#' @section Perceptual space:
#' Five methods are available through the \code{method} argument of
#' \code{\link{PrefMap}}:
#' \itemize{
#'   \item \code{"acp"}: Principal Component Analysis of the sensory profiles;
#'   \item \code{"qpls"}: Quadratic PLS regression without interactions
#'     (\code{\link{QPLSniPond}}; Cariou et al., 2014);
#'   \item \code{"plss"}: PLS regression on B-splines of degree 1 with one
#'     knot at the median of each descriptor (Durand, 2001);
#'   \item \code{"x2pls"}: PLS regression on the descriptors and their squares;
#'   \item \code{"ox2pls"}: PLS regression on the descriptors and their squares
#'     orthogonalised against them.
#' }
#'
#' @section Regression models:
#' The scores of each consumer are regressed on the two components with a
#' vectorial, circular, elliptic or full quadratic model (\code{regmod} = 1 for vectorial to
#' 4 for full quadratic in \code{\link{PrefMap}}).
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{PrefMap}}: for one method and one regression model, the
#'     preference map, the variable importance in the projection (VIP), the
#'     sensory profile of the ideal product and its gap to the closest tested
#'     product;
#'   \item \code{\link{reg_choice}}: the weighted adjusted R-squared of the 20
#'     method x model pairs;
#'   \item \code{\link{bootstrap}}: the bootstrap standard deviation of these
#'     R-squared, to assess the stability of each pair;
#'   \item \code{\link{QPLSniPond}}: the quadratic PLS regression itself.
#' }
#'
#' @references
#' Cariou, V., Verdun, S., & Qannari, E. M. (2014). Quadratic PLS regression
#' applied to external preference mapping. \emph{Food Quality and Preference},
#' 32, 28-34. \doi{10.1016/j.foodqual.2013.07.003}
#'
#' Danzart, M., Sieffermann, J.-M., & Delarue, J. (2004). New developments in
#' preference mapping techniques: finding out a consumer optimal product, its
#' sensory profile and the key sensory attributes. Communication at
#' Sensometrics 2004, Davis (CA).
#'
#' Durand, J.-F. (2001). Local polynomial additive regression through PLS and
#' splines: PLSS. \emph{Chemometrics and Intelligent Laboratory Systems},
#' 58(2), 235-246. \doi{10.1016/S0169-7439(01)00162-9}
#'
#' @seealso \code{\link{beer}}, an example dataset.
#' @export
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

## mockable bindings: start
## mockable bindings: end
NULL
