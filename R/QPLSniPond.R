#' Quadratic PLS regression without interactions
#'
#' Fits a quadratic PLS2 model (QPLS2) without interaction terms. Each
#' component t = Xw maximises a criterion combining a linear term (covariance
#' with t) and a quadratic term (covariance with t^2), weighted by a parameter
#' alpha chosen on a grid to minimise the residual variance of Y after
#' regression on (1, t, t^2).
#'
#' @param X the sensory characterisation dataset (matrix or data frame,
#'   products in rows). Centred internally, not scaled.
#' @param Y the hedonic scores dataset (matrix, data frame or vector), with the
#'   same products in the same order as X. Centred internally, not scaled.
#' @param ncomp the number of components to extract.
#' @param eps the convergence tolerance of the iterative loops on u and w.
#' @param def how Y is deflated after each component: '1' (default) for a
#'   linear deflation (regression on t), any other value (e.g. '2') for a
#'   quadratic deflation (regression on t and t^2).
#' @param critere the criterion maximised to build t: '1' (default) summarises
#'   Y by a single latent variable u; any other value (e.g. '2') uses the Y
#'   variables directly, without a latent variable.
#' @param alpha the grid of candidate values, between 0 and 1, for the weight of the
#'   linear term in the criterion (1: linear term only, 0: quadratic term
#'   only). The best value is selected for each component (see \code{aopt}).
#' @param maxit.u the maximum number of iterations of the loop on u (only used
#'   when \code{critere = '1'}).
#' @param maxit.w the maximum number of iterations of the loop on w.
#'
#' @returns A list containing:
#' \itemize{
#'   \item \code{X_score}: the X scores t (n x ncomp).
#'   \item \code{X_weight}: the X weights w (p x ncomp), of norm 1, applied to
#'     the deflated X.
#'   \item \code{X_weightstar}: the X weights w* (p x ncomp), applied to the
#'     centred original X; use them to project new samples.
#'   \item \code{RegVar}: the regressors of the quadratic model, t1, t1^2, t2,
#'     t2^2, ... (n x 2ncomp), e.g. for \code{lm(Y ~ RegVar)}.
#'   \item \code{VarExpl}: the cumulative proportion of variance of X and Y
#'     explained by the first k components (ncomp x 2). For Y, only the linear
#'     part is counted when \code{def = '1'}.
#'   \item \code{method}: the name of the method ("Quadratic PLS2").
#'   \item \code{ncomp}: the number of components.
#'   \item \code{Xcol}: the number of columns of X.
#'   \item \code{Ycol}: the number of columns of Y.
#'   \item \code{nbInd}: the number of products (rows).
#'   \item \code{critere}: the value of the argument \code{critere}.
#'   \item \code{deflation}: the value of the argument \code{def}.
#'   \item \code{alpha}: a matrix whose first column is the alpha grid and
#'     column k + 1 the residual variance of Y obtained with each alpha value
#'     for component k.
#'   \item \code{aopt}: the alpha value selected for each component.
#'   \item \code{RegCoeff}: an array (ncomp x 2 x q) of the standardised
#'     coefficients of t and t^2 for each Y variable; \code{NA} when
#'     \code{def = '1'}.
#' }
#' When \code{critere = '1'}, the list also contains:
#' \itemize{
#'   \item \code{Y_score}: the Y scores u (n x ncomp).
#'   \item \code{Y_weight}: the Y weights c (q x ncomp), of norm 1.
#'   \item \code{Ydef}, \code{Xdef}: lists of the deflated Y and X matrices
#'     after each component.
#' }
#' When \code{critere = '2'}, the list also contains:
#' \itemize{
#'   \item \code{Y}: the Y matrix deflated by all the components (final
#'     residuals).
#' }
#'
#' @export QPLSniPond
#' @importFrom stats sd
#'
#' @references
#' Cariou, V., Verdun, S., & Qannari, E. M. (2014). Quadratic PLS regression
#' applied to external preference mapping. \emph{Food Quality and Preference},
#' 32, 28-34. \doi{10.1016/j.foodqual.2013.07.003}
#'
#' @examples
#'
#' data(beer)
#'
#' resqpls <- QPLSniPond(senso.beer, hedo.beer, ncomp = 2)
#'
#' resqpls$aopt
#' resqpls$VarExpl
#'
QPLSniPond <- function(X, Y, ncomp, eps = 0.00001, def = '1', critere = '1',
                       alpha = seq(0, 1, by = 0.05),
                       maxit.u = 150, maxit.w = 1000) {

  X <- scale(as.matrix(X), scale = FALSE)
  if (is.vector(Y)) Y <- matrix(Y, length(Y), 1)
  Y <- scale(as.matrix(Y), scale = FALSE)

  n  <- nrow(X); p <- ncol(X); q <- ncol(Y)
  na <- length(alpha)
  crit1 <- (critere == 1)

  ## Stockage
  Tsc <- matrix(0, n, ncomp)                       # scores t
  W   <- matrix(0, p, ncomp)                       # poids w
  R   <- matrix(0, p, ncomp)                       # poids w*
  P   <- matrix(0, p, ncomp)                       # loadings X (pour w*)
  U   <- matrix(0, n, ncomp)                       # scores u   (critere 1)
  C   <- matrix(0, q, ncomp)                       # poids de Y (critere 1)
  H   <- matrix(0, n, 2 * ncomp)
  colnames(H) <- as.vector(rbind(paste0('t', 1:ncomp),
                                 paste0('t', 1:ncomp, '.', 1:ncomp)))
  CoefReg <- array(dim = c(ncomp, 2, q),
                   dimnames = list(paste('Comp.', 1:ncomp), c('t', 't^2'),
                                   paste('Y-var.', 1:q)))
  VarYalpha <- matrix(0, na, ncomp + 1)
  VarYalpha[, 1] <- alpha
  VarResid <- matrix(0, ncomp, 2)
  VarTot   <- c(.qpls_totvar(X), .qpls_totvar(Y))
  aopt <- numeric(ncomp)
  Xdef <- Ydef <- vector('list', ncomp)

  ## Tolerance de rang, relative a l'echelle de X de depart
  tol <- max(n, p) * .Machine$double.eps * svd(X, nu = 0, nv = 0)$d[1]

  for (k in seq_len(ncomp)) {

    sp <- .qpls_reduce(X, tol)
    if (sp$r == 0)
      stop("X est de rang insuffisant pour extraire la composante ", k)
    nul <- sp$r < p   # G a alors des valeurs propres nulles hors de l'espace des lignes

    if (crit1) {
      ## ---- Critere 1 : recherche du alpha optimal ----------------------------
      u0 <- .qpls_init_u(X, Y)          # premiere composante PLS2 (ex-mvr)
      for (i in seq_len(na)) {
        a  <- alpha[i]
        cp <- .qpls_comp1(sp, Y, u0, a, c(a, 1 - a), nul, eps, maxit.u, maxit.w)
        VarYalpha[i, k + 1] <- .qpls_totvar(qr.resid(qr(cbind(1, cp$t, cp$t^2)), Y))
      }
      ## ---- Composante pour le alpha optimal (comme l'original : u initialise
      ## ---- par Y[,1] et poids sqrt(a), sqrt(1-a) pour le pas sur u) ----------
      a  <- alpha[which.min(VarYalpha[, k + 1])]
      cp <- .qpls_comp1(sp, Y, Y[, 1], a, sqrt(c(a, 1 - a)), nul, eps,
                        maxit.u, maxit.w)
      if (cp$iter >= maxit.u || cp$niter >= maxit.w)
        warning("Nombre d'iterations maximal atteint. Axe ", k)
      U[, k] <- cp$u
      C[, k] <- cp$c

    } else {
      ## ---- Critere 2 : B, G, wb, wg ne dependent pas d'alpha -----------------
      M  <- crossprod(sp$Xr, Y)                      # V'X'Y (r x q)
      Bz <- tcrossprod(M)                            # V'X'YY'XV
      Gz <- crossprod(sp$Xr, sp$Xr * rowSums(Y))     # V'X' diag(sum_j y_j) XV
      wb <- eigen(Bz, symmetric = TRUE)$vectors[, 1]
      wg <- .qpls_wg(Gz, nul)
      res <- vector('list', na)
      for (i in seq_len(na)) {
        ws <- .qpls_wstep(Gz, wb, wg, alpha[i], sp$V, eps, maxit.w, Bz = Bz)
        ws$t <- drop(sp$Xr %*% ws$z)
        VarYalpha[i, k + 1] <- .qpls_totvar(qr.resid(qr(cbind(1, ws$t, ws$t^2)), Y))
        res[[i]] <- ws
      }
      ## La composante du alpha optimal est exactement celle de la recherche
      iopt <- which.min(VarYalpha[, k + 1])
      a  <- alpha[iopt]
      cp <- res[[iopt]]
      if (cp$niter >= maxit.w)
        warning("Nombre d'iterations maximal atteint. Axe ", k)
    }
    aopt[k] <- a

    ## ---- Stockage et deflation ---------------------------------------------
    w  <- cp$w
    tt <- cp$t
    W[, k]   <- w
    Tsc[, k] <- tt
    H[, 2 * k - 1] <- tt
    H[, 2 * k]     <- tt^2

    tsq  <- sum(tt^2)
    proj <- drop(crossprod(X, tt)) / tsq
    X <- X - tt %o% proj

    if (def == '1') {
      Y <- Y - tt %o% (drop(crossprod(Y, tt)) / tsq)
    } else {
      qrH <- qr(cbind(1, tt, tt^2))
      cf  <- qr.coef(qrH, Y)
      Y[] <- qr.resid(qrH, Y)
      CoefReg[k, 1, ] <- cf[2, ] * sd(tt)
      CoefReg[k, 2, ] <- cf[3, ] * sd(tt^2)
    }
    Xdef[[k]] <- X
    Ydef[[k]] <- Y

    VarResid[k, ] <- c(.qpls_totvar(X), .qpls_totvar(Y))

    ## w* = w - sum_{j<k} w*_j (p_j' w)   (equivalent a la recurrence sur E)
    R[, k] <- if (k == 1) w else
      w - R[, 1:(k - 1), drop = FALSE] %*% crossprod(P[, 1:(k - 1), drop = FALSE], w)
    P[, k] <- proj
  }

  VarExpl <- 1 - VarResid / rep(VarTot, each = ncomp)
  colnames(VarExpl) <- c('X', 'Y')

  if (crit1) {
    list(X_score = Tsc, Y_score = U, X_weight = W, X_weightstar = R,
         Y_weight = C, RegVar = H, VarExpl = VarExpl,
         method = 'Quadratic PLS2', ncomp = ncomp, Xcol = p, Ycol = q,
         nbInd = n, critere = critere, deflation = def,
         alpha = VarYalpha, aopt = aopt, RegCoeff = CoefReg,
         Ydef = Ydef, Xdef = Xdef)
  } else {
    list(X_score = Tsc, X_weight = W, X_weightstar = R,
         RegVar = H, VarExpl = VarExpl,
         method = 'Quadratic PLS2', ncomp = ncomp, Xcol = p, Ycol = q,
         nbInd = n, critere = critere, deflation = def,
         alpha = VarYalpha, aopt = aopt, RegCoeff = CoefReg, Y = Y)
  }
}


########
#' @noRd

.qpls_totvar <- function(M) {
  M <- M - rep(colMeans(M), each = nrow(M))
  sum(M^2) / (nrow(M) - 1)
}

.qpls_reduce <- function(X, tol) {
  s <- svd(X)
  r <- sum(s$d > tol)
  idx <- seq_len(r)
  list(V  = s$v[, idx, drop = FALSE],
       Xr = s$u[, idx, drop = FALSE] * rep(s$d[idx], each = nrow(X)),
       r  = r)
}

.qpls_init_u <- function(X, Y) {
  X <- X - rep(colMeans(X), each = nrow(X))
  Y <- Y - rep(colMeans(Y), each = nrow(Y))
  XtY <- crossprod(X, Y)
  if (ncol(Y) == 1) {
    w <- XtY / sqrt(c(crossprod(XtY)))
  } else if (ncol(Y) < ncol(X)) {
    qv <- eigen(crossprod(XtY), symmetric = TRUE)$vectors[, 1]
    w  <- XtY %*% qv
    w  <- w / sqrt(c(crossprod(w)))
  } else {
    w <- eigen(XtY %*% t(XtY), symmetric = TRUE)$vectors[, 1]
  }
  drop(X %*% w)
}

.qpls_wg <- function(Gz, nul) {
  e <- eigen(Gz, symmetric = TRUE)
  if (nul && e$values[1] < 0) NULL else e$vectors[, 1]
}

.qpls_wstep <- function(Gz, wb, wg, a, V, eps, maxit, vz = NULL, Bz = NULL) {
  Bmul <- if (is.null(Bz)) function(z) vz * sum(vz * z) else function(z) drop(Bz %*% z)
  crit <- function(z) a * sum(z * Bmul(z)) + (1 - a) * sum(z * (Gz %*% z))^2

  z <- if (is.null(wg) || crit(wb) > crit(wg)) wb else wg
  w <- drop(V %*% z)
  zold <- z
  wold <- rep(100000, length(w))
  niter <- 0
  ## Le test de convergence reste celui de l'original (norme L1 sur w dans R^p)
  while (sum(abs(wold - w)) > eps && niter < maxit) {
    niter <- niter + 1
    zold <- z
    wold <- w
    g   <- drop(Gz %*% z)
    num <- a * Bmul(z) + (1 - a) * sum(z * g) * g
    z   <- num / sqrt(sum(num^2))
    w   <- drop(V %*% z)
  }
  if (niter >= maxit && crit(zold) > crit(z)) { z <- zold; w <- wold }

  ## Convention de signe : plus grand |w| positif
  if (w[which.max(abs(w))] < 0) { z <- -z; w <- -w }
  list(z = z, w = w, niter = niter)
}

.qpls_topvec2 <- function(c1, c2) {
  a11 <- sum(c1^2); a22 <- sum(c2^2); a12 <- sum(c1 * c2)
  lam <- (a11 + a22) / 2 + sqrt(((a11 - a22) / 2)^2 + a12^2)
  y1 <- c(a12, lam - a11)
  y2 <- c(lam - a22, a12)
  y  <- if (sum(y1^2) >= sum(y2^2)) y1 else y2
  v  <- c1 * y[1] + c2 * y[2]
  nv <- sqrt(sum(v^2))
  if (!is.finite(nv) || nv == 0) NULL else v / nv
}

.qpls_comp1 <- function(sp, Y, u, a, wy, nul, eps, maxit.u, maxit.w) {
  Xr <- sp$Xr
  uold <- rep(1, length(u))
  iter <- 0
  cc <- NULL
  while (sum((u - uold)^2) > eps && iter < maxit.u) {
    iter <- iter + 1
    uold <- u
    vz <- drop(crossprod(Xr, u))                 # V'X'u
    Gz <- crossprod(Xr, Xr * u)                  # V'X' diag(u) X V
    ws <- .qpls_wstep(Gz, vz / sqrt(sum(vz^2)), .qpls_wg(Gz, nul), a,
                      sp$V, eps, maxit.w, vz = vz)
    tt <- drop(Xr %*% ws$z)                      # t = X w
    c1 <- wy[1] * drop(crossprod(Y, tt))
    c2 <- wy[2] * drop(crossprod(Y, tt^2))
    cn <- if (is.null(cc)) NULL else .qpls_topvec2(c1, c2)
    if (is.null(cn)) {
      cn <- eigen(c1 %o% c1 + c2 %o% c2, symmetric = TRUE)$vectors[, 1]
    } else if (sum(cn * cc) < 0) {
      cn <- -cn
    }
    cc <- cn
    u  <- drop(Y %*% cc)
  }
  list(w = ws$w, t = tt, u = u, c = cc, iter = iter, niter = ws$niter)
}

