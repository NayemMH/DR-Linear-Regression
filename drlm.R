library(stats)

`%||%` <- function(a, b) if (!is.null(a)) a else b

.safe_chol <- function(M, ctx = "") {
  kap <- tryCatch(kappa(M, exact = FALSE), error = function(e) Inf)
  if (kap > 1e12)
    warning(sprintf("drlm%s: near-singular matrix (kappa ~ %.2e). Check collinearity.",
                    if (nchar(ctx)) paste0(" [", ctx, "]") else "", kap), call. = FALSE)
  tryCatch(chol(M), error = function(e)
    stop("drlm: Cholesky failed. Check for perfect multicollinearity.\n",
         conditionMessage(e), call. = FALSE))
}
.chol_solve <- function(ch, rhs) as.vector(backsolve(ch, forwardsolve(t(ch), rhs)))

.chunk_ranges <- function(n, k) {
  rpc <- ceiling(n / k)
  Filter(Negate(is.null), lapply(seq_len(k), function(i) {
    lo <- (i - 1L) * rpc + 1L; hi <- min(i * rpc, n)
    if (lo > hi) NULL else c(lo, hi)
  }))
}
.check_chunk_size <- function(n, k, p) {
  if (floor(n / k) <= p)
    stop(sprintf("drlm: chunk size (%d) <= parameters (%d). Reduce k.", floor(n / k), p),
         call. = FALSE)
}

.simple_numeric <- function(formula, data) {
  tm <- terms(formula); labs <- attr(tm, "term.labels")
  yname <- if (attr(tm, "response")) all.vars(formula)[1] else NA_character_
  if (length(labs) == 0L || any(grepl("[():*^|/]", labs)) || !all(labs %in% names(data)) ||
      (!is.na(yname) && !yname %in% names(data))) return(list(ok = FALSE))
  if (!all(vapply(labs, function(v) is.numeric(data[[v]]) && !is.factor(data[[v]]), logical(1))))
    return(list(ok = FALSE))
  list(ok = TRUE, yname = yname, xvars = labs, intercept = attr(tm, "intercept") == 1L)
}
.fast_X <- function(block, xvars, intercept) {
  M <- as.matrix(block[, xvars, drop = FALSE]); storage.mode(M) <- "double"
  if (intercept) { X <- cbind(1, M); colnames(X) <- c("(Intercept)", xvars); X }
  else { colnames(M) <- xvars; M }
}
.design_of <- function(formula, block, contrasts = NULL) {
  sn <- .simple_numeric(formula, block)
  if (sn$ok) {
    list(X = .fast_X(block, sn$xvars, sn$intercept), y = as.numeric(block[[sn$yname]]))
  } else {
    mf <- model.frame(formula, block)
    mm <- model.matrix(attr(mf, "terms"), mf, contrasts.arg = contrasts)
    attr(mm, "contrasts") <- NULL; attr(mm, "assign") <- NULL
    list(X = mm, y = as.numeric(model.response(mf)))
  }
}

recombine_suffstats <- function(XtX, XtY, YtY, n, sy, pnames = NULL, alpha = 0.05) {
  p  <- NROW(XtX)
  ch <- .safe_chol(XtX, "recombine")
  B  <- .chol_solve(ch, XtY)
  RSS <- max(as.numeric(YtY - crossprod(B, XtY)), 0)
  TSS <- YtY - sy^2 / n
  sigma2 <- RSS / (n - p)
  XtXinv <- chol2inv(ch)
  se <- sqrt(diag(XtXinv) * sigma2)
  pnames <- pnames %||% rownames(XtX) %||% paste0("b", seq_len(p) - 1L)
  names(B) <- names(se) <- pnames
  zc <- qnorm(1 - alpha / 2); tval <- B / se; pv <- 2 * (1 - pnorm(abs(tval)))
  ct <- data.frame("Estimate" = B, "Std. Error" = se, "t value" = tval, "Pr(>|t|)" = pv,
                   "95% CI" = paste0("[", round(B - zc * se, 4), ", ", round(B + zc * se, 4), "]"),
                   check.names = FALSE, row.names = pnames)
  aic <- n * (log(RSS / n * 2 * pi) + 1) + 2 + 2 * p
  list(coefficients = ct, B = B, se = se, vcov = XtXinv * sigma2, XtXinv = XtXinv,
       deviance = RSS, null.deviance = TSS, r.squared = 1 - RSS / TSS,
       df.residual = n - p, df.null = n - 1L, sigma = sqrt(sigma2), aic = aic,
       n = n, p = p, XtX = XtX, XtY = XtY, YtY = YtY, sy = sy, alpha = alpha)
}

drlm <- function(formula, data, k = 1L, weights = NULL, offset = NULL,
                 subset = NULL, contrasts = NULL, na.action = na.omit,
                 alpha = 0.05, keep.fit = FALSE, x = FALSE, y = TRUE) {
  cl <- match.call()
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0,1).", call. = FALSE)

  if (is.function(data)) {
    XtX <- XtY <- NULL; YtY <- 0; n <- 0; sy <- 0; pnames <- NULL; nblk <- 0L
    blk <- data(reset = TRUE)
    while (!is.null(blk) && NROW(blk) > 0L) {
      d <- .design_of(formula, blk, contrasts)
      if (is.null(XtX)) { pnames <- colnames(d$X); XtX <- matrix(0, ncol(d$X), ncol(d$X)); XtY <- numeric(ncol(d$X)) }
      XtX <- XtX + crossprod(d$X); XtY <- XtY + as.vector(crossprod(d$X, d$y))
      YtY <- YtY + sum(d$y * d$y); n <- n + length(d$y); sy <- sy + sum(d$y); nblk <- nblk + 1L
      blk <- data(reset = FALSE)
    }
    if (is.null(XtX)) stop("drlm: the block reader returned no data.", call. = FALSE)
    res <- recombine_suffstats(XtX, XtY, YtY, n, sy, pnames, alpha)
    return(.assemble(res, cl, formula, engine = "streaming", k = nblk, yvec = NULL, X = NULL, w = NULL))
  }

  if (!is.data.frame(data)) stop("'data' must be a data.frame or a function(reset).", call. = FALSE)

  eval_in_data <- function(expr) if (is.null(expr)) NULL else eval(expr, envir = data, enclos = parent.frame(3))
  w_arg <- eval_in_data(cl$weights); off_arg <- eval_in_data(cl$offset); sub_arg <- eval_in_data(cl$subset)
  sn <- .simple_numeric(formula, data)
  use_fast <- sn$ok && is.null(sub_arg) && identical(na.action, na.omit) &&
              !anyNA(data[, c(sn$yname, sn$xvars), drop = FALSE])
  if (use_fast) {
    X <- .fast_X(data, sn$xvars, sn$intercept); yv <- as.numeric(data[[sn$yname]]); mt <- terms(formula)
  } else {
    mf_call <- match.call(expand.dots = FALSE)
    mf_call <- mf_call[c(1L, match(c("formula", "data", "subset", "weights", "offset", "na.action"), names(mf_call), 0L))]
    mf_call[[1L]] <- quote(stats::model.frame); mf_call$drop.unused.levels <- TRUE
    mf <- eval(mf_call, parent.frame()); mt <- attr(mf, "terms")
    yv <- as.numeric(model.response(mf, "any"))
    w_arg <- as.vector(model.weights(mf)) %||% w_arg; off_arg <- as.vector(model.offset(mf)) %||% off_arg
    X <- model.matrix(mt, mf, contrasts.arg = contrasts); attr(X, "contrasts") <- NULL; attr(X, "assign") <- NULL
  }
  n <- NROW(X); p <- NCOL(X); pnames <- colnames(X)
  w <- if (is.null(w_arg)) NULL else as.numeric(w_arg)
  off <- if (is.null(off_arg)) NULL else as.numeric(off_arg)
  ya <- if (is.null(off)) yv else yv - off
  k <- as.integer(k %||% 1L); if (k < 1L) k <- 1L
  if (k > 1L) .check_chunk_size(n, k, p)

  XtX <- matrix(0, p, p); XtY <- numeric(p); YtY <- 0; sy <- 0
  for (rg in .chunk_ranges(n, k)) {
    idx <- rg[1]:rg[2]; Xs <- X[idx, , drop = FALSE]; ys <- ya[idx]
    if (is.null(w)) {
      XtX <- XtX + crossprod(Xs); XtY <- XtY + as.vector(crossprod(Xs, ys)); YtY <- YtY + sum(ys * ys); sy <- sy + sum(ys)
    } else {
      ws <- w[idx]; rw <- sqrt(ws)
      XtX <- XtX + crossprod(Xs * rw); XtY <- XtY + as.vector(crossprod(Xs, ws * ys))
      YtY <- YtY + sum(ws * ys * ys); sy <- sy + sum(ws * ys)
    }
  }
  n_eff <- if (is.null(w)) n else sum(w)
  res <- recombine_suffstats(XtX, XtY, YtY, n_eff, sy, pnames, alpha)
  res$df.residual <- n - p; res$df.null <- n - 1L; res$n <- n
  if (keep.fit) {
    eta <- as.vector(X %*% res$B) + (off %||% 0)
    res$fitted.values <- eta; res$residuals <- yv - eta
  }
  .assemble(res, cl, formula, engine = "in-memory", k = k, yvec = if (y) yv else NULL,
            X = if (x) X else NULL, w = w, mt = mt)
}

.assemble <- function(res, cl, formula, engine, k, yvec, X, w, mt = NULL) {
  res$call <- cl; res$formula <- formula; res$terms <- mt %||% terms(formula)
  res$engine <- engine; res$k <- k; res$y <- yvec; res$x <- X; res$prior.weights <- w
  res$family <- "gaussian"; res$fkey <- "gaussian"
  class(res) <- "drlm"; res
}

drglm <- function(formula, family = gaussian, data, k = NULL, parallel = FALSE, ...) {
  fname <- if (is.character(family)) family else if (is.function(family)) family()$family else family$family
  if (!grepl("gaussian", tolower(fname))) stop("drlm.R implements the Gaussian (linear) model only; use the drglm package for other families.", call. = FALSE)
  drlm(formula, data = data, k = k %||% 1L, ...)
}

add_data <- function(object, newdata) {
  d <- .design_of(object$formula, newdata)
  res <- recombine_suffstats(object$XtX + crossprod(d$X), object$XtY + as.vector(crossprod(d$X, d$y)),
                             object$YtY + sum(d$y * d$y), object$n + length(d$y), object$sy + sum(d$y),
                             names(object$B), object$alpha)
  .assemble(res, object$call, object$formula, object$engine, object$k, NULL, NULL, NULL, object$terms)
}

robust_se <- function(object, X, y = NULL, type = c("HC1", "HAC"), lag = NULL) {
  type <- match.arg(type); B <- object$B; XtXinv <- object$XtXinv; p <- length(B)
  if (is.function(X)) {
    if (type == "HAC") stop("HAC needs the full in-memory design (row order across blocks).", call. = FALSE)
    meat <- matrix(0, p, p); n <- 0; blk <- X(reset = TRUE)
    while (!is.null(blk) && NROW(blk) > 0L) {
      d <- .design_of(object$formula, blk); e <- as.numeric(d$y - d$X %*% B)
      meat <- meat + crossprod(d$X * e); n <- n + length(e); blk <- X(reset = FALSE)
    }
    V <- (n / (n - p)) * XtXinv %*% meat %*% XtXinv
  } else {
    n <- nrow(X); e <- as.numeric(y - X %*% B); U <- X * e; meat <- crossprod(U)
    if (type == "HAC") {
      L <- lag %||% floor(4 * (n / 100)^(2 / 9)); w <- 1 - seq_len(L) / (L + 1)
      Upad <- rbind(matrix(0, L, p), U)
      M <- as.matrix(stats::filter(Upad, c(0, w), sides = 1))[-(1:L), , drop = FALSE]
      A <- crossprod(U, M); meat <- meat + A + t(A)
      V <- XtXinv %*% meat %*% XtXinv
    } else V <- (n / (n - p)) * XtXinv %*% meat %*% XtXinv
  }
  se <- sqrt(diag(V)); names(se) <- names(B); attr(se, "vcov") <- V; attr(se, "type") <- type; se
}

print.drlm <- function(x, digits = 4, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat(sprintf("Divide-and-recombine linear regression (%s, S = %s subsets)\n", x$engine, x$k))
  cat(sprintf("Obs: %s | Predictors: %d\n", formatC(x$n, format = "d", big.mark = ","), x$p))
  cat(strrep("-", 72), "\n\nCoefficients:\n")
  ct <- x$coefficients; ct2 <- ct
  for (cn in names(ct)) if (is.numeric(ct[[cn]])) ct2[[cn]] <- formatC(ct[[cn]], digits = digits, format = "g")
  ct2[["Pr(>|t|)"]] <- ifelse(ct[["Pr(>|t|)"]] < 2e-16, "< 2e-16", formatC(ct[["Pr(>|t|)"]], format = "g", digits = 3))
  print(ct2)
  cat(sprintf("\nResidual std error: %s on %s df;  R-squared: %s\n",
              formatC(x$sigma, digits = digits, format = "g"),
              formatC(x$df.residual, format = "d", big.mark = ","),
              formatC(x$r.squared, digits = 4, format = "f")))
  invisible(x)
}
summary.drlm      <- function(object, ...) object
coef.drlm         <- function(object, ...) object$B
vcov.drlm         <- function(object, ...) { V <- object$vcov; dimnames(V) <- list(names(object$B), names(object$B)); V }
confint.drlm      <- function(object, parm, level = 0.95, ...) {
  zc <- qnorm((1 + level) / 2); ci <- cbind(object$B - zc * object$se, object$B + zc * object$se)
  colnames(ci) <- paste0(formatC(c((1 - level) / 2, (1 + level) / 2) * 100, format = "g"), " %")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]; ci
}
deviance.drlm     <- function(object, ...) object$deviance
nobs.drlm         <- function(object, ...) object$n
df.residual.drlm  <- function(object, ...) object$df.residual
formula.drlm      <- function(x, ...) x$formula
fitted.drlm       <- function(object, ...) object$fitted.values %||% stop("Refit with keep.fit = TRUE.", call. = FALSE)
residuals.drlm    <- function(object, ...) object$residuals %||% stop("Refit with keep.fit = TRUE.", call. = FALSE)
predict.drlm      <- function(object, newdata, ...) {
  rhs <- as.formula(paste("~", deparse(object$formula[[3]])))
  X <- model.matrix(rhs, data = newdata)
  if (ncol(X) != length(object$B)) stop("'newdata' column mismatch with model formula.", call. = FALSE)
  as.vector(X %*% object$B)
}
logLik.drlm <- function(object, ...) structure(object$p + 1L - object$aic / 2, df = object$p + 1L, nobs = object$n, class = "logLik")
AIC.drlm    <- function(object, ..., k = 2) -2 * as.numeric(logLik(object)) + k * (object$p + 1L)
BIC.drlm    <- function(object, ...) AIC(object, k = log(object$n))
