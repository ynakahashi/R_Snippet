### gam::s.wam
function (x, y, w, s, which, smooth.frame, maxit = 30, tol = 1e-07, 
    trace = FALSE, se = TRUE, ...) 
{
    ### smooth.frame は data.frame なのでここが評価される
    if (is.data.frame(smooth.frame)) {
        first <- TRUE
        data <- smooth.frame[, names(which), drop = FALSE]

        ### 入力 x について、指定の digits で丸めた上でユニークなデータ数を得る
        ### デフォルトでは digits は 6
        smooth.frame <- gam.match(data)
        dx <- as.integer(dim(x))
        smooth.frame$n <- dx[1] # number of row（ユニーク）
        smooth.frame$p <- dx[2] # number of column
        oldClass(data) <- NULL
        smooth.frame$spar <- unlist(lapply(data, attr, "spar"))
        smooth.frame$df <- unlist(lapply(data, attr, "df"))
    }
    else first <- FALSE

    ### 後でFortranに渡すために必要な指定
    storage.mode(tol) <- "double"
    storage.mode(maxit) <- "integer"
    which <- unlist(which)
    storage.mode(which) <- "integer"
    storage.mode(y) <- "double"
    storage.mode(w) <- "double"
    p <- smooth.frame$p
    n <- smooth.frame$n

    ### ここから平滑化対象の変数を順番に処理する
    for (ich in which) x[, ich] = signif(x[, ich], 6) # 6桁に丸める

    ### ここが本体。Fortran で書かれた bakfit を呼び出す
    fit <- .Fortran("bakfit", x, npetc = as.integer(c(n, p, length(which), 
        se, 0, maxit, 0)), y = y, w = w, which, spar = as.double(smooth.frame$spar), 
        df = as.double(smooth.frame$df), as.integer(smooth.frame$o), 
        as.integer(smooth.frame$nef), etal = double(n), s = s, 
        eta = double(n), beta = double(p), var = s, tol, qr = x, 
        qraux = double(p), qpivot = as.integer(1:p), effects = double(n), 
        double((10 + 2 * 4 + 5) * (max(smooth.frame$nef) + 2) + 
            15 * n + 15 + length(which)), PACKAGE = "gam")
    
    
    nit <- fit$npetc[5]
    qrank <- fit$npetc[7]
    if ((nit == maxit) & maxit > 1) 
        warning(paste("s.wam convergence not obtained in ", maxit, 
            " iterations"))
    if (first) {
        smooth.frame$spar <- fit$spar
        first <- FALSE
    }
    names(fit$df) <- dimnames(s)[[2]]
    names(fit$beta) <- labels(x)[[2]]
    qrx <- structure(list(qr = fit$qr, qraux = fit$qraux, rank = qrank, 
        pivot = fit$qpivot, tol = 1e-07), class = "qr")
    effects <- fit$effects
    r1 <- seq(len = qrx$rank)
    dn <- colnames(x)
    if (is.null(dn)) 
        dn <- paste("x", 1:p, sep = "")
    names(effects) <- c(dn[qrx$pivot[r1]], rep.int("", n - qrx$rank))
    rl <- list(coefficients = fit$beta, residuals = fit$y - fit$eta, 
        fitted.values = fit$eta, effects = effects, weights = w, 
        rank = qrank, assign = attr(x, "assign"), qr = qrx, smooth = fit$s, 
        nl.df = fit$df - 1)
    rl$df.residual <- n - qrank - sum(rl$nl.df) - sum(fit$w == 
        0)
    if (se) 
        rl <- c(rl, list(var = fit$var))
    c(list(smooth.frame = smooth.frame), rl)
}