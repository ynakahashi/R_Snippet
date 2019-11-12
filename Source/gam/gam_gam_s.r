### gam::gam.s
function (x, y, w = rep(1, length(x)), df = 4, spar = 1, xeval) 
{
    ### 後から Fortran に渡すので型の指定
    storage.mode(x) <- storage.mode(y) <- storage.mode(w) <- storage.mode(spar) <- storage.mode(df) <- "double"

    n <- as.integer(length(x))
    x <- signif(x, 6)

    ### 入力 x について、指定の digits で丸めた上でユニークなデータ数を得る
    ### デフォルトでは digits は 6
    mat <- gam.match(x)
    omat <- mat$o
    nef <- mat$nef

    ### 何これ
    work.len <- max(3 * nef + 2 * n + 10, (10 + 2 * 4) * (nef + 
        2) + 5 * nef + n + 15)
    
    ### 
    fit <- .Fortran("splsm", x, y, w, n, omat, nef, spar = spar, 
        df = df, s = double(n), s0 = double(1), var = double(nef), 
        FALSE, work = double(work.len), PACKAGE = "gam")
    if (missing(xeval)) 
        list(residuals = y - fit$s, nl.df = fit$df - 1, var = fit$var[omat])
    else {
        skn <- .Fortran("sknotl", fit$work[seq(nef)], nef, knot = double(nef + 
            6), k = integer(1), PACKAGE = "gam")
        smallest <- x[omat == 1][1]
        largest <- x[omat == nef][1]
        k <- skn$k
        gam.sp(xeval, skn$knot[seq(k)], k - 4, fit$work[seq(3 * 
            nef + n + 10, length = k - 4)], smallest, largest - 
            smallest)
    }
}