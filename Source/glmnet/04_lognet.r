function (x, is.sparse, ix, jx, y, weights, offset, alpha, nobs, 
    nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd, 
    intr, vnames, maxit, kopt, family) 
{
    nc = dim(y)
    maxit = as.integer(maxit)
    if (is.null(nc)) {
        y = as.factor(y)
        ntab = table(y)
        minclass = min(ntab)
        if (minclass <= 1) 
            stop("one multinomial or binomial class has 1 or 0 observations; not allowed")
        if (minclass < 8) 
            warning("one multinomial or binomial class has fewer than 8  observations; dangerous ground")
        classnames = names(ntab)
        nc = as.integer(length(ntab))
        y = diag(nc)[as.numeric(y), ]
    }
    else {
        noo = nc[1]
        if (noo != nobs) 
            stop("x and y have different number of rows in call to glmnet", 
                call. = FALSE)
        nc = as.integer(nc[2])
        classnames = colnames(y)
    }
    maxvars = .Machine$integer.max/(nlam * nc)
    if (nx > maxvars) 
        stop(paste("Integer overflow; num_classes*num_lambda*pmax should not exceed .Machine$integer.max. Reduce pmax to be below", 
            trunc(maxvars)))
    if (!missing(weights)) 
        y = y * weights
    weights = drop(y %*% rep(1, nc))
    o = weights > 0
    if (!all(o)) {
        y = y[o, ]
        if (is.sparse) {
            x = sparseMatrix(i = jx, p = ix - 1, x = x, dims = c(nobs, 
                nvars))[o, , drop = FALSE]
            ix = as.integer(x@p + 1)
            jx = as.integer(x@i + 1)
            x = as.double(x@x)
        }
        else x = x[o, , drop = FALSE]
        nobs = sum(o)
    }
    else o = NULL
    if (family == "binomial") {
        if (nc > 2) 
            stop("More than two classes; use multinomial family instead in call to glmnet", 
                call. = FALSE)
        nc = as.integer(1)
        y = y[, c(2, 1)]
    }
    storage.mode(y) = "double"
    if (is.null(offset)) {
        offset = y * 0
        is.offset = FALSE
    }
    else {
        offset = as.matrix(offset)
        if (!is.null(o)) 
            offset = offset[o, , drop = FALSE]
        do = dim(offset)
        if (do[[1]] != nobs) 
            stop("offset should have the same number of values as observations in binomial/multinomial call to glmnet", 
                call. = FALSE)
        if ((do[[2]] == 1) & (nc == 1)) 
            offset = cbind(offset, -offset)
        if ((family == "multinomial") & (do[[2]] != nc)) 
            stop("offset should have same shape as y in multinomial call to glmnet", 
                call. = FALSE)
        storage.mode(offset) = "double"
        is.offset = TRUE
    }
    fit = if (is.sparse) 
        .Fortran("splognet", parm = alpha, nobs, nvars, nc, x, 
            ix, jx, y, offset, jd, vp, cl, ne = ne, nx, nlam, 
            flmin, ulam, thresh, isd, intr, maxit, kopt, lmu = integer(1), 
            a0 = double(nlam * nc), ca = double(nx * nlam * nc), 
            ia = integer(nx), nin = integer(nlam), nulldev = double(1), 
            dev = double(nlam), alm = double(nlam), nlp = integer(1), 
            jerr = integer(1), PACKAGE = "glmnet")
    else .Fortran("lognet", parm = alpha, nobs, nvars, nc, as.double(x), 
        y, offset, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, 
        isd, intr, maxit, kopt, lmu = integer(1), a0 = double(nlam * 
            nc), ca = double(nx * nlam * nc), ia = integer(nx), 
        nin = integer(nlam), nulldev = double(1), dev = double(nlam), 
        alm = double(nlam), nlp = integer(1), jerr = integer(1), 
        PACKAGE = "glmnet")
    if (fit$jerr != 0) {
        errmsg = jerr(fit$jerr, maxit, pmax = nx, family)
        if (errmsg$fatal) 
            stop(errmsg$msg, call. = FALSE)
        else warning(errmsg$msg, call. = FALSE)
    }
    if (family == "binomial") {
        outlist = getcoef(fit, nvars, nx, vnames)
    }
    else outlist = getcoef.multinomial(fit, nvars, nx, vnames, 
        nc, classnames)
    dev = fit$dev[seq(fit$lmu)]
    outlist = c(outlist, list(dev.ratio = dev, nulldev = fit$nulldev, 
        npasses = fit$nlp, jerr = fit$jerr, offset = is.offset, 
        classnames = classnames))
    if (family == "multinomial") {
        if (kopt == 2) 
            grouped = TRUE
        else grouped = FALSE
        outlist$grouped = grouped
    }
    class(outlist) = switch(family, binomial = "lognet", multinomial = "multnet")
    outlist
}