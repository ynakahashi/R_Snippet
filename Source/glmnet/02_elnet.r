function (x, is.sparse, ix, jx, y, weights, offset, type.gaussian = c("covariance", 
    "naive"), alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, 
    ulam, thresh, isd, intr, vnames, maxit) 
{
    # 1. パラメータの受け取り
    ### maxit
    maxit = as.integer(maxit)
    ### weights
    weights = as.double(weights)
    ### type.gaussian
    type.gaussian = match.arg(type.gaussian)
    ka = as.integer(switch(type.gaussian, covariance = 1, naive = 2, 
        ))
    ### y の storage.mode
    storage.mode(y) = "double"
    ### offset
    if (is.null(offset)) {
        is.offset = FALSE
    }
    else {
        storage.mode(offset) = "double"
        is.offset = TRUE
        y = y - offset
    }
    ### 重み付き平均
    ybar = weighted.mean(y, weights)
    ### Null Deviance（帰無モデルの残差逸脱度）
    nulldev = sum(weights * (y - ybar)^2)

    if (nulldev == 0) 
        stop("y is constant; gaussian glmnet fails at standardization step")
    
    # 2. フィッティング
    ## 疎行列であるかで関数を変える
    fit = if (is.sparse) 
        .Fortran("spelnet", ka, parm = alpha, nobs, nvars, x, 
            # ix, jx は疎行列における非ゼロの要素の累積個数と行番号
            ix, jx, 
            y, weights, jd, vp, cl, ne, nx, nlam, flmin, 
            ulam, thresh, isd, intr, maxit, lmu = integer(1), 
            a0 = double(nlam), ca = double(nx * nlam), ia = integer(nx), 
            nin = integer(nlam), rsq = double(nlam), alm = double(nlam), 
            nlp = integer(1), jerr = integer(1), PACKAGE = "glmnet")
    else .Fortran("elnet", ka, parm = alpha, nobs, nvars, as.double(x), 
        y, weights, jd, vp, cl, ne, nx, nlam, flmin, 
        ulam, thresh, isd, intr, maxit, lmu = integer(1), 
        a0 = double(nlam), ca = double(nx * nlam), ia = integer(nx), 
        nin = integer(nlam),  rsq = double(nlam), alm = double(nlam), 
        nlp = integer(1), jerr = integer(1), PACKAGE = "glmnet")
        # nx は 非ゼロの変数の個数
        # nlam は検証する lambda の個数
        # なので ca は変数の数 * lambda の数

    # 3. 後処理
    ## エラーチェック
    if (fit$jerr != 0) {
        errmsg = jerr(fit$jerr, maxit, pmax = nx, family = "gaussian")
        if (errmsg$fatal) 
            stop(errmsg$msg, call. = FALSE)
        else warning(errmsg$msg, call. = FALSE)
    }
    ## パラメータ（切片、回帰係数、自由度、次元、lambda）を取ってくる
    outlist = getcoef(fit, nvars, nx, vnames)
    ## パラメータ（xxxxxxxxxx）を取ってきて outlist に結合する
    dev = fit$rsq[seq(fit$lmu)]
    outlist = c(outlist, list(dev.ratio = dev, nulldev = nulldev, 
        npasses = fit$nlp, jerr = fit$jerr, offset = is.offset))
    ## elnet クラスを付与する
    class(outlist) = "elnet"
    outlist
}