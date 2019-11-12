### gam::gam.fit
### gam.fit は x, y に加えて smooth.frame を受け取る。これは gam で作った mf で、中身は平滑化に関する情報を持った data.frame 
function (x, y, smooth.frame, weights = rep(1, nobs), start = NULL, 
    etastart = NULL, mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = gam.control()) 
{
    ### 列名
    ynames <- if (is.matrix(y)) 
        dimnames(y)[[1]]
    else names(y)
    xnames <- dimnames(x)[[2]]

    ### データの行数、列数
    nobs <- NROW(y)
    nvars <- ncol(x)

    ## その他の gam.control オプション
    maxit <- control$maxit
    bf.maxit <- control$bf.maxit
    epsilon <- control$epsilon
    bf.epsilon <- control$bf.epsilon
    trace <- control$trace

    ### digits, weigths, offset
    digits <- -log10(epsilon) + 1
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    
    ### family に関するパラメータ
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("illegal `family' argument")
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    eval(family$initialize)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    ### eta の初期値
    eta <- if (!is.null(etastart)) 
        etastart
    
    ### エラーチェック
    else if (!is.null(start)) 
        if (length(start) != nvars) 
            stop("Length of start should equal ", nvars, " and correspond to initial coefs for ", 
                deparse(xnames))
        else {
            coefold <- start
            offset + as.vector(if (NCOL(x) == 1) 
                x * start
            else x %*% start)
        }
    else family$linkfun(mustart)

    ### mu の初期値
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
        stop("Can't find valid starting values: please specify some")

    ### デビアンス（残差平方和）
    new.dev <- sum(dev.resids(y, mu, weights))

    ### smoothers が指定されている場合に smoother を抽出する
    ### 今回のケースでは s が入る
    a <- attributes(attr(smooth.frame, "terms"))
    smoothers <- a$specials
    if (length(smoothers) > 0) {
        smoothers <- smoothers[sapply(smoothers, length) > 0]
        for (i in seq(along = smoothers)) {
            tt <- smoothers[[i]]
            ff <- apply(a$factors[tt, , drop = FALSE], 2, any)
            smoothers[[i]] <- if (any(ff)) 
                seq(along = ff)[a$order == 1 & ff]
            else NULL
        }
    }

    ### smoother が指定されている場合、ここが処理される
    if (length(smoothers) > 0) {
        gam.wlist = gam.smoothers()$wlist
        smooth.labels <- a$term.labels[unlist(smoothers)]
        assignx <- attr(x, "assign")
        assignx <- assign.list(assignx, a$term.labels)
        which <- assignx[smooth.labels]

        ### ２つ以上の smoothers が指定されている場合は backfit として general.wam が指定される
        ### wam は weighted additive model
        if (length(smoothers) > 1) 
            bf <- "general.wam"

        ### 今回のケース（平滑化の関数として s を指定）はこちらで、 s.wam が指定される。
        else {
            sbf <- match(names(smoothers), gam.wlist, FALSE)
            bf <- if (sbf) 
                paste(gam.wlist[sbf], "wam", sep = ".")
            else "general.wam"
        }

        ### bf にオプション部分を文字列で結合
        bf.call <- parse(text = paste(bf, "(x, z, wz, fit$smooth, which, fit$smooth.frame,bf.maxit,bf.epsilon, trace)", 
            sep = ""))[[1]]
        s <- matrix(0, length(y), length(which))
        dimnames(s) <- list(names(y), names(which))
        fit <- list(smooth = s, smooth.frame = smooth.frame)
    }

    ### 平滑化しない場合。 通常の lm.wfit に渡す。なお general.wam でも lm.wfit が使われる。
    ### ここの深掘りは以前のGLMの記事を紹介
    else {
        bf.call <- expression(lm.wfit(x, z, wz, method = "qr", 
            singular.ok = TRUE))
        bf <- "lm.wfit"
    }

    ### ここから反復に入る
    old.dev <- 10 * new.dev + 10
    for (iter in 1:maxit) {

        ### weight が 0 のデータは除外する
        good <- weights > 0
        varmu <- variance(mu)
        if (any(is.na(varmu[good]))) 
            stop("NAs in V(mu)")
        if (any(varmu[good] == 0)) 
            stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good]))) 
            stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)

        ### z を生成。ただし bf.call で二番目の引数は y なので、この z は目的変数の意味
        z <- eta - offset
        z[good] <- z[good] + (y - mu)[good]/mu.eta.val[good]

        ### wz を生成。重み。今回のケースでは無視して良い
        wz <- weights
        wz[!good] <- 0
        wz[good] <- wz[good] * mu.eta.val[good]^2/varmu[good]

        ### ここで bf.call が評価される。 s.wam の場合、bakfit が呼ばれる。
        ### bf.call で指定されている smooth.frame はこの時点では単なる data.frame 
        fit <- eval(bf.call)

        ### 予測値にオフセットを加算する
        eta <- fit$fitted.values + offset

        ### eta から mu に変換する
        mu <- linkinv(eta)

        ### デビアンスを更新
        old.dev <- new.dev
        new.dev <- sum(dev.resids(y, mu, weights))

        ### ここの trace は対角和ではなく、 gam のオプションでループごとのデビアンスをモニターできる
        if (trace) 
            cat("GAM ", bf, " loop ", iter, ": deviance = ", 
                format(round(new.dev, digits)), " \n", sep = "")

        ### ループの打ち切り判定
        ### デビアンスが NA となった場合、警告を出して打ち切る
        if (is.na(new.dev)) {
            one.more <- FALSE
            warning("iterations terminated prematurely because of singularities")
        }
        ### 差が十分に小さければ終了
        else one.more <- abs(old.dev - new.dev)/(old.dev + 0.1) > 
            epsilon
        if (!one.more) 
            break
    }

    fitqr <- fit$qr
    xxnames <- xnames[fitqr$pivot]
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
        Rmat <- diag(nvars)
        Rmat[1:nr, 1:nvars] <- fitqr$qr[1:nr, 1:nvars]
    }
    else Rmat <- fitqr$qr[1:nvars, 1:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    dimnames(Rmat) <- list(xxnames, xxnames)
    names(fit$residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames

    ### eta and mu
    fit$additive.predictors <- eta
    fit$fitted.values <- mu
    names(fit$weights) <- ynames
    names(fit$effects) <- c(xxnames[seq(len = fitqr$rank)], rep.int("", 
        sum(good) - fitqr$rank))
    if (length(fit$smooth) > 0) 
        fit$smooth.frame <- smooth.frame[smooth.labels]
    wtdmu <- if (a$intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(a$intercept)
    rank <- n.ok - fit$df.residual
    aic.model <- aic(y, nobs, mu, weights, new.dev) + 2 * rank
    if (!is.null(fit$smooth)) {
        nonzeroWt <- (wz > 0)
        nl.chisq <- gam.nlchisq(fit$qr, fit$residuals, wz, fit$smooth)
    }
    else nl.chisq <- NULL
    fit <- c(fit, list(R = Rmat, rank = fitqr$rank, family = family, 
        deviance = new.dev, aic = aic.model, null.deviance = nulldev, 
        iter = iter, prior.weights = weights, y = y, df.null = nulldf, 
        nl.chisq = nl.chisq))
    fit
}