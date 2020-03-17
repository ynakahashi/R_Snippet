function (x, y, family = c("gaussian", "binomial", "poisson", 
    "multinomial", "cox", "mgaussian"), weights, offset = NULL, 
    alpha = 1, nlambda = 100, lambda.min.ratio = ifelse(nobs < 
        nvars, 0.01, 1e-04), lambda = NULL, standardize = TRUE, 
    intercept = TRUE, thresh = 1e-07, dfmax = nvars + 1, pmax = min(dfmax * 
        2 + 20, nvars), exclude, penalty.factor = rep(1, nvars), 
    lower.limits = -Inf, upper.limits = Inf, maxit = 1e+05, type.gaussian = ifelse(nvars < 
        500, "covariance", "naive"), type.logistic = c("Newton", 
        "modified.Newton"), standardize.response = FALSE, type.multinomial = c("ungrouped", 
        "grouped"), relax = FALSE, trace.it = 0, ...) 
{
    # パラメータの設定、前処理、エラーチェック
    ## 指定したfamilyが引数としてOKかチェック
    ## https://hoxo-m.hatenablog.com/entry/20111025/p1
    family = match.arg(family)

    ## alpha
    ## LassoとRidgeそれぞれに対するペナルティの配分を決めるパラメータ
    ## glmnetにおける罰則項は以下で定義
    ## (1 − α)/2||β||^2_2 + α||β||_1
    ## alphaは0~1で、1ならLasso、0ならRidgeに対応
    if (alpha > 1) {
        warning("alpha >1; set to 1")
        alpha = 1
    }
    if (alpha < 0) {
        warning("alpha<0; set to 0")
        alpha = 0
    }
    alpha = as.double(alpha)
    
    ## match.call
    ## 引数を順番通り、正式な名前に直してくれる
    ## match.call returns a call in which all of the specified arguments are specified by their full names.
    ## https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/match.call
    ## myfun <- function(x = 100, pos = 3, dim = 100) { return(x + pos + dim) }
    ## myfun(1, 2, 3)
    ## match.call(myfun, call("myfun", 1, d = 3, p = 5))
    this.call = match.call()
    
    ### lambda
    nlam = as.integer(nlambda)
    y = drop(y)
    np = dim(x)
    if (is.null(np) | (np[2] <= 1)) 
        stop("x should be a matrix with 2 or more columns")
    nobs = as.integer(np[1])
    if (missing(weights)) 
        weights = rep(1, nobs)
    else if (length(weights) != nobs) 
        stop(paste("number of elements in weights (", length(weights), 
            ") not equal to the number of rows of x (", nobs, 
            ")", sep = ""))
    nvars = as.integer(np[2])
    dimy = dim(y)
    nrowy = ifelse(is.null(dimy), length(y), dimy[1])
    if (nrowy != nobs) 
        stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (", 
            nobs, ")", sep = ""))
    vnames = colnames(x)
    if (is.null(vnames)) 
        vnames = paste("V", seq(nvars), sep = "")
    ### dfmax = nvars + 1
    ne = as.integer(dfmax)
    ### pmax = min(dfmax * 2 + 20, nvars)
    nx = as.integer(pmax)
    if (missing(exclude)) 
        exclude = integer(0)
    
    ###
    if (any(penalty.factor == Inf)) {
        exclude = c(exclude, seq(nvars)[penalty.factor == Inf])
        exclude = sort(unique(exclude))
    }
    if (length(exclude) > 0) {
        jd = match(exclude, seq(nvars), 0)
        if (!all(jd > 0)) 
            stop("Some excluded variables out of range")
        penalty.factor[jd] = 1
        jd = as.integer(c(length(jd), jd))
    }
    else jd = as.integer(0)
    vp = as.double(penalty.factor)
    internal.parms = glmnet.control()
    if (internal.parms$itrace) 
        trace.it = 1
    else {
        if (trace.it) {
            glmnet.control(itrace = 1)
            on.exit(glmnet.control(itrace = 0))
        }
    }

    ### 上限・下限
    if (any(lower.limits > 0)) {
        stop("Lower limits should be non-positive")
    }
    if (any(upper.limits < 0)) {
        stop("Upper limits should be non-negative")
    }
    lower.limits[lower.limits == -Inf] = -internal.parms$big
    upper.limits[upper.limits == Inf] = internal.parms$big
    if (length(lower.limits) < nvars) {
        if (length(lower.limits) == 1) 
            lower.limits = rep(lower.limits, nvars)
        else stop("Require length 1 or nvars lower.limits")
    }
    else lower.limits = lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars) {
        if (length(upper.limits) == 1) 
            upper.limits = rep(upper.limits, nvars)
        else stop("Require length 1 or nvars upper.limits")
    }
    else upper.limits = upper.limits[seq(nvars)]
    cl = rbind(lower.limits, upper.limits)
    if (any(cl == 0)) {
        fdev = glmnet.control()$fdev
        if (fdev != 0) {
            glmnet.control(fdev = 0)
            on.exit(glmnet.control(fdev = fdev))
        }
    }
    storage.mode(cl) = "double"
    isd = as.integer(standardize)
    intr = as.integer(intercept)
    if (!missing(intercept) && family == "cox") 
        warning("Cox model has no intercept")
    jsd = as.integer(standardize.response)
    thresh = as.double(thresh)
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1) 
            stop("lambda.min.ratio should be less than 1")
        flmin = as.double(lambda.min.ratio)
        ulam = double(1)
    }
    else {
        flmin = as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam = as.double(rev(sort(lambda)))
        nlam = as.integer(length(lambda))
    }
    is.sparse = FALSE
    ix = jx = NULL
    if (inherits(x, "sparseMatrix")) {
        is.sparse = TRUE
        x = as(x, "CsparseMatrix")
        x = as(x, "dgCMatrix")
        ix = as.integer(x@p + 1)
        jx = as.integer(x@i + 1)
        x = as.double(x@x)
    }
    if (trace.it) {
        if (relax) 
            cat("Training Fit\n")
        pb <- createPB(min = 0, max = nlam, initial = 0, style = 3)
    }

    ### 最適化の手法
    kopt = switch(match.arg(type.logistic), Newton = 0, modified.Newton = 1)
    if (family == "multinomial") {
        type.multinomial = match.arg(type.multinomial)
        if (type.multinomial == "grouped") 
            kopt = 2
    }
    kopt = as.integer(kopt)

    # フィッティング
    fit = switch(family, 
                    gaussian = elnet(x, is.sparse, ix, jx, 
        y, weights, offset, type.gaussian, alpha, nobs, nvars, 
        jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd, intr, 
        vnames, maxit), 
                    poisson = fishnet(x, is.sparse, ix, jx, 
        y, weights, offset, alpha, nobs, nvars, jd, vp, cl, ne, 
        nx, nlam, flmin, ulam, thresh, isd, intr, vnames, maxit),
                    binomial = lognet(x, is.sparse, ix, jx, y, weights, offset, 
            alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, 
            ulam, thresh, isd, intr, vnames, maxit, kopt, family), 
                    multinomial = lognet(x, is.sparse, ix, jx, y, weights, 
            offset, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, 
            flmin, ulam, thresh, isd, intr, vnames, maxit, kopt, 
            family), 
                    cox = coxnet(x, is.sparse, ix, jx, y, weights, 
            offset, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, 
            flmin, ulam, thresh, isd, vnames, maxit), 
                    mgaussian = mrelnet(x, 
            is.sparse, ix, jx, y, weights, offset, alpha, nobs, 
            nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, 
            isd, jsd, intr, vnames, maxit))
    if (trace.it) {
        utils::setTxtProgressBar(pb, nlam)
        close(pb)
    }
    if (is.null(lambda)) 
        fit$lambda = fix.lam(fit$lambda)
    fit$call = this.call
    fit$nobs = nobs
    class(fit) = c(class(fit), "glmnet")
    
    # リターン
    ## relax???
    if (relax) 
        relax.glmnet(fit, x = x, y = y, weights = weights, offset = offset, 
            lower.limits = lower.limits, upper.limits = upper.limits, 
            check.args = FALSE, ...)
    else fit
}