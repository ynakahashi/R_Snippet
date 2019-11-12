### gam::gam
function (formula, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, control = gam.control(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, ...) 
{
    ### 関数の引数を名前付きで確定。gam(wage ~ s(year, 4) + education, Wage)として与えた場合、
    ### formula = と data = がそれぞれ保持される。
    ### match.call returns a call in which all of the specified arguments are specified by their full names.
    call <- match.call()

    ### familyの判定
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }

    ### データが指定されていない場合？
    if (missing(data)) 
        data <- environment(formula)
    
    ### 指定されている引数の取り出し
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "etastart", 
        "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]

    ### 指定されていない引数の指定し、stats::model.frame()の形式に仕立てる
    mf$na.action = quote(na.pass)
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)

    ### 平滑化の関数を取ってくる（s, lo, random)
    gam.slist <- gam.smoothers()$slist

    ### term を mf$formula に渡す
    mt <- if (missing(data)) 
        terms(formula, gam.slist)
    else terms(formula, gam.slist, data = data)
    mf$formula <- mt

    ### ここで mf 、つまり model.frame が実行されて data.frame になる
    ### ただし平滑化は実行されず、平滑化のパラメータは attribute として持っている
    mf <- eval(mf, parent.frame())
    if (missing(na.action)) {
        naa = getOption("na.action", "na.fail")
        na.action = get(naa)
    }
    mf = na.action(mf)
    mt = attributes(mf)[["terms"]]

    ### method の指定によって処理を分ける。 glm.fit または glm.fit.null 以外の場合はエラー
    switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1, 
        stop("invalid `method': ", method))
    
    ### Y
    Y <- model.response(mf, "any")

    ### X を matrix で。
    ### gamを実行したときのエラーメッセージ（ `non-list contrasts argument ignored` ）はここで出ている。 contrasts の指定が良くない様子。
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)

    ### その他パラメータ（weights, offset, mustart, etastart）
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0)) 
        stop("Negative wts not allowed")
    if (!is.null(offset) && length(offset) != NROW(Y)) 
        stop("Number of offsets is ", length(offset), ", should equal ", 
            NROW(Y), " (number of observations)")
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

    ### ここが本体。 gam.fit を呼び出している
    fit <- gam.fit(x = X, y = Y, smooth.frame = mf, weights = weights, 
        start = start, etastart = etastart, mustart = mustart, 
        offset = offset, family = family, control = control)
    
    ### offset が指定されており intercept 項がある場合
    if (length(offset) && attr(mt, "intercept") > 0) {
        fit$null.dev <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, weights = weights, offset = offset, family = family, 
            control = control[c("epsilon", "maxit", "trace")], 
            intercept = TRUE)$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c("Gam", "glm", "lm")
    if (!is.null(fit$df.residual) && !(fit$df.residual > 0)) 
        warning("Residual degrees of freedom are negative or zero.  This occurs when the sum of the parametric and nonparametric degrees of freedom exceeds the number of observations.  The model is probably too complex for the amount of data available.")
    fit
}