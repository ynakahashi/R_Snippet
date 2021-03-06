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
    ### https://hoxo-m.hatenablog.com/entry/20111025/p1
    family = match.arg(family)

    ## alpha
    ### LassoとRidgeそれぞれに対するペナルティの配分を決めるパラメータ
    ### glmnetにおける罰則項は以下で定義
    ### (1 − α)/2||β||^2_2 + α||β||_1
    ### alphaは0~1で、1ならLasso、0ならRidgeに対応
    ### なんで1/2なんだろ？？
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
    ### 引数を順番通り、正式な名前に直してくれる
    ### match.call returns a call in which all of the specified arguments are specified by their full names.
    ### https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/match.call
    ### myfun <- function(x = 100, pos = 3, dim = 100) { return(x + pos + dim) }
    ### myfun(1, 2, 3)
    ### match.call(myfun, call("myfun", 1, d = 3, p = 5))
    this.call = match.call()
    
    ## nlambda
    ### number of lambda
    ### 検証する lambda の数を指定する。デフォルトは 100
    nlam = as.integer(nlambda)

    ## drop は length が 1 である冗長な次元を落とす
    y = drop(y)

    ## x
    ### x は２列以上持つ必要があるので、単回帰はできない
    np = dim(x)
    if (is.null(np) | (np[2] <= 1)) 
        stop("x should be a matrix with 2 or more columns")
    ### x のレコード数
    nobs = as.integer(np[1])
    ### weights
    ### 未入力のときは 1 を与え、weights と nobs が一致しないときはエラー
    if (missing(weights)) 
        weights = rep(1, nobs)
    else if (length(weights) != nobs) 
        stop(paste("number of elements in weights (", length(weights), 
            ") not equal to the number of rows of x (", nobs, 
            ")", sep = ""))
    ### 変数の数
    nvars = as.integer(np[2])

    ## y
    dimy = dim(y)
    ### y のレコード数 
    nrowy = ifelse(is.null(dimy), length(y), dimy[1])
    ### y と x でレコード数が合わないときはエラー
    if (nrowy != nobs) 
        stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (", 
            nobs, ")", sep = ""))

    ## 変数名
    vnames = colnames(x)
    if (is.null(vnames)) 
        vnames = paste("V", seq(nvars), sep = "")

    ## 自由度
    ### モデルに含まれる変数の上限を指定
    ### dfmax = nvars + 1
    ne = as.integer(dfmax)

    ### 非ゼロとする変数の数の上限？
    ### Limit the maximum number of variables ever to be nonzero
    ### pmax = min(dfmax * 2 + 20, nvars)
    nx = as.integer(pmax)

    ### 除外対象となる変数の指定
    if (missing(exclude)) 
        exclude = integer(0)
    
    ## 変数ごとに異なるペナルティを与える
    ### デフォルトは 1 が入る
    ### Inf が指定されている変数は exclude として扱われる
    ### Default is 1 for all variables (and implicitly infinity for variables listed in exclude)
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

    ## 内部でデフォルトで持っているパラメータ 
    internal.parms = glmnet.control()
    ### プログレスバーを表示する！
    if (internal.parms$itrace) 
        trace.it = 1
    else {
        if (trace.it) {
            glmnet.control(itrace = 1)
            on.exit(glmnet.control(itrace = 0))
        }
    }

    ## 上限・下限
    ### lower.limit としては非正の値のみ指定できる
    if (any(lower.limits > 0)) {
        stop("Lower limits should be non-positive")
    }
    ### upper.limtit は逆に非負の値のみ指定できる
    if (any(upper.limits < 0)) {
        stop("Upper limits should be non-negative")
    }
    ### Inf になっているものについてはデフォルトの値(9.9e35)に差し替え  
    lower.limits[lower.limits == -Inf] = -internal.parms$big
    upper.limits[upper.limits == Inf] = internal.parms$big
    ### nvars との整合性チェック
    if (length(lower.limits) < nvars) {
        ### lower.limits としてスカラが指定されている場合は nvars 全てに適用
        if (length(lower.limits) == 1) 
            lower.limits = rep(lower.limits, nvars)
        else stop("Require length 1 or nvars lower.limits")
    }
    ### lower.limits が nvars よりも長い場合は前から利用する
    else lower.limits = lower.limits[seq(nvars)]
    ### nvars との整合性チェック（lower.limits と同様）
    if (length(upper.limits) < nvars) {
        if (length(upper.limits) == 1) 
            upper.limits = rep(upper.limits, nvars)
        else stop("Require length 1 or nvars upper.limits")
    }
    else upper.limits = upper.limits[seq(nvars)]
    cl = rbind(lower.limits, upper.limits)
    ### lower または upper に 0 を含む場合
    ### ここ何をやっているのかよくわからない
    if (any(cl == 0)) {
        ### fdev は最小となるデビアンスの変化量(割合)
        ### minimum fractional change in deviance for stopping path; factory default = 1.0e5
        fdev = glmnet.control()$fdev
        if (fdev != 0) {
            glmnet.control(fdev = 0)
            on.exit(glmnet.control(fdev = fdev))
        }
    }
    storage.mode(cl) = "double"

    ## 標準化
    ### standardize と intercept はデフォルトは TRUE なので 1 になる
    isd = as.integer(standardize)
    intr = as.integer(intercept)
    ### Cox回帰における警告
    if (!missing(intercept) && family == "cox") 
        warning("Cox model has no intercept")
    ### standardize.response は family="mgaussian" のときに目的変数を標準化するかの指定
    jsd = as.integer(standardize.response)

    ## 収束判定
    ### coordinate descent における収束の閾値
    thresh = as.double(thresh)

    ## lambda
    ### ペナルティの大きさ
    ### 指定がない場合、flmin(下限？)とulam(上限？)は lambda.min.ratio および 1 に指定される
    ### lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04)
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1) 
            stop("lambda.min.ratio should be less than 1")
        flmin = as.double(lambda.min.ratio)
        ulam = double(1)
    }
    ### 指定がある場合、flmin(下限？)とulam(上限？)は 1 および lambdaの降順 に指定される
    else {
        flmin = as.double(1)
        if (any(lambda < 0)) 
            stop("lambdas should be non-negative")
        ulam = as.double(rev(sort(lambda)))
        nlam = as.integer(length(lambda))
    }

    ## sparse matrix 
    ### x が Matrix::sparseMatrix の場合は Matrix::dgCMatrix に変換する
    ### 違いはわからない
    ### dgCMatrix: csc順に並び替えて(csc形式)の疎行列圧縮保管
    ### https://ryamada.hatenadiary.jp/entry/20150310/1425944605
    is.sparse = FALSE
    ix = jx = NULL
    if (inherits(x, "sparseMatrix")) {
        is.sparse = TRUE
        x = as(x, "CsparseMatrix")
        x = as(x, "dgCMatrix")
        ### x@p は各列の非ゼロの値の個数を積み上げたものが格納されている（列数 + 1）
        ### diff(x@p + 1) すれば各列の非ゼロの値の個数がわかる
        ix = as.integer(x@p + 1)
        ### x@i は各列の非ゼロの値の行番号が格納されている（なので length(x@i) が非ゼロの値の個数と一致する）
        ### 0-index なので R のスタイルと合わせるために +1 しているのでしょう
        jx = as.integer(x@i + 1)
        ### x@x は非ゼロである値そのもののベクトル
        x = as.double(x@x)
    }

    ## プログレスバー
    if (trace.it) {
        if (relax) 
            cat("Training Fit\n")
        pb <- createPB(min = 0, max = nlam, initial = 0, style = 3)
    }

    ## 最適化の手法（ロジスティックおよび多項ロジスティックの時）
    ### type.logistic = c("Newton", "modified.Newton")
    ### Newton を指定なら 0、modified.Newton を指定なら 1 を返す
    ### If "Newton" then the exact hessian is used (default), while "modified.Newton" uses an upper-bound on the hessian, and can be faster.
    kopt = switch(match.arg(type.logistic), Newton = 0, modified.Newton = 1)
    ### type.multinomial = c("ungrouped", "grouped")
    ### 多項ロジスティックで更にgroupedの場合は kopt は 2 となる
    ### If "grouped" then a grouped lasso penalty is used on the multinomial coefficients for a variable. This ensures they are all in our out together. 
    ### The default is "ungrouped"
    if (family == "multinomial") {
        type.multinomial = match.arg(type.multinomial)
        if (type.multinomial == "grouped") 
            kopt = 2
    }
    kopt = as.integer(kopt)


    # フィッティング
    ## family に応じてその後に呼び出す関数を変える
    fit = switch(family,
                    ### gaussian のときは elnet 
                    gaussian = elnet(x, is.sparse, ix, jx, 
        y, weights, offset, type.gaussian, alpha, nobs, nvars, 
        jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, isd, intr, 
        vnames, maxit), 
                    ### poisson のときは fishnet
                    poisson = fishnet(x, is.sparse, ix, jx, 
        y, weights, offset, alpha, nobs, nvars, jd, vp, cl, ne, 
        nx, nlam, flmin, ulam, thresh, isd, intr, vnames, maxit),
                    ### binomial のときは lognet
                    binomial = lognet(x, is.sparse, ix, jx, y, weights, offset, 
            alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin, 
            ulam, thresh, isd, intr, vnames, maxit, kopt, family), 
                    ### multinomial のときも lognet
                    multinomial = lognet(x, is.sparse, ix, jx, y, weights, 
            offset, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, 
            flmin, ulam, thresh, isd, intr, vnames, maxit, kopt, 
            family), 
                    ### cox のときは coxnet
                    cox = coxnet(x, is.sparse, ix, jx, y, weights, 
            offset, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, 
            flmin, ulam, thresh, isd, vnames, maxit), 
                    ### mgaussian のときは mrelnet
                    mgaussian = mrelnet(x, 
            is.sparse, ix, jx, y, weights, offset, alpha, nobs, 
            nvars, jd, vp, cl, ne, nx, nlam, flmin, ulam, thresh, 
            isd, jsd, intr, vnames, maxit))
    ## プログレスバー
    if (trace.it) {
        utils::setTxtProgressBar(pb, nlam)
        close(pb)
    }
    
    # 後処理
    ## lambda が指定されておらず fit$lambda が 3 パターン以上検証されている場合、先頭を差し替える
    ## glmnet::fix.lam
    ## function (lam) {
    ## if (length(lam) > 2) {
    ##     llam = log(lam)
    ##     lam[1] = exp(2 * llam[2] - llam[3])
    ## }
    ## lam
    ## }
    if (is.null(lambda)) 
        fit$lambda = fix.lam(fit$lambda)
    ## call
    fit$call = this.call
    ## レコード数
    fit$nobs = nobs
    ## class に glmnet を追加
    class(fit) = c(class(fit), "glmnet")
    
    # リターン
    ## relax が TRUE の場合、解パスの各セットについて罰則なしでモデルをフィッティングする   
    ## If TRUE then for each active set in the path of solutions, the model is refit without any regularization. See details for more information. 
    ## This argument is new, and users may experience convergence issues with small datasets, especially with non-gaussian families. 
    ## Limiting the value of ’maxp’ can alleviate these issues in some cases.
    if (relax) 
        relax.glmnet(fit, x = x, y = y, weights = weights, offset = offset, 
            lower.limits = lower.limits, upper.limits = upper.limits, 
            check.args = FALSE, ...)
    else fit
}