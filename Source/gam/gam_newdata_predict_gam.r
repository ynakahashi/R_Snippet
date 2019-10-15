### newdata.predict.gam
### missing(newdata)がFALSE（新規データに対する予測）の場合はこちらが呼ばれる
function (object, newdata, type = c("link", "response", "terms"), 
    dispersion = NULL, se.fit = FALSE, na.action = na.pass, terms = labels(object), 
    ...) 
{
    out.attrs <- attr(newdata, "out.attrs")
    ### Gam クラスを継承し、かつ平滑化が行われている場合（Gam使うなら基本はこちらとなるはず）
    is.Gam <- inherits(object, "Gam") && !is.null(object$smooth)
    if (is.Gam) {

        if (se.fit) {
            se.fit <- FALSE
            ### 新規データに対するSEは出せない
            warning("No standard errors (currently) for gam predictions with newdata")
        }

        type <- match.arg(type)
        local.type <- type
        if (type == "response") 
            local.type <- "link"

        ### まずは predcit.glm で予測値を求める
        pred <- predict.glm(object, newdata, type = local.type, 
            dispersion = dispersion, se.fit = FALSE, terms = terms)

        ### モデルオブジェクトの term から目的変数を消去する
        tt <- terms(object)
        Terms <- delete.response(tt)

        ### 使用している変数で model.frame を作成する
        smooth.frame <- model.frame(Terms, newdata, na.action = na.action, 
            xlev = object$xlevels)
        

        nrows <- nrow(smooth.frame)

        ### モデルオブジェクトの smooth を取り出す（平滑化された説明変数の列を持つmatrix）
        old.smooth <- object$smooth
        ### モデルオブジェクトの smooth.frame（平滑化された説明変数の列を持つdata.frame）
        data <- object$smooth.frame
        smooth.labels <- names(data)

        ### 平滑化された説明変数の数
        n.smooths <- length(smooth.labels)

        ### 説明変数の型の検査
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, smooth.frame)
        out.attrs <- attr(newdata, "out.attrs")

        ### データの重み
        w <- object$weights

        ### newdataの行数、平滑化された説明変数の列数の0ベクトルを用意し、行と列それぞれにラベルを与える
        pred.s <- array(0, c(nrows, n.smooths), list(row.names(smooth.frame), 
            smooth.labels))
        
        ### 平滑化された説明変数だけを取り出す
        smooth.wanted <- smooth.labels[match(smooth.labels, terms, 
            0) > 0]
        pred.s <- pred.s[, smooth.wanted, drop = FALSE]

        ### モデルオブジェクトの残差
        residuals <- object$residuals

        ### 平滑化された説明変数の数だけループを回す
        for (TT in smooth.wanted) {

            ### object$smooth.frameの call attributeを取り出す
            Call <- attr(data[[TT]], "call")

            ### Call にオプションを追加する
            Call$xeval <- substitute(smooth.frame[[TT]], list(TT = TT))

            ### 残差に対象の平滑化された説明変数を加算する
            z <- residuals + object$smooth[, TT]

            ### データの各点の重み w と、上で作った z を渡して gam.s
            pred.s[, TT] <- eval(Call)
        }
        if (type == "terms") 
            pred[, smooth.wanted] <- pred[, smooth.wanted] + 
                pred.s[, smooth.wanted]
        else pred <- pred + rowSums(pred.s)
        if (type == "response") {
            famob <- family(object)
            pred <- famob$linkinv(pred)
        }
    }
    else {
        pred <- predict.glm(object, newdata, type = type, dispersion = dispersion, 
            se.fit = se.fit, terms = terms)
    }
    if (type != "terms" && !is.null(out.attrs)) {
        if (!is.null(out.attrs)) {
            if (se.fit) {
                attributes(pred$fit) <- out.attrs
                attributes(pred$se.fit) <- out.attrs
            }
            else attributes(pred) <- out.attrs
        }
    }
    pred
}