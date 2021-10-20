function (formula, data = NULL, family = gaussian, control = glmerControl(), 
    start = NULL, verbose = 0L, nAGQ = 1L, subset, weights, na.action, 
    offset, contrasts = NULL, mustart, etastart, devFunOnly = FALSE) 
{
    if (!inherits(control, "glmerControl")) {
        if (!is.list(control)) 
            stop("'control' is not a list; use glmerControl()")
        if (class(control)[1] == "lmerControl") {
            warning("please use glmerControl() instead of lmerControl()", 
                immediate. = TRUE)
            control <- c(control[!names(control) %in% c("checkConv", 
                "checkControl")], control$checkControl, control$checkConv)
            control["restart_edge"] <- NULL
        }
        else {
            msg <- "Use control=glmerControl(..) instead of passing a list"
            if (length(cl <- class(control))) {
                msg <- paste(msg, "of class", dQuote(cl[1]))
            }
            warning(msg, immediate. = TRUE)
        }
        control <- do.call(glmerControl, control)
    }
    mc <- mcout <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame(2))
    if (is.function(family)) 
        family <- family()
    if (isTRUE(all.equal(family, gaussian()))) {
        warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;", 
            " please call lmer() directly")
        mc[[1]] <- quote(lme4::lmer)
        mc["family"] <- NULL
        return(eval(mc, parent.frame()))
    }
    mc[[1]] <- quote(lme4::glFormula)
    glmod <- eval(mc, parent.frame(1L))
    mcout$formula <- glmod$formula
    glmod$formula <- NULL
    nAGQinit <- if (control$nAGQ0initStep) 
        0L
    else 1L
    devfun <- do.call(mkGlmerDevfun, c(glmod, list(verbose = verbose, 
        control = control, nAGQ = nAGQinit)))
    if (nAGQ == 0 && devFunOnly) 
        return(devfun)
    if (is.list(start)) {
        start.bad <- setdiff(names(start), c("theta", "fixef"))
        if (length(start.bad) > 0) {
            stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s", 
                paste(start.bad, collapse = ", "), shQuote("theta"), 
                shQuote("fixef")), call. = FALSE)
        }
        if (!is.null(start$fixef) && nAGQ == 0) 
            stop("should not specify both start$fixef and nAGQ==0")
    }
    if (control$nAGQ0initStep) {
        opt <- optimizeGlmer(devfun, optimizer = control$optimizer[[1]], 
            restart_edge = if (nAGQ == 0) 
                control$restart_edge
            else FALSE, boundary.tol = if (nAGQ == 0) 
                control$boundary.tol
            else 0, control = control$optCtrl, start = start, 
            nAGQ = 0, verbose = verbose, calc.derivs = FALSE)
    }
    if (nAGQ > 0L) {
        devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
        if (control$nAGQ0initStep) {
            start <- updateStart(start, theta = opt$par)
        }
        if (devFunOnly) 
            return(devfun)
        opt <- optimizeGlmer(devfun, optimizer = control$optimizer[[2]], 
            restart_edge = control$restart_edge, boundary.tol = control$boundary.tol, 
            control = control$optCtrl, start = start, nAGQ = nAGQ, 
            verbose = verbose, stage = 2, calc.derivs = control$calc.derivs, 
            use.last.params = control$use.last.params)
    }
    cc <- if (!control$calc.derivs) 
        NULL
    else {
        if (verbose > 10) 
            cat("checking convergence\n")
        checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, 
            lbound = environment(devfun)$lower)
    }
    mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr, 
        mc = mcout, lme4conv = cc)
}