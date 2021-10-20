function (devfun, optimizer = formals(lmerControl)$optimizer, 
    restart_edge = formals(lmerControl)$restart_edge, boundary.tol = formals(lmerControl)$boundary.tol, 
    start = NULL, verbose = 0L, control = list(), ...) 
{
    verbose <- as.integer(verbose)
    rho <- environment(devfun)
    opt <- optwrap(optimizer, devfun, getStart(start, rho$pp), 
        lower = rho$lower, control = control, adj = FALSE, verbose = verbose, 
        ...)
    if (restart_edge) {
        if (length(bvals <- which(rho$pp$theta == rho$lower)) > 
            0) {
            theta0 <- new("numeric", rho$pp$theta)
            d0 <- devfun(theta0)
            btol <- 1e-05
            bgrad <- sapply(bvals, function(i) {
                bndval <- rho$lower[i]
                theta <- theta0
                theta[i] <- bndval + btol
                (devfun(theta) - d0)/btol
            })
            devfun(theta0)
            if (any(bgrad < 0)) {
                if (verbose) 
                  message("some theta parameters on the boundary, restarting")
                opt <- optwrap(optimizer, devfun, opt$par, lower = rho$lower, 
                  control = control, adj = FALSE, verbose = verbose, 
                  ...)
            }
        }
    }
    if (boundary.tol > 0) 
        check.boundary(rho, opt, devfun, boundary.tol)
    else opt
}