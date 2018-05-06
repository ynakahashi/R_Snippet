
n <- 100
q <- 1

set.seed(123)
x <- runif(100)
e <- rnorm(100, 0, 1)
y <- 3 + x * 1.5 + e

res <- lm(y ~ x)

RSS <- sum((y - fitted(res))^2)
AIC_cal <- n * log(2 * pi) + n * log(RSS / n) + n + 2 * q + 4

AIC(res)
AIC_cal

sig1 <- tmp$sigma^2
sig2 <- sd(resid(tmp))
loglik1 <- (-n/2) * log(2 * pi) - (n/2) * log(sig1) - (n/2)
loglik2 <- (-n/2) * log(2 * pi) - (n/2) * log(sig2) - (n/2)
val <- 0.5 * (n * (log(2 * pi) + 1 - log(n) + log(sum(resi^2))))
val <- 0.5 * (n * (log(2 * pi) + 1 - log(n) + log(sum(resi^2))))

logLik(res)
loglik1
loglik2
logLik(res, REML = TRUE)
val

stats:::logLik.lm

w <- rep(1, n)
N <- n
resi <- res$residuals
# if (REML) 
#    N <- N - p
val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * resi^2))))
# if (REML) 
#    val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))