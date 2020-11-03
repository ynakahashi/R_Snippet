## Sample data
set.seed(123)
n <- 1000 # レコード数
ni <- 15  # 説明変数の数
nx <- ni  # nx は本来非ゼロとする変数の数だけど、ここを実装するには nin, mm も必要になってくるのでパス
alpha <- 1 # L1 / L2 それぞれに対する罰則の配分を決定するパラメータ

b0 <- rnorm(ni) # 回帰係数
b <- ifelse(abs(b0) < 0.5, 0, b0) # Lasso っぽくするために不要な説明変数を作る
x <- matrix(rnorm(n*ni), nrow = n, ncol = ni)
# x <- matrix(0, nrow = n, ncol = ni) # 説明変数ごとに N(0, 1) にしたかったので列ごとに生成
# for (i in 1:ni) {
#   x[, i] <- rnorm(n)
# }
y <- x %*% b + rnorm(n)


## Set lambda
lambda <- c(9.9E+2, 0, rev(seq(0, 5, length.out = 100))) # Warm start
nlam <- length(lambda)

a <- rep(0, ni)
rsqo <- rep(0, nlam)
almo <- rep(0, nlam)

## --- ここは Standard ----
xv <- colSums(x^2) / n
g <- g0 <- t(y) %*% x / n
vp <- rep(1, ni)

## --- ここから elnet1 ----
bta <- alpha
omb <- 1-bta
alm <- 0.0
alf <- 1.0

eqs <- max(0, 1e-04)
alf <- eqs^(1.0/(nlam-1))

rsq <- 0.0
nlp <- 0
mm <- rep(0, ni)
ia <- rep(0, ni)

c <- matrix(0, nrow = ni, ncol = nx)

# 1番目のループ
for (m in 1:nlam) {
  
  alm <- lambda[m] # このあと goto 10291
  
  # if(m == 1) {
  #   alm <- 9E+10
  # } else if(m == 2) {
  #   alm <- 0
  # } else {
  #   alm <- alf*alm
  # }
  # 
  # # 2番目のループ
  # for (j in 1:ni) {
  #   # ここで alm を update
  #   alm <- max(alm, abs(g[j])/vp[j])
  # }
  
  dem <- alm * omb # ここから 10291
  ab <- alm * bta
  rsq0 <- rsq
  
  # while(dlx >= thr & nin <= nx) {
  nlp <- nlp + 1
  dlx <- 0.0
  
  # 3番目のループ
  # 全変数に対する回帰係数の推定
  for (k in 1:ni) {
    ak <- a[k]
    u <- g[k] + ak * xv[k]
    v <- abs(u[1]) - ab
    # print(v)
    
    a[k] <- 0
    if(v > 0) a[k] <- (sign(u) * abs(v))/(xv[k] + vp[k]*dem)
    
    nin <- k # 本当は if(mm[k] == 0) nin <- nin + 1
    
    # 4番目のループ
    # 分散共分散行列を作る
    for (j in 1:ni) {
      if(j != k) {
        c[j, nin] <- (t(x[, j]) %*% x[, k]) / n
      } else {
        c[j, nin] <- xv[j]
      }
    }
    
    mm[k] <- nin
    ia[nin] <- k
    del <- a[k] - ak
    rsq <- rsq + del*(2.0 * g[k] - del*xv[k])
    
    # 5番目のループ
    for (j in 1:ni) {
      g[j] <- g[j] - c[j, mm[k]]*del
    }
  }
  
  
  # while(dlx >= thr & nlp <= maxit) {
  # 6番目のループ
  for (l in 1:nin) {
    k <- ia[l]
    ak <- a[k]
    u <- g[k] + ak*xv[k]
    v <- abs(u) - vp[k]*ab
    a[k] <- 0.0
    if(v > 0.0) a[k] <- (sign(u) * abs(v))/(xv[k]+vp[k]*dem)
    if(a[k] == ak) next
    del <- a[k] - ak
    rsq <- rsq + del*(2.0*g[k] - del*xv[k])
    dlx <- max(xv[k]*del^2, dlx)
    
    # 7番目のループ
    for (j in 1:nin) {
      g[ia[j]] <- g[ia[j]] - c[ia[j], mm[k]]*del
    }
  }
  # }
  
  # 8番目のループ
  # for (j in 1:ni) {
  #   g[j] <- g[j] - dot_product(da(1:nin),c(j,1:nin))       
  # }
  # }
  
  rsqo[m] <- rsq
  almo[m] <- alm
  
  # 9番目のループ
  # for (j in 1:nin) {
  #   if(ao[j,m] != 0.0) me <- me+1
  # }
  
}

oldpar <- par(no.readonly = T)
par(mfrow = c(1, 3))
plot(rsqo, xlab = "Times", ylab = "Residual Squared", type = "b")
plot(log(lambda), rsqo, xlab = "Log-Lambda", ylab = "Residual Squared.", type = "b")
plot(a, b, xlab = "Estimated Coef.", ylab = "True Coef.")
par(oldpar)
cbind(a, b)
