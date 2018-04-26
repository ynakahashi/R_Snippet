
x1 <- rnorm(100)
x2 <- rnorm(100)
y  <- ifelse(x1 + x2 > rnorm(100, 0, 0.3), 1, 0)
dat <- data.frame(y, x1, x2)

plot(dat$x1, dat$x2, col = c(1, 6)[as.factor(dat$y)], bg = c(1, 6)[as.factor(dat$y)])

tmp1 <- glm(y ~ ., data = dat, family = binomial(link = "logit"))
summary(tmp1)

tmp2 <- glm(y ~ x1:x2, data = dat, family = binomial(link = "logit"))
summary(tmp2)

dat$x3 <- dat$x1 + dat$x2
tmp3 <- glm(y ~ x3, data = dat, family = binomial(link = "logit"))
summary(tmp3)

table(
   dat$y,
   ifelse(fitted(tmp1) > 0.5, 1, 0))
table(
   dat$y,
   ifelse(fitted(tmp2) > 0.5, 1, 0))
table(
   dat$y,
   ifelse(fitted(tmp3) > 0.5, 1, 0))



