

library(partykit)
set.seed(290875)
airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq,
               controls = ctree_control(maxsurrogate = 3))
### distribution of responses in the terminal nodes
plot(airq$Ozone ~ as.factor(where(airct)))
plot(airct)

res_ct <- ctree(Sepal.Length ~ Species, iris)
plot(res_ct)

dat <- epitools::expand.table(Titanic)

res_ct <- ctree(Survived ~ ., dat, control = ctree_control(multiway = TRUE))
plot(res_ct)



set.seed(1234)
x <- rnorm(20)
df <- data.frame(x = x,
                 y = x + rnorm(20))

plot(y ~ x, data = df)

# model
mod <- lm(y ~ x, data = df)

# predicts + interval
newx <- seq(min(df$x), max(df$x), length.out=100)
preds <- predict(mod, newdata = data.frame(x=newx), 
                 interval = 'confidence')

# plot
plot(y ~ x, data = df, type = 'n')
# add fill
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
# model
abline(mod)
# intervals
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')