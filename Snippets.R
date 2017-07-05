################################################################################
##
## R snippets
## written by Y.Nakahashi 
##
################################################################################
library(tidyverse)

################################################################################
### Machine Learning Algorithm, Analytics functions
################################################################################
## glmnet
library(glmnet)
resGLM <- glmnet(trn_x, iris$Species, "multinomial")
table(predict = predict(resGLM, newx = trn_x, type = "class", s = 0.01)[, 1],
      true = iris$Species)

## glmnet 02
y <- as.matrix(datLrn$TARGET)
x <- as.matrix(datLrn[, unlist(predVars)])
resgn  <- glmnet(x = x, y = y, alpha = 1)
rescv  <- cv.glmnet(x = x, y = y, alpha = 0, nfolds = 10)
resfix <- glmnet(x, y, alpha = 0, lambda = rescv$lambda.min)

## random forest by ranger
library(ranger)
resRF <- ranger(Species ~ ., data = iris, mtry = 2, num.trees = 500, 
                write.forest = TRUE)
table(predict = predict(resRF, data = iris)$predictions, 
      true    = iris$Species)

## random forest by randomForest
library(randomForest)
tuneRF(dat[,-8], dat[,8], doBest = TRUE)
resRF <- randomForest(cv ~., dat, mtry = 2)
importance(resRF)


################################################################################
### Time-Series specific
################################################################################
## Seasonal decomposition by STL
datSTL <- stl(ts(dat$TARGET, frequency = 52),
              s.window = "periodic")$time.series %>% data.frame()



################################################################################
### Plot
################################################################################
## plot time series in single view
plot(as.zoo(as.ts(
   cbind(datLrn$TARGET, fitted(mod)))), plot.type = "single", col=1:2)

## Scatter plot by ggplot2
library(ggplot2)
ggplot(datasaurus, aes(x = x, y = y)) +
   facet_wrap(~ dataset, nrow = 3) +
   geom_point()


################################################################################
### Exploratory Analysis
################################################################################
## Explolatory Data Analysis via Random Forest
library(edarf)
library(randomForest)
datRF <- datLrn[, c("TARGET", unlist(predVars))]
fit   <- randomForest(TARGET ~ ., datRF)
imp   <- variable_importance(fit, data = datRF, vars = names(datRF)[-1])
plot_imp(imp)


