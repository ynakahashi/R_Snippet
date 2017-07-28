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

## Random forest by ranger
library(ranger)
resRF <- ranger(Species ~ ., data = iris, mtry = 2, num.trees = 500, 
                write.forest = TRUE)
table(predict = predict(resRF, data = iris)$predictions, 
      true    = iris$Species)

## Random forest by randomForest
library(randomForest)
tuneRF(dat[,-8], dat[,8], doBest = TRUE)
resRF <- randomForest(cv ~., dat, mtry = 2)
importance(resRF)

## Decision tree
library(rpart)
Titanic.rpart <- rpart(Survived ~ Sex + Age + Class, data = tat)
library(rpart.plot) 
prp(Titanic.rpart, type = 2, extra = 101,
    nn = TRUE, fallen.leaves = TRUE, faclen = 0, varlen = 0,
    shadow.col = "grey", branch.lty = 3, cex = 1.2, split.cex = 1.2,
    under.cex = 1.2)



## Cox Propotional-Hazard
library(MASS)
library(survival)
res.cox <- coxph(Surv(time, cens) ~ treat, data = gehan)
res.fit <- survfit(res.cox)
res.zph <- cox.zph(res.cox)
scatter.smooth(residuals(res.cox, type="deviance"))
scatter.smooth(residuals(res.cox))
colnames(res.out) <- c("time", "n.risk", "n.event", "n.censor", 
                       "surv", "cumhaz", "std_err", "upper", "lower")

res.stat <- rbind(summary(res.cox)$rsq[1], summary(res.cox)$concordance[1])
res.stat.names <- c("R_square", "Concordance_AUC")
res.stat       <- cbind(res.stat.names, res.stat)


################################################################################
### Time-Series specific
################################################################################
## Seasonal decomposition by STL
datSTL <- stl(ts(dat$TARGET, frequency = 52),
              s.window = "periodic")$time.series %>% data.frame()



################################################################################
### Data visualization
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
fit   <- randomForest(Species ~ ., iris)
imp   <- variable_importance(fit, data = iris, vars = names(iris)[-5])
plot_imp(imp)

## MIC
library(minerva)
mine(iris[, -5])$MIC

## Pivot table
library(rpivotTable)
rpivotTable(iris, rows="Species", "100%", "40pt")

## glance at data
glimpse(iris)

################################################################################
### Data wrangling
################################################################################
## Data expansion
library(epitools)
tat <- expand.table(Titanic) %>% tbl_df %>% 
   mutate("Survived" = if_else(Survived == "Yes", 1, 0))
rpivotTable(tat, rows = "Sex", "100%", "40pt")


## Select, Rename, Filter, etc...




################################################################################
### Make tidy data
################################################################################
## lm 01
library(broom)
library(gapminder)
gapminder %>%
   group_by(country) %>%
   do(data = lm(lifeExp ~ year, data = .) %>% tidy()) %>%
   unnest()

## lm 02
gapminder %>%
   split(.$country) %>%
   map(., ~ lm(lifeExp ~ year, data = .) %>% tidy()) %>%
   bind_rows(.id = "country") %>%
   tbl_df()

