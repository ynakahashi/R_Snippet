# library(discretization)
# library(ggplot2)
library(tidyverse)  # for data wrangling
library(infotheo)   # for discretization
library(ROCR)       # for model evaluation
library(pROC)       # for model evaluation
library(caret)      # for machine learning algorithm

Work_Dir <- "C:/Users/ynakahashi/Desktop/Work/"
File_NM  <- "credit_data_track2_part_A.CSV" 

Dat_Sample <- read_csv(paste0(Work_Dir, File_NM))
# glimpse(Dat_Sample)

Var_Char <- colnames(Dat_Sample)[sapply(Dat_Sample, is.character)]
Var_Num <- colnames(Dat_Sample)[sapply(Dat_Sample, is.numeric)]
# ncol(Dat_Sample)
# length(Var_Char) + length(Var_Num)

# for (i in Var_Char) {
#    print(i)
#    print(table(Dat_Sample[, i], useNA = "always"))
# }
# 
# for (i in Var_Num) {
#    print(i)
#    Dat_Sample %>% select_(i) %>% collect() %>% .[[1]] %>% hist()
#    readline()
# }
# 
# colSums(apply(Dat_Sample, c(1, 2), is.na))


Obj_Var <- "accepted"

Var_Cont <- colnames(Dat_Sample[, Var_Num])[
   lapply(Dat_Sample[, Var_Num], function(x) {length(table(x))}) > 4]
Dat_Cont  <- as.data.frame(Dat_Sample[, c(Var_Cont, Obj_Var)])
# Dat_Cont <- Dat_Sample %>% 
#    select_(.dots = c(Var_Cont, Obj_Var)) %>% 
#    mutate_each(funs("log" = log1p)) %>% 
#    select(-starts_with("accepted"), accepted)

# Res_ChiM <- chiM(Dat_Num, 0.05)

Dat_Disc <- discretize(Dat_Cont, disc="equalfreq", nbins = 5)
colnames(Dat_Disc) <- paste0(colnames(Dat_Cont), "_Cat")
# for (i in colnames(Dat_Disc)) {
#    tmp <- table(Dat_Disc[, c(i, "accepted_Cat")])
#    plot(tmp)
#    readline()
# }
# 
# for (i in Var_Char) {
#    tmp <- table(Dat_Sample[, c(i, Obj_Var)])
#    plot(tmp)
#    readline()
# }


## Remove NAs
## Remove anormaly on credit_amount
## Categorize less-recorded levels into Others on purpose


Dat_Ana <- Dat_Sample %>% 
   bind_cols(Dat_Disc) %>% 
   mutate("purpose_02" = if_else(purpose %in% c("domestic_appliances", 
                                                "repairs", "retraining"),
                                 "other", purpose)) %>% 
   mutate_each(funs("Char" = as.character), ends_with("Cat"),
               installment_commitment) %>% 
   select(-purpose) %>% 
   select(-c(duration, credit_amount, age, asnm, ends_with("Cat"), accepted,
             installment_commitment, accepted_Cat_Char), accepted) %>% 
   drop_na() 

# glimpse(Dat_Ana)

Dat_Mat <- model.matrix(~.-1, data = Dat_Ana)
Dat_Cor <- cor(Dat_Mat, use = "pairwise.complete")
# diag(solve(Dat_Cor))
# max(diag(MASS::ginv(Dat_Cor)))
# sapply(as.data.frame(lower.tri(Dat_Cor) * Dat_Cor), max) %>% 
#    as_data_frame() %>% arrange(desc(value))


set.seed(123)
p <-0.8
Train     <- sample(1:nrow(Dat_Mat), size = nrow(Dat_Ana) * p, replace = FALSE)
Dat_Train <- Dat_Ana[rownames(Dat_Mat) %in% Train, ] 
Dat_Test  <- Dat_Ana[!rownames(Dat_Mat) %in% Train, ] 


Res_GLM      <- glm(accepted ~ ., data = Dat_Train, family = binomial(link = logit))
Res_Step_GLM <- step(Res_GLM, direction = "both")
Dat_Test$Prd <- predict(Res_Step_GLM, type = "response", newdata = Dat_Test)
Dat_Test$Pred_Class <- ifelse(Dat_Test$Prd < 0.3266702, 0, 1)


## ROC
pred <- prediction(Dat_Test$Prd, Dat_Test$accepted)
perf <- performance(pred, "tpr", "fpr")
plot(perf)

Get_Optimal_CutOff <- function(perf, pred){
   cut.ind <- mapply(FUN = function(x, y, p){
      d <- (x - 0)^2 + (y - 1)^2
      ind <- which(d == min(d))
      c(sensitivity <- y[[ind]], 
        specificity <- 1 - x[[ind]], 
        cutoff      <- p[[ind]])
   }, perf@x.values, perf@y.values, pred@cutoffs)
}

print(Get_Optimal_CutOff(perf, pred))


## AUC
performance(pred, "auc")

## KS
max(attr(perf, 'y.values')[[1]] - attr(perf, 'x.values')[[1]])

## Get important variable
Preds <- predict(Res_Step_GLM, type = "term", newdata = Dat_Test)
Get_Important_Variable <- function(x, top = 3) {
   res <- names(x)[order(x, decreasing = TRUE)][1:top]
   paste(res, collapse = " / ", sep = "")
}

Var_Imp <- apply(Preds, 1, Get_Important_Variable, top = 3)
head(Var_Imp)


## evaluate using pROC
ROC_GLM <- roc(Dat_Test$accepted, Dat_Test$Prd)
plot.roc(ROC_GLM)
coords(ROC_GLM,
       x = "best", 
       ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"), 
       best.method = "closest.topleft")

AUC_GLM <- auc(ROC_GLM)
print(AUC_GLM)


## Try machine learning algorithm 
Res_GLM <- train(
   as.factor(accepted) ~ ., 
   data = Dat_Train, 
   method = "glm", 
   tuneLength = 4,
   preProcess = c('center', 'scale'),
   trControl = trainControl(method = "cv")
)
GLM_PredictClass <- predict(Res_GLM, Dat_Test)

Res_DT <- train(
   as.factor(accepted) ~ ., 
   data = Dat_Train, 
   method = "rpart", 
   tuneLength = 4,
   preProcess = c('center', 'scale'),
   trControl = trainControl(method = "cv")
)
DT_PredictClass <- predict(Res_DT, Dat_Test)

Res_RF <- train(
   as.factor(accepted) ~ ., 
   data = Dat_Train, 
   method = "rf", 
   tuneLength = 4,
   preProcess = c('center', 'scale'),
   trControl = trainControl(method = "cv")
)
RF_PredictClass <- predict(Res_RF, Dat_Test)

Res_NN <- train(
   as.factor(accepted) ~ ., 
   data = Dat_Train, 
   method = "nnet", 
   tuneLength = 4,
   preProcess = c('center', 'scale'),
   trControl = trainControl(method = "cv")
)
NN_PredictClass <- predict(Res_NN, Dat_Test)

Res_EN <- train(
   as.factor(accepted) ~ ., 
   data = Dat_Train, 
   method = "glmnet", 
   tuneLength = 4,
   preProcess = c('center', 'scale'),
   trControl = trainControl(method = "cv")
)
EN_PredictClass <- predict(Res_EN, Dat_Test)

## Model evaluation
data.frame(
   "GLM_Step" = confusionMatrix(Dat_Test$Pred_Class, Dat_Test$accepted)$byClass,
   "GLM" = confusionMatrix(GLM_PredictClass, Dat_Test$accepted)$byClass,
   "DT"  = confusionMatrix(DT_PredictClass, Dat_Test$accepted)$byClass,
   "RF"  = confusionMatrix(RF_PredictClass, Dat_Test$accepted)$byClass,
   "NN"  = confusionMatrix(NN_PredictClass, Dat_Test$accepted)$byClass,
   "EN"  = confusionMatrix(EN_PredictClass, Dat_Test$accepted)$byClass)



## Random forest
Dat_Mat_Train <- Dat_Mat[rownames(Dat_Mat) %in% Train, ] 
Dat_Mat_Test  <- Dat_Mat[!rownames(Dat_Mat) %in% Train, ] 
res_rf_tune <- tuneRF(Dat_Mat_Train[, -60], as.factor(Dat_Mat_Train[, 60]), 
                      doBest = T)
prd_rf <- predict(res_rf_tune, Dat_Mat_Test)
confusionMatrix(prd_rf, Dat_Test$accepted)


