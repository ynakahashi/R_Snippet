## plot time series
plot(as.zoo(as.ts(
   cbind(datLrn$TARGET, fitted(mod)))),
   plot.type = "single",
   col=1:2)

## Seasonal decomposition of time series panel data
datTmpAll <- c()
for(Car in model_name){
   for(Reg in region_name){
      for (Cust in customer) {

         ## select target carmodel, region, variable
         obj_name <- paste0("NUM_ORDERS_", Cust)
         datSel   <- datOrd_Target %>%
            filter(CAR.MODEL == Car, Region == Reg) %>%
            select_(obj_name, "YEARWEEK") %>% # use "select_" instead of "select"
            rename_("TARGET" = obj_name)      # use "rename_" instead of "rename"
         
         ## run stl
         datSTL    <- stl(ts(datSel$TARGET, frequency = 52),
                          s.window = "periodic")$time.series %>% data.frame()
         
         ## rbind
         datTmpAll <- datSTL %>%
            mutate("CarModel" = Car, "Region" = Reg, "Customer" = Cust) %>%
            mutate("YEARWEEK" = datSel$YEARWEEK) %>%
            bind_rows(datTmpAll)
      }
   }
}

## calculate Xb in MMM
n   <- 7
NM  <- tail(unlist(predVars), n)
X   <- as.matrix(datLrn[, NM])
b   <- summary(mod)$coef[which(rownames(summary(mod)$coef) %in% NM)]
B   <- matrix(rep(b, each=nrow(datLrn)), nrow(datLrn))
Mat <- X * B
head(exp(Mat), 20)


## calculate correlations of 1 against all
resCor <- c()
myCor <- function(vec) {
   if (is.numeric(vec)) {
      return(cor(datLrn$TARGET, vec))
   } else {
      return()
   }
}
resCor <- as.data.frame(do.call("rbind", lapply(datLrn, myCor)))
resCor$Var <- row.names(resCor)
write.csv(resCor %>% arrange(desc(V1)), file="cor.csv")


## glmnet
y <- as.matrix(datLrn$TARGET)
x <- as.matrix(datLrn[, unlist(predVars)])
resgn <- glmnet(x=x, y=y, alpha=1)
plot(resgn, label = T, xvar="norm")
plot(resgn, xvar="lambda", label=T)

resgn <- glmnet(x=x, y=y, alpha=-0)
plot(resgn, label = T, xvar="norm")
plot(resgn, xvar="lambda", label=T)

rescv <- cv.glmnet(x=x, y=y, alpha=0, nfolds = 10)
plot(rescv)
resfix <- glmnet(x, y, alpha=0, lambda = rescv$lambda.min)


## Explolatory Data Analysis via Random Forest
library(edarf)
library(randomForest)
datRF <- datLrn[, c("TARGET", unlist(predVars))]
fit   <- randomForest(TARGET ~ ., datRF)
imp   <- variable_importance(fit, data = datRF, vars = names(datRF)[-1])
plot_imp(imp)
do.call("rbind", imp)


## Extract main columns
out <- datLrn[, c("YEARWEEK", "MONTH1", 
                  "NUM_ORDERS_ALL", "LOGSALES", 
                  "TVGRP_CarModel.LOG", "TVGRP_CarModel.ADSTOCK",
                  "DIGITAL_CARS_TOTAL.LOG", "DIGITAL_CARS_TOTAL.ADSTOCK",
                  "DIGITAL_CARS_TOTAL",
                  "VME_CCL_Intrade_Op_MileageUP_NSXJE")]
Decomp <- as.data.frame(stl(ts(datLrn$LOGSALES, frequency = 52), 
                            "periodic")$time.series)
out$AdjLogOrder <- Decomp$trend + Decomp$remainder
write.csv(out, paste0(progdir, "MajorColumns.csv"), row.names=F)


## Plot in single view using zoo
plot(zoo::as.zoo(as.ts(
   cbind(datLrn$TARGET, fitted(mod)))), plot.type = "single", col=1:2)

plot(zoo::as.zoo(as.ts(
   cbind(datLrn$TARGET, datLrn$CarLife_Mix))), plot.type = "single", col=1:2)

plot(zoo::as.zoo(as.ts(
   cbind(datLrn$TaxRateUp, datLrn$CarLife_Mix_02))), plot.type = "single", col=1:2)

plot(zoo::as.zoo(as.ts(
   cbind(datLrn$TARGET, datLrn$MASS3_Carmodel.LOG))), plot.type = "single", col=1:2)

plot(zoo::as.zoo(as.ts(
   cbind(datLrn$TARGET, datLrn$DIGITAL_CARS_TOTAL.ADSTOCK))), plot.type = "single", col=1:2)


## Pick up decay rate
datSpd %>% 
   filter(Region == model_region, CarModel == model_car, TacticType == "FMI",
          TacticName_ROICalc_eng == "Digital")



