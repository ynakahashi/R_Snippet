## plot time series
plot(as.zoo(as.ts(
   cbind(datLrn$TARGET, fitted(mod)))),
   plot.type = "single",
   col=1:2)