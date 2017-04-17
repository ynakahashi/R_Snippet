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
