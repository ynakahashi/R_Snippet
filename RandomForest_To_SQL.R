library(randomForest)
set.seed(123)
res <- randomForest(Species ~ ., data = iris, ntree = 2)

# for (i in 1:(res$ntree)) {
#    Rule_Table <- getTree(res, k = i, labelVar = TRUE)
#    cat(paste0("INSERT INTO tbl_rn \nSELECT"))
#    Make_SQL <- function(res, Rule_Table, cnt, ind = 1) {
#       Rule   <- Rule_Table[cnt, ]
#       Indent <- paste(rep("   ", ind), collapse = "")
#       var    <- as.character(Rule[, "split var"])
#       val    <- Rule[, "split point"]
#       
#       if(Rule[, "status"] != -1) {
#          cat(paste0("\n", Indent, "CASE\n", Indent, "   WHEN ", var, " <= ", val, " THEN"))
#          Make_SQL(res, Rule_Table, Rule[, "left daughter"], ind = (ind + 2))
#          cat(paste0("\n", Indent, "   ELSE"))
#          Make_SQL(res, Rule_Table, Rule[, "right daughter"], ind = (ind + 2))
#          cat(paste0("\n", Indent, "END"))
#       } else { 
#          cat(paste0(" '", Rule[, "prediction"], "'"))
#       }
#    }
#    Make_SQL(res, Rule_Table, 1, 1)
#    cat(paste0(" as tree", i, " \nFROM \n", "   input_data", ";\n\n"))
# }


## Create recursive-function 
Make_SQL <- function(Rule_Table, cnt, ind = 1) {
   Rule   <- Rule_Table[cnt, ]
   Indent <- paste(rep("   ", ind), collapse = "")
   var    <- as.character(Rule[, "split var"])
   val    <- Rule[, "split point"]
   
   if(Rule[, "status"] != -1) {
      cat(paste0("\n", Indent, "CASE\n", Indent, "   WHEN ", var, " <= ", val, " THEN"))
      Make_SQL(Rule_Table, Rule[, "left daughter"], ind = (ind + 2))
      cat(paste0("\n", Indent, "   ELSE"))
      Make_SQL(Rule_Table, Rule[, "right daughter"], ind = (ind + 2))
      cat(paste0("\n", Indent, "END"))
   } else { 
      cat(paste0(" '", Rule[, "prediction"], "'"))
   }
}

for (i in 1:(res$ntree)) {
   Rule_Table <- getTree(res, k = i, labelVar = TRUE)
   cat(paste0("INSERT INTO tbl_rf \nSELECT\n   id, "))
   Make_SQL(Rule_Table, 1, 1)
   cat(paste0(" as tree", i, "\nFROM\n", "   input_data", ";\n\n"))
}


cat(paste0("INSERT INTO rf_predictions\n",
           "SELECT\n   a.id,\n   a.pred\n",
           "FROM\n   (\n      SELECT id as id, pred, COUNT(*) as cnt \n",
           "      FROM tbl_rf\n      GROUP BY id, pred\n   ) a\n",
           "   INNER JOIN\n   (\n      ",
           "SELECT id, MAX(cnt) as cnt\n",
           "      FROM\n",
           "         (\n", 
           "            SELECT id as id, pred, COUNT(*) as cnt\n", 
           "            FROM tbl_rf\n",
           "            GROUP BY id, pred\n", 
           "         )\n",
           "      GROUP BY id\n   ) b\n",
           "      ON a.id = b.id AND a.cnt = b.cnt;\n\n"))






