### load library
library(tidyverse)
library(topicmodels)

### load trained model
load("/Users/yoshinobu.nakahashi/EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/topicmodel_result_quant/2021-12-11-00-08-59_Obj.Rdata")

### function を読み込む
source("/Users/yoshinobu.nakahashi/EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/script/function_for_lda.R")

### 分析対象4,000件の df_mcb を作る
df_mcb <- create_mecab_df("./Desktop/Aflac_data/result/quan_data/量重視_テキスト_1210_300sec_OP.csv")

### df_lda を作る
df_lda <- create_df(df_mcb, limit_num = 50, 
                    term_type = c("名詞", "代名詞", "動詞", "形容詞", "形容動詞",
                                  "副詞", "連体詞", "接続詞", "感動詞",
                                  "記号"))

### 行動変容したデータに絞る
cntr_file <- "./Desktop/Aflac_data/result/quan_data/量重視_行動変容_1210_300sec_OP_201906.csv"
cntr_id <- read_csv(cntr_file)$r_call_id

df_lda_cntr <- df_lda$df[rownames(df_lda$df) %in% cntr_id, ]

## 
st <- Sys.time()
res_mat <- estimate_topic_in_call_02(out$res_lda, df_lda_cntr, df_mcb, 
                                     window = 3, sep = 10, num_call = nrow(df_lda_cntr))
ed <- Sys.time()
ed - st

# save(list = c("df_lda", "df_mcb"), 
#      file = "./EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/topicmodel_result_quant/2021-1217_TM_Div.Rdata")


df_topic_in_call <- as_tibble(res_mat) %>% 
   inner_join(df_cntr, by = c("V1" = "r_call_id")) %>% 
   rename(Call_ID = V1,
          Topic_at_1 = V2,
          Topic_at_2 = V3,
          Topic_at_3 = V4,
          Topic_at_4 = V5,
          Topic_at_5 = V6,
          Topic_at_6 = V7,
          Topic_at_7 = V8,
          Topic_at_8 = V9) %>% 
   select(Call_ID, starts_with("is_"), everything())
save(list = c("df_lda", "df_mcb", "df_topic_in_call"),
     file = "./EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/topicmodel_result_quant/2021-1217_TM_Div.Rdata")


load("./EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/topicmodel_result_quant/2021-1217_TM_Div.Rdata")
df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_3 == 32, Topic_at_4 == 32)

df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_3 == 32, Topic_at_4 == 1)


df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_1 == 9, Topic_at_2 == 9)


df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_1 == 6, Topic_at_2 == 19)

df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_1 == 33, Topic_at_4 == 33)


df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_1 == 33, Topic_at_4 == 33)

df_topic_in_call %>% 
   select(-2) %>% 
   filter(Topic_at_1 == 6)



trans_mat <- create_trans_mat(res_mat)
trans_mat <- as.data.frame(trans_mat)



write_excel_csv(trans_mat, file = "./EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/topicmodel_result_quant/2021-1217_TransMat.csv")


# write_excel_csv(trans_mat, paste0(out_file_path, "trans_mat_15_100.csv"))
