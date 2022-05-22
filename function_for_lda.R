hankana2zenkana <- function(x) {
   # character変換 
   if (!is.character(x)){ x <- as.character(x) }
   
   # 濁点、半濁点文字の置換
   dh <- c("ｶﾞ","ｷﾞ","ｸﾞ","ｹﾞ","ｺﾞ","ｻﾞ","ｼﾞ","ｽﾞ","ｾﾞ","ｿﾞ","ﾀﾞ","ﾁﾞ","ﾂﾞ","ﾃﾞ",
           "ﾄﾞ","ﾊﾞ","ﾋﾞ","ﾌﾞ","ﾍﾞ","ﾎﾞ","ﾊﾟ","ﾋﾟ","ﾌﾟ","ﾍﾟ","ﾎﾟ")
   dz <- c("ガ","ギ","グ","ゲ","ゴ","ザ","ジ","ズ","ゼ","ゾ","ダ","ヂ","ヅ","デ","ド",
           "バ","ビ","ブ","ベ","ボ","パ","ピ","プ","ペ","ポ")
   for( i in 1:length(dz) ){ x <- gsub(dh[i],　dz[i],　x) }
   
   # 1bite文字の置換
   x <- chartr("ｱｲｳｴｵｶｷｸｹｺｻｼｽｾｿﾀﾁﾂﾃﾄﾅﾆﾇﾈﾉﾊﾋﾌﾍﾎﾏﾐﾑﾒﾓﾔﾕﾖﾗﾘﾙﾚﾛﾜｦﾝ｡｢｣､･ｦｧｨｩｪｫｬｭｮｯｰ"
               , "アイウエオカキクケコサシスセソタチツテトナニヌネノハヒフヘホマミムメモヤユヨラリルレロワヲン。「」、・ヲァィゥェォャュョッー"
               , x)
   x
}

kanji_num_2_num <- function(x) {
   x <- chartr("一二三四五六七八九", "123456789", x)
   # x <- chartr("１２３４５６７８９", "123456789", x)
   x
}

hiragana_1 <- c("ぁ","あ","ぃ","い","ぅ","う","ぇ","え","ぉ","お",
                "か","が","き","ぎ","く","ぐ","け","げ","こ","ご",
                "さ","ざ","し","じ","す","ず","せ","ぜ","そ","ぞ",
                "た","だ","ち","ぢ","っ","つ","づ","て","で","と","ど",
                "な","に","ぬ","ね","の",
                "は","ば","ぱ","ひ","び","ぴ","ふ","ぶ","ぷ","へ","べ","ぺ","ほ","ぼ","ぽ",
                "ま","み","む","め","も","ゃ","や","ゅ","ゆ","ょ","よ",
                "ら","ら゚","り","り゚","る","る゚","れ","れ゚","ろ","ろ゚",
                "ゎ","わ","わ゙","ゐ","ゐ゙","ゑ","ゑ゙","を","を゙","ん",
                "ゔ","ゕ","ゖ")
katakana_1 <- c("ァ","ア","ィ","イ","ゥ","ウ","ェ","エ","ォ","オ",
                "カ","ガ","キ","ギ","ク","グ","ケ","ゲ","コ","ゴ",
                "サ","ザ","シ","ジ","ス","ズ","セ","ゼ","ソ","ゾ",
                "タ","ダ","チ","ヂ","ッ","ツ","ヅ","テ","デ","ト","ド",
                "ナ","ニ","ヌ","ネ","ノ",
                "ハ","バ","パ","ヒ","ビ","ピ","フ","ブ","プ","ヘ","ベ","ペ","ホ","ボ","ポ",
                "マ","ミ","ム","メ","モ","ャ","ヤ","ュ","ユ","ョ","ヨ",
                "ラ","ラ゚","リ","リ゚","ル","ル゚","レ","レ゚","ロ","ロ゚",
                "ヮ","ワ","ヰ","ヱ","ヲ","ン","ヴ","ヵ","ヶ","ヷ","ヸ","ヹ","ヺ",
                "・","ー")

stop_words <- as.vector(unlist(stopwords::data_stopwords_marimo$ja))
stop_words <- unique(c(stop_words,
                       as.vector(read.csv("http://svn.sourceforge.jp/svnroot/slothlib/CSharp/Version1/SlothLib/NLP/Filter/StopWord/word/Japanese.txt")[,1])))
stop_words <- stop_words[-which(stop_words %in% c("都", "道", "府", "県", "市", "区", "町", "村"))]
stop_words <- c(stop_words, "．")
stop_words <- unique(c(stop_words,
                       as.vector(read.csv("/Users/yoshinobu.nakahashi/EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/script/stop_words.txt")[,1])))


emotional_trigger_path <- "EY/【Aflac】アウトバウンドコール(データサイエンスxBX) - ドキュメント/分析用チャネル/work/Nakahashi/Emotional_Trigger_v001.csv"
emotional_trigger <- read.csv(emotional_trigger_path)

create_mecab_df <- function(txt_file) {

   ### Load libraries
   library(tidyverse)
   library(RMeCab)
   
   ### Read data
   df <- read_csv(txt_file)
   colnames(df) <- c("Call_ID", "Text")
   
   ### 半角を全角に変換、漢数字を数字に変換、形態素解析を行う
   txt <- df %>% 
      mutate(Text_01 = gsub(" ", "", Text)) %>% # これを入れるかは要検討
      mutate(Text_02 = kanji_num_2_num(Text_01)) %>%
      mutate(Text_03 = stringi::stri_trans_nfkc(Text_02)) %>% 
      mutate(Text_04 = stringi::stri_trans_toupper(Text_03)) %>% 
      as.data.frame() %>% 
      RMeCabDF("Text_04", 1)
   
   ### dataframeに変換
   df1 <- purrr::pmap_df(list(nv = txt, id = df$Call_ID),
                         function(nv, id) {
                            tibble(id = id, word = nv, type = names(nv))
                         })
   
   return(df1)
}

create_df <- function(df1, limit_num = 10, 
                      term_type = c("名詞", "代名詞", "動詞", "形容詞", "形容動詞",
                                    "副詞", "連体詞", "接続詞", "感動詞", "フィラー",
                                    "記号", "助詞", "助動詞")) {
   
   ### Load libraries
   library(tidyverse)
   library(RMeCab)

   ### 必要な品詞を残し、頻度の合計値の列を追加
   df2 <- df1 %>% 
      filter(!grepl("^[ぁ-ゖァ-ー]{1}$", word)) %>% 
      filter(!word %in% stop_words) %>% 
      filter(type %in% term_type) %>% 
      group_by(word) %>% 
      mutate(num_ttl = n()) %>% 
      group_by(id, word) %>% 
      mutate(num_in_doc = n()) %>%
      # mutate(flag = 1) %>% 
      # group_by(word) %>% 
      # summarise(num_doc = sum(flag)) %>% 
      ungroup() %>% 
      left_join(df1 %>% select(id, word) %>% unique() %>% 
                   group_by(word) %>% summarize(num_doc = n()) %>% ungroup(),
                by = "word")
   
   ### 頻度を指定して単語を抽出する
   df3 <- df2 %>% 
      filter(num_ttl >= limit_num)
   
   ### 単語分布用
   df_tf <- df3 %>% 
      select(word, num_ttl, num_doc) %>% 
      unique() %>% 
      arrange(desc(num_ttl)) %>% 
      left_join(emotional_trigger %>% mutate(Flag = 1), 
                by = c("word" = "キーワード")) 
   
   ### 語数の確認
   print(paste0("Number of Terms is: ", nrow(df3)))
   print(paste0("Number of Unique-terms is: ", length(unique(df3$word))))
   
   ### データフレームを整形
   df4 <- df3 %>% 
      select(id, word, num_in_doc) %>%
      unique() %>% 
      pivot_wider(names_from = word, values_from = num_in_doc)
   
   
   df5 <- df4 %>% 
      select(-id) %>% 
      as.data.frame()
   df5[is.na(df5)] <- 0
   rownames(df5) <- df4$id
   
   out <- list()
   out$df <- df5
   out$df_tf <- df_tf
   out$limit_num <- limit_num
   out$term_type <- term_type
   return(out)
}   

create_df_divided <- function(df1, div_num = 10, limit_num = 10, 
                              term_type = c("名詞", "代名詞", "動詞", "形容詞", "形容動詞",
                                            "副詞", "連体詞", "接続詞", "感動詞", "フィラー",
                                            "記号", "助詞", "助動詞")) {
   
   ### Load libraries
   library(tidyverse)
   library(RMeCab)
   
   ### 必要な品詞を残し、頻度の合計値の列を追加
   df2 <- df1 %>% 
      filter(!grepl("^[ぁ-ゖァ-ー]{1}$", word)) %>% 
      filter(!word %in% stop_words) %>% 
      filter(type %in% term_type) %>% 
      group_by(id) %>% # 以下4行は同一コールを div_num で分ける
      mutate(num_term = n(), num_order = row_number()) %>% 
      ungroup() %>% 
      mutate(idx = ceiling(num_order / num_term * div_num)) %>% 
      mutate(id_ord = paste(id, idx, sep = "_")) %>% 
      group_by(word) %>% 
      mutate(num_ttl = n()) %>% 
      # group_by(id, word) %>% 
      group_by(id_ord, word) %>% 
      mutate(num_in_doc = n()) %>%
      # mutate(flag = 1) %>% 
      # group_by(word) %>% 
      # summarise(num_doc = sum(flag)) %>% 
      ungroup() %>% 
      left_join(df1 %>% select(id, word) %>% unique() %>% 
                   group_by(word) %>% summarize(num_doc = n()) %>% ungroup(),
                by = "word")
   
   ### 頻度を指定して単語を抽出する
   df3 <- df2 %>% 
      filter(num_ttl >= limit_num)
   
   ### 単語分布用
   df_tf <- df3 %>% 
      select(word, num_ttl, num_doc) %>% 
      unique() %>% 
      arrange(desc(num_ttl)) %>% 
      left_join(emotional_trigger %>% mutate(Flag = 1), 
                by = c("word" = "キーワード")) 
   
   ### 語数の確認
   print(paste0("Number of Terms is: ", nrow(df3)))
   print(paste0("Number of Unique-terms is: ", length(unique(df3$word))))
   
   ### データフレームを整形
   df4 <- df3 %>% 
      # select(id, word, num_in_doc) %>%
      select(id_ord, word, num_in_doc) %>%
      unique() %>% 
      pivot_wider(names_from = word, values_from = num_in_doc)
   
   
   df5 <- df4 %>% 
      # select(-id) %>% 
      select(-id_ord) %>% 
      as.data.frame()
   df5[is.na(df5)] <- 0
   rownames(df5) <- df4$id_ord
   
   out <- list()
   out$df <- df5
   out$df_tf <- df_tf
   out$limit_num <- limit_num
   out$term_type <- term_type
   return(out)
}   


# create_df <- function(txt_file, limit_num = 10, 
#                       term_type = c("名詞", "代名詞", "動詞", "形容詞", "形容動詞",
#                                     "副詞", "連体詞", "接続詞", "感動詞", "フィラー",
#                                     "記号", "助詞", "助動詞")) {
#    
#    ### Load libraries
#    library(tidyverse)
#    library(RMeCab)
#    library(topicmodels)
# 
#    ### Read data
#    df <- read_csv(txt_file, header = T)
#    colnames(df) <- c("Call_ID", "Text")
#    
#    ### 半角を全角に変換、漢数字を数字に変換、形態素解析を行う
#    txt <- df %>% 
#       mutate(Text_01 = gsub(" ", "", Text)) %>% # これを入れるかは要検討
#       mutate(Text_02 = kanji_num_2_num(Text_01)) %>%
#       mutate(Text_03 = stringi::stri_trans_nfkc(Text_02)) %>% 
#       mutate(Text_04 = stringi::stri_trans_toupper(Text_03)) %>% 
#       as.data.frame() %>% 
#       RMeCabDF("Text_04", 1)
#    
#    ### dataframeに変換
#    df1 <- purrr::pmap_df(list(nv = txt, id = df$Call_ID),
#                              function(nv, id) {
#                                 tibble(id = id, word = nv, type = names(nv))
#                              })
#    
#    
#    ### 必要な品詞を残し、頻度の合計値の列を追加
#    df2 <- df1 %>% 
#       filter(!grepl("^[ぁ-ゖァ-ー]{1}$", word)) %>% 
#       filter(!word %in% stop_words) %>% 
#       filter(type %in% term_type) %>% 
#       group_by(word) %>% 
#       mutate(num_ttl = n()) %>% 
#       group_by(id, word) %>% 
#       mutate(num_in_doc = n()) %>%
#       # mutate(flag = 1) %>% 
#       # group_by(word) %>% 
#       # summarise(num_doc = sum(flag)) %>% 
#       ungroup() %>% 
#       left_join(df1 %>% select(id, word) %>% unique() %>% 
#                    group_by(word) %>% summarize(num_doc = n()) %>% ungroup(),
#                 by = "word")
#    
#    ### 頻度を指定して単語を抽出する
#    df3 <- df2 %>% 
#       filter(num_ttl >= limit_num)
#    
#    ### 単語分布用
#    df_tf <- df3 %>% 
#       select(word, num_ttl, num_doc) %>% 
#       unique() %>% 
#       arrange(desc(num_ttl))
#    
#    ### 語数の確認
#    print(paste0("Number of Terms is: ", nrow(df3)))
#    print(paste0("Number of Unique-terms is: ", length(unique(df3$word))))
#    
#    ### データフレームを整形
#    df4 <- df3 %>% 
#       select(id, word, num_in_doc) %>%
#       unique() %>% 
#       pivot_wider(names_from = word, values_from = num_in_doc)
#    
#    
#    df5 <- df4 %>% 
#       select(-id) %>% 
#       as.data.frame()
#    df5[is.na(df5)] <- 0
#    rownames(df5) <- df4$id
#    
#    out <- list()
#    out$df <- df5
#    out$df_tf <- df_tf
#    out$limit_num <- limit_num
#    out$term_type <- term_type
#    return(out)
# }   

get_topic_terms <- function(obj) {
   n <- obj@Dim[2]
   topic_terms <- terms(obj, n)
   topic_num <- obj@k
   prob <- c()
   for (i in 1:topic_num) {
      prob <- cbind(prob,
                    sort(t(exp(obj@beta))[1:n, i], decreasing = T))
   }
   out <- c()
   for (i in 1:topic_num) {
      out <- cbind(out,
                   topic_terms[, i],
                   prob[, i])
   }
   colnames(out) <- paste0(rep(c("Topic_", "Prob_"), topic_num),
                           rep(1:topic_num, each = 2))
   out
}

# do_lda <- function(obj, k, out_file_path, ...){
#    
#    library(topicmodels)
#    
#    ### LDAの実行
#    res_lda <- LDA(obj$df, method = "Gibbs", k = k, 
#                   control = list(burnin = 1000, verbose = T, seed = 123))
#    
#    ### 各文章のトピック（確率が最大となるトピック）
#    print(topics(res_lda))
#    
#    ### 各トピックに所属する確率の高い単語
#    out <- as.data.frame(get_topic_terms(res_lda))
#    
#    ### 時刻を取得
#    t <- paste0(gsub(":", "-", gsub(" ", "-", as.character(Sys.time()))))
#    
#    ### パラメータ
#    # params <- data.frame(
#    #    Variable = c("num_doc", "term_type", "num_term", "minimum_term_freq", "num_topic"),
#    #    Value = c(nrow(df), paste(term_type, collapse = "/"), ncol(df), limit_num, k))
#    
#    ### ファイル名
#    term_prob_file <- paste0(t, "_TermProb.csv")
#    term_freq_file <- paste0(t, "_TermFreq.csv")
#    # param_file <- paste0(t, "_Params.csv")
#    obj_file <- paste0(t, "_Obj.Rdata")
#    
#    ### 書き出し
#    readr::write_excel_csv(out,       file = paste0(out_file_path, term_prob_file))
#    readr::write_excel_csv(obj$df_tf, file = paste0(out_file_path, term_freq_file))
#    # readr::write_excel_csv(params, file = paste0(out_file_path, param_file))
#    save(res_lda, file = paste0(out_file_path, obj_file))
#    out
# }

write_result <- function(df_lda, res_lda) {
   
   # LDAの結果を受け取り、timestampを同じとして以下を出力する：
   # ・トピックごとの各単語の所属確率
   # ・単語の出現頻度
   # ・行動変容の有無による各トピックの構成比率
   # ・結果オブジェクト
   
   ### Load library
   library(tidyverse)
   library(topicmodels)
   
   ### 時刻を取得
   t <- paste0(gsub(":", "-", gsub(" ", "-", as.character(Sys.time()))))
   
   ### ファイル名を作成
   term_prob_file <- paste0(t, "_TermProb.csv")
   term_freq_file <- paste0(t, "_TermFreq.csv")
   obj_file <- paste0(t, "_Obj.Rdata")
   topic_comp_file <- paste0(t, "_Topic_Composition_by_contract.csv")
   
   ### オブジェクト作成
   term_prob <- as_tibble(get_topic_terms(res_lda))
   term_freq <- df_lda$df_tf
   df_pst <- posterior(res_lda)$topics
   df_topic_cntr <- df_pst %>% 
      as_tibble() %>% 
      mutate(Call_ID = as.numeric(dimnames(df_pst)[[1]])) %>% 
      inner_join(df_lda$df_cntr, by = c("Call_ID" = "r_call_id")) %>% 
      select(Call_ID, is_info_material, is_contract, everything())
   
   out <- list(
      "term_prob" = term_prob, 
      "term_freq" = term_freq,
      "df_topic_cntr" = df_topic_cntr,
      "res_lda" = res_lda,
      "t" = t)
   
   ### 書き出し
   readr::write_excel_csv(term_prob, file = paste0(out_file_path, term_prob_file))
   readr::write_excel_csv(term_freq, file = paste0(out_file_path, term_freq_file))
   readr::write_excel_csv(df_topic_cntr, file = paste0(out_file_path, topic_comp_file))
   save(out, file = paste0(out_file_path, obj_file))
   return(out)
}


write_result_div <- function(df_lda, res_lda) {
   
   # LDAの結果を受け取り、timestampを同じとして以下を出力する：
   # ・トピックごとの各単語の所属確率
   # ・単語の出現頻度
   # ・行動変容の有無による各トピックの構成比率
   # ・結果オブジェクト
   
   ### Load library
   library(tidyverse)
   library(topicmodels)
   
   ### 時刻を取得
   t <- paste0(gsub(":", "-", gsub(" ", "-", as.character(Sys.time()))))
   
   ### ファイル名を作成
   term_prob_file <- paste0(t, "_TermProb.csv")
   term_freq_file <- paste0(t, "_TermFreq.csv")
   obj_file <- paste0(t, "_Obj.Rdata")
   topic_comp_file <- paste0(t, "_Topic_Composition_by_contract.csv")
   
   ### オブジェクト作成
   term_prob <- as_tibble(get_topic_terms(res_lda))
   term_freq <- df_lda$df_tf
   df_pst <- posterior(res_lda)$topics
   df_topic_cntr <- df_pst %>% 
      as_tibble() %>% 
      mutate(Call_ID = dimnames(df_pst)[[1]]) %>% 
      mutate(Call_ID2 = str_split(.$Call_ID, "_", simplify = T)[, 1]) %>% 
      mutate(Call_ID3 = as.numeric(Call_ID2)) %>% 
      inner_join(df_lda$df_cntr, by = c("Call_ID3" = "r_call_id")) %>% 
      select(Call_ID, Call_ID2, Call_ID3, is_info_material, is_contract, everything())
   
   out <- list(
      "term_prob" = term_prob, 
      "term_freq" = term_freq,
      "df_topic_cntr" = df_topic_cntr,
      "res_lda" = res_lda,
      "t" = t)
   
   ### 書き出し
   readr::write_excel_csv(term_prob, file = paste0(out_file_path, term_prob_file))
   readr::write_excel_csv(term_freq, file = paste0(out_file_path, term_freq_file))
   readr::write_excel_csv(df_topic_cntr, file = paste0(out_file_path, topic_comp_file))
   save(out, file = paste0(out_file_path, obj_file))
   return(out)
}



# calc_and_plot_cntr <- function(out) {
#    
#    ### Load libraries
#    library(ggplot2)
#    
#    df_plt <- out$df_topic_cntr
#    colnames(df_plt)[-c(1:3)] <- paste0("Topic_", 
#                                        sprintf("%02s", colnames(df_plt)[-c(1:3)]))
#    k <- ncol(df_plt) - 3
#    
#    ### Calc odds ratio
#    df_plt2 <- df_plt %>% 
#       mutate(Contract = if_else(is_contract == 1, "Yes", "No")) %>% 
#       select(-Call_ID, -is_info_material, -is_contract) %>% 
#       pivot_longer(!Contract,
#                    names_to = "Topic", values_to = "Perc") %>% 
#       group_by(Contract, Topic) %>% 
#       summarise(Perc_AVG = mean(Perc, na.rm = T)) %>% 
#       ungroup() %>%
#       ungroup() %>% 
#       pivot_wider(names_from = Contract, values_from = Perc_AVG) %>% 
#       mutate(Yes_Odds = Yes / (1-Yes),
#              No_Odds  = No  / (1-No)) %>% 
#       mutate(Odds_Ratio = Yes_Odds / No_Odds) %>% 
#       select(Topic, Yes, No, Odds_Ratio)
#    print(as.data.frame(df_plt2))
#    
#    ### Plot
#    ggplot(df_plt %>% 
#              select(-is_info_material) %>% 
#              pivot_longer(!c(Call_ID, is_contract),
#                           names_to = "Topic", values_to = "Perc"),
#           aes(x = Topic, y = Perc, fill = factor(is_contract))) +
#       geom_boxplot(outlier.shape = NA) +
#       # geom_point(position = position_jitter(width=0.05)) +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       theme(legend.position = "top")  +
#       labs(title = paste0("# of Topic : ", k), x = "")
#    
# }

calc_and_plot <- function(out) {
   
   ### Load libraries
   library(ggplot2)
   
   ### Create file name
   out_file_name <- paste(out$t, file_name, sep = "_")
   
   df_plt <- out$df_topic_cntr
   colnames(df_plt)[-c(1:3)] <- paste0("Topic_", 
                                       sprintf("%02s", colnames(df_plt)[-c(1:3)]))
   k <- ncol(df_plt) - 3
   
   ### Calc odds ratio
   df_plt2 <- df_plt %>% 
      mutate(Info_Req = if_else(is_info_material == 1, "Info_Yes", "Info_No")) %>% 
      select(-Call_ID, -is_info_material, -is_contract) %>% 
      pivot_longer(!Info_Req,
                   names_to = "Topic", values_to = "Perc") %>% 
      group_by(Info_Req, Topic) %>% 
      summarise(Info_Perc_AVG = mean(Perc, na.rm = T), .groups = "drop") %>% 
      pivot_wider(names_from = Info_Req, values_from = Info_Perc_AVG) %>% 
      mutate(Yes_Odds = Info_Yes / (1-Info_Yes),
             No_Odds  = Info_No  / (1-Info_No)) %>% 
      mutate(Info_Odds_Ratio = Yes_Odds / No_Odds) %>% 
      select(Topic, Info_Odds_Ratio)
   df_plt3 <- df_plt %>% 
      mutate(Cntr = if_else(is_contract == 1, "Cntr_Yes", "Cntr_No")) %>% 
      select(-Call_ID, -is_info_material, -is_contract) %>% 
      pivot_longer(!Cntr,
                   names_to = "Topic", values_to = "Perc") %>% 
      group_by(Cntr, Topic) %>% 
      summarise(Cntr_Perc_AVG = mean(Perc, na.rm = T), .groups = "drop") %>% 
      ungroup() %>%
      ungroup() %>% 
      pivot_wider(names_from = Cntr, values_from = Cntr_Perc_AVG) %>% 
      mutate(Yes_Odds = Cntr_Yes / (1-Cntr_Yes),
             No_Odds  = Cntr_No  / (1-Cntr_No)) %>% 
      mutate(Cntr_Odds_Ratio = Yes_Odds / No_Odds) %>% 
      select(Topic, Cntr_Odds_Ratio)
   df_out <- as.data.frame(df_plt2 %>% left_join(df_plt3, by = "Topic"))
   print(as.data.frame(df_plt2 %>% left_join(df_plt3, by = "Topic")))
   write_excel_csv(df_out, paste0(out_file_path, out_file_name, "_OddsRatio.csv"))
   
   ### Plot
   tmp1 <- df_plt %>% 
      mutate(Outcome = "01_Info_Req",
             Result  = if_else(is_info_material == 1, "01_Yes", "02_No")) %>% 
      select(Call_ID, Outcome, Result, everything(), -is_info_material, -is_contract)
   tmp2 <- df_plt %>% 
      mutate(Outcome = "02_Contract",
             Result  = if_else(is_contract == 1, "01_Yes", "02_No")) %>% 
      select(Call_ID, Outcome, Result, everything(), -is_info_material, -is_contract)
   tmp3 <- bind_rows(tmp1, tmp2)
   
   ggplot(tmp3 %>% 
             pivot_longer(!c(Call_ID, Outcome, Result),
                          names_to = "Topic", values_to = "Perc"),
          aes(x = Topic, y = Perc, fill = factor(Result))) +
      geom_boxplot(outlier.shape = NA) +
      facet_wrap(~Outcome, nrow = 2) +
      # geom_point(position = position_jitter(width=0.05)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "top")  +
      labs(title = paste0("# of Topic : ", k), x = "")
   ggsave(paste0(out_file_path, out_file_name, "_TopicCompMean.png"))
   
}


calc_and_plot_div <- function(out) {
   
   ### Load libraries
   library(ggplot2)
   
   ### Create file name
   out_file_name <- paste(out$t, file_name, sep = "_")
   
   df_plt <- out$df_topic_cntr
   colnames(df_plt)[-c(1:5)] <- paste0("Topic_", 
                                       sprintf("%02s", colnames(df_plt)[-c(1:5)]))
   k <- ncol(df_plt) - 5
   
   ### Calc odds ratio
   df_plt2 <- df_plt %>% 
      mutate(Info_Req = if_else(is_info_material == 1, "Info_Yes", "Info_No")) %>% 
      select(-Call_ID, -Call_ID2, -Call_ID3, -is_info_material, -is_contract) %>% 
      pivot_longer(!Info_Req,
                   names_to = "Topic", values_to = "Perc") %>% 
      group_by(Info_Req, Topic) %>% 
      summarise(Info_Perc_AVG = mean(Perc, na.rm = T), .groups = "drop") %>% 
      pivot_wider(names_from = Info_Req, values_from = Info_Perc_AVG) %>% 
      mutate(Yes_Odds = Info_Yes / (1-Info_Yes),
             No_Odds  = Info_No  / (1-Info_No)) %>% 
      mutate(Info_Odds_Ratio = Yes_Odds / No_Odds) %>% 
      select(Topic, Info_Odds_Ratio)
   df_plt3 <- df_plt %>% 
      mutate(Cntr = if_else(is_contract == 1, "Cntr_Yes", "Cntr_No")) %>% 
      select(-Call_ID, -Call_ID2, -Call_ID3, -is_info_material, -is_contract) %>% 
      pivot_longer(!Cntr,
                   names_to = "Topic", values_to = "Perc") %>% 
      group_by(Cntr, Topic) %>% 
      summarise(Cntr_Perc_AVG = mean(Perc, na.rm = T), .groups = "drop") %>% 
      ungroup() %>%
      ungroup() %>% 
      pivot_wider(names_from = Cntr, values_from = Cntr_Perc_AVG) %>% 
      mutate(Yes_Odds = Cntr_Yes / (1-Cntr_Yes),
             No_Odds  = Cntr_No  / (1-Cntr_No)) %>% 
      mutate(Cntr_Odds_Ratio = Yes_Odds / No_Odds) %>% 
      select(Topic, Cntr_Odds_Ratio)
   df_out <- as.data.frame(df_plt2 %>% left_join(df_plt3, by = "Topic"))
   print(as.data.frame(df_plt2 %>% left_join(df_plt3, by = "Topic")))
   write_excel_csv(df_out, paste0(out_file_path, out_file_name, "_OddsRatio.csv"))
   
   ### Plot
   tmp1 <- df_plt %>% 
      mutate(Outcome = "01_Info_Req",
             Result  = if_else(is_info_material == 1, "01_Yes", "02_No")) %>% 
      select(Call_ID, Outcome, Result, everything(), -is_info_material, -is_contract, -Call_ID2, -Call_ID3)
   tmp2 <- df_plt %>% 
      mutate(Outcome = "02_Contract",
             Result  = if_else(is_contract == 1, "01_Yes", "02_No")) %>% 
      select(Call_ID, Outcome, Result, everything(), -is_info_material, -is_contract, -Call_ID2, -Call_ID3)
   tmp3 <- bind_rows(tmp1, tmp2)
   
   ggplot(tmp3 %>% 
             pivot_longer(!c(Call_ID, Outcome, Result),
                          names_to = "Topic", values_to = "Perc"),
          aes(x = Topic, y = Perc, fill = factor(Result))) +
      geom_boxplot(outlier.shape = NA) +
      facet_wrap(~Outcome, nrow = 2) +
      # geom_point(position = position_jitter(width=0.05)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "top")  +
      labs(title = paste0("# of Topic : ", k), x = "")
   ggsave(paste0(out_file_path, out_file_name, "_TopicCompMean.png"))
   
}

# calc_and_plot_2016 <- function(out) {
#    
#    ### Load libraries
#    library(ggplot2)
#    
#    ### Create file name
#    out_file_name <- paste(out$t, file_name, sep = "_")
#    
#    df_plt <- out$df_topic_cntr
#    colnames(df_plt)[-c(1:3)] <- paste0("Topic_", 
#                                        sprintf("%02s", colnames(df_plt)[-c(1:3)]))
#    k <- ncol(df_plt) - 3
#    
#    ### Calc odds ratio
#    df_plt2 <- df_plt %>% 
#       mutate(Info_Req = if_else(is_info_material == 1, "Info_Yes", "Info_No")) %>% 
#       select(-Call_ID, -is_info_material, -is_contract) %>% 
#       pivot_longer(!Info_Req,
#                    names_to = "Topic", values_to = "Perc") %>% 
#       group_by(Info_Req, Topic) %>% 
#       summarise(Info_Perc_AVG = mean(Perc, na.rm = T), .groups = "drop") %>% 
#       pivot_wider(names_from = Info_Req, values_from = Info_Perc_AVG) %>% 
#       mutate(Yes_Odds = Info_Yes / (1-Info_Yes),
#              No_Odds  = Info_No  / (1-Info_No)) %>% 
#       mutate(Info_Odds_Ratio = Yes_Odds / No_Odds) %>% 
#       select(Topic, Info_Odds_Ratio)
#    df_plt3 <- df_plt %>% 
#       mutate(Cntr = if_else(is_contract == 1, "Cntr_Yes", "Cntr_No")) %>% 
#       select(-Call_ID, -is_info_material, -is_contract) %>% 
#       pivot_longer(!Cntr,
#                    names_to = "Topic", values_to = "Perc") %>% 
#       group_by(Cntr, Topic) %>% 
#       summarise(Cntr_Perc_AVG = mean(Perc, na.rm = T), .groups = "drop") %>% 
#       ungroup() %>%
#       ungroup() %>% 
#       pivot_wider(names_from = Cntr, values_from = Cntr_Perc_AVG) %>% 
#       mutate(Yes_Odds = Cntr_Yes / (1-Cntr_Yes),
#              No_Odds  = Cntr_No  / (1-Cntr_No)) %>% 
#       mutate(Cntr_Odds_Ratio = Yes_Odds / No_Odds) %>% 
#       select(Topic, Cntr_Odds_Ratio)
#    df_out <- as.data.frame(df_plt2 %>% left_join(df_plt3, by = "Topic"))
#    print(as.data.frame(df_plt2 %>% left_join(df_plt3, by = "Topic")))
#    write_excel_csv(df_out, paste0(out_file_path, out_file_name, "_OddsRatio.csv"))
#    
#    ### Plot
#    tmp1 <- df_plt %>% 
#       mutate(Outcome = "01_Info_Req",
#              Result  = if_else(is_info_material == 1, "01_Yes", "02_No")) %>% 
#       select(Call_ID, Outcome, Result, everything(), -is_info_material, -is_contract)
#    tmp2 <- df_plt %>% 
#       mutate(Outcome = "02_Contract",
#              Result  = if_else(is_contract == 1, "01_Yes", "02_No")) %>% 
#       select(Call_ID, Outcome, Result, everything(), -is_info_material, -is_contract)
#    tmp3 <- bind_rows(tmp1, tmp2)
#    
#    ggplot(tmp3 %>% 
#              pivot_longer(!c(Call_ID, Outcome, Result),
#                           names_to = "Topic", values_to = "Perc"),
#           aes(x = Topic, y = Perc, fill = factor(Result))) +
#       geom_boxplot(outlier.shape = NA) +
#       facet_wrap(~Outcome, nrow = 2) +
#       # geom_point(position = position_jitter(width=0.05)) +
#       theme_classic() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       theme(legend.position = "top")  +
#       labs(title = paste0("# of Topic : ", k), x = "")
#    ggsave(paste0(out_file_path, out_file_name, "_TopicCompMean.png"))
#    
# }

estimate_topic_in_call <- function(res_lda, df_lda, df_mcb, num = 1000) {
   
   library(tidyverse)
   library(topicmodels)
   
   ids <- rownames(df_lda$df)
   ids <- sample(x = ids, size = num, replace = F)
   num_doc <- length(ids)
   res_mat <- matrix(0, nrow = num_doc, ncol = 11) # 11 は 10 分割に id 列
   
   for (i in 1:num_doc) {
      tmp_df    <- df_mcb %>% filter(id == ids[i])
      cut_point <- round(quantile(1:nrow(tmp_df), seq(0, 1.0, by = 0.1)))
      id <- tmp_df$id[1]
      
      for (j in 1:10) {
         row_start <- cut_point[j]
         row_end   <- cut_point[j+1]
         tmp <- tmp_df[row_start:row_end, ]
         tmp <- tmp %>% 
            group_by(id, word) %>% 
            summarize(num = n(), .groups = "drop") %>% 
            pivot_wider(!id, names_from = word, values_from = num)
         est_topic <- which.max(posterior(res_lda, tmp, 
                                          control = list(seed = 123))$topic)
         res_mat[i, 1] <- id
         res_mat[i, j+1] <- est_topic
      }
      if(i %% 5 == 0) {
         print(paste("Doc", i, "th Complete at :", Sys.time()))
      }
   }
   
   return(res_mat)
}

estimate_topic_in_call_02 <- function(res_lda, df, df_mcb, window = 3, sep = 10, 
                                      seed = 123, num_call = 1000) {
   
   library(tidyverse)
   library(topicmodels)
   
   ## 対象とするコールをサンプルする   
   set.seed(seed)
   ids <- rownames(df)
   ids <- sample(x = ids, size = num_call, replace = F)
   num_doc <- length(ids)
   if (num_doc != num_call) print("Not Match")
   
   ## df_mcb を絞る
   df_mcb_tmp <- df_mcb %>% filter(id %in% ids)
   
   ## テキストの分割数と window に応じて確保するカラム数を決める
   ncol_to_est <- sep - window + 1 # 分割数 - window + 1
   res_mat <- matrix(0, nrow = num_doc, ncol = ncol_to_est + 1) # id 列を確保 
   
   ## コールごとにループ   
   for (i in 1:num_doc) {
      
      tmp_df    <- df_mcb_tmp %>% filter(id == ids[i])
      # cut_point <- round(quantile(1:nrow(tmp_df), seq(0, 1.0, by = 0.1)))
      cut_point <- round(seq(1, nrow(tmp_df), length = sep))
      id <- tmp_df$id[1]
      
      ## 分割されたテキストごとにループ
      for (j in 1:ncol_to_est) {
         row_start <- cut_point[j]
         row_end   <- cut_point[j+window-1] # 微妙に重複する行が生じるけど仕方ない
         tmp <- tmp_df[row_start:row_end, ]
         tmp <- tmp %>% 
            group_by(id, word) %>% 
            summarize(num = n(), .groups = "drop") %>% 
            pivot_wider(!id, names_from = word, values_from = num)
         est_topic <- which.max(posterior(res_lda, tmp, 
                                          control = list(seed = 123))$topic)
         res_mat[i, 1] <- id
         res_mat[i, j+1] <- est_topic
      }
      if(i %% 5 == 0) {
         print(paste("Doc", i, "th Complete at :", Sys.time()))
      }
   }
   
   return(res_mat)
}


create_trans_mat <- function(res_mat) {
   
   mat <- res_mat[, -1]
   num_slice <- ncol(mat)
   num_topic <- max(mat)
   trans_prob_mat <- matrix(0, num_topic, num_topic)
   
   for (i in 1:nrow(mat)) {
      for (j in 2:num_slice) {
         loc1 <- mat[i, j]
         loc2 <- mat[i, j-1]
         trans_prob_mat[loc2, loc1] <- trans_prob_mat[loc2, loc1] + 1
      }
   }
   
   ttl_cnt <- rep(0, num_topic)
   ttl_cnt[as.integer(names(table(mat[, -num_slice])))] <- table(mat[, -num_slice])
   names(ttl_cnt) <- 1:num_topic
   
   return(apply(trans_prob_mat, 2, function(x) x / ttl_cnt))
   
}


draw_hclust <- function(res){
   n_topic <- ncol(res) / 2
   df <- data.frame()
   for (i in 1:n_topic){
      tmp <- res[, c(i*2-1, i*2)]
      colnames(tmp) <- c("Topic", "Prob")
      tmp$id <- i
      df <- rbind(df, tmp)
   }
   df_cls <- 
      df %>% 
      pivot_wider(Topic, names_from = id, names_prefix = "Topic_",
                  values_from = Prob) %>% 
      select(-Topic) %>% 
      t()
   topic_dist <- dist(df_cls)
   plot(hclust(topic_dist))
}
