

## 祖先となる個体が何人いれば、遺伝的多様性を保ったまま個体数を安定して増やしていいけるかシミュレーションしたい

## 条件
## 各世代の個体は世代ごとの男女間で繁殖し、世代をまたがない
## 各個体は繁殖をしないことはあるが、２人以上の個体とは繁殖しない
## 生まれる個体の男女比は 1:1 とする
## 子供の数はパラメータ Lambda_n のポアソン分布にしたがう
## 祖先個体の遺伝的多様性は極めて大であるとする（近交係数はいずれも０）

## 変化させたい変数
## N_m_0, N_f_0 : 祖先個体(Generation = 0)となる男女の個体数。N_m_0 = N_f_0 とする

## 評価値
## N_j : j 番目の世代の個体数
## D_j : j 番目の世代の遺伝的多様性（集団の有効サイズ）

## その他
## Lambda_n : 子供の数を決めるパラメータ
## J : シミュレーションを実施する世代数
## tv_eff : 集団の有効サイズが小さくなりすぎて以降の繁殖ができないと判断され、シミュレーションを打ち切る閾値
## tv_brd_m : 男性個体が繁殖行為を行うと判断するときの閾値
## tv_brd_f : 女性個体が繁殖行為を行うと判断するときの閾値

## 血統ファイル
## ID, Father_ID, Mother_ID, Inbreeding_Coef, Generation, Sex, Will

## イメージ
## 結果を格納しておくマトリックス
result <- matrix(NA, J+1, 2)

## まずは０世代目の繁殖意思決定を行う
pedigree_file <- generate_ancestors(N_m_0, N_f_0)
N_col <- 7
for (i in 1:N_m_0) {
    pedigree_file$Will[i] <- ifelse(runif(1) > th_brd_m, 1, 0)
    pedigree_file$Will[(i + N_f_0)] <- ifelse(runif(1) > th_brd_f, 1, 0)
}

## ０世代目の評価値を得る
result[1, 1] <- N_j <- N_m_0 + N_f_0
result[1, 2] <- D_j <- calculate_inbreed_coef(pedigree_file)

## 繁殖意思を持つ個体について順番に子を作る
N_m_j_brd <- sum(pedigree_file$Will[pedigree_file$Generation == 0 && pedigree_file$Sex == 1])
N_f_j_brd <- sum(pedigree_file$Will[pedigree_file$Generation == 0 && pedigree_file$Sex == 0])
N_couple_j <- min(N_m_j_brd, N_f_j_brd)
count <- 1
for (i in 1:N_couple_j) {
    N_child <- rpois(1, lambda = Lambda_n)
    pedigree_file <- rbind(pedigree_file, matrix(0, N_child, N_col))

    for (j in 1:N_child) {
        pedigree_file$ID[N_j + count] <- N_j + count
        pedigree_file$Father_ID <- 
        pedigree_file$Father_ID <- 
#        pedigree_file$Inbreeding_Coef <- calculate_inbreed_coef()
        pedigree_file$Generation <- 1
        pedigree_file$Sex <- ifelse(runif(1) > 0.5, 1, 0)
        pedigree_file$Will <- 0
        count <- count + 1
    }
}


## 続いて１世代目から
for (j in 1:J) {
    N_m_j <- nrow(Pedigree_file[, Pedigree_file$Sex == "Male"])
    N_f_j <- nrow(Pedigree_file[, Pedigree_file$Sex == "Female"])
    N_j <- N_m_j + N_f_j

    ## 男性個体が繁殖行為を行うかの判断
    ## ここは将来的にはベクトル化できる
    N_m_j_brd <- matrix(0, N_m_j, 1)
    for (k in 1:N_m_i) {
        N_m_j_brd[k] <- ifelse(runif(1) > th_brd_m, 1, 0)
    }    

    ## 女性個体が繁殖行為を行うかの判断
    ## ここも将来的にはベクトル化できる
    N_f_j_brd <- matrix(0, N_f_j, 1)
    for (l in 1:N_f_i) {
        N_f_j_brd[l] <- ifelse(runif(1) > th_brd_f, 1, 0)
    }    




}



