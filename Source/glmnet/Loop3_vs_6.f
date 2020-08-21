      do 10371 k=1,ni                                                   
      if(ju(k).eq.0)goto 10371 
      ! ju は各説明変数列のバラツキを示す 1/0 のベクトル。バラツキがなければ（ju(k) == 0 であれば）評価をスキップ                                         
      ak=a(k)  
      ! k 番目の変数の a の値を ak に代入。 a = 0.0 で初期化されているので ak も 0。
      ! ただし lambda のループ（一番外側のループ） の二回目以降は 縮小された回帰係数が入っている                                                     

      ! ここは大きなポイント!!!!
      ! 回帰係数を縮小
      u=g(k)+ak*xv(k)  
      ! g(k) に ak を加算する。g(k) は y と k 番目の x の内積。 標準化前提なら回帰係数に相当（beta = cov(y, x) / var(x)なので）
      ! ak は lamda のループ（一番外側）１回目は 0、２回目以降は前回の lambda で得られた推定値。 xv で重みをつけたものを加算する
      ! xv は weight を乗じた x の二乗和                                                
      v=abs(u)-vp(k)*ab 
      ! x の回帰係数から vp * ab を減じる、つまり回帰係数を縮小している。 b(alpha) が 1（lasso）なら lambda の分だけ縮小される。
      ! vp(k) は k 番目の変数に対する罰則の重みでデフォルトは 1 なので ab を減じる。 ab は alm * bta                                               
      
      ! a(k) を更新
      ! a は回帰係数なので、 v が 0.0 より大きくないと 0 として推定される!!!
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem))) 
      ! v が 0 より大きい（前回ループの推定値などで調整した回帰係数の絶対値が lambda * alpha より大きい）場合、 a を更新する。そうでなければ 0のまま。
      ! sign(a, b) は a の絶対値に b の符号を乗じる関数。 u は y と x の内積なので回帰係数の符号を取ってきて v に乗じる
      ! xv は weight を乗じた x の二乗和
      ! vp(k) は k 番目の変数に対する罰則の重みでデフォルトは 1                                               
      ! d*em になってるけど dem じゃないの？？？？
      ! cl は glmnet.r で cl = rbind(lower.limits, upper.limits) として定義されているので、 1行目 は 下限、2行目は上限が入っている
      ! max(下限, min(上限, 縮小された回帰係数)) という形なので、上限・下限の範囲内で回帰係数を縮小する
      ! d*em は dem のことなら alm * (1-bta)
      ! x の weight 調整済み二乗和に罰則 * 配分 を加算したもので 縮小された回帰係数を調整したものを a に格納している

      ! a(k) が ak と同じなら（このループで更新がなければ）ループを抜けて次の変数へ
      if(a(k).eq.ak)goto 10371                                          
      ! mm は lambda ループの１回目では 0 なので１回目だけ処理を行う？                                        
      if(mm(k) .ne. 0)goto 10391 ! 10391 は４番目のループの先なので、 mm が 0 でなければ ４番目のループをスキップ
      nin=nin+1  ! mm(k) が 0 なら nin を +1 する。 おそらく、パラメータが 0 でないときに mm は 0 となる。                                                      
      if(nin.gt.nx)goto 10372 ! nx は非ゼロとする変数の上限なので、推定したパラメータ数がそれを越えると２番目のループを抜ける
      
      ! ４番目のループ
      ! 分散共分散行列のようなものを作っている
      ! ここでもループの対象は説明変数（ただしインデックスは k ではなく j）
      do 10401 j=1,ni                                                   
      ! バラツキがなければ以降の処理をスキップ                                          
      if(ju(j).eq.0)goto 10401
      ! mm が 0（パラメータが 0 でない）なら 以降をスキップ。                                        
      if(mm(j) .eq. 0)goto 10421
      c(j,nin)=c(k,mm(j))                                               
      goto 10401                                                  
10421 continue                                                          
      if(j .ne. k)goto 10441  ! 変数が同一でなければ 10441 に飛ぶ                                          
      c(j,nin)=xv(j) ! 同一だったらここ                                                   
      goto 10401                                                        
10441 continue                                                          
      c(j,nin)=dot_product(x(:,j),x(:,k)) ! 同一でなかったら j と k の内積をとる                               
10401 continue ! ４番目のループはここまで

      continue 
      ! mm に nin を入れる                                                         
      mm(k)=nin
      ! ia に k を格納                                                         
      ia(nin)=k ! 0 でないパラメータが推定された変数の位置                                                         
10391 continue   
      ! a(k) の差分をとる。 a(k)、 ak は推定された回帰係数。                                                       
      del=a(k)-ak
      ! 残差平方和を更新する
      rsq=rsq+del*(2.0*g(k)-del*xv(k))
      ! rsq = rsq + del * (2.0 * g(k) - del * xv(k))
      ! rsq は残差平方和
      ! del は a(k)-ak
      ! g(k) は縮小前の回帰係数（y と x(k) の内積）、そこから weight 調整済みの x の二乗和 を減じる                                                      
      dlx=max(xv(k)*del**2,dlx)

      ! ５番目のループ          
      ! 探索範囲は三度説明変数                               
      do 10451 j=1,ni ! インデックスは再度 j を使う                                                  
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del            