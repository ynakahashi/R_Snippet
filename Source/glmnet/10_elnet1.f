      subroutine elnet1(beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,th
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam)     
      double precision rsqo(nlam),almo(nlam),xv(ni)                     
      double precision cl(2,ni)                                         
      integer ju(ni),ia(nx),kin(nlam)                                   
      double precision, dimension (:), allocatable :: a,da              
      integer, dimension (:), allocatable :: mm                         
      double precision, dimension (:,:), allocatable :: c               
      allocate(c(1:ni,1:nx),stat=jerr)                                  
      if(jerr.ne.0) return;                                             
      
      ! 初期パラメータを取得
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      
      ! a, mm, da を allocate
      allocate(a(1:ni),stat=jerr)  ! a は説明変数の数の次元をもつベクトル                                     
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr) ! mm は説明変数の数の次元をもつベクトル                                     
      if(jerr.ne.0) return                                              
      allocate(da(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      
      bta=beta ! beta はもともと parm として受け渡されたもので、これは elnet.r で parm = alpha として渡している。この alpha は glmnet.r で定義されたもの。
               ! L1 と L2 それぞれに対する罰則の配分                                                         
      omb=1.0-bta                                                       
      alm=0.0                                                           
      alf=1.0                                                           
      
      ! flmin が 1.0 以上の場合（lambda が指定されている場合）は以下をスキップして 10271 に移動
      if(flmin .ge. 1.0)goto 10271                                      
      eqs=max(eps,flmin) ! eps と flmin の大きい方を eqs とする。 eps は get_int_parms で eps0 を受け取っており、 1.0d-6 なので 0.000001                                 
      alf=eqs**(1.0/(nlam-1))  ! alf を eqs の (1/(nlam-1)) で定義する                                         

10271 continue
      ! パラメータの初期化                                                          
      rsq=0.0 ! 残差平方和                                                          
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp ! nin はここでは 0                                                          
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      
      ! ここからがループの本体
      ! 一番外側のループは lambda の個数（nlam）
      !! 以下のまとまりは： 
      !!! flmin が 1.0 より小さい場合（glmnet において lambda に指定がない場合）は goto 10301 なので alm = ulam(m) をスキップ
      !!! flmin が 1.0 以上の場合は alm = ulam(m) とし、10301 ~ 10331 の処理をスキップ
      do 10281 m=1,nlam ! nlambda なので lambda の個数（デフォルトは 100）だけループ。 m はループのインデックス                                                
      if(itrace.ne.0) call setpb(m-1)  ! プログレスバー                                 
      if(flmin .lt. 1.0)goto 10301                                      
      alm=ulam(m) ! flmin が 1.0 以上の場合は alm = ulam(m) とする。 glmnet.r で lambda の指定がない場合は ulam は 1 が入る。基本指定しないと思うので 1 。指定があれば 1 より大きい値もありうる。                                                     
      goto 10291                                                        
10301 if(m .le. 2)goto 10311 ! ループの１回目と２回目はここをスキップ                                           
      alm=alm*alf ! ループの３回目からは alm を alf を乗じて 10311、10321、10341、10331 の処理をスキップ                                                      
      goto 10291                                                        
10311 if(m .ne. 1)goto 10321 ! ループの２回目はここをスキップ                                           
      alm=big     ! ループの１回目は alm = 0.0 を big = 9.9d35（アホみたいに大きな値） にして 10321、10341の処理をスキップ                                                      
      goto 10331                                                        
10321 continue                                                          
      alm=0.0     ! ループの２回目は alm を いったん 0 にする                                                      
      
      ! ２番目のループ
      ! alm の更新
      do 10341 j=1,ni  ! ni は変数の数                                                 
      if(ju(j).eq.0)goto 10341  ! ju は 各変数列において全く同じ値が入っていないかを確認したもの（chkvars で行われる）で、全く同じなら 0 となりスキップされる                                   
      if(vp(j).le.0.0)goto 10341   ! vp は各変数に対する罰則の重み（デフォルトは 1） が入ったベクトル（vp = as.double(penalty.factor) in glmnet.r）                                     
      alm=max(alm,abs(g(j))/vp(j)) ! ここで alm と abs(g(j)) / vp(j) の max をとる。罰則を考慮した内積の最大値を取りたい。 g は j 番目の変数にバラツキがあるとき y と x の内積（共分散）を格納したもの。内積の絶対値を罰則で割ったものを alm とする。でも alm=0.0 としていて、 vp は [0,1) なので、要は内積が負なら 0、正なら内積をとるので [0, 内積] にしている。
10341 continue  ! ２番目のループここまで                                                        

      continue                                                          
      alm=alf*alm/max(bta,1.0d-3)  ! さらにここで bta（つまり alpha、 L1 と　L2 への重みの配分パラメータ）と 0.001 の max で alm を割る                              
10331 continue                                                          
10291 continue                                                          
      dem=alm*omb ! omb は 1-bta                                                      
      ab=alm*bta  ! さらに alm に bta を乗じたものを ab とする                                                      
      rsq0=rsq  ! rsq はこの時点で 0.0 なので rsq0 も 0                                                        
      jz=1                                                              
      continue                                                          
10351 continue                                                          
      if(iz*jz.ne.0) goto 10360   ! iz = 0, jz = 1                                     
      nlp=nlp+1  ! nlp は 0 だったので一番外側のループ（lambda）のカウンタになってる？                                                      
      dlx=0.0                                                           
      
      ! ３番目のループは本体
      ! ここで回帰係数を推定している
      ! ３番目のループは ni なので説明変数。k をインデックスとして各説明変数をさらう。
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
      if(nin.gt.nx)goto 10372 ! nx は非ゼロとする変数の上限なので、推定したパラメータ数がそれを越えると３番目のループを抜ける
      
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
10451 continue ! ５番目のループはここまで                                                         
      continue                                                          
10371 continue ! ３番目のループはここまで

10372 continue                                                          
      if(dlx.lt.thr)goto 10352                                          
      if(nin.gt.nx)goto 10352                                           
      if(nlp .le. maxit)goto 10471                                      
      jerr=-m                                                           
      return                                                            
10471 continue                                                          
10360 continue                                                          
      iz=1                                                              
      da(1:nin)=a(ia(1:nin))                                            
      continue                                                          
10481 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0

      ! ６番目のループ
      ! ３番目のループと同じことを ni ではなく nin に対して再度実行                                                           
      do 10491 l=1,nin                                                  
      k=ia(l) ! k を取り出す（ ia には 0 ではないパラメータが推定された変数の列が格納されてる）                                                         
      ak=a(k) ! a を取り出す                                                            
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))
      if(a(k).eq.ak)goto 10491                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)

      ! ７番目のループ
      ! 上と同様、 nin に対して g を更新                                         
      do 10501 j=1,nin                                                  
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                              
10501 continue ! ７番目のループここまで

      continue                                                          
10491 continue ! ６番目のループここまで

      continue                                                          
      if(dlx.lt.thr)goto 10482                                          
      if(nlp .le. maxit)goto 10521                                      
      jerr=-m                                                           
      return                                                            
10521 continue                                                          
      goto 10481  ! えっ！！！                                                      
10482 continue                                                          
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                  
      
      ! ８番目のループ
      do 10531 j=1,ni                                                   
      if(mm(j).ne.0)goto 10531                                          
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))        
10531 continue ! ８番目のループここまで
                                                         
      continue                                                          
      jz=0                                                              
      goto 10351  ! えっ！！ ３番目のループの開始まで戻すの！                                                     
10352 continue                                                          
      if(nin .le. nx)goto 10551  ! nin が nx を超えた場合はここにくる                                       
      jerr=-10000-m                                                     
      goto 10282 ! jerr を 更新して elnet1 を抜ける                                                      
10551 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin   ! m 回目のループの nin を kin[m] に格納する                                                     
      rsqo(m)=rsq  ! m 回目のループの rsq を rsqo[m] に格納する                                                     
      almo(m)=alm  ! m 回目のループの alm を almo[m] に格納する                                                     
      lmu=m                                                             
      if(m.lt.mnl)goto 10281                                            
      if(flmin.ge.1.0)goto 10281                                        
      me=0 

      ! ９番目のループ                                                             
      do 10561 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
10561 continue ! ９番目のループここまで 

      continue                                                          
      if(me.gt.ne)goto 10282                                            
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                 
      if(rsq.gt.rsqmax)goto 10282                                       
10281 continue ! lambda のループはここまで

10282 continue                                                          
      deallocate(a,mm,c,da)                                             
      return                                                            
      end                                                               
