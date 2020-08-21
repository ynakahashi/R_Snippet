      subroutine elnetu(parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,  flmin,ulam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)  
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam   
      integer, dimension (:), allocatable :: ju                         
      allocate(g(1:ni),stat=jerr)  ! g は変数の数を次元にもつベクトル                                     
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr) ! ju は変数の数だけの次元をもつベクトル                                     
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr) ! xv は変数の数だけの次元を持つベクトル                                     
      if(jerr.ne.0) return                                              
      allocate(vlam(1:nlam),stat=jerr)                                  
      if(jerr.ne.0) return

      ! 1. 前処理
      ! 説明変数のバラツキのチェック  
      call chkvars(no,ni,x,ju)

      ! jd は glmnet の中で以下のように定義されている
      ! if (length(exclude) > 0) { jd = match(exclude, seq(nvars), 0)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0
      
      ! ju の maxval が 0 より大きい場合はスキップ
      ! そうでない場合は jerr = 7777 で return
      ! jerr 7777 は　msg = "All used predictors have zero variance"                              
      if(maxval(ju) .gt. 0)goto 10071                                   
      jerr=7777                                                         
      return                                                            
10071 continue                                                          
      
      ! 2. 標準化
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)      
      
      ! jerr に 0 でない値が入っていると return
      if(jerr.ne.0) return                                              

      ! ys は y の重み調整済み標準偏差のようなやつ
      ! cl は lower.limit と upper.limit が格納されている
      ! cl = rbind(lower.limits, upper.limits) in glmnet
      ! 目的変数を標準化しているから回帰係数も標準化するのか
      cl=cl/ys

      ! 標準化の指定が 0 であれば以下はスキップ                                                          
      if(isd .le. 0)goto 10091                                          
      ! 説明変数の標準偏差を乗じる
      ! ni は nvars
      do 10101 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
10101 continue                                                          
      continue                                                          
10091 continue                                                          
      
      ! flmin は glmnet のなかで flmin = as.double(lambda.min.ratio) で定義される
      ! ここで lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04)
      ! また if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
      ! なので基本的にここは引っかからないはず
      if(flmin.ge.1.0) vlam=ulam/ys 

      ! 3. フィッティング
      ! 本体である elnet1 の呼び出し
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)

      ! jerr が 0 でなければ return
      if(jerr.gt.0) return

      ! 4. 後処理
      ! lmu はデフォルトで 1 
      ! lmu = integer(1)                                             
      do 10111 k=1,lmu
      ! alm は nlambda、lambdaの長さ（ただし double）
      ! alm = double(nlam)                                                  
      ! ys は y の重み調整済み標準偏差
      alm(k)=ys*alm(k)                                                  
      
      ! nin は nlambda（ただし integer）
      ! nin = integer(nlam)
      nk=nin(k)

      do 10121 l=1,nk                                                   
      ! 回帰係数に y の重み調整済み標準偏差を乗じ、説明変数の重み調整済み標準偏差で除す
      ! ここで lambda の個数だけ係数が格納されている
      ! ca は変数の数 * lambda の数
      ! ys は重み調整済み標準偏差
      ! v は sqrt(w)
      ! w は w = w/sum(w)
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
10121 continue

      continue
      ! a0 は elnet.r で以下のように定義： a0 = double(nlam)                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk))) ! y の重み付き平均 - y_hat の平均   
10111 continue

      continue                                                          
      deallocate(xm,xs,g,ju,xv,vlam)                                    
      return                                                            
      end                                                               
