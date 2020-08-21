      subroutine standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)  
      integer ju(ni)                                                    
      double precision, dimension (:), allocatable :: v                 
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=w/sum(w)  ! ここで重みを 重みの総和あたりの重みに変換している。                                                      
      v=sqrt(w)  ! 二乗和を求めるときに各観測値の重みを二乗することになるため、事前に平方根を取っておくのだと思う（平方和に対して重みを乗じるイメージ）  
      
      ! intr は intercept なので切片が 0 であるかで判定
      ! 切片が 0 でない場合は 10141 に飛ばされる
      if(intr .ne. 0)goto 10141                                         

      ! 以下のセクションでは y と x それぞれについて観測値の重み weights を使って色々と調整する
      ! まずは y      
      ym=0.0                                                            
      y=v*y ! v は sqrt(w), w = w/sum(w), w は weights なので観測値に対する重み   
      ! dot_product は内積                                                         
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2) ! y の重み付き二乗和から、y の重み付き平均（v を使って、その後に二乗しているので重み付き平均ぽいもの）を減じ（分散）、平方根をとる（SD）
      y=y/ys ! さらに 元の値と重み調整後の値を比率を y として返す -> 標準偏差で割っている！

      ! 次に x
      do 10151 j=1,ni ! ni は nvars

      ! ju は各変数に 2 つ以上の異なる値があるか（バラツキがあるか）を示すフラグで chkvars で作られる
      ! 0 なら全ての値が同じなのでスキップ                                        
      if(ju(j).eq.0)goto 10151
      ! y と同様の操作を各変数に対して実施    
      xm(j)=0.0                                                         
      x(:,j)=v*x(:,j) ! x にも重みを乗じる                                                  
      xv(j)=dot_product(x(:,j),x(:,j)) ! xv に重みを乗じた x の二乗和を格納                                 
      if(isd .eq. 0)goto 10171 ! isd は標準化するかの指定なので、標準化する場合は 1 が入っておりここはスキップされる（10171 に飛ばされない）                                         
      xbq=dot_product(v,x(:,j))**2 ! 重み付き平均ぽいものの二乗                                      
      vc=xv(j)-xbq ! 重み付き二乗和から重み付き平均ぽいものを引く（分散）                                                     
      xs(j)=sqrt(vc) ! 平方根をとる（SD）     ys と対応している。                                              
      x(:,j)=x(:,j)/xs(j) ! 標準偏差で割っている    ! y/ys と対応している                                          
      xv(j)=1.0+xbq/vc ! 標準偏差 / 分散 に 1 を加算している                                                  
      goto 10181                                                        
10171 continue                                                          
      xs(j)=1.0                                                         
10181 continue                                                          
      continue                                                          
10151 continue                                                          
      continue                                                          
      goto 10191                                                        

      ! 切片が 0 でない場合
10141 continue                                                          
      do 10201 j=1,ni           

      ! ju は各変数に 2 つ以上の異なる値があるか（バラツキがあるか）を示すフラグで chkvars で作られる 
      ! 0 なら全ての値が同じなのでスキップ                                        
      if(ju(j).eq.0)goto 10201                                        
      
      xm(j)=dot_product(w,x(:,j)) ! x の重み付き平均                                 
      x(:,j)=v*(x(:,j)-xm(j))     ! 重み付き平均を引いたあとに重みを乗じる                                      
      xv(j)=dot_product(x(:,j),x(:,j)) ! その二乗和（分散）                                 
      if(isd.gt.0) xs(j)=sqrt(xv(j)) ! その平方根（SD）                                   
10201 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 10221                                          
      xs=1.0                                                            
      goto 10231                                                        
10221 continue                                                          
      do 10241 j=1,ni                                                   
      if(ju(j).eq.0)goto 10241                                          
      x(:,j)=x(:,j)/xs(j) ! 標準化はここで実行                                              
10241 continue                                                          
      continue                                                          
      xv=1.0                                                            
10231 continue                                                          
      continue                                                          
      ym=dot_product(w,y) ! y の重み付き平均                                              
      y=v*(y-ym)          ! y から重み付き平均を引いたものに重みを乗じる                                              
      ys=sqrt(dot_product(y,y)) ! 二乗和（分散）の平方根（SD）                                        
      y=y/ys ! 標準化

      ! 切片の有無それぞれについてのケースがここで合流                                                            
10191 continue                                                          
      continue                                                          
      g=0.0                                                             
      do 10251 j=1,ni 
      ! j 番目の変数にバラツキがあるなら g に　y と x の内積（共分散）を格納する                                                  
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                          
10251 continue                                                          
      continue                                                          
      deallocate(v)                                                     
      return                                                            
      end                                                               
