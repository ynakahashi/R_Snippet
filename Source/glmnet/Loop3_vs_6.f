! Loop 3 
      do 10371 k=1,ni                                                   
      if(ju(k).eq.0)goto 10371 
      ak=a(k)  
      u=g(k)+ak*xv(k)  
      v=abs(u)-vp(k)*ab 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem))) 
      if(a(k).eq.ak)goto 10371                                          
      if(mm(k) .ne. 0)goto 10391 ! 10391 は４番目のループの先なので、 mm が 0 でなければ ４番目のループをスキップ
      nin=nin+1  ! mm(k) が 0 なら nin を +1 する。 おそらく、パラメータが 0 でないときに mm は 0 となる。                                                      
      if(nin.gt.nx)goto 10372 ! nx は非ゼロとする変数の上限なので、推定したパラメータ数がそれを越えると２番目のループを抜ける
      continue 
      mm(k)=nin
      ia(nin)=k ! 0 でないパラメータが推定された変数の位置                                                         
10391 continue   
      del=a(k)-ak
      rsq=rsq+del*(2.0*g(k)-del*xv(k))
      dlx=max(xv(k)*del**2,dlx)
      do 10451 j=1,ni ! インデックスは再度 j を使う                                                  
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del        


! Loop 6
      do 10491 l=1,nin                                                  
      k=ia(l)                                                           
      ak=a(k)                                                        
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))
      if(a(k).eq.ak)goto 10491                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 10501 j=1,nin                                                  
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                              
    