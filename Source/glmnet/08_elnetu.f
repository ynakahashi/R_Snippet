      subroutine elnetu(parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,  flmin,ula
     *m,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)  
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam   
      integer, dimension (:), allocatable :: ju                         
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vlam(1:nlam),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 10071                                   
      jerr=7777                                                         
      return                                                            
10071 continue                                                          
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)      
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 10091                                          
      do 10101 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
10101 continue                                                          
      continue                                                          
10091 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return
      ! lmu はデフォルトで 1 
      ! lmu = integer(1)                                             
      do 10111 k=1,lmu
      ! alm は nlambda、lambdaの長さ（ただし double）
      ! alm = double(nlam)                                                  
      alm(k)=ys*alm(k)                                                  
      ! nin は nlambda（ただし integer）
      ! nin = integer(nlam)
      nk=nin(k)                                                         
      do 10121 l=1,nk                                                   
      ! ここで lambda の個数だけ係数が格納されている
      ! ca は変数の数 * lambda の数
      ! ys は重み調整した y で、２パターンある（切片が 0 であるかで分かれる）
      ! ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)   
      ! ys=sqrt(dot_product(y,y))
      ! 違いは dot_product(v,y)^2 を減じるか
      ! v は sqrt(w)
      ! w は w = w/sum(w)
      ! weight でしょう
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
10121 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
10111 continue                                                          
      continue                                                          
      deallocate(xm,xs,g,ju,xv,vlam)                                    
      return                                                            
      end                                                               
