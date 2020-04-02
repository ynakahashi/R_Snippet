      subroutine lognet(parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,ul
     *am,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nla
     *m)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv    
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 12221                                 
      jerr=10000                                                        
      return                                                            
12221 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(kopt .ne. 2)goto 12241      ! kopt が 2 の時は xv も allocate （メモリの動的割付け）する                                    
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
12241 continue                                                          
      if(isd .le. 0)goto 12261       ! isd が 1 の時は xs を allocate する                                   
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
12261 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 12281                                   
      jerr=7777                                                         
      return                                                            
12281 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      do 12291 i=1,no                                                   
      ww(i)=sum(y(i,:))                                                 
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                              
12291 continue                                                          
      continue                                                          
      sw=sum(ww)                                                        
      ww=ww/sw                                                          
      if(nc .ne. 1)goto 12311  ! nc .ne. 1 なら lognet2n 、そうでないなら lognetn or mullognetn                                       
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 12331 ! ここで isd .le. 0(標準化の指定が0) なら直接 *netn に飛ぶ、そうでないなら cl() の処理を挟む                                      
      do 12341 j=1,ni  ! 列の数？                                                 
      cl(:,j)=cl(:,j)*xs(j) ! 元の値に xs を乗じる（多分標準化してる）                                            
12341 continue                                                          
      continue                                                          
12331 continue                                                          
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl　 ! lognet2n
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12301    ! lognet2n を実行したら multlognetn、lognetnを飛ばす（それはそう） 
12311 if(kopt .ne. 2)goto 12351    ! kopt .ne. 2 なら lognetn、 そうでないなら multlognetn                                     
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)              
      if(isd .le. 0)goto 12371                                          
      do 12381 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12381 continue                                                          
      continue                                                          
12371 continue                                                          
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin, 
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12391                                                        
12351 continue                                                          
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 12411                                          
      do 12421 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12421 continue                                                          
      continue                                                          
12411 continue                                                          
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam 
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12391 continue                                                          
12301 continue                                                          
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      do 12431 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 12441 ic=1,nc                                                  
      if(isd .le. 0)goto 12461                                          
      do 12471 l=1,nk                                                   
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                   
12471 continue                                                          
      continue                                                          
12461 continue                                                          
      if(intr .ne. 0)goto 12491                                         
      a0(ic,k)=0.0                                                      
      goto 12501                                                        
12491 continue                                                          
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))         
12501 continue                                                          
      continue                                                          
12441 continue                                                          
      continue                                                          
12431 continue                                                          
      continue                                                          
      deallocate(ww,ju,vq,xm)                                           
      if(isd.gt.0) deallocate(xs)                                       
      if(kopt.eq.2) deallocate(xv)                                      
      return                                                            
      end       