      subroutine elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,  flmin,u
     *lam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni) 
      double precision ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)          
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: vq;               
      if(maxval(vp) .gt. 0.0)goto 10021                                 
      jerr=10000                                                        
      return                                                            
10021 continue                                                          
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ! elnetu か elnetn のどちらを呼ぶかは ka .ne. 1 であるかで判断している
      ! 1 でなければ elnetn 、 1 なら elnetu
      ! ka は elnet の第一引数であり、ka = as.integer(switch(type.gaussian, covariance = 1, naive = 2, ))
      ! ということなので、 type.gaussian を指定している
      if(ka .ne. 1)goto 10041                                           
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,
     *isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                        
10041 continue                                                          
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i
     *sd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                          
      continue                                                          
      deallocate(vq)                                                    
      return                                                            
      end                                                               
