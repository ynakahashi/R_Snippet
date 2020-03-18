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
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(da(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      alm=0.0                                                           
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 10271                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
10271 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 10281 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      if(flmin .lt. 1.0)goto 10301                                      
      alm=ulam(m)                                                       
      goto 10291                                                        
10301 if(m .le. 2)goto 10311                                            
      alm=alm*alf                                                       
      goto 10291                                                        
10311 if(m .ne. 1)goto 10321                                            
      alm=big                                                           
      goto 10331                                                        
10321 continue                                                          
      alm=0.0                                                           
      do 10341 j=1,ni                                                   
      if(ju(j).eq.0)goto 10341                                          
      if(vp(j).le.0.0)goto 10341                                        
      alm=max(alm,abs(g(j))/vp(j))                                      
10341 continue                                                          
      continue                                                          
      alm=alf*alm/max(bta,1.0d-3)                                       
10331 continue                                                          
10291 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      continue                                                          
10351 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10371 k=1,ni                                                   
      if(ju(k).eq.0)goto 10371                                          
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10371                                          
      if(mm(k) .ne. 0)goto 10391                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10372                                           
      do 10401 j=1,ni                                                   
      if(ju(j).eq.0)goto 10401                                          
      if(mm(j) .eq. 0)goto 10421                                        
      c(j,nin)=c(k,mm(j))                                               
      goto 10401                                                        
10421 continue                                                          
      if(j .ne. k)goto 10441                                            
      c(j,nin)=xv(j)                                                    
      goto 10401                                                        
10441 continue                                                          
      c(j,nin)=dot_product(x(:,j),x(:,k))                               
10401 continue                                                          
      continue                                                          
      mm(k)=nin                                                         
      ia(nin)=k                                                         
10391 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 10451 j=1,ni                                                   
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                           
10451 continue                                                          
      continue                                                          
10371 continue                                                          
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
      do 10491 l=1,nin                                                  
      k=ia(l)                                                           
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10491                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 10501 j=1,nin                                                  
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                              
10501 continue                                                          
      continue                                                          
10491 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10482                                          
      if(nlp .le. maxit)goto 10521                                      
      jerr=-m                                                           
      return                                                            
10521 continue                                                          
      goto 10481                                                        
10482 continue                                                          
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                  
      do 10531 j=1,ni                                                   
      if(mm(j).ne.0)goto 10531                                          
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))        
10531 continue                                                          
      continue                                                          
      jz=0                                                              
      goto 10351                                                        
10352 continue                                                          
      if(nin .le. nx)goto 10551                                         
      jerr=-10000-m                                                     
      goto 10282                                                        
10551 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 10281                                            
      if(flmin.ge.1.0)goto 10281                                        
      me=0                                                              
      do 10561 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
10561 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10282                                            
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                 
      if(rsq.gt.rsqmax)goto 10282                                       
10281 continue                                                          
10282 continue                                                          
      deallocate(a,mm,c,da)                                             
      return                                                            
      end                                                               
