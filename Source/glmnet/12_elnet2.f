      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam)     
      double precision rsqo(nlam),almo(nlam),xv(ni)                     
      double precision cl(2,ni)                                         
      integer ju(ni),ia(nx),kin(nlam)                                   
      double precision, dimension (:), allocatable :: a,g               
      integer, dimension (:), allocatable :: mm,ix                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(ix(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      ix=0                                                              
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 10771                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
10771 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      alm=0.0                                                           
      do 10781 j=1,ni                                                   
      if(ju(j).eq.0)goto 10781                                          
      g(j)=abs(dot_product(y,x(:,j)))                                   
10781 continue                                                          
      continue                                                          
      do 10791 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 10811                                      
      alm=ulam(m)                                                       
      goto 10801                                                        
10811 if(m .le. 2)goto 10821                                            
      alm=alm*alf                                                       
      goto 10801                                                        
10821 if(m .ne. 1)goto 10831                                            
      alm=big                                                           
      goto 10841                                                        
10831 continue                                                          
      alm0=0.0                                                          
      do 10851 j=1,ni                                                   
      if(ju(j).eq.0)goto 10851                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
10851 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
10841 continue                                                          
10801 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 10861 k=1,ni                                                   
      if(ix(k).eq.1)goto 10861                                          
      if(ju(k).eq.0)goto 10861                                          
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                    
10861 continue                                                          
      continue                                                          
      continue                                                          
10871 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10891 k=1,ni                                                   
      if(ix(k).eq.0)goto 10891                                          
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10891                                          
      if(mm(k) .ne. 0)goto 10911                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10892                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
10911 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y=y-del*x(:,k)                                                    
      dlx=max(xv(k)*del**2,dlx)                                         
10891 continue                                                          
10892 continue                                                          
      if(nin.gt.nx)goto 10872                                           
      if(dlx .ge. thr)goto 10931                                        
      ixx=0                                                             
      do 10941 k=1,ni                                                   
      if(ix(k).eq.1)goto 10941                                          
      if(ju(k).eq.0)goto 10941                                          
      g(k)=abs(dot_product(y,x(:,k)))                                   
      if(g(k) .le. ab*vp(k))goto 10961                                  
      ix(k)=1                                                           
      ixx=1                                                             
10961 continue                                                          
10941 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 10872                                                        
10931 continue                                                          
      if(nlp .le. maxit)goto 10981                                      
      jerr=-m                                                           
      return                                                            
10981 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
10991 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11001 l=1,nin                                                  
      k=ia(l)                                                           
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11001                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y=y-del*x(:,k)                                                    
      dlx=max(xv(k)*del**2,dlx)                                         
11001 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10992                                          
      if(nlp .le. maxit)goto 11021                                      
      jerr=-m                                                           
      return                                                            
11021 continue                                                          
      goto 10991                                                        
10992 continue                                                          
      jz=0                                                              
      goto 10871                                                        
10872 continue                                                          
      if(nin .le. nx)goto 11041                                         
      jerr=-10000-m                                                     
      goto 10792                                                        
11041 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 10791                                            
      if(flmin.ge.1.0)goto 10791                                        
      me=0                                                              
      do 11051 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
11051 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10792                                            
      if(rsq-rsq0.lt.sml*rsq)goto 10792                                 
      if(rsq.gt.rsqmax)goto 10792                                       
10791 continue                                                          
10792 continue                                                          
      deallocate(a,mm,g,ix)                                             
      return                                                            
      end                                                               
