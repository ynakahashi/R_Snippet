      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:,:), allocatable :: q               
      double precision, dimension (:), allocatable :: sxp,sxpl          
      double precision, dimension (:), allocatable :: di,v,r,ga         
      double precision, dimension (:,:), allocatable :: b,bs,xv         
      integer, dimension (:), allocatable :: mm,is,ixx                  
      allocate(b(0:ni,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(q(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      exmn=-exmx                                                        
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(is(1:max(nc,ni)),stat=jerr)                              
      if(jerr.ne.0) return                                              
      allocate(sxp(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sxpl(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(di(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      pmax=1.0-pmin                                                     
      emin=pmin/pmax                                                    
      emax=1.0/emin                                                     
      pfm=(1.0+pmin)*pmin                                               
      pfx=(1.0-pmin)*pmax                                               
      vmin=pfm*pmax                                                     
      bta=parm                                                          
      omb=1.0-bta                                                       
      dev1=0.0                                                          
      dev0=0.0                                                          
      do 13361 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 13381                                        
      jerr =8000+ic                                                     
      return                                                            
13381 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 13401                                    
      jerr =9000+ic                                                     
      return                                                            
13401 continue                                                          
      if(intr .ne. 0)goto 13421                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 13431                                                        
13421 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
13431 continue                                                          
      continue                                                          
      b(1:ni,ic)=0.0                                                    
13361 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      ixx=0                                                             
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 13451                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 13461 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
13461 continue                                                          
      continue                                                          
      goto 13471                                                        
13451 continue                                                          
      do 13481 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
13481 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 13501                                         
      b(0,:)=0.0                                                        
      goto 13511                                                        
13501 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
13511 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 13521 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
13521 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 13531 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
13531 continue                                                          
      continue                                                          
13471 continue                                                          
      continue                                                          
      do 13541 ic=1,nc                                                  
      do 13551 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
13551 continue                                                          
      continue                                                          
13541 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      if(kopt .le. 0)goto 13571                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 13591                         
      xv=0.25                                                           
      goto 13601                                                        
13591 continue                                                          
      do 13611 j=1,ni                                                   
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)              
13611 continue                                                          
      continue                                                          
13601 continue                                                          
      continue                                                          
13571 continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 13631                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
13631 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 13641 ic=1,nc                                                  
      r=w*(y(:,ic)-q(:,ic)/sxp)                                         
      do 13651 j=1,ni                                                   
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))        
13651 continue                                                          
      continue                                                          
13641 continue                                                          
      continue                                                          
      do 13661 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 13681                                      
      al=ulam(ilm)                                                      
      goto 13671                                                        
13681 if(ilm .le. 2)goto 13691                                          
      al=al*alf                                                         
      goto 13671                                                        
13691 if(ilm .ne. 1)goto 13701                                          
      al=big                                                            
      goto 13711                                                        
13701 continue                                                          
      al0=0.0                                                           
      do 13721 j=1,ni                                                   
      if(ju(j).eq.0)goto 13721                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
13721 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
13711 continue                                                          
13671 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 13731 k=1,ni                                                   
      if(ixx(k).eq.1)goto 13731                                         
      if(ju(k).eq.0)goto 13731                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
13731 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
13741 continue                                                          
      ix=0                                                              
      jx=ix                                                             
      ig=0                                                              
      do 13751 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      xmz=0.0                                                           
      do 13761 i=1,no                                                   
      pic=q(i,ic)/sxp(i)                                                
      if(pic .ge. pfm)goto 13781                                        
      pic=0.0                                                           
      v(i)=0.0                                                          
      goto 13771                                                        
13781 if(pic .le. pfx)goto 13791                                        
      pic=1.0                                                           
      v(i)=0.0                                                          
      goto 13801                                                        
13791 continue                                                          
      v(i)=w(i)*pic*(1.0-pic)                                           
      xmz=xmz+v(i)                                                      
13801 continue                                                          
13771 continue                                                          
      r(i)=w(i)*(y(i,ic)-pic)                                           
13761 continue                                                          
      continue                                                          
      if(xmz.le.vmin)goto 13751                                         
      ig=1                                                              
      if(kopt .ne. 0)goto 13821                                         
      do 13831 j=1,ni                                                   
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                 
13831 continue                                                          
      continue                                                          
13821 continue                                                          
      continue                                                          
13841 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13851 k=1,ni                                                   
      if(ixx(k).eq.0)goto 13851                                         
      bk=b(k,ic)                                                        
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13871                                         
      b(k,ic)=0.0                                                       
      goto 13881                                                        
13871 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
13881 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 13851                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r=r-d*v*x(:,k)                                                    
      if(mm(k) .ne. 0)goto 13901                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 13921                                         
      jx=1                                                              
      goto 13852                                                        
13921 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
13901 continue                                                          
13851 continue                                                          
13852 continue                                                          
      if(jx.gt.0)goto 13842                                             
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13941                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13941 continue                                                          
      if(dlx.lt.shr)goto 13842                                          
      if(nlp .le. maxit)goto 13961                                      
      jerr=-ilm                                                         
      return                                                            
13961 continue                                                          
      continue                                                          
13971 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13981 l=1,nin                                                  
      k=m(l)                                                            
      bk=b(k,ic)                                                        
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 14001                                         
      b(k,ic)=0.0                                                       
      goto 14011                                                        
14001 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
14011 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 13981                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r=r-d*v*x(:,k)                                                    
13981 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 14031                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
14031 continue                                                          
      if(dlx.lt.shr)goto 13972                                          
      if(nlp .le. maxit)goto 14051                                      
      jerr=-ilm                                                         
      return                                                            
14051 continue                                                          
      goto 13971                                                        
13972 continue                                                          
      goto 13841                                                        
13842 continue                                                          
      if(jx.gt.0)goto 13752                                             
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                         
      if(ix .ne. 0)goto 14071                                           
      do 14081 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14101             
      ix=1                                                              
      goto 14082                                                        
14101 continue                                                          
14081 continue                                                          
14082 continue                                                          
14071 continue                                                          
      do 14111 i=1,no                                                   
      fi=b(0,ic)+g(i,ic)                                                
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))      
      fi=min(max(exmn,fi),exmx)                                         
      sxp(i)=sxp(i)-q(i,ic)                                             
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                 
      sxp(i)=sxp(i)+q(i,ic)                                             
14111 continue                                                          
      continue                                                          
13751 continue                                                          
13752 continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      di=s                                                              
      do 14121 j=1,nin                                                  
      l=m(j)                                                            
      if(vp(l) .gt. 0.0)goto 14141                                      
      s=sum(b(l,:))/nc                                                  
      goto 14151                                                        
14141 continue                                                          
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                  
14151 continue                                                          
      continue                                                          
      b(l,:)=b(l,:)-s                                                   
      di=di-s*x(:,l)                                                    
14121 continue                                                          
      continue                                                          
      di=exp(di)                                                        
      sxp=sxp*di                                                        
      do 14161 ic=1,nc                                                  
      q(:,ic)=q(:,ic)*di                                                
14161 continue                                                          
      continue                                                          
      if(jx.gt.0)goto 13742                                             
      if(ig.eq.0)goto 13742                                             
      if(ix .ne. 0)goto 14181                                           
      do 14191 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14191                                         
      if(ju(k).eq.0)goto 14191                                          
      ga(k)=0.0                                                         
14191 continue                                                          
      continue                                                          
      do 14201 ic=1,nc                                                  
      r=w*(y(:,ic)-q(:,ic)/sxp)                                         
      do 14211 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14211                                         
      if(ju(k).eq.0)goto 14211                                          
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                       
14211 continue                                                          
      continue                                                          
14201 continue                                                          
      continue                                                          
      do 14221 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14221                                         
      if(ju(k).eq.0)goto 14221                                          
      if(ga(k) .le. al1*vp(k))goto 14241                                
      ixx(k)=1                                                          
      ix=1                                                              
14241 continue                                                          
14221 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 13742                                                        
14181 continue                                                          
      goto 13741                                                        
13742 continue                                                          
      if(jx .le. 0)goto 14261                                           
      jerr=-10000-ilm                                                   
      goto 13662                                                        
14261 continue                                                          
      devi=0.0                                                          
      do 14271 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 14281 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 14281                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
14281 continue                                                          
      continue                                                          
14271 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ig.eq.0)goto 13662                                             
      if(ilm.lt.mnl)goto 13661                                          
      if(flmin.ge.1.0)goto 13661                                        
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13662          
      if(dev(ilm).gt.devmax)goto 13662                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13662                          
13661 continue                                                          
13662 continue                                                          
      g=log(q)                                                          
      do 14291 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
14291 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                        
      return                                                            
      end                                 