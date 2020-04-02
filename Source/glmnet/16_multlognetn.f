      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam
     *),cl(2,ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(
     *ni)
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:,:), allocatable :: q,r,b,bs        
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk,del
      integer, dimension (:), allocatable :: mm,is,ixx,isc              
      allocate(b(0:ni,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(q(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(r(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      exmn=-exmx                                                        
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(is(1:max(nc,ni)),stat=jerr)                              
      if(jerr.ne.0) return                                              
      allocate(sxp(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sxpl(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(gk(1:nc),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(del(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(isc(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      pmax=1.0-pmin                                                     
      emin=pmin/pmax                                                    
      emax=1.0/emin                                                     
      bta=parm                                                          
      omb=1.0-bta                                                       
      dev1=0.0                                                          
      dev0=0.0                                                          
      do 21441 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 21461                                        
      jerr =8000+ic                                                     
      return                                                            
21461 continue                                                          
      if(q0 .lt. pmax)goto 21481                                        
      jerr =9000+ic                                                     
      return                                                            
21481 continue                                                          
      if(intr .ne. 0)goto 21501                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 21511                                                        
21501 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
21511 continue                                                          
      continue                                                          
      b(1:ni,ic)=0.0                                                    
21441 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      ixx=0                                                             
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 21531                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 21541 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
21541 continue                                                          
      continue                                                          
      goto 21551                                                        
21531 continue                                                          
      do 21561 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
21561 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 21581                                         
      b(0,:)=0.0                                                        
      goto 21591                                                        
21581 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
21591 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 21601 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
21601 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 21611 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
21611 continue                                                          
      continue                                                          
21551 continue                                                          
      continue                                                          
      do 21621 ic=1,nc                                                  
      do 21631 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
21631 continue                                                          
      continue                                                          
21621 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 21651                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
21651 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 21661 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 21671 j=1,ni                                                   
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2         
21671 continue                                                          
      continue                                                          
21661 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 21681 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 21701                                      
      al=ulam(ilm)                                                      
      goto 21691                                                        
21701 if(ilm .le. 2)goto 21711                                          
      al=al*alf                                                         
      goto 21691                                                        
21711 if(ilm .ne. 1)goto 21721                                          
      al=big                                                            
      goto 21731                                                        
21721 continue                                                          
      al0=0.0                                                           
      do 21741 j=1,ni                                                   
      if(ju(j).eq.0)goto 21741                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
21741 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
21731 continue                                                          
21691 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 21751 k=1,ni                                                   
      if(ixx(k).eq.1)goto 21751                                         
      if(ju(k).eq.0)goto 21751                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
21751 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
21761 continue                                                          
      ix=0                                                              
      jx=ix                                                             
      kx=jx                                                             
      t=0.0                                                             
      do 21771 ic=1,nc                                                  
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                    
21771 continue                                                          
      continue                                                          
      if(t .ge. eps)goto 21791                                          
      kx=1                                                              
      goto 21762                                                        
21791 continue                                                          
      t=2.0*t                                                           
      alt=al1/t                                                         
      al2t=al2/t                                                        
      do 21801 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                 
      d=0.0                                                             
      if(intr.ne.0) d=sum(r(:,ic))                                      
      if(d .eq. 0.0)goto 21821                                          
      b(0,ic)=b(0,ic)+d                                                 
      r(:,ic)=r(:,ic)-d*w                                               
      dlx=max(dlx,d**2)                                                 
21821 continue                                                          
21801 continue                                                          
      continue                                                          
      continue                                                          
21831 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 21841 k=1,ni                                                   
      if(ixx(k).eq.0)goto 21841                                         
      gkn=0.0                                                           
      do 21851 ic=1,nc                                                  
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                  
      gkn=gkn+gk(ic)**2                                                 
21851 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      del=b(k,:)                                                        
      if(u .gt. 0.0)goto 21871                                          
      b(k,:)=0.0                                                        
      goto 21881                                                        
21871 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
21881 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21841                             
      do 21891 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                  
21891 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 21911                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 21931                                         
      jx=1                                                              
      goto 21842                                                        
21931 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
21911 continue                                                          
21841 continue                                                          
21842 continue                                                          
      if(jx.gt.0)goto 21832                                             
      if(dlx.lt.shr)goto 21832                                          
      if(nlp .le. maxit)goto 21951                                      
      jerr=-ilm                                                         
      return                                                            
21951 continue                                                          
      continue                                                          
21961 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 21971 l=1,nin                                                  
      k=m(l)                                                            
      gkn=0.0                                                           
      do 21981 ic=1,nc                                                  
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                  
      gkn=gkn+gk(ic)**2                                                 
21981 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      del=b(k,:)                                                        
      if(u .gt. 0.0)goto 22001                                          
      b(k,:)=0.0                                                        
      goto 22011                                                        
22001 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22011 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21971                             
      do 22021 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                  
22021 continue                                                          
      continue                                                          
21971 continue                                                          
      continue                                                          
      if(dlx.lt.shr)goto 21962                                          
      if(nlp .le. maxit)goto 22041                                      
      jerr=-ilm                                                         
      return                                                            
22041 continue                                                          
      goto 21961                                                        
21962 continue                                                          
      goto 21831                                                        
21832 continue                                                          
      if(jx.gt.0)goto 21762                                             
      do 22051 ic=1,nc                                                  
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                             
      if(ix .ne. 0)goto 22071                                           
      do 22081 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22101                
      ix=1                                                              
      goto 22082                                                        
22101 continue                                                          
22081 continue                                                          
22082 continue                                                          
22071 continue                                                          
      do 22111 i=1,no                                                   
      fi=b(0,ic)+g(i,ic)                                                
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))      
      fi=min(max(exmn,fi),exmx)                                         
      sxp(i)=sxp(i)-q(i,ic)                                             
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                 
      sxp(i)=sxp(i)+q(i,ic)                                             
22111 continue                                                          
      continue                                                          
22051 continue                                                          
      continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      if(jx.gt.0)goto 21762                                             
      if(ix .ne. 0)goto 22131                                           
      do 22141 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22141                                         
      if(ju(k).eq.0)goto 22141                                          
      ga(k)=0.0                                                         
22141 continue                                                          
      continue                                                          
      do 22151 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 22161 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22161                                         
      if(ju(k).eq.0)goto 22161                                          
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                        
22161 continue                                                          
      continue                                                          
22151 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 22171 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22171                                         
      if(ju(k).eq.0)goto 22171                                          
      if(ga(k) .le. al1*vp(k))goto 22191                                
      ixx(k)=1                                                          
      ix=1                                                              
22191 continue                                                          
22171 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 21762                                                        
22131 continue                                                          
      goto 21761                                                        
21762 continue                                                          
      if(kx .le. 0)goto 22211                                           
      jerr=-20000-ilm                                                   
      goto 21682                                                        
22211 continue                                                          
      if(jx .le. 0)goto 22231                                           
      jerr=-10000-ilm                                                   
      goto 21682                                                        
22231 continue                                                          
      devi=0.0                                                          
      do 22241 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 22251 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 22251                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
22251 continue                                                          
      continue                                                          
22241 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 21681                                          
      if(flmin.ge.1.0)goto 21681                                        
      me=0                                                              
      do 22261 j=1,nin                                                  
      if(a(j,1,ilm).ne.0.0) me=me+1                                     
22261 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 21682                                            
      if(dev(ilm).gt.devmax)goto 21682                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21682                          
21681 continue                                                          
21682 continue                                                          
      g=log(q)                                                          
      do 22271 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
22271 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                 
      return                                                            
      end