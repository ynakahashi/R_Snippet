      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)          
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga  
      integer, dimension (:), allocatable :: mm,ixx                     
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      allocate(b(0:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(q(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      fmax=log(1.0/pmin-1.0)                                            
      fmin=-fmax                                                        
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                   
      bta=parm                                                          
      omb=1.0-bta                                                       
      q0=dot_product(w,y)                                               
      if(q0 .gt. pmin)goto 12681                                        
      jerr=8001                                                         
      return                                                            
12681 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 12701                                    
      jerr=9001                                                         
      return                                                            
12701 continue                                                          
      if(intr.eq.0.0) q0=0.5                                            
      ixx=0                                                             
      al=0.0                                                            
      bz=0.0                                                            
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                 
      if(nonzero(no,g) .ne. 0)goto 12721                                
      vi=q0*(1.0-q0)                                                    
      b(0)=bz                                                           
      v=vi*w                                                            
      r=w*(y-q0)                                                        
      q=q0                                                              
      xmz=vi                                                            
      dev1=-(bz*q0+log(1.0-q0))                                         
      goto 12731                                                        
12721 continue                                                          
      b(0)=0.0                                                          
      if(intr .eq. 0)goto 12751                                         
      b(0)=azero(no,y,g,w,jerr)                                         
      if(jerr.ne.0) return                                              
12751 continue                                                          
      q=1.0/(1.0+exp(-b(0)-g))                                          
      v=w*q*(1.0-q)                                                     
      r=w*(y-q)                                                         
      xmz=sum(v)                                                        
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                     
12731 continue                                                          
      continue                                                          
      if(kopt .le. 0)goto 12771                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 12791                         
      xv=0.25                                                           
      goto 12801                                                        
12791 continue                                                          
      do 12811 j=1,ni                                                   
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                
12811 continue                                                          
      continue                                                          
12801 continue                                                          
      continue                                                          
12771 continue                                                          
      dev0=dev1                                                         
      do 12821 i=1,no                                                   
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                     
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))           
12821 continue                                                          
      continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 12841                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
12841 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      b(1:ni)=0.0                                                       
      shr=shri*dev0                                                     
      do 12851 j=1,ni                                                   
      if(ju(j).eq.0)goto 12851                                          
      ga(j)=abs(dot_product(r,x(:,j)))                                  
12851 continue                                                          
      continue                                                          
      do 12861 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 12881                                      
      al=ulam(ilm)                                                      
      goto 12871                                                        
12881 if(ilm .le. 2)goto 12891                                          
      al=al*alf                                                         
      goto 12871                                                        
12891 if(ilm .ne. 1)goto 12901                                          
      al=big                                                            
      goto 12911                                                        
12901 continue                                                          
      al0=0.0                                                           
      do 12921 j=1,ni                                                   
      if(ju(j).eq.0)goto 12921                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
12921 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
12911 continue                                                          
12871 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 12931 k=1,ni                                                   
      if(ixx(k).eq.1)goto 12931                                         
      if(ju(k).eq.0)goto 12931                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
12931 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
12941 continue                                                          
      bs(0)=b(0)                                                        
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                             
      if(kopt .ne. 0)goto 12961                                         
      do 12971 j=1,ni                                                   
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                    
12971 continue                                                          
      continue                                                          
12961 continue                                                          
      continue                                                          
12981 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 12991 k=1,ni                                                   
      if(ixx(k).eq.0)goto 12991                                         
      bk=b(k)                                                           
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13011                                         
      b(k)=0.0                                                          
      goto 13021                                                        
13011 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
13021 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 12991                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r=r-d*v*x(:,k)                                                    
      if(mm(k) .ne. 0)goto 13041                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 12992                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
13041 continue                                                          
12991 continue                                                          
12992 continue                                                          
      if(nin.gt.nx)goto 12982                                           
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13061                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13061 continue                                                          
      if(dlx.lt.shr)goto 12982                                          
      if(nlp .le. maxit)goto 13081                                      
      jerr=-ilm                                                         
      return                                                            
13081 continue                                                          
      continue                                                          
13091 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13101 l=1,nin                                                  
      k=m(l)                                                            
      bk=b(k)                                                           
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13121                                         
      b(k)=0.0                                                          
      goto 13131                                                        
13121 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
13131 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 13101                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r=r-d*v*x(:,k)                                                    
13101 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13151                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13151 continue                                                          
      if(dlx.lt.shr)goto 13092                                          
      if(nlp .le. maxit)goto 13171                                      
      jerr=-ilm                                                         
      return                                                            
13171 continue                                                          
      goto 13091                                                        
13092 continue                                                          
      goto 12981                                                        
12982 continue                                                          
      if(nin.gt.nx)goto 12942                                           
      do 13181 i=1,no                                                   
      fi=b(0)+g(i)                                                      
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))         
      if(fi .ge. fmin)goto 13201                                        
      q(i)=0.0                                                          
      goto 13191                                                        
13201 if(fi .le. fmax)goto 13211                                        
      q(i)=1.0                                                          
      goto 13221                                                        
13211 continue                                                          
      q(i)=1.0/(1.0+exp(-fi))                                           
13221 continue                                                          
13191 continue                                                          
13181 continue                                                          
      continue                                                          
      v=w*q*(1.0-q)                                                     
      xmz=sum(v)                                                        
      if(xmz.le.vmin)goto 12942                                         
      r=w*(y-q)                                                         
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                        
      ix=0                                                              
      do 13251 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                        
      ix=1                                                              
      goto 13252                                                        
13251 continue                                                          
13252 continue                                                          
      if(ix .ne. 0)goto 13271                                           
      do 13281 k=1,ni                                                   
      if(ixx(k).eq.1)goto 13281                                         
      if(ju(k).eq.0)goto 13281                                          
      ga(k)=abs(dot_product(r,x(:,k)))                                  
      if(ga(k) .le. al1*vp(k))goto 13301                                
      ixx(k)=1                                                          
      ix=1                                                              
13301 continue                                                          
13281 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 12942                                                        
13271 continue                                                          
13241 continue                                                          
      goto 12941                                                        
12942 continue                                                          
      if(nin .le. nx)goto 13321                                         
      jerr=-10000-ilm                                                   
      goto 12862                                                        
13321 continue                                                          
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                             
      kin(ilm)=nin                                                      
      a0(ilm)=b(0)                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      devi=dev2(no,w,y,q,pmin)                                          
      dev(ilm)=(dev1-devi)/dev0                                         
      if(xmz.le.vmin)goto 12862                                         
      if(ilm.lt.mnl)goto 12861                                          
      if(flmin.ge.1.0)goto 12861                                        
      me=0                                                              
      do 13331 j=1,nin                                                  
      if(a(j,ilm).ne.0.0) me=me+1                                       
13331 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 12862                                            
      if(dev(ilm).gt.devmax)goto 12862                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                          
12861 continue                                                          
12862 continue                                                          
      g=log(q/(1.0-q))                                                  
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                               
      return                                                            
      end