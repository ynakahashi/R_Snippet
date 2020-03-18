      subroutine standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)  
      integer ju(ni)                                                    
      double precision, dimension (:), allocatable :: v                 
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=w/sum(w)                                                        
      v=sqrt(w)                                                         
      ! intr は intercept なので切片が 0 であるかで判定
      if(intr .ne. 0)goto 10141                                         
      ym=0.0                                                            
      y=v*y                                                             
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                     
      y=y/ys                                                            
      do 10151 j=1,ni                                                   
      if(ju(j).eq.0)goto 10151                                          
      xm(j)=0.0                                                         
      x(:,j)=v*x(:,j)                                                   
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd .eq. 0)goto 10171                                          
      xbq=dot_product(v,x(:,j))**2                                      
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
      goto 10181                                                        
10171 continue                                                          
      xs(j)=1.0                                                         
10181 continue                                                          
      continue                                                          
10151 continue                                                          
      continue                                                          
      goto 10191                                                        
10141 continue                                                          
      do 10201 j=1,ni                                                   
      if(ju(j).eq.0)goto 10201                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=v*(x(:,j)-xm(j))                                           
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
10201 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 10221                                          
      xs=1.0                                                            
      goto 10231                                                        
10221 continue                                                          
      do 10241 j=1,ni                                                   
      if(ju(j).eq.0)goto 10241                                          
      x(:,j)=x(:,j)/xs(j)                                               
10241 continue                                                          
      continue                                                          
      xv=1.0                                                            
10231 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=v*(y-ym)                                                        
      ys=sqrt(dot_product(y,y))                                         
      y=y/ys                                                            
10191 continue                                                          
      continue                                                          
      g=0.0                                                             
      do 10251 j=1,ni                                                   
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                         
10251 continue                                                          
      continue                                                          
      deallocate(v)                                                     
      return                                                            
      end                                                               
