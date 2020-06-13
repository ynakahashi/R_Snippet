      ! chkvars は x の各変数について、 1 行目の値と異なる値が 2 行目以降にあるかを確認している
      ! 異なる値がなければ、全ての値は同じということなので情報がない
      ! 結果は ju に格納されるので、 ju は変数ごとのバラツキを確認している

      subroutine chkvars(no,ni,x,ju)                                    
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni)                                         
      integer ju(ni)                                                    
      do 11061 j=1,ni                                                   
      ju(j)=0                                                           
      t=x(1,j)                                                          
      do 11071 i=2,no                                                   
      if(x(i,j).eq.t)goto 11071                                         
      ju(j)=1                                                           
      goto 11072                                                        
11071 continue                                                          
11072 continue                                                          
11061 continue                                                          
      continue                                                          
      return                                                            
      end              