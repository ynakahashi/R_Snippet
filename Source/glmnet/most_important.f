! ab
u=g(k)+ak*xv(k)
v=abs(u)-vp(k)*ab

! dem
a(k)=0.0                                                          
if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)))