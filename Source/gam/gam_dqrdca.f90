https://github.com/cran/gam/blob/master/inst/ratfor/linear.r

# qr decomposition, modified from linpack routines to give stable
# ordering and rank estimation
subroutine dqrdca(x,ldx,n,p,qraux,jpvt,work,rank,eps)
integer ldx,n,p,rank
integer jpvt(*)
double precision x(ldx,*),qraux(*),work(*),eps
integer j,jj,jp,l,lup,curpvt
double precision dnrm2,tt
double precision ddot,nrmxl,t,ww
do j=1,p {
	qraux(j) = dnrm2(n,x(1,j),1)
	work(j) = qraux(j); work(j+p) =  qraux(j)
}
l=1; lup = min0(n,p); curpvt = p
while(l<=lup) {
	qraux(l) = 0.0d0
	nrmxl = dnrm2(n-l+1,x(l,l),1)
	t = work(l+p); if(t > 0.)t = nrmxl/t
	if(t < eps){
		call dshift(x,ldx,n,l,curpvt)
		jp = jpvt(l); t=qraux(l); tt=work(l); ww = work(l+p)
		for(j=l+1; j<=curpvt; j=j+1){
			jj=j-1
			jpvt(jj)=jpvt(j); qraux(jj)=qraux(j)
			work(jj)=work(j); work(jj+p) = work(j+p)
		}
		jpvt(curpvt)=jp; qraux(curpvt)=t;
		work(curpvt)=tt; work(curpvt+p) = ww
		curpvt=curpvt-1; if(lup>curpvt)lup=curpvt
	}
	else {
		if(l==n)break
		if (x(l,l)!=0.0d0)
			nrmxl = dsign(nrmxl,x(l,l))
		call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
		x(l,l) = 1.0d0+x(l,l)
		for(j=l+1; j<=curpvt; j=j+1) {
			t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
			call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
			if (qraux(j)!=0.0d0) {
				tt = 1.0d0-(dabs(x(l,j))/qraux(j))**2
				tt = dmax1(tt,0.0d0)
				t = tt
				tt = 1.0d0+0.05d0*tt*(qraux(j)/work(j))**2
				if (tt!=1.0d0)
					qraux(j) = qraux(j)*dsqrt(t)
				else {
					qraux(j) = dnrm2(n-l,x(l+1,j),1)
					work(j) = qraux(j)
				}
			}
		}
		qraux(l) = x(l,l)
		x(l,l) = -nrmxl
		l=l+1
	}
}
rank = lup
return
end