#### bakfit
subroutine bakfit(x,npetc,y,w,which,spar,dof,match,nef,
			etal,s,eta,beta,var,tol,
			qr,qraux,qpivot,effect,work)
#integer npetc(7)
#1:n
#2:p
#3:q
#4:ifvar
#5:nit
#6:maxit
#7:qrank
#subroutine bakfit(x,n,p,y,w,q,which,spar,dof,match,nef,
#			etal,s,eta,beta,var,ifvar,tol,nit,maxit,
#			qr,qraux,qrank,qpivot,work)
#This subroutine fits an additive spline fit to y
#All arguments are either double precision or integer


# bakfit uses the modified backfitting algorithm described in Buja, Hastie 
# and Tibshirani, Annals of Statistics, 1989. It calls splsm, and some
# linpack based routines


# This was written by Trevor Hastie in 1990
# It has been modified from the S3 version by Trevor Hastie
# in March 2005, to accommodate the modified sbart routine in R
# Note that spar has changed, and we change it here to conform with
# the smooth.spline routine in R







#OUTPUT
#
#x,y,w,n,p,which,q,maxit,match,nef  are untouched
#spar	for each element of spar:
# 	if spar was 0 and dof was 0, then spar is that spar 
#			that minimized gcv
#	if spar was 0 and dof > 0, then spar is that which achieves dof
#dof	the dof of the fitted smooth. Note: even if dof was given
#		as 4, it will be returned as say 3.995 which is what
#		spar produces
#etal	double length n linear component of the fit
#s	double n by q nonlinear part of the smooth functions
#eta	double length n fitted values
#beta	double length p linear coefficients
#       So, the centered fitted functions are:
#			 b(j)*(x(i,j)-mean(x(.,j)) +s(i,j)
#		where j is an element of which
#var	double n by q 
#	if ifvar was .true.
#		the unscaled variance elements for the NONLINEAR
#		and UNIQUE part of s, in the order of sort(unique(x))
#		var is lev(i)/w(i) -h(i)/w where h(i) is the hat element from
#		the simple weighted least squares fit. This is  used in gamcov
#
#nit	number of iterations used
#qr etc the qr is returned

# 変数宣言
implicit double precision(a-h,o-z) # i~n 以外は double
logical ifvar
integer npetc(7),iter
integer n,p,q,which(*),match(*),nef(*),nit,maxit,qrank,qpivot(*)
double precision x(*),y(*),w(*),spar(*),dof(*),
			etal(*),s(*),eta(*),beta(*),var(*),tol,
			qr(*),qraux(*),effect(*),work(*)

# 代入
#INPUT
#
#x	double dim  n by p ; x variables, includes constant
#n	integer number of rows in x
#p	integer number of columns of x
#y	double length n	; y variable for smoothing
#w	double length n ; prior weights for smoothing, > 0
#q	integer number of nonlinear terms
#which 	integer length q indices of columns of x for nonlinear fits
#spar	double length q spars for smoothing; see below
#dof	double length q dof for/from smoothing; see below
#match 	integer n by q  matrix of match'es; see below
#nef	integer q vector of nef's; see below
#s	double n by q nonlinear part of the smooth functions
#		used as starting values. the linear part is
#		irrelevant
#ifvar	logical should the variance information be computed
#tol	double tolerance for backfitting convergence; 0.0005 is good
#maxit	integer maximum number of iterations; 15 is good
#qr	double n by p weighted qr decomposition of x
#qraux 	double p belongs with qr
#qrank	integer rank of x ; if qrank=0, then bakfit computes qr and qraux
#qpivot	integer p the columns of qr are rearranged according to pivot
#effec  double n effect vector
#work	double 
#	Let nk=max(nef)+2, then 
# 	work should be (10+2*4)*nk+5*nef+5*n+15 +q double

n=npetc(1) # num row
p=npetc(2) # num col
q=npetc(3) # num nonlinear term
ifvar=.false.
if(npetc(4)==1)ifvar=.true. # npetc の四番目は se (npetc = as.integer(c(n, p, length(which), se, 0, maxit, 0)))
maxit=npetc(6) # max iteration
qrank=npetc(7) # q rank

## 平滑化対象の変数の自由度を work に入れる
do i=1,q{work(i)=dof(i)}

## backf1 を呼び出す
#BELOW
#
#the following comments come from documentation for splsm
#	they apply to each element of spar,dof match etc
#spar	double smoothing parameter -1.5 <spar<1.5; default is 1
#dof	double  equivalent degrees of freedom
#		if dof is 0, spar is used
#               if 0< dof <1, dof = 1                         
#		if dof >=1, dof is used
#		note: dof does not use the constant term
#match 	integer length n -- in S language x[i] == sort(unique(x)[match[i]]
#		match is produced by subroutine namat
#nef	number of unique elements in x; so match has values between 1 and nef+1
#		missing data are given the match number nef+1
#work	double workspace of length (10+2*4)*(nef+2)+5*nef+n+15
call backf1(x,n,p,y,w,q,which,spar,dof,
			match, # match がわからない。 namat モジュールで作成されるよう。
			nef, # num unique element
            etal,s,eta,beta,var,ifvar,tol,nit,maxit,
            qr,qraux,qrank,qpivot,effect,work(q+1),work(q+n+1),
            work(q+2*n+1),work(q+3*n+1),work(q+4*n+1))

npetc(7)=qrank
return
end



#### backf1
subroutine backf1(x,n,p,y,w,q,which,spar,dof,match,nef,
			etal,s,eta,beta,var,ifvar,tol,nit,maxit,
			qr,qraux,qrank,qpivot,effect,z,old,sqwt,sqwti,work)
implicit double precision(a-h,o-z)
logical ifvar
integer n,p,q,which(q),match(n,q),nef(q),nit,maxit,qrank,qpivot(p)
double precision x(n,p),y(n),w(n),spar(q),dof(q),
			etal(n),s(n,q),eta(n),beta(p),var(n,q),tol,
			qr(n,p),qraux(p),effect(n),work(*)
double precision z(*),old(*),dwrss,ratio
double precision sqwt(n),sqwti(n)
logical anyzwt
double precision deltaf, normf,onedm7
integer job,info

onedm7=1d-7
job=1101;info=1
if(q==0)maxit=1 # nonlinear termの個数なので、0なら反復せず終了
ratio=1d0 # ループの終了基準。1.0

### weight の平方根をとって逆数にする
# fix up sqy's for weighted problems.
anyzwt=.false.
do i=1,n{
	if(w(i)>0d0){
		sqwt(i)=dsqrt(w(i)) ## Fortranの組み込み関数。sqrt()
		sqwti(i)=1d0/sqwt(i) ## weight の逆数
	}
	else{
		sqwt(i)=0d0
		sqwti(i)=0d0
		anyzwt=.true.
		}
	}

### qrank が 0 のときは、x に重みを加えた上で QR 分解する
# if qrank > 0 then qr etc contain the qr decomposition
# else bakfit computes it. 
if(qrank==0){
	do i=1,n{ # 全データ
		do j=1,p{ # 全変数
			### 重みを乗じた x を生成
			qr(i,j)=x(i,j)*sqwt(i) # qr は n*p の行列
			}
		}
	do j=1,p{qpivot(j)=j} ### qvivot は長さ p の integer。1~p の数字を格納

	### QR分解
	### どこで定義されてる？
	### https://github.com/cran/gam/blob/master/inst/ratfor/linear.r
	call dqrdca(qr,n,n,p,qraux,qpivot,work,qrank,onedm7)
	}

### eta に 非線形項 s(i,j) を加算する
### ただしここで eta は fitted value ではなく 0 スタート
do i=1,n{ # 全データ
	#s	double n by q nonlinear part of the smooth functions
	#		used as starting values. the linear part is
	#		irrelevant
	eta(i)=0d0 # 0.0 にリセット 
	for(j=1;j<=q;j=j+1){ # 1~q まで
		eta(i)=eta(i)+s(i,j) # 平滑化部分を順次加算していく
		}
	}


### ループ開始
nit=0
while ((ratio > tol )&(nit < maxit)){ # ratio のデフォルトは 1.0 。 tol は多分 0.0005
	# first the linear fit
	deltaf=0d0
	nit=nit+1 # イテレータに +1

	### まずは Linear part + error 部分に最小二乗法を当てはめる

	### y から非線形項を減じたもの（Linear part + error）に重みを乗じて z を更新(これを y として dqrsl に渡す)
	# etal を old に格納。ループ１回目の時点では etal は Null？
	do i=1,n{
		z(i)=(y(i)-eta(i))*sqwt(i) 
		old(i)=etal(i)
	}

	### https://github.com/wch/r-source/blob/trunk/src/appl/dqrsl.f
	### dqrsl の6個目の引数として z を渡しているが、これは y になる
	### etal は11番目の引数として渡しているが、これは xb 。
	#	call dqrsl1(qr,dq,qraux,qrank,sqz,one,work(1),etal,two,three)
	#job=1101 -- computes fits, effects and beta
	call dqrsl(qr,n,n,qrank,qraux,z,work(1),effect(1),beta,
		work(1),etal,job,info)

# now unsqrt the fits
#Note: we dont have to fix up the zero weights till the end, since their fits
#are always immaterial to the computation

	### 重みの逆数で戻し、 etal を更新
	### etal は y の予測値の線形部分
	#etal	double length n linear component of the fit
	do i=1,n{
		etal(i)=etal(i)*sqwti(i)
		}

	### ここで平滑化対象の変数についてループ
	# now a single non-linear backfitting loop 
	for(k=1;k<=q;k=k+1){ # 平滑化対象の変数ごとに
		j=which(k)
		do i=1,n{
			### old に k 番目の平滑化対象変数の値を入れる
			old(i)=s(i,k)
			### y から etal (Linear part)と eta (Nonlinear part) を引き、old (k 番目の Nonlinear part) を足している
			### 残差 + k 番目の平滑化変数の値になっている
			z(i)=y(i)-etal(i)-eta(i)+old(i)
		}
                ### df は 0 にリセットされてしまう
				# this uses spar to set smoothing after iteration 1
                if(nit>1){dof(k)=0d0}
		
		### splsm を呼びだす
		### 切片と k 番目の平滑化変数を x として渡す
		### 9番目の引数が平滑化変数（splsm 内では smo）なので s が更新される
		call splsm(x(1,j),z,w,n,match(1,k),nef(k),spar(k),
			dof(k),s(1,k),s0,var(1,k),ifvar,work)
		### 	
		do i=1,n{ # データごとに
			### eta を更新
			### 古い eta にk番目の非線形部分を加算して、古い非線形部分を減じる)
			eta(i)=eta(i)+s(i,k)-old(i)
			### etal を更新(古い etal に s0 を乗じる。)
			### s0 は weighted mean of y
			etal(i)=etal(i)+s0
			}
		### deltaf を更新（while ループの打ち切り判定）
		deltaf=deltaf+dwrss(n,old,s(1,k),w)
		} # ここまで for ループ
	
	### normf の更新（while ループの打ち切り判定）
	normf=0d0
	do i=1,n{
		normf=normf+w(i)*eta(i)*eta(i)
		}
	### ratio（while ループの判定に使われる）を計算
	if(normf>0d0){
		ratio=dsqrt(deltaf/normf)
		}
	 else {ratio = 0d0}
#         call DBLEPR("ratio",-1,ratio,1)
	} # ここまで while ループ
#now package up the results
do j=1,p {work(j)=beta(j)}
do j=1,p {beta(qpivot(j))=work(j)}
if(anyzwt){
	do i=1,n {
		if(w(i) <= 0d0){
			etal(i)=0d0
			do j=1,p{
				etal(i)=etal(i)+beta(j)*x(i,j)
				}
			}
		}
	}

### 非線形項 + 線形項として eta を再定義
do i=1,n
	eta(i)=eta(i)+etal(i)
	
do j=1,q {
	call unpck(n,nef(j),match(1,j),var(1,j),old)
	do i=1,n {var(i,j)=old(i)}
	}

return
end
