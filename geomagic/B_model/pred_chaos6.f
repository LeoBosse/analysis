      program plotobs

c     plotting program for a sv models and obsy data
c
c		CAREFULL FORMAT OF MODEL :
c		(NOT A SPLINE MODEL BUT A SUCCESSION OF SV
c		MODELS AT DIFFERENT EPOCHS)
c
c     modified August 1988 & July 1990 & August 5 1991 (AJ) to use plotxy
c     3 January 2001 change to plot model only up to tmax (AJ)
c
      implicit none

      integer lmax, n, nl, nspl, nepochs, jorder, ltrunc, ntrunc
      character*50 path

      real*8 ts,te
      parameter(jorder	= 6,
     &		lmax	= 20,
     &		nspl	= 48,
     &		n       = lmax*(lmax+2),
     &		ltrunc	= 14,
     &		ntrunc  = ltrunc*(ltrunc+2),
     &		nl	= (ltrunc+1)*(ltrunc+2)/2,
     &		path	=
     &'/home/bossel/These/Analysis/src/geomagic/B_model/')
      integer Nangle
      parameter (Nangle=10000)
c      parameter (Nangle=10000)
c      compile with "gfortran pred_chaos6.f -o chaos_auto[N]"

      real*8 pi, r_a,r_c
      parameter(pi	= 4.0*ATAN(1.0),
     &		r_a	= 6371.2,
     &		r_c	= 3485.0)

      real*8 theta_phi_H(Nangle,2)
      real*8 tknts(nspl+jorder),spl(nspl)
      real*8 GT(N,NSPL)

      real*8 g(ntrunc), dg(ntrunc), dx(ntrunc), dy(ntrunc), dz(ntrunc)
      real*8 p(nl),dp(nl)
      real*8 dtheta, dphi

      real*8 ai, alt, dtmod, d, sd, cd, costh, sinth, f,h, t, x, y, z
      real*8 tt1, tt2, o1,o2, splj, theta, phi, rad, tmp

      integer i,j,k,l,m, ii, jj, kk, iplot, nels,
     &	nleft, JJJ, lm, ns, jo, id, nderivs, nob, it, ip
      integer dig1, dig2, dig3

      character*80 file, head
      character*60 name, filename, mod_in
      character*8 lunits
      character*5 obsc,obscode
      character*1 chel(3), resp
      character*4 year
      character*3 strnumb

      data chel/'Z','Y','X'/
      DATA JJJ/1/


001   continue

!      write(6,*) 'magnetic model filename ?'
!      read(5,*) mod_in

      write(6,*) 'epoch ?'
C      read(5,*) t
      t = 2018.5
      write(6,*) 'altitude (km) ?'
C      read(5,*) alt
      alt = 220

c-----
c     read input mfsv model:
c
      filename='CHAOS-6-x6_spline-coefficients.dat.txt'

      OPEN(1,FILE=trim(path)//trim(filename))
      READ(1,1000) HEAD
      READ(1,*) lm,ns,jo
      read(1,*) (TKNTS(K),K=1,ns+jo)
      write(6,"(4f12.4)")tknts
      if(lm.ne.lmax) stop 'bad lmax !'
      if(ns.ne.nspl) stop 'bad nspl !'
      if(jo.ne.jorder) stop 'bad jorder !'

      write(6,"('lmax   = ',i5)") lm
      write(6,"('nspl   = ',i5)") ns
      write(6,"('jorder = ',i5)") jo

      write(6,"('start at  :',f15.4)") tknts(jo)
      write(6,"('finish at :',f15.4)") tknts(ns+1)
      do j=1,ns
         read(1,*) (gt(k,j),k=1,lm*(lm+2))
      enddo
      CLOSE(1)
1000  FORMAT(A80)
1100  FORMAT(4E20.12)

      print*,(gt(k,35),k=1,10)

      write(6,"('mfsv model read')")

      call interv_general(tknts,t,nspl,nleft,jorder)
      print*,nleft

      do j=1,nspl
         spl(j)=0.
      enddo
      call bspline(tknts,t,nspl,jorder,
     &		nleft,spl(nleft-jorder+1))

      do k=1,ntrunc
         g(k)=0.0
         do j=1,jorder
            jj=j+nleft-jorder
            splj=spl(jj)
            g(k) = g(k) + splj*gt(K,JJ)
         enddo
      enddo

      do 450 j=1,nspl
         spl(j)=0.0
450   continue
      do 460 k=1,ntrunc
         dg(k)=0.0
460   continue
      nderivs=1
      call bspline_deriv(tknts,t,nspl,nleft,
     &		spl, jorder, nderivs)

      do 470 j=1,nspl
         do 470 k=1,ntrunc
            dg(k) = dg(k) + spl(j)*gt(K,J)
470   continue

      print*,(g(k),k=1,10)

c theta,phi,alt,nob,id

      open(1100,file = trim(path)//'theta_phi_H.dat', form='formatted')
      do jj=1,Nangle
         read(1100,*) (theta_phi_H(jj,kk),kk=1,2)
      enddo
      close(1100)

      open(33,file=trim(path)//'/Bxyz_H_chaos6.dat')
      do jj=1,Nangle
         theta = theta_phi_H(jj,1)
         phi = theta_phi_H(jj,2)

         call coords(alt,theta,rad,sd,cd)

         sinth=sin(theta)
         costh=cos(theta)
         CALL PLMBAR(P,DP,COSTH,ltrunc,JJJ)

         CALL MAGFDZ(P,DP,THETA,PHI,RAD,Ltrunc,G,DX,DY,DZ,
     &	                x,y,z,h,f,ai,d,
     &			1,SD,CD)
         print*,jj,theta*180/pi, phi*180/pi,x,y,z
         write(33,*) theta,phi,x,y,z
      enddo
      close(33)

      stop
      end
***************************************************************
      SUBROUTINE COORDS(H,THETA,R,SD,CD)
*
      REAL*8 H, THETA, R, SD, CD, PI, B1, B2,CLAT, SLAT
     &     , ONE, TWO, THREE, FOUR, SINTH, COSTH
*
      PI=4.0*ATAN(1.0)
      B1=40680925.0
      B2=40408585.0
      THETA=PI/2-THETA
      CLAT=COS(THETA)
      SLAT=SIN(THETA)
      ONE=B1*CLAT*CLAT
      TWO=B2*SLAT*SLAT
      THREE=ONE+TWO
      FOUR=SQRT(THREE)
      R=SQRT(H*(H+2.0*FOUR)+(B1*ONE+B2*TWO)/THREE)
      CD=(H+FOUR)/R
      SD=(B1-B2)/FOUR*SLAT*CLAT/R
      SINTH=SLAT*CD-CLAT*SD
      COSTH=CLAT*CD+SLAT*SD
      THETA=PI/2.0-ATAN2(SINTH,COSTH)
*
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MAGFDZ(P,DP,THETA,PHI,R,LMAX,G,DX,DY,DZ
     &                             ,X,Y,Z,H,F,I,D,IGEO,SD,CD)
*
** gives field components at radius r
** if igeo=1 then the field elements are rotated to
** the local geodetic frame at the point specified by
** (r,theta,phi) in the spherical coordinate system.
*

      INTEGER IGEO, LMAX, L, K, K1, L1, M
      REAL*8 G(LMAX*(LMAX+2)), DX(LMAX*(LMAX+2)), DY(LMAX*(LMAX+2))
     &     , DZ(LMAX*(LMAX+2)), P((LMAX+1)*(LMAX+2)/2)
     &     , DP((LMAX+1)*(LMAX+2)/2)
      REAL*8 B, THETA, PHI, R, X, Y, Z, H, F
     &     , I, D, SD, CD, SINTH, BB, T, SINT, COST, DXD, DXY
     &     , DZD, DXK, DZK, XS

      B = 6371.2d0/R
      X = 0.0d0
      Y = 0.0d0
      Z = 0.0d0

      SINTH = SIN(THETA)
      IF(ABS(SINTH).LT.1E-10) SINTH=1E-10

      DO 20 L=1,LMAX
       L1=L+1
       BB=B**(L+2)
       K=L*L
       K1=(L*L1)/2+1
       DX(K)=DP(K1)*BB
       DY(K)=0.0d0
       DZ(K)=-P(K1)*L1*BB
       X=X+G(K)*DX(K)
       Z=Z+G(K)*DZ(K)

       DO 20 M=1,L
        T=FLOAT(M)*PHI
        K=L*L+2*M-1
        K1=(L*L1)/2+M+1
        SINT=SIN(T)
        COST=COS(T)
        DXD = DP(K1)*BB
        DX(K) = DXD*COST
        DX(K+1) = DXD*SINT
        X = X + (G(K)*DX(K)) + (G(K+1)*DX(K+1))
        DXY = M*P(K1)*BB/SINTH
        DY(K) = DXY*SINT
        DY(K+1) = -DXY*COST
        Y = Y + (G(K)*DY(K)) + (G(K+1)*DY(K+1))

        DZD = -L1*P(K1)*BB
        DZ(K) = DZD*COST
        DZ(K+1) = DZD*SINT
        Z = Z + (G(K)*DZ(K)) + (G(K+1)*DZ(K+1))
20    CONTINUE
*
** Rotate to geocentric coordinates?
*
      IF(IGEO.EQ.1) THEN
       XS = X
       X = X*CD + Z*SD
       Z = Z*CD - XS*SD
       DO 50 K=1,LMAX*(LMAX+2)
        DXK = DX(K)
        DZK = DZ(K)
        DX(K) = DXK*CD + DZK*SD
        DZ(K) = DZK*CD - DXK*SD
50     CONTINUE
      ENDIF
*
      H=SQRT(X*X+Y*Y)
      F=SQRT(H*H+Z*Z)
      I=ASIN(Z/F)
      D=ATAN2(Y,X)
*
      RETURN
      END
c==================================================================
c
       subroutine bspline(tknts,t,nspl,jorder,nleft,spl)
       implicit real*8(a-h,o-z)

c calculate splines of order jorder where 1 <= jorder <= jmax
       parameter(jmax=8)
       dimension tknts(nspl+jmax)
       dimension spl(jmax)
       dimension deltal(jmax),deltar(jmax)

       if(jorder.gt.jmax)then
          write(6,*)' increase jmax in bspline'
          stop
       endif

       spl(1)=1.0

       do 200 j=1,jorder-1

       deltar(j) = tknts(nleft+j) - t
       deltal(j) = t - tknts(nleft+1-j)
       saved=0.0

       do 100 i=1,j
        term = spl(i)/(deltar(i)+deltal(j+1-i))
        spl(i) = saved + deltar(i)*term
        saved = deltal(j+1-i)*term
100    continue

       spl(j+1) = saved

200    continue

       return
       end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine bspline_deriv(tknts,t,nspl,nleft,spl, jorder, nderivs)
      implicit real*8(a-h,o-z)

c     get derivative of general order bspline by calling deboor's routines
c     this is after abandoning specific routines like bspline1
c     see AJ's beige book for derivations
      parameter(max=500)
      dimension tknts(nspl+jorder)
      dimension spl(nspl),coefs(max)
      if(nspl.gt.max)then
      print*,' increase max in bspline_deriv'
      stop
      endif

      do i=nleft-jorder+1,nleft

c zero everything except one coefficient
      do kk=1,nspl
      coefs(kk)=0.
      enddo
      coefs(i)=1

c now get the spline for that coefficient
      spl(i)=bvalue(tknts,coefs,nspl,jorder,t,nderivs)

      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function bvalue ( t, bcoef, n, k, x, jderiv )
      implicit double precision (a-h,o-z)
c  from  * a practical guide to splines *  by c. de boor
calls  interv
c
calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1
     *                     ,mflag,nmi,jdrvp1
      parameter (kmax = 20)
      real*8 bcoef(n),t(n+k),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
      bvalue = 0.
      if (jderiv .ge. k)                go to 99
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t(n+k) where it is leftcontinuous.)
      call cdb_interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do 5 j=1,i
    5    dl(j) = x - t(i+1-j)
      do 6 j=i,km1
         aj(k-j) = 0.
    6    dl(j) = dl(i)
                                        go to 10
    8 do 9 j=1,km1
    9    dl(j) = x - t(i+1-j)
c
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do 15 j=1,jcmax
   15    dr(j) = t(i+j) - x
      do 16 j=jcmax,km1
         aj(j+1) = 0.
   16    dr(j) = dr(jcmax)
                                        go to 20
   18 do 19 j=1,km1
   19    dr(j) = t(i+j) - x
c
   20 do 21 jc=jcmin,jcmax
   21    aj(jc) = bcoef(imk + jc)
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do 23 j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do 23 jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
   23       ilo = ilo - 1
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1
      do 33 j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do 33 jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
   33       ilo = ilo - 1
   39 bvalue = aj(1)
c
   99                                   return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine interv_general(tknts,time,nspl,nleft, jorder)
      implicit real*8(a-h,o-z)
      dimension tknts(nspl+1)
c     in the calling programme tknts should be set up as an array of dimension nspl+jorder

c  amendment of interv.f to treat general order splines (AJ 5/6/08)
c  new last argument jorder is introduced
c  interval starts in general at tknts(jorder) and ends at tknts(nspl+1)

c old comments appropriate to interv.f [CUBIC SPLINES, jorder=4]
c  interv3 altered from subroutine interv to treat cubic splines where
c  knots equally spaced (no stripping of repeated knots needed)
c  interval starts at tknts(4) and ends at tknts(nspl+1)

c-----
c    calculate nleft:
c                  tknts(nleft) < tknts(nleft+1)
c                  tknts(nleft) <= time <= tknts(nleft+1)
c

c  check we are in-range

      if(time.lt.tknts(jorder).or.time.gt.tknts(nspl+1)) return

      do 200 n=jorder+1,nspl+1
       if(time.le.tknts(n)) then
        nleft=n-1
        goto 210
       endif
200   continue
210   continue

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cdb_interv ( xt, lxt, x, left, mflag )
c  from  * a practical guide to splines *  by C. de Boor
c	computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      double precision x,xt(lxt)
      data ilo /1/
      save ilo
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end
*
******************************************************************
*
      SUBROUTINE PLMBAR(P,DP,Z,LMAX,INORM)
*
** Evaluates normalized associated Legendre function P(l,m) as
** function of z=cos(colatitude) using recurrence relation starting
** with P(l,l) and then increasing l keeping m fixed.  Normalization
** is: Integral(Y(l,m)*Y(l,m))=4.*pi, where Y(l,m) = P(l,m)*
** exp(i*m*longitude), which is incorporated into the recurrence
** relation. p(k) contains p(l,m) with k=(l+1)*l/2+m+1; i.e. m
** increments through range 0 to l before incrementing l. Routine
** is stable in single and double precision to l,m = 511 at least;
** timing proportional to lmax**2. R.J.O'Connell 7 Sept. 1989
*
** (1) Choice of normalisation added: If inorm = 1 Schmidt
** normalisation is chosen, where P[Schmidt](l,m) = 1/sqrt(2*l+1)*
** P[Normalised](l,m). Note the 2**0.5 missing in the normalisation.
*
** Derivatives added and stored in dp(k) using same arrangement
** as for p(k)
*
      INTEGER  LMAX, INORM, K, L, KSTART, M
*
      REAL*8  P((LMAX+1)*(LMAX+2)/2), DP((LMAX+1)*(LMAX+2)/2)
      REAL*8  PM2, PM1, PLM, PMM, F1, F2, FDEN, FNUM, Z
*
** The P(L,0) terms.
*
      PM2   = 1.0D0
      P(1)  = 1.0D0
      DP(1) = 0.0D0
      IF(LMAX.EQ.0) RETURN
      PM1   = Z
      P(2)  = PM1*3.0D0**0.5
      K     = 2
      DO 4 L = 2,LMAX
       K = K + L
       PLM = (DFLOAT(2*L-1)*Z*PM1 - DFLOAT(L-1)*PM2)/DFLOAT(L)
       P(K) = PLM*(DFLOAT(2*L+1))**0.5
       PM2 = PM1
       PM1 = PLM
4     CONTINUE
*
** The P(L,M) terms (M>0)
*
      PMM  =  1.0D0
      FNUM = -1.0D0
      FDEN =  0.0D0
      KSTART = 1
*
      DO 20 M = 1, LMAX
*
** The P(m,m) terms
*
       KSTART = KSTART + M + 1
       FNUM = FNUM + 2.0D0
       FDEN = FDEN + 2.0D0
       PMM = (1.0D0 - Z)*(1.0D0 + Z)*PMM*FNUM/FDEN
       PM2 = (DFLOAT(4*M + 2)*PMM)**0.5
       P(KSTART) = PM2
       IF(M.EQ.LMAX) GOTO 100
*
** The P(m+1,m) terms
*
       PM1=Z*(DFLOAT(2*M+3))**0.5*PM2
       K = KSTART + M + 1
       P(K) = PM1
*
** The P(l,m) terms, with l> m+1
*
       IF(M.LT.LMAX-1) THEN
        DO 10 L = M+2,LMAX
         K = K + L
         F1 = (DFLOAT((2*L+1)*(2*L-1))/DFLOAT((L+M)*(L-M)))**0.5
         F2 = (DFLOAT((2*L+1)*(L-M-1)*(L+M-1))/DFLOAT((2*L-3)
     &      * (L+M)*(L-M)))**0.5
         PLM = Z*F1*PM1 - F2*PM2
	 P(K) = PLM
	 PM2 = PM1
	 PM1 = PLM
10      CONTINUE
       ENDIF
20    CONTINUE
*
100   CONTINUE
*
** Schmitt normalise? (1=yes, 0 = no)
*
      IF(INORM.EQ.1) THEN
      K = 1
      DO 30 L = 1,LMAX
       F1 = 1.0D0/DFLOAT(2*L+1)**0.5
       DO 30 M = 0,L
        K = K + 1
        P(K) = P(K)*F1
30    CONTINUE
      ENDIF
*
** Evaluate dP/dz where z = cos(theta)
*

      DP(2) = -P(3)
      DP(3) =  P(2)
      K = 3
      DO 200 L =2,LMAX
       K = K + 1
*
** The P(l,0) terms
*
       DP(K)   = - (DFLOAT(L*(L+1))/2.D0)**0.5*P(K+1)
       DP(K+L) =   (DFLOAT(L)/2.D0)**0.5*P(K+L-1)
*
** The P(l,m) terms (m not equal to 0)
*
       DO 300 M = 1,L-1
        K = K + 1
        F1 = DFLOAT((L-M)*(L+M+1))**0.5
        F2 = DFLOAT((L+M)*(L-M+1))**0.5
        IF(M.EQ.1) F2 = F2*2.0D0**0.5
        DP(K) = (F2*P(K-1) - F1*P(K+1))/2.0D0
300    CONTINUE
       K = K + 1
*
200   CONTINUE
*
      RETURN
      END
