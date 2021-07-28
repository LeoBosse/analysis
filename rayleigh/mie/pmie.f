c *******************************************************************
c
c              calculs de mie  sur 387 mugauss  (C. Devaux, LOA)
c              -------------------------------------------------
c
c     ici:  granulometrie log-normale :
c  rn = rn0*exp(-(log(r/rm))**2/(2*(log(sig))**2))/r/log(sig)/(2*pi)**0.5 pour r<rmax
c         (voir dans la routine pp2).
c
c     entrees sur fichier lu en read(5) :
c      -  lambda, rn, in  <0 ;
c                        ||||
c      -  les parametres de la granulometrie (rn0, rm, rlnsig).
c      -  les rayons limitant la granulometrie rmin, rmax.
c
c     parametres figurant dans le programme mais pouvant etre changes :
c      -  le alpha initial (0.05) par ex
c      -  les pas en alpha.
c
c compile with f77 -o pmie pmie.f (Leo)
c Fichier d'output: i=fct de phase, q=fct de phase polarisee, u=?, -q/i = DoLP
c **********************************************************************
c
      program pri
      parameter (ialfmx=10000,ialfmx2=2*ialfmx,
     &           nmumax=193)
c
      implicit double precision (a-h,o-z)
      double precision in,ia,ib,k1,k2,k3,int,q,u,igna,idnb
      double precision kmat1,kmat2,kmat3,np,it,i1
      dimension rmu(-nmumax:nmumax),chr(-nmumax:nmumax)
      dimension int(-nmumax:nmumax),q(-nmumax:nmumax),
     &  u(-nmumax:nmumax)
      dimension ra(0:ialfmx2),rb(0:ialfmx2),ia(0:ialfmx2),
     &  ib(0:ialfmx2),cna(-1:ialfmx2),rgna(-1:ialfmx2),
     &  igna(-1:ialfmx2),rdna(0:ialfmx2), rdnb(0:ialfmx2),
     &  idnb(0:ialfmx2),sna(-1:ialfmx2)
c
      dimension i1(-nmumax:nmumax),q1(-nmumax:nmumax),
     &  it(-nmumax:nmumax),
     &  u1(-nmumax:nmumax),qt(-nmumax:nmumax),ut(-nmumax:nmumax),
     &  pl(-1:2*nmumax+10),pol(0:2*nmumax+10),beta(0:2*nmumax+10),
     &  gamma(0:2*nmumax+10),
     &  alp(0:2*nmumax+10),zeta(0:2*nmumax+10),delta(0:2*nmumax+10)
c
c     parametres de la granulometrie (gam.std.) :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     common /gr/ rn0,aa,bb,rmin,rmax
      common /gr/ rn0,rmg,rlnsig,rmin,rmax
c
      common/a/ra,ia,rb,ib,alpha,k1,k2,k3,n2
      common /b/ cna,rgna,igna,rdna,rdnb,idnb,sna
      common /c/ int, q, u
      common /g/ nbmu,rmu,chr
      common /l/ rlam,rn,in
      common /cpi/pi,cons
      common /p/pas
      common /d/ pl,pol,beta,gamma,alp,zeta,delta,it,
     &  qt,ut
      common /r/ kmat1,kmat2,kmat3,i1,q1,u1,np
c
      pi=dacos(-1.d+00)
      cons=dsqrt(2*pi)
c
      open(4,file='mie/mu387')
      nbmu=nmumax
      read(4,222)(rmu(j),chr(j),j=nbmu,-nbmu,-1)
      lll=2*nbmu+1
C
      read(5,*)rlam,rn,in
      read(5,*)rn0,rmg,rlnsig,rmin,rmax           ! param LND : rmg (r median) et ln(sigmag)
c
      alphao=0.001
      alphaf = 2*pi*rmax/rlam
c     alphaf augmente de 10 (petite securite...)
      alphaf=dint(alphaf+10.)
      ralfmx = ialfmx*1.0
      if(alphaf.gt.ralfmx)then
	  write(6,*)'alpha max trop grand'
	  stop
      endif
      write(6,*)'alpha final = ',alphaf
      write(6,*)
      call zero(ra,ialfmx2+1)
      call zero(rb,ialfmx2+1)
      call zero(ia,ialfmx2+1)
      call zero(ib,ialfmx2+1)
      call zero(cna,ialfmx+2)
      call zero(rgna,ialfmx+2)
      call zero(igna,ialfmx+2)
      call zero(rdnb,ialfmx2+1)
      call zero(sna,ialfmx2+1)
      call zero(idnb,ialfmx2+1)
      call zero(rdna,ialfmx2+1)
c
      call zero(beta,2*nmumax+11)
      call zero(gamma,2*nmumax+11)
      call zero(alp,2*nmumax+11)
      call zero(zeta,2*nmumax+11)
      call zero(delta,2*nmumax+11)
      call zero(i1,2*nmumax+1)
      call zero(q1,2*nmumax+1)
      call zero(u1,2*nmumax+1)
      call zero(it,2*nmumax+1)
      call zero(qt,2*nmumax+1)
      call zero(ut,2*nmumax+1)
c
      np=0.0
      kmat1=0.0
      kmat2=0.0
      kmat3=0.0
c
      if(alphao.gt.alphaf)goto301
      alpha=alphao
      pas=0.001
 7777 alpha1=alpha
c
c-----------------------------
c     definition des pas :
c-----------------------------
      if(alpha.ge.0.1)pas=0.01
      if(alpha.ge.10)pas=0.05
      if(alpha.gt.100)pas=0.1
      if(alpha.gt.500)pas=1.0
      if(alpha.gt.1000)pas=5.0
c-----------------------------
c
c     write(6,*)'alpha :',alpha
c
      call pp1
      call pp2
c
      alpha=alpha+pas
      if (alpha.le.alphaf+pas/100.0) goto 7777
  301 continue
c
      do 55 j=-nbmu,nbmu
      i1(j)=i1(j)/kmat2
      q1(j)=q1(j)/kmat2
      u1(j)=u1(j)/kmat2
  55  continue
c
c     facteur d'assymetrie non tronque, sections ext/dif, pizero :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sext=kmat1/np
      scatt=kmat2/np
      gasym=kmat3/kmat2
      piz=kmat2/kmat1
c
      write(6,100) rn0,rmg,rlnsig,rmax
      write(6,156)rlam
      write(6,158)rn,in
      write(6,22)np
c
      write(6,21)sext,scatt,piz,gasym
      write(6,155)
      write(6,154)
      write(6,155)
      do 444 j=-nbmu,nbmu
      write(6,157)j,dacos(rmu(j))*180/pi,rmu(j),i1(j),q1(j),u1(j),
     1-q1(j)/i1(j)*100.
  444 continue
      write(6,*)
c
  888 continue
      stop
c
c -------------------------------------------------------
c
  222 format(2d21.14)
  203 format(5(2x,e14.7))
  204 format(8(2x,f8.5))
  100 format(5x,'granulo :  rn = ',f4.1,'  rm = ',f4.2,'  ln sig = ',
     & f6.3,'  pour r < ',f4.0,' mic.',/,5x,70(1h=),/)
  154 format(2x,' j',6x,'theta',7x,'mu',9x,'i',13x,'q',13x,'u',
     s9x,'-q/i*100')
  155 format(1x,78(1h-))
  156 format(20x,'longueur d onde :',f7.3,' microns'/,20x,32(1h=),/)
  157 format(1x,i4,3x,f7.2,3x,f8.5,1x,3(e12.5,2x),2x,f7.2)
  158 format(5x,'indice reel :',f10.5,5x,'indice imaginaire : ',
     1e12.5,/,/)
  21  format(7x,'Sections d extinction et de diffusion ramenees a',
     & ' une particule :',/,7x,64(1h-),/,5x,'extinction : ',
     & e12.5,' mic**2',5x,'diffusion : ',e12.5,' mic**2',/,/,
     & 12x,'Pizero : ',f8.6,17x,'g.assymetrie :',f8.5,/)
  22  format(20x,'nombre de particules :',e13.5,/)
      end
c
c*********************************************************************
c
      subroutine pp1
c     ~~~~~~~~~~~~~~
c
c     calcul des series de Mie, facteurs de diffusion et d'extinction,
c     et facteur d'assymetrie pour un parametre de Mie alpha donne :
c
c     ----------------------------------------------------------------
c
      parameter (ialfmx=10000,ialfmx2=2*ialfmx,
     &           nmumax=193)
      implicit double precision (a-h,o-z)
      integer un,r1,test
      double precision k1,k2,k3,ia,ib,in,igna,idnb,ibeta
c
      dimension cna(-1:ialfmx2),rgna(-1:ialfmx2),igna(-1:ialfmx2),
     srdna(0:ialfmx2),ra(0:ialfmx2),rb(0:ialfmx2),ia(0:ialfmx2),
     &  ib(0:ialfmx2),rdnb(0:ialfmx2),idnb(0:ialfmx2),sna(-1:ialfmx2),
     s  rmu(-nmumax:nmumax),chr(-nmumax:nmumax)
      common/a/ra,ia,rb,ib,alpha,k1,k2,k3,n2
      common /b/ cna,rgna,igna,rdna,rdnb,idnb,sna
      common /g/ nbmu,rmu,chr
      common /l/ rlam,rn,in
      common /cpi/pi,cons
c
      r1=20
      rbeta=rn*alpha
      ibeta=in*alpha
      n1=dint(alpha+alpha+r1)
      n2=dint(alpha+alpha+5)
      n2p1=n2+1
      rsurl=alpha/2.0/pi+1.0d-10
      irl=dint(rsurl)
      frc=rsurl-irl*1.0d00
      if(frc.lt.1.0d-08) alpha=alpha*(1.+1.0d-08)
      cna(-1)=-dsin(alpha)
      rgna(-1)=0
      rgna(0)=0
      cna(0)=dcos(alpha)
      igna(-1)=0
      igna(0)=-1
      do 25 i=1,n2
      cna(i)=(2*i-1)*cna(i-1)/alpha-cna(i-2)
      x=rgna(i-1)
      z=i/alpha
      y=igna(i-1)
      w=((z-x)*(z-x)+(y*y))
      rgna(i)=(z-x)/w-z
      igna(i)=y/w
      if(cna(i).lt.1d+100) goto 25
      n2=i
      n2p1=i+1
      n1=i+15
      goto 26
   25 continue
  26  continue
      rdna(n1)=0
      rdnb(n1)=0
      idnb(n1)=0
      x1=rbeta*rbeta+ibeta*ibeta
      x2=rbeta/x1
      x3=ibeta/x1
      sna(n1)=0
      sna(n1-1)=1
      do 30 kk=1,n1
      i=n1-kk
      x=rdnb(i+1)
      y=idnb(i+1)
      z=x+(i+1)*x2
      w=y-(i+1)*x3
      x4=z*z+w*w
      rdnb(i)=(i+1)*x2-z/x4
      idnb(i)=-(i+1)*x3+w/x4
      z=(i+1)/alpha
      x=rdna(i+1)
      rdna(i)=z-1/(x+z)
      sna(i-1)=(2*i+1)*sna(i)/alpha-sna(i+1)
      if (sna(i-1).le.1d+60) goto 40
      test=i-1
      x=sna(test)
      do 35 j=test,n2
   35 sna(j)=sna(j)/x
   40 continue
   30 continue
      q=-sna(0)/cna(-1)
      do 45 i1=1,n2p1
   45 sna(i1-1)=sna(i1-1)/q
      test=0
      un=1
      do 50 i=1,n2
      x1=sna(i)
      x2=cna(i)
      x3=rdnb(i)
      x4=idnb(i)
      x5=rdna(i)
      x6=rgna(i)
      x7=igna(i)
      y1=x3-rn*x5
      y2=x4-in*x5
      y3=x3-rn*x6+in*x7
      y4=x4-rn*x7-in*x6
      y5=rn*x3-in*x4-x5
      y6=in*x3+rn*x4
      y7=rn*x3-in*x4-x6
      y8=in*x3+rn*x4-x7
      x4=y2*y3-y1*y4
      x3=y1*y3+y2*y4
      x5=x1*x1+x2*x2
      x6=y3*y3+y4*y4
      x7=y5*y7+y6*y8
      x8=y6*y7-y5*y8
      x9=y7*y7+y8*y8
      q=(i+i+1.)/i/(i+1.)*un
      y1=x1*(x1*x3+x2*x4)/x5/x6
      y2=x1*(x1*x4-x2*x3)/x5/x6
      y3=x1*(x1*x7+x2*x8)/x5/x9
      y4=x1*(x1*x8-x2*x7)/x5/x9
      ra(i)=y2*q
      ib(i)=y3*q
      q=-q
      rb(i)=y4*q
      ia(i)=y1*q
      un=-un
   50 continue
      ra(0)=0
      ia(0)=0
      rb(0)=0
      ib(0)=0
      ra(n2p1)=0
      ia(n2p1)=0
      rb(n2p1)=0
      ib(n2p1)=0
c
      k1=0.0
      k2=0.0
      k3=0.0
      j=-1
      do 55 n=1,n2
      l1=n*n
      a2=(n+1.)*(n+1.)
      x=ra(n)
      y=ia(n)
      z=rb(n)
      t=ib(n)
      xp=ra(n+1)
      yp=ia(n+1)
      zp=rb(n+1)
      tp=ib(n+1)
      k1=k1+n*(n+1)*j*(y-t)
      k2=k2+l1*a2/(n+n+1.)*(x*x+y*y+z*z+t*t)
      y10=l1*(n+2.)*(n+2.)*(n+1.)/(2*n+1.)/(2*n+3.)
      k3=k3-y10*(x*xp+y*yp+z*zp+t*tp)
      k3=k3-n*(n+1.)/(2*n+1.)*(x*z+y*t)
      j=-j
   55 continue
      w6=2/alpha/alpha
c
c     k1 = Qext
c     k2 = Qscat
c     k3 donne g
c
      k1=w6*k1
      k2=w6*k2
      k3=k3*4./alpha/alpha/k2
      call calcul
 9999 return
      end
c
c******************************************************
c
      subroutine zero(it,n)
c
      parameter (ialfmx=10000,ialfmx2=2*ialfmx,
     &           nmumax=193)
      double precision it(n)
      do 1 i=1,n
    1 it(i)=0
      return
      end
c
c*******************************************************
c
      subroutine calcul
c     ~~~~~~~~~~~~~~~~~
c
c     calcul de i, q, u pour un parametre de Mie alpha :
c
c     --------------------------------------------------
c
      parameter (ialfmx=10000,ialfmx2=2*ialfmx,
     &           nmumax=193)
      implicit double precision (a-h,o-z)
c     implicit integer*4 (i-n)
      double precision kma1,kma2,kma3,ia,ib,in,ims1,ims2
      double precision int,q,u
c     real int,ay1,ay2,ay3,ay4,q,u
      dimension rmu(-nmumax:nmumax),chr(-nmumax:nmumax),
     &  ra(0:ialfmx2),rb(0:ialfmx2),ib(0:ialfmx2), ia(0:ialfmx2),
     &  int(-nmumax:nmumax),q(-nmumax:nmumax),u(-nmumax:nmumax)
      common/a/ra,ia,rb,ib,alpha,kma1,kma2,kma3,n2
      common /c/ int, q, u
      common /g/ nbmu,rmu,chr
      common /l/ rlam,rn,in
      common /cpi/p,cons
      ay1=alpha
      ay2=kma1
      ay3=kma2
      ay4=kma3
      coef=2./kma2/alpha**2
      do 1 j=-nbmu,nbmu
      x=-rmu(j)
      pim=0.
      pi=1.
      tau=x
      res1=0
      res2=0
      ims1=0
      ims2=0
      do 2 n=1,n2
      ai=ia(n)
      bi=ib(n)
      ar=ra(n)
      br=rb(n)
      res1=res1-ai*pi-bi*tau
      res2=res2+ai*tau+bi*pi
      ims1=ims1+ar*pi+br*tau
      ims2=ims2-ar*tau-br*pi
      pip=((2*n+1)*x*pi-(n+1)*pim)/n
      pim=pi
      pi=pip
      tau=(n+1)*x*pi-(n+2)*pim
    2 continue
      y1=res1**2+ims1**2
      y2=res2**2+ims2**2
      y3=2*res2*res1
      y4=2*ims2*ims1
      int(j)=coef*(y1+y2)
      q(j)=coef*(y2-y1)
      u(j)=coef*(y3+y4)
    1 continue
c     write(1)ay1,ay2,ay3,ay4,int
  900 format(5d15.8)
c     write(6,1234)alpha,kma1,kma2,kma3,int(-nbmu),
c    1 int(nbmu)
 1234 format(1x,f7.2,5x,3e14.7,3x,2e13.6)
      return
      end
c
c*******************************************************************
c
      subroutine pp2
c     ~~~~~~~~~~~~~`
c
c     integrations sur la granulometrie
c
c     --------------------------------------------------------------
c
      parameter (ialfmx=10000,ialfmx2=2*ialfmx,
     &           nmumax=193)
c
      implicit double precision (a-h,o-z)
      double precision in,ia,ib,kma1,kma2,kma3,i2,q2,u2,igna,idnb
      double precision kmat1,kmat2,kmat3,np,nr,it,i1
      dimension rmu(-nmumax:nmumax),chr(-nmumax:nmumax),
     & i2(-nmumax:nmumax),q2(-nmumax:nmumax),u2(-nmumax:nmumax)
      dimension ra(0:ialfmx2),rb(0:ialfmx2),ia(0:ialfmx2),
     &  ib(0:ialfmx2),cna(-1:ialfmx2),rgna(-1:ialfmx2),
     &  igna(-1:ialfmx2),rdna(0:ialfmx2), rdnb(0:ialfmx2),
     &  idnb(0:ialfmx2),sna(-1:ialfmx2)
c
      dimension i1(-nmumax:nmumax),q1(-nmumax:nmumax),
     &  pl(-1:2*nmumax+10),pol(0:2*nmumax+10),beta(0:2*nmumax+10),
     &  gamma(0:2*nmumax+10),it(-nmumax:nmumax),
     &  alp(0:2*nmumax+10),zeta(0:2*nmumax+10),delta(0:2*nmumax+10),
     &  u1(-nmumax:nmumax),qt(-nmumax:nmumax),ut(-nmumax:nmumax)
c
      common/a/ra,ia,rb,ib,alpha,kma1,kma2,kma3,n2
      common /c/ i2, q2, u2
      common /g/ nbmu,rmu,chr
      common /l/ wa,rn,in
      common /cpi/pi,cons
C     common /gr/ rn0,aa,bb,rmin,rmax
      common /gr/ rn0,rmg,rlnsig,rmin,rmax
      common /d/ pl,pol,beta,gamma,alp,zeta,delta,it,
     &  qt,ut
      common /r/ kmat1,kmat2,kmat3,i1,q1,u1,np
      common /p/ pas
c
      r=alpha*wa/2/pi
      nr=0.0
c     write(6,1234)alpha,kma1,kma2,kma3,i2(-nbmu),
c    1i2(nbmu)
 1234 format(1x,f7.2,5x,3e14.7,3x,2e13.6)
c
c---------------------------------------
c     la granulometrie est definie ici :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(r.lt.rmin) then
	  nr=0.0
	  goto 121
	  endif
      if(r.gt.rmax) then
	  nr=0.0
	  goto 121
	  endif
C     nr=rn0*r**aa*exp(-bb*r)           ! gamma
      XX=LOG(R/RMG)
      XX=XX*XX/(2.*RLNSIG*RLNSIG)
      YY=DSQRT(2.*PI)*RLNSIG
      nr=exp(-XX)/r                     ! LND en log neperien
      nr=nr/YY
  121 continue
c---------------------------------------
c
c     write(6,151)alpha,r,nr
  151 format(e10.4,2x,e10.4,3x,i3,3x,e10.4)
c
c
      pr=wa*pas/2/pi
      x1=nr*pr*pi*r**2
      kmat1=kmat1+x1*kma1
      kmat2=kmat2+kma2*x1
      kmat3=kmat3+x1*kma2*kma3
      np=np+nr*pr
      nr=nr*pr*wa**2/4/pi
      x1=kma2*alpha*alpha*nr
      do 1 j=-nbmu,nbmu
      i1(j)=i1(j)+i2(j)*x1
      q1(j)=q1(j)+q2(j)*x1
      u1(j)=u1(j)+u2(j)*x1
  1   continue
c     write(6,1234)alpha,kmat1,kmat2,kmat3,i1(-nbmu),i1(nbmu)
      return
      end
