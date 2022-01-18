
 	subroutine polarisation
c
 	implicit none
 	integer nbrrot,nbrang,nang
 	parameter(nbrrot = 90000)
 	parameter(nbrang = 80,nang = 80)

c 	Entrees
c 	-------
!!!     La calibration en angle a été prise en compte dans lect.
 	real pi,calib_angle,calib_dark,trans_filtre
 	integer nmoy,nrot,entergreg,dess_polar
 	real dec,chi(nbrrot),chideg(nbrrot)
 	character inputfile*70,labelplot*50,entete*80
        character*50 dirin
        real an,numjour,alat,along,elevation,azimut
	real angdeg
 	data dess_polar/1/
c
c 	Interne
c 	-------
 	real bid,V_inst,angrad
 	integer iang,irot,jrot,ibid,lchaine
 	integer nang_reel
c
c 	Sur une rotation
c 	----------------
 	real hrloc(nbrrot),UT(nbrrot),V_rot(nbrrot),Vcos_rot(nbrrot),
     .       Vsin_rot(nbrrot),DoLp_rot(nbrrot),AoLP_rot(nbrrot)
c	Moyenne sur nmoy rotations
c 	--------------------------
 	real V_glis(nbrrot),Vcos_glis(nbrrot),
     .       Vsin_glis(nbrrot),DoLp_glis(nbrrot),AoLP_glis(nbrrot)
c 	Moyenne sur toute la manip
c 	--------------------------
 	real hrloc_tot(2),V_tot(2),Vcos_tot(2),
     .       Vsin_tot(2), DoLp_tot(2),AoLP_tot(2)
c
 	real tabmoy
 	real tmin,tmax,ymin,ymax,yymin,yymax
c
 	pi = 4.*atan(1.)

 	open(10,file='polar.input')
        read(10,*)dirin
 	call xline(2,10)
 	read(10,*)calib_dark
 	trans_filtre = 0.95	! a mettre en lecture plus tard
 	call xline(1,10)
 	read(10,*)nmoy
 	call xline(14,10)
 	read(10,*)calib_angle
 	read(10,*)an,numjour
 	read(10,*)alat,along
 	read(10,*)elevation,azimut
 	read(10,1000)labelplot
 	read(10,*)entergreg
 	close(10)
c
        open(10,file=
     .   '../'//dirin(1:lchaine(dirin))//'/lect_spp_modsynch.output')
 	call xline(1,10)
 	read(10,*)nrot
 	call xline(3,10)
	do irot = 1,nrot
 	  read(10,*) jrot,hrloc(irot),UT(irot),chideg(irot),
     .               bid,bid,bid,bid,bid,V_rot(irot)
 	enddo
c
	do irot = 1,nrot
 	  call xline(3,10)
 	  Vsin_rot(irot) = 0.
 	  Vcos_rot(irot) = 0.
 	  nang_reel=1
 	  do iang = 1,nang
 	    read(10,*)ibid,bid,angdeg,bid,bid,bid,bid,V_inst
 	    if(V_inst.ne.0.)then
 	      angrad = 2.*angdeg*pi/180.
 	      Vsin_rot(irot)=Vsin_rot(irot)+V_inst*sin(angrad)
 	      Vcos_rot(irot)=Vcos_rot(irot)+V_inst*cos(angrad)
 	      nang_reel = nang_reel+1
 	    endif
 	  enddo
 	  Vsin_rot(irot)=Vsin_rot(irot)/float(nang_reel)
 	  Vcos_rot(irot)=Vcos_rot(irot)/float(nang_reel)
 	enddo
c
c 	Calcul de AoLP et DoLP à chaque rotation
c 	----------------------------------------
 	do irot = 1,nrot
 	  V_rot(irot) = max(0.,V_rot(irot)-calib_dark)
          DoLP_rot(irot) = 0.
          AoLP_rot(irot) = 0.
 	  if (V_rot(irot).ne.0.)then
            DoLP_rot(irot)=2.*sqrt(Vcos_rot(irot)**2+Vsin_rot(irot)**2)
     .                    /V_rot(irot)
            DoLP_rot(irot) = 100.*DoLP_rot(irot)/trans_filtre
            AoLP_rot(irot)=atan2(Vsin_rot(irot),Vcos_rot(irot))/2.
            AoLP_rot(irot)=180./pi* AoLP_rot(irot)
 	  endif
 	enddo
c
c 	Moyenne sur toute la manip
c	--------------------------
 	hrloc_tot(1) = hrloc(1)
 	hrloc_tot(2) = hrloc(nrot)
c
 	V_tot(1)= max(0.,tabmoy(V_rot,nrot))
 	V_tot(2) = V_tot(1)
 	if (V_tot(1).ne.0.)then
 	  Vcos_tot(1)= tabmoy(Vcos_rot,nrot)
 	  Vcos_tot(2)= Vcos_tot(1)
 	  Vsin_tot(1)= tabmoy(Vsin_rot,nrot)
 	  Vsin_tot(2)= Vsin_tot(1)
c
 	  DoLP_tot(1)=2.*sqrt(Vcos_tot(1)**2+Vsin_tot(1)**2)/V_tot(1)
          DoLP_tot(1) = 100.*DoLP_tot(1)/trans_filtre
 	  DoLP_tot(2)= DoLP_tot(1)
 	  AoLP_tot(1)=180./pi* atan2(Vsin_tot(1),Vcos_tot(1))/2.
 	  AoLP_tot(2)= AoLP_tot(1)
 	else
 	  DoLP_tot(1)=0.
 	  DoLP_tot(2)=0.
 	  AoLP_tot(1)=0.
 	  AoLP_tot(2)=0.
 	endif
c
c 	Moyenne sur nmoy points
c	-----------------------
        call moyglis(V_rot,nrot,nmoy,V_glis)
c
c 	Retire les éclairs sur une moyenne de 20 rotations
        call moyglis(V_rot,nrot,20,V_glis)
 	do irot=2,nrot
  	  if (V_rot(irot).gt.V_glis(irot)+0.5)then
  	    V_rot(irot)=V_rot(irot-1)
  	    Vcos_rot(irot)=Vcos_rot(irot-1)
  	    Vsin_rot(irot)=Vsin_rot(irot-1)
  	  endif
  	enddo

c 	Retourne à la moyenne demandée
        call moyglis(V_rot,nrot,nmoy,V_glis)
        call moyglis(Vcos_rot,nrot,nmoy,Vcos_glis)
        call moyglis(Vsin_rot,nrot,nmoy,Vsin_glis)
c
 	do irot=1,nrot
 	  DoLP_glis(irot)=2.*sqrt(Vcos_glis(irot)**2+Vsin_glis(irot)**2)
     .                    /V_glis(irot)
 	  DoLP_glis(irot)=100.* DoLP_glis(irot)/trans_filtre
 	  AoLP_glis(irot)=atan2(Vsin_glis(irot),Vcos_glis(irot))/2.
 	  AoLP_glis(irot)=180./pi* AoLP_glis(irot)
c 	  On va essayer de supprimer les sauts de phase...
c	  if(AoLP_glis(irot).lt.AoLP_tot(1)-100.)
c    .       AoLP_glis(irot)=AoLP_glis(irot)+180.
 	enddo
c
c       Prints
c  	------
c
 	write(6,*)'DoLP et AoLP moyens',DoLP_tot(1),AoLP_tot(1)
        open(20,file='../'//dirin(1:lchaine(dirin))//
     .       '/polar_spp_modsynch.output')
        write(20,1000)labelplot
 	write(20,*)nrot,'  Number of rotations'
 	write(20,1030)
 	do irot = 1,nrot
          write(20,1040)irot, hrloc(irot),UT(irot),V_rot(irot),
     .       DoLP_rot(irot),DoLP_glis(irot),DoLP_tot(1),
     .       AoLP_rot(irot),AoLP_glis(irot),AoLP_tot(1),chideg(irot)
 	enddo
c
c       Plots
c  	-----
        call GR_EXEC('set /default')
c       call GR_EXEC('SET PLOT_PAGE PORTRAIT')
        call GR_EXEC('clear ')
        call GR_EXEC('set char .8')
c
 	tmin = hrloc(1)
 	tmax = hrloc(nrot)
c
c 	Intensity (mV)
c 	--------------
 	if (dess_polar.eq.1)then
        call mnmx(V_rot,nrot,ymin,ymax,0)
 	ymax = ymax*1.02
        call GR_EXEC('set view .2 .9 .65 .9 ')
        call GR_EXEC('pen 0 /weight 3')
        call gr_limi(4,tmin,tmax,ymin,ymax)
        call GR_EXEC('box N')
        call GR_EXEC('label "Intensity" /Y')
        call GR_EXEC('set view .24 .9 .65 .9 ')
        call GR_EXEC('label "[mV]" /Y')
        call GR_EXEC('set view .2 .9 .65 .9 ')
c
        call GR_EXEC('pen 0 /weight 3')
        call gr4_give('X',2,hrloc_tot)
        call gr4_give('Y',2,V_tot)
c       call GR_EXEC('connect')
        call gr4_give('X',nrot,hrloc)
        call gr4_give('Y',nrot,V_rot)
        call GR_EXEC('connect')
        call GR_EXEC('pen 1 /weight 3')
        call gr4_give('X',nrot,hrloc)
        call gr4_give('Y',nrot,V_glis)
c       call GR_EXEC('connect')
c
c       DoLP (%)
c       --------
        call mnmx(DoLP_glis,nrot,ymin,ymax,0)
c       call mnmx(DoLP_rot,nrot,yymin,yymax,0)
c	if(ymax.lt.yymax*1.5)ymax = yymax*1.5
 	ymax = 8.
        call GR_EXEC('set view .2 .9 .40 .65')
        call GR_EXEC('pen 0 /weight 3')
        call gr_limi(4,tmin,tmax,ymin,ymax)
        call GR_EXEC('box N')
        call GR_EXEC('label "DoLP" /Y')
        call GR_EXEC('set view .24 .9 .40 .65')
        call GR_EXEC('label "[%]" /Y')
        call GR_EXEC('set view .2 .9 .40 .65')
c
        call GR_EXEC('pen 0 /weight 3 /dash 2')
        call gr4_give('X',2,hrloc_tot)
        call gr4_give('Y',2,DoLP_tot)
        call GR_EXEC('connect')
        call GR_EXEC('pen 1 /weight 3')
        call gr4_give('X',nrot,hrloc)
        call gr4_give('Y',nrot,DoLP_rot)
c       call GR_EXEC('connect')
        call GR_EXEC('pen 0 /weight 3 /dash 1')
        call gr4_give('X',nrot,hrloc)
        call gr4_give('Y',nrot,DoLP_glis)
        call GR_EXEC('connect')
c
c       AoLP (%)
c       --------
        call mnmx(AoLP_glis,nrot,ymin,ymax,0)
        call mnmx(AoLP_rot,nrot,yymin,yymax,0)
 	ymin = min(ymin,yymin)
 	ymax = max(ymax,yymax)
 	ymin = -20.
        call GR_EXEC('set view .2 .9 .15 .40')
        call GR_EXEC('pen 0 /weight 3')
        call gr_limi(4,tmin,tmax,ymin,ymax)
        call GR_EXEC('box')
        call GR_EXEC('label "AoLP" /Y')
        call GR_EXEC('set view .24 .9 .15 .40')
        call GR_EXEC('label "[deg]" /Y')
        call GR_EXEC('set view .2 .9 .15 .40')
c
        call GR_EXEC('pen 0 /weight 3')
        call gr4_give('X',2,hrloc_tot)
        call gr4_give('Y',2,AoLP_tot)
        call GR_EXEC('connect')
        call GR_EXEC('pen 1 /weight 3 /dashe 1')
        call gr4_give('X',nrot,hrloc)
        call gr4_give('Y',nrot,AoLP_rot)
c       call GR_EXEC('connect')
        call GR_EXEC('pen 0 /weight 3 /dashe 1')
        call gr4_give('X',nrot,hrloc)
        call gr4_give('Y',nrot,AoLP_glis)
        call GR_EXEC('connect')

        call GR_EXEC('pen 0 /weight 3')
        call GR_EXEC('set view .2 .9 .15 .90')
        call GR_EXEC('label "Local Time" /X')
        call GR_EXEC('DRAW TEXT .0 0.1 "'//labelplot//'" 8 0. /BOX 8')
        if(entergreg.eq.1)call gmaster_enter_loop('')
        write(entete,1053)dirin
        call GR_EXEC(entete)
1053    format
     .   ('hardcopy "../',a12,'/polar_spp_modsynch" /dev eps col /over')
 	endif

1000    format(a50)
1030 	format(/,14x,'Loc Time',7x,'UT',6x,'Intensity',4x,
     .        'DoLP_rot',3x,'DoLP_glis',
     .    3x,'DoLP tot',5x,'AoLP_rot    AoLP_glis   AoLP_tot',
     .    6x,'chideg')
1040 	format(i10,2f12.6,1f11.3,3x,3(f10.7,2x),3(f10.5,2x),f10.2)
c
 	stop
 	end
c
c --------
c
 	subroutine zeroit(tab,ntab)
c
  	implicit none
c
 	integer ntab,itab
 	real tab(ntab)
c
 	do itab = 1,ntab
 	  tab(itab) = 0.
 	enddo
 	return
 	end

c
c-------------------------- mnmxplt -----------------------------------
c
        subroutine mnmxplt(tmin,tmax,linlog)
c
c       find best min and max for nice plot.
c       linlog = 0 if linear axis
c              = 1 if logarithmic axis
        implicit none
        integer linlog
        real tmin, tmax,ttmin,ttmax,ttminint,ttmaxint,deltmin,deltmax,
     .       delta
c
        if(linlog.eq.1)then
c         Quand on est en axe log et qu'il y a moins d'une decade entre
c         le min et le max, GREG ne marque pas les unites sur l'axe >
c         On se premuni de cela.
          if(tmin.ne.0.)then
            ttmin = alog10(tmin)
            ttmax = alog10(tmax)
            delta = ttmax-ttmin
            if(delta .lt. 1.)then
              ttminint = float(ifix(ttmin))
              ttmaxint = float(ifix(ttmax)+1)
c             Quel est le plus proche d'un tick?
              deltmin = ttmin-ttminint
              deltmax = ttmaxint-ttmax
              if(deltmin.le.deltmax)then
c               C'est tmin!
                tmin = 10**ttminint
                tmax = tmax*1.2
              else
c               C'est tmax!
                tmax = 10**ttmaxint
                tmin = tmin/1.2
              endif
              return
            endif
          endif
c
          tmin=tmin/1.2
          if(tmin.eq.0.)tmin=1.e-05
          tmax=tmax*1.2
        else
          if (tmin.lt.200.) then
            tmin = float(ifix(tmin)/10 - 1) * 10.
          elseif (tmin.lt.1000.) then
            tmin = float(ifix(tmin)/10 - 5) * 10.
          else
            tmin = float(ifix(tmin)/10 -10) * 10.
          endif
          if (tmax.lt.200.) then
            tmax = float(ifix(tmax)/10 + 2) * 10.
          elseif (tmax.lt.1000.) then
            tmax = float(ifix(tmax)/10 + 5) * 10.
          else
            tmax = float(ifix(tmax)/10 +10) * 10.
          endif
        endif
c
        return
        end

c
c -------
c
 	function tabmoy(tab,ntab)
c
 	implicit none
 	integer itab,ntab
 	real tabmoy, tab(ntab)
c
 	tabmoy = 0.
 	do itab = 1,ntab
 	  tabmoy = tabmoy + tab(itab)
 	enddo
 	tabmoy = tabmoy/float(ntab)
 	return
 	end
