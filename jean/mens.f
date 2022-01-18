c
 	subroutine mens
c
 	implicit none
 	integer nbrrot
 	parameter (nbrrot = 9000)
 	character inputfile*70,labelplot*50,entete*50
 	real time_tot(2),Vavg_tot(2),Vcos_tot(2),
     .       Vsin_tot(2), DoLp_tot(2),AoLP_tot(2)
 	real t_inst(nbrrot),V_inst(nbrrot),Vcos_inst(nbrrot),
     .       Vsin_inst(nbrrot),DoLp_inst(nbrrot),AoLP_inst(nbrrot)
 	real V_inst_corr(nbrrot),time_inst(nbrrot)
 	real V_glis(nbrrot),Vcos_glis(nbrrot),
     .       Vsin_glis(nbrrot),DoLp_glis(nbrrot),AoLP_glis(nbrrot)
 	real DoLp_moy(nbrrot),AoLP_moy(nbrrot)
 	real pi,calib_angle,calib_dark
 	real tabmoy,bidV,bidVcos,bidVsin
 	integer nmoy,irot,nrot,entergreg
 	real tmin,tmax,ymin,ymax,yymin,yymax
 	real hdeb,mndeb,tdeb
c
1000 	format(a)
 	data calib_angle/58.1/,calib_dark/0.14/
c
 	pi = 4.*atan(1.)
	
 	open(10,file='mens.input')
 	read(10,1000)inputfile
 	write(6,1000)inputfile
 	read(10,*)nmoy
 	read(10,1000)labelplot
 	read(10,*)hdeb,mndeb
 	read(10,*)entergreg
 	close(10)
 	tdeb = hdeb+mndeb/60.
c
 	open(10,file=inputfile)
 	open(20,file='mens.dat')
 	call xline(7,10)
 	read(10,*)nrot
 	write(20,*)nrot,'  Number of rotations'
 	call xline(1,10)
	do irot=1,nrot
 	  read(10,*) t_inst(irot),V_inst(irot),Vcos_inst(irot),
     .       Vsin_inst(irot),DoLp_inst(irot),AoLP_inst(irot)
 	  V_inst_corr(irot)=V_inst(irot)-calib_dark
 	  time_inst(irot)=tdeb+t_inst(irot)/3600.
 	enddo
c
c 	Calcul des paramètres sur toute la manip
 	
 	time_tot(1) = time_inst(1)
 	time_tot(2) = time_inst(nrot)
c
 	bidV= tabmoy(V_inst_corr,nrot)
 	bidVcos= tabmoy(Vcos_inst,nrot)
 	bidVsin= tabmoy(Vsin_inst,nrot)
c
 	Vavg_tot(1) = bidV
 	Vavg_tot(2) = Vavg_tot(1)
 	DoLP_tot(1)=200.*sqrt(bidVcos**2+bidVsin**2)/bidV
 	DoLP_tot(2)= DoLP_tot(1)
 	AoLP_tot(1)=180./pi* atan2(bidVsin,bidVcos)/2.- calib_angle
 	AoLP_tot(2)= AoLP_tot(1)
c
        call moyglis(V_inst_corr,nrot,nmoy,V_glis)
c
c 	Retire les éclairs sur une moyenne de 20 rotations
        call moyglis(V_inst_corr,nrot,20,V_glis)
	do irot=2,nrot
 	  if (V_inst_corr(irot).gt.V_glis(irot)+0.5)then
 	    V_inst_corr(irot)=V_inst_corr(irot-1) 
 	    Vcos_inst(irot)=Vcos_inst(irot-1) 
 	    Vsin_inst(irot)=Vsin_inst(irot-1) 
 	  endif
 	enddo
 	
c 	Retourne à la moyenne demandée
        call moyglis(V_inst_corr,nrot,nmoy,V_glis)
c
        call moyglis(Vcos_inst,nrot,nmoy,Vcos_glis)
        call moyglis(Vsin_inst,nrot,nmoy,Vsin_glis)
 	do irot=1,nrot
 	  DoLP_glis(irot)=2.*sqrt(Vcos_glis(irot)**2+Vsin_glis(irot)**2)
     .                    /V_glis(irot)
 	  DoLP_glis(irot)=100.* DoLP_glis(irot)
 	  AoLP_glis(irot)=atan2(Vsin_glis(irot),Vcos_glis(irot))/2.
 	  AoLP_glis(irot)=180./pi* AoLP_glis(irot) - calib_angle
c 	  On va essayer de supprimer les sauts de phase...
 	  if(AoLP_glis(irot).lt.AoLP_tot(1)-100.)
     .       AoLP_glis(irot)=AoLP_glis(irot)+180.
 	enddo
        call moyglis(DoLP_inst,nrot,nmoy,DoLP_moy)
        call moyglis(AoLP_inst,nrot,nmoy,AoLP_moy)
c
c       Prints
c  	------
 	write(20,1010)
 	do irot = 1,nrot
 	  write(20,1020)t_inst(irot),time_inst(irot),
     .       V_inst_corr(irot),V_glis(irot),
     .       Vcos_inst(irot),Vcos_glis(irot),
     .       Vsin_inst(irot),Vsin_glis(irot)
 	enddo
c
 	write(20,1030)
 	do irot = 1,nrot
 	  write(20,1040)t_inst(irot),time_inst(irot),
     .       DoLP_inst(irot),DoLP_glis(irot),
     .       DoLP_moy(irot),DoLP_tot(1),
     .       AoLP_inst(irot),AoLP_glis(irot),
     .       AoLP_moy(irot),AoLP_tot(1)
 	enddo
c
c       Plots
c  	-----
        call GR_EXEC('set /default')
        call GR_EXEC('clear ')
        call GR_EXEC('set char .8')
c
        call mnmx(time_inst,nrot,tmin,tmax,0)
c
c 	Intensity (mV)
c 	--------------
        call mnmx(V_inst_corr,nrot,ymin,ymax,0)
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
        call gr4_give('X',2,time_tot)
        call gr4_give('Y',2,Vavg_tot)
        call GR_EXEC('connect')
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,V_inst_corr)
        call GR_EXEC('connect')
        call GR_EXEC('pen 1 /weight 3')
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,V_glis)
        call GR_EXEC('connect')
c
c       DoLP (%)
c       --------
        call mnmx(DoLP_inst,nrot,ymin,ymax,0)
        call mnmx(DoLP_glis,nrot,yymin,yymax,0)
 	if(ymax.gt.yymax*1.5)ymax = yymax*1.5
        call GR_EXEC('set view .2 .9 .40 .65')
        call GR_EXEC('pen 0 /weight 3') 
        call gr_limi(4,tmin,tmax,ymin,ymax)
        call GR_EXEC('box N') 
        call GR_EXEC('label "DoLP" /Y')
        call GR_EXEC('set view .24 .9 .40 .65')
        call GR_EXEC('label "[%]" /Y')
        call GR_EXEC('set view .2 .9 .40 .65')
c
        call GR_EXEC('pen 0 /weight 3') 
        call gr4_give('X',2,time_tot)
        call gr4_give('Y',2,DoLP_tot)
        call GR_EXEC('connect')
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,DoLP_inst)
        call GR_EXEC('connect')
        call GR_EXEC('pen 3 /weight 3') 
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,DoLP_moy)
c       call GR_EXEC('connect')
        call GR_EXEC('pen 1 /weight 3') 
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,DoLP_glis)
        call GR_EXEC('connect')
c
c       AoLP (%)
c       --------
        call mnmx(AoLP_inst,nrot,ymin,ymax,0)
        call mnmx(AoLP_glis,nrot,yymin,yymax,0)
 	ymin = min(ymin,yymin)
 	ymax = max(ymax,yymax)
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
        call gr4_give('X',2,time_tot)
        call gr4_give('Y',2,AoLP_tot)
        call GR_EXEC('connect')
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,AoLP_inst)
        call GR_EXEC('connect')
        call GR_EXEC('pen 3 /weight 3')
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,AoLP_moy)
c       call GR_EXEC('connect')
        call GR_EXEC('pen 1 /weight 3')
        call gr4_give('X',nrot,time_inst)
        call gr4_give('Y',nrot,AoLP_glis)
        call GR_EXEC('connect')

        call GR_EXEC('pen 0 /weight 3')
        call GR_EXEC('set view .2 .9 .15 .90')
        call GR_EXEC('label "Local Time" /X')
        call GR_EXEC('DRAW TEXT .0 0.1 "'//labelplot//'" 8 0. /BOX 8')
        if(entergreg.eq.1)call gmaster_enter_loop('')
        call GR_EXEC('hardcopy "mens" /dev eps col /over')

 	  
1010 	format(5x,'T since   Local Time  V_inst_cor  V_glis',
     .    '     Vcos_inst   Vcos_glis   Vsin_inst    Vsin_glis',
     .       /,7x,'start        [s]        [V]        [V]',
     .    '         [V]       [V]           [V]        [V]')
1020 	format(2f12.4,3x,6(0pf10.7,2x))
1030 	format(/,5x,'  T/0      Loc Time   DoLP_inst   DoLP_glis',
     .          '   DoLP_moy    DoLP tot',
     .          '   AoLP_inst   AoLP_glis',
     .          '   AoLP_moy   AoLP_tot')
1040 	format(2f12.4,3x,4(0pf10.7,2x),4(0pf10.5,2x))
c
 	stop
 	end
c
c -----------------------------------------------------------------
c
        subroutine xline(nbline,ijfile)
c
        character nc*30
        do i=1,nbline
           read(ijfile,1000) nc
        enddo
c
1000    format(a1)
c
        return
        end
c
c-------------------------------------------------------------------
c
        subroutine moyglis(tab,ntab,nmoy,tabmoy)
c
        implicit none
        integer ntab,itab,nmoy,imoy
        integer ndeb,nfin
        real tab(ntab),tabmoy(ntab)
        real parite
c
        if(nmoy.eq.1)then
          do itab = 1,ntab
            tabmoy(itab)=tab(itab)
          enddo
          return
        endif
        if (nmoy.ge.ntab)then
         do itab=1,ntab
           tabmoy(itab)=tab(itab)
         enddo
         return
        endif
        call zeroit(tabmoy,ntab)
c
        parite = float(nmoy/2) - float(nmoy)/2.
        if (parite.eq.0.)then           !nombre pair
          ndeb = -(nmoy/2-1)
          nfin = nmoy/2
        else
          ndeb = -nmoy/2
          nfin = nmoy/2
        endif
c
        do itab = 1-ndeb,ntab-nfin
          tabmoy(itab)=0.
          do imoy=ndeb,nfin
            tabmoy(itab)=tabmoy(itab)+tab(itab+imoy)
          enddo
          tabmoy(itab)=tabmoy(itab)/float(nmoy)
        enddo
c       premiers points
        tabmoy(1)=tab(1)
        if (abs(ndeb).ge.2)then
          do imoy=2,abs(ndeb)
            tabmoy(imoy)=tabmoy(imoy-1)+tab(imoy)
          enddo
          do imoy=2,abs(ndeb)
            tabmoy(imoy)=tabmoy(imoy)/float(imoy)
          enddo
        endif
c       derniers points
        tabmoy(ntab) = tab(ntab)
        do itab = ntab-1,ntab-nfin+1,-1
          tabmoy(itab)=tabmoy(itab+1)+tab(itab)
        enddo
        do itab = ntab-1,ntab-nfin+1,-1
          tabmoy(itab)=tabmoy(itab)/float(ntab-itab+1)
        enddo

        return
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
c------------------------- mnmx --------------------------------------
c
        subroutine mnmx(tab,ntab,tmin,tmax,linlog)
c
        implicit none
        integer ntab,linlog,i
        real tab(ntab),tmin,tmax
c
        if(linlog.eq.0)then
          tmin=tab(1)
          tmax=tab(1)
          do 1 i=2,ntab
            tmax=max(tmax,tab(i))
            tmin=min(tmin,tab(i))
1         continue
        else
c         finds the first min non equal to zero, and the max of tab
          tmax=tab(1)
          do 20 i=2,ntab
            tmax=max(tmax,tab(i))
20         continue
          if(tmax.le.0.)then
            write(6,*)'Max <= 0, dessin log impossible'
            go to  60
          endif
          do 30 i=1,ntab
            if (tab(i).gt.0.)then
              tmin=tab(i)
              go to 40
            endif
 30       continue
 40       continue
          do 50 i=1,ntab
            if (tab(i).lt.tmin.and.tab(i).gt.0.) tmin=tab(i)
 50       continue
 60       continue
        endif
c
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
 	
