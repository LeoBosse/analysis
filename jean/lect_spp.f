        subroutine lect
c
        implicit none
c
c       Un cycle dure 4 secondes environ, soit 21600 cycles max par jour
c----   Pendant la première campagne, 248 parametres sont sauvegardés
c       correspondant à :
c       an mois jour hr Mn Sec 
c       pan_pos tilt_pos  = azimuth, elevation, geog. coord.
c       pol_flt_arr(80),canalP(80),canalT(80)
c----   A partir de la seconde campagne,254 parametres sont sauvegardés 
c       correspondant à :
c       an mois jour hr Mn Sec 
c       voltage (of main supply)
c       PosV = positive voltage, NegV = negative voltage
c       MotV = motor voltage 
c       Btemperature, Ftemperature = temperature of Base (photomutip), Filter
c       tilt_angle, pan_angle = elevation, azimuth (geog. coord.)
c       pol_flt_ang(iang = 1,80),canalP(iang = 1,80),canalT(iang = 1,80)
c       Chaque paramètre  est sauvegardé sur 1 mesure
c       On ne peut pas lire un fichier binaire en plusieurs fois: il
c       faut lire tout d'un coup.
c
c       canalP = canal polarisé (1 dans la 1ère camp., 2 dans la  2ème)
c       canalT = canal Int. Tot (2 dans la 1ère camp., 1 dans la  2ème)
c
        integer nbrmax_points_par_cycle,nbrdump,nbrang
        integer nbr_oct_par_point,lchaine
        parameter (nbrmax_points_par_cycle=254,nbrdump=45000)
        parameter (nbrang=80,nbr_oct_par_point = 2)
c
c       Donnnees brutes en binaire
        integer*2 tab(2*nbrmax_points_par_cycle*nbrdump)
c       nmoy est le pas pour la moyenne glissante amont
        integer  taille_en_bits,nmoy,nmoyplot
        integer nbr_points_total,ipoint
        character*50 ficin
        character*50 dirin
c
c       Lus dans le fichier spp
        real an,mois,jour,hr,Mn,Sec
 	real numjour_hrloc(nbrdump),numjour_ut(nbrdump)
        integer voltage,PosV,NegV, MotV
        real Btemper(nbrdump),Ftemper(nbrdump)
        real azym(nbrdump),elevation(nbrdump)
        real pol_flt_ang(nbrdump,nbrang)
        real canalP_brut(nbrdump,nbrang),canalT_brut(nbrdump,nbrang)
c
c       calculé dans le programme
        real ut(nbrdump),hrloc(nbrdump)
        real dec,chi(nbrdump),chideg(nbrdump)
        real angrad(nbrang),angdeg(nbrang),utdump(nbrdump,nbrang)
        real canalP_corr(nbrdump,nbrang)
        real canalT_corr(nbrdump,nbrang)
        real canalP_norm(nbrdump,nbrang)
        real canalP_moy(nbrdump,nbrang),canalT_moy(nbrdump,nbrang)
        real rap_TsurP(nbrdump,nbrang),rap_TsurP_moy(nbrdump,nbrang)
c       secondes represente le nombre de secondes depuis le début de 
c       la manip
        real deltaT(nbrdump),secondes(nbrdump)
 	real DarkP(nbrdump),DarkT(nbrdump)
 	real DarkPmeas,DarkTmeas ! negative if estimation requested
c
 	real Ang_Calib_rad,Ang_Calib_deg
 	integer errstat
c
c 	intbrut : données spp
c 	intcorr : corrigees du Dark etc.
c 	intnorm : P est normalisée à T sur l'intensité pendant 4s
c 	          et cette valeur est donc égale à intcorrT
c 	          (mais pas bien sûr les valeurs angle par angle)
c 	intmoy : moyenne du précédent sur norm pas (programmé dans input)
        real intbrutP(nbrdump),intbrutT(nbrdump)
        real intmoy(nbrdump),rapmoy(nbrdump)
        real intnormP(nbrdump)
        real intcorrP(nbrdump),intcorrT(nbrdump)
c
c       internes:
        real utdeb,utfin
        real y0,ym1,yp1,ym2,yp2,ym3,yp3
 	integer test
        integer ndump,idump,jdump,nang,iang,nbr_points_par_cycle
 	integer irot
 	integer pointeur,pointeurP,pointeurT,idess
        integer kelanalyse,ndepart
        real tampon,pi,moyP,moyT,facteur
        real workin(nbrdump),workout(nbrdump)
 	integer n_bad_dumpsT,n_bad_dumpsP
        integer num_bad_dumpT(nbrdump),num_bad_dumpP(nbrdump)
 	integer na2v
 	data na2v/20/		! numero angle pour les 2 voies
 	logical mode_dark	! Estimation du dark si pas connu
c
 	real ymin,ymax,yymin,yymax,utmin,utmax
        character*8 yyyymmdd
        integer quelle_campagne 
c 	 1 = 1ère campagne Advantalen 
c 	 2 = 2ème campagne KHO
c 	 3 = 3ème campagne Hornsund
c 	 4 = 4ème campagne ESR Longyear / KHO
c 	 5 = 5ème campagne ESR Longyear / EISCAT Skybotn
c 	 6 = à partir de 2016, les caneaux remis dans le bon ordre
 	character entete*80
 	character entete_date*60
c
1000    format(a13)
1005    format(1a12)
1010    format(3i5,'   POLARLIS',/,f10.3,'  declinaison')
1015    format(i8,' points. ',
     .          ' Moyenne pour la phase faite sur ',i5,' points')
1020    format(10x,
     .   'hrloc       UT        chideg      delta T',
     .   '   intbrutP    intbrutT   intcorrP    intcorrT  ',
     .   '  P Norm. a T    intmoy    ',
     .   ' rapp T/P  Azym     Elev    T Base   T filt')
1030    format(i5,1x,2(1f10.5,1x),1f11.5,1x,1f10.3,1x,
     .         7(1f11.5,1x),2(1f7.2,1x),2(1f9.2,1x),i5)
1040    format(/,'dump ',i5,' sza=',f12.6,'°. UT=',f12.6,' hrloc=',
     .     f12.6,'   DarkP 'f12.6,/,
     .     '           UT        angle      P brut      T brut',
     .     '      P corr       T corr    Pnorm à T  ',
     .     ' canalP_moy  canalT_moy   moy. T/P')
1050    format(i5,1x,f11.7,10(1f11.5,1x))
1060    format(i4,'0',i1,'0',i1)
1070    format(i4,'0',i1,i2)
1080    format(i4,i2,'0',i1)
1090    format(i4,i2,i2)
1100    format(' Il y a ',i5,' cycles sur ',i5,
     .          ' avec un angle solaire > 108°')
1105    format('Il y a ',i5,' cycles sur ',i5)
1110    format(' Analyse du fichier ',a12,'/',a18)
1120    format('Il y a ',i5,' cycles sur ',i5,
     .          '  entre ',f5.2,' et ',f5.2,' UT')
1130 	format('Estimated Dark current for channel P: ',1f10.3,/,
     .         '                       for channel T: ',1f10.3)
c
        pi = 4.*atan(1.)
c
        open(10,file='polar.input',status='old')
        rewind(10)
        read(10,*)dirin
        read(10,*)ficin
        read(10,*)taille_en_bits
        read(10,*)DarkPmeas,DarkTmeas  ! ajout fevrier 2014
        read(10,*)nmoy
 	if (nmoy.eq.0)nmoy=1
        read(10,*)nmoyplot
        read (10,*)kelanalyse
 	utdeb = 0.
 	utfin = 240.    	! 10 jours...
        if (kelanalyse.eq.3)then
          read (10,*)utdeb,utfin
          if (utdeb.ge.utfin)then
            write(6,*)'Heure de debut d''impression >  heure de fin!!!'
            write(6,*)'Programme arrete'
            stop
          endif
 	else
          call xline(1,10)
        endif
c ---   Ajout 2012
 	read(10,*)idess
        call xline(9,10)
  	read(10,*)n_bad_dumpsT
 	if (n_bad_dumpsT.ne.0)then
          read(10,*)(num_bad_dumpT(idump),idump=1,n_bad_dumpsT)
 	else
 	  call xline(1,10)
 	endif
 	n_bad_dumpsP = 0
c ---   Ajout 2015
 	errstat = 0
 	read (10,*,iostat=errstat)Ang_Calib_deg
c	write(6,*)'Ang_Calib_deg?'
c	read (5,*)Ang_Calib_deg
 	if(errstat.ne.0)then
 	  write(6,*)'Attention ! pas de calibration d''angle en fin du'
 	  write(6,*)'fichier polar.input.'
 	  write(6,*)'Mise a 0 artificiellement'
 	  Ang_Calib_deg = 0.
 	endif
 	Ang_Calib_rad = Ang_Calib_deg*pi/180.
        close(10)
c
        write(6,1110)dirin, ficin
c
c       Determine la campagne de mesure (dont dépend le format de
c       fichier).
        open(10,file='../'//dirin(1:lchaine(dirin))//
     .       '/'//ficin(1:lchaine(ficin)),
     .       form='unformatted',access='direct',recl=4)

        read(10,rec=1)(tab(idump),idump=1,2)        ! an, mois
        close(10)
c
        nang=nbrang
        quelle_campagne = 5
        nbr_points_par_cycle  = 254
        if(tab(1).eq.2006)quelle_campagne=1
        if(tab(1).eq.2007 .and. tab(2).eq.1)quelle_campagne=1
        if(quelle_campagne.eq.1) nbr_points_par_cycle  = 248
        if(tab(1).eq.2010)quelle_campagne=3	! Hornsund
        if(tab(1).eq.2012)quelle_campagne=4	! ESR/KHO 2012
        if(tab(1).eq.2014)quelle_campagne=5	! ESR/KHO 2014 / Skybotn
        if(tab(1).ge.2016)quelle_campagne=6	! Suivantes
 	
c
c       On saute le dernier cycle qui peut être mal stocké
        nbr_points_total= taille_en_bits/nbr_oct_par_point
c    .          - nbr_points_par_cycle
        ndump=nbr_points_total/nbr_points_par_cycle 
        write(6,*) nbr_points_par_cycle,'points par cycle et',
     .        ndump,' cycles de mesure'
c
        write(6,*) 'Lecture'
        open(10,file='../'//dirin(1:lchaine(dirin))//
     .       '/'//ficin(1:lchaine(ficin)),
     .       form='unformatted',access='direct',
     .       recl=nbr_oct_par_point*nbr_points_total)
        read(10,rec=1)(tab(ipoint),ipoint=1,nbr_points_total)
        close(10)
c
c       an mois jour hr Mn Sec 
c       voltage, PosV, NegV, MotV, Btemper, Ftemper
c       tilt_angle, pan_angle
c       pol_flt_ang(iang = 1,80),canalP(iang = 1,80),canalT(iang = 1,80)

c       ndepart permet de sauter les premier dumps si nécessaire
        ndepart = 0
        ndump = ndump - ndepart
        hr  = 0.
        mn  = 0.
        sec = 0.
        do idump = 1,ndump
          pointeur = nbr_points_par_cycle*(ndepart+idump-1)
          an       = float(int(tab( 1+pointeur)))
          mois     = float(int(tab( 2+pointeur)))
          jour     = float(int(tab( 3+pointeur)))
          hr       = float(int(tab( 4+pointeur)))
          mn       = float(int(tab( 5+pointeur)))
          sec      = float(int(tab( 6+pointeur)))
          if(mois.lt.10 .and. jour.lt.10)then
            write(yyyymmdd,1060)int(an),int(mois),int(jour)
          else if (mois .lt.10 .and. jour.ge.10)then
            write(yyyymmdd,1070)int(an),int(mois),int(jour)
          else if (mois .ge.10 .and. jour.lt.10)then
            write(yyyymmdd,1080)int(an),int(mois),int(jour)
          else 
            write(yyyymmdd,1090)int(an),int(mois),int(jour)
          endif
c
          ut(idump) = hr+(mn*60.+sec)/3600.
          secondes(idump)=hr*3600.+mn*60.+sec
          call asz(jour,UT(idump),78.20,15.83,hrloc(idump),dec,
     .          chi(idump),chideg(idump))
c
          if(quelle_campagne.eq.1) then        ! campagne 1
            azym(idump)   =float(int(tab( 7+pointeur)))   ! pan angle
            elevation(idump) =float(int(tab( 8+pointeur)))   ! tilt angle
            do iang=1,nbrang
c             angdeg(idump,iang) = 
c    .          float(tab(iang+ 8+pointeur))
c             angrad(idump,iang)=angdeg(idump,iang)*pi/180.
              canalP_brut(idump,iang)= 
     .          float(int(tab(iang+88+ pointeur)))
              canalT_brut(idump,iang)= 
     .          float(int(tab(iang+168+ pointeur)))
            enddo
          elseif(quelle_campagne.gt.2 .and. quelle_campagne.lt.6)then
c           campagne 2 et plus
            voltage   = int(tab( 7+pointeur))
            PosV      = int(tab( 8+pointeur))
            NegV      = int(tab( 9+pointeur))
            MotV      = int(tab(10+pointeur))
            Btemper(idump)   = float(int(tab(11+pointeur)))
            Ftemper(idump)   = float(int(tab(12+pointeur)))
            elevation(idump) =float(int(tab(13+pointeur)))   ! tilt angle
            azym(idump)   =float(int(tab(14+pointeur)))   ! pan angle
            do iang=1,nbrang
              canalT_brut(idump,iang)= 
     .          abs(float(int(tab(iang+94+ pointeur))))
              canalP_brut(idump,iang)= 
     .          abs(float(int(tab(iang+174+ pointeur))))
            enddo
          elseif(quelle_campagne.eq.6)then
c           Avec Alain, on a remis les voies P et T dans l'ordre
            voltage   = int(tab( 7+pointeur))
            PosV      = int(tab( 8+pointeur))
            NegV      = int(tab( 9+pointeur))
            MotV      = int(tab(10+pointeur))
            Btemper(idump)   = float(int(tab(11+pointeur)))
            Ftemper(idump)   = float(int(tab(12+pointeur)))
            elevation(idump) =float(int(tab(13+pointeur)))   ! tilt angle
            azym(idump)   =float(int(tab(14+pointeur)))   ! pan angle
            do iang=1,nbrang
              canalP_brut(idump,iang)= 
     .          abs(float(int(tab(iang+94+ pointeur))))
              canalT_brut(idump,iang)= 
     .          abs(float(int(tab(iang+174+ pointeur))))
            enddo
          endif
        enddo           ! sur les dumps (= cycles)
c 	
c       write(6,*)'Lecture des donnees corrigees par Pierre-Olivier'
 	if(an.eq.2010. .and. mois.eq.11. .and. jour.eq.7.)then
          open(37,file='../dir.20101107/20101107_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.11. .and. jour.eq.8.)then
          open(37,file='../dir.20101108/20101108_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.11. .and. jour.eq.9.)then
          open(37,file='../dir.20101109/20101109_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.11. .and. jour.eq.10.)then
          open(37,file='../dir.20101110/20101110_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.2.)then
          open(37,file='../dir.20101202/20101202_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.3.)then
          open(37,file='../dir.20101203/20101203_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.5.)then
          open(37,file='../dir.20101205/20101205_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.6.)then
          open(37,file='../dir.20101206/20101206_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.7.)then
          open(37,file='../dir.20101207/20101207_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.8.)then
          open(37,file='../dir.20101208/20101208_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.9.)then
          open(37,file='../dir.20101209/20101209_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.10.)then
          open(37,file='../dir.20101210/20101210_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	elseif(an.eq.2010. .and. mois.eq.12. .and. jour.eq.11.)then
          open(37,file='../dir.20101211/20101211_denoise.txt')
 	  write(6,*)'lecture val. debruitees POA pour ',an,mois,jour
 	else
 	  go to 38
 	endif
 	do idump = 1,ndump
 	  do iang = 1,nbrang
            read (37,*)canalP_brut(idump,iang)
 	  enddo
 	enddo
 	close(37)
38 	continue
c
        do idump=1,ndump-1
          secondes(idump)= secondes(idump)-secondes(1)
          if (hrloc(idump+1).ge.hrloc(idump)) then
            deltaT(idump)=(hrloc(idump+1)-hrloc(idump))*3600.
          else          ! On passe minuit
            deltaT(idump)=(hrloc(idump+1)+24.-hrloc(idump))*3600.
          endif
        enddo
        deltaT(ndump)=deltaT(ndump-1)
        secondes(ndump)= secondes(ndump)-secondes(1)
c
 	numjour_hrloc(1) = 0.
 	numjour_ut(1) = 0.
        do idump=2,ndump
          if (hrloc(idump-1).gt.hrloc(idump))then
             numjour_hrloc(idump) = numjour_hrloc(idump-1)+1.
 	  else
             numjour_hrloc(idump) = numjour_hrloc(idump-1)
 	  endif

          if (ut(idump-1).gt. ut(idump))then
             numjour_ut(idump) = numjour_ut(idump-1)+1.
 	  else
             numjour_ut(idump) = numjour_ut(idump-1)
 	  endif
 	enddo
c
        do idump=2,ndump
          hrloc(idump) = hrloc(idump)+24.*numjour_hrloc(idump)
          ut(idump) = ut(idump)+24.*numjour_ut(idump)
 	enddo
c
c       Calcul des heures de mesures des 80 pts de cycles
        do idump=1,ndump
          utdump(idump,1)=ut(idump)
          do iang=2,nbrang
            utdump(idump,iang)=utdump(idump,iang-1)
     .        + deltaT(idump)/(3600.*float(nbrang))
          enddo
        enddo
c 
c       Calcul des angles corrects
        iang=1
c       angdeg(iang)=0.
c       angrad(iang)=0.
c 	Attention ! Il faut compter les angles a partir du premier
c 	span, sinon on n'obtient pas les bonnes calibrations!
        angdeg(iang)=4.0 + Ang_Calib_deg
        angrad(iang)=angdeg(iang)*pi/180.
 	
        do iang= 2,nbrang
          angdeg(iang)= angdeg(iang-1)+360./float(nbrang)
          angrad(iang)=angdeg(iang)*pi/180.
        enddo
c
c       write(6,*) 'Soustraction du bruit et normalisation des données'
c
c       Correction du courant d'obscurité
        do idump=1,ndump
          if(quelle_campagne.eq.1)then         ! campagne 1
            DarkP(idump)=1.+ 20.*exp(-secondes(idump)/1650.)
     .           + 35.*exp(-secondes(idump)/12000.)
 	    DarkT(idump)=DarkP(idump)
          else if(quelle_campagne.eq.2)then         ! campagnes KHO
            DarkP(idump)=4.      ! On retire 4 coups
 	    DarkT(idump)=DarkP(idump)
          else if(quelle_campagne.eq.3)then         ! campagne Hornsund 
            DarkP(idump)=0.
 	    DarkT(idump)=DarkP(idump)
          else if(quelle_campagne.eq.4)then         ! ESR / KHO 2012
            DarkP(idump)=2.
 	    DarkT(idump)=DarkP(idump)
          else if(quelle_campagne.gt.4)then         ! Suivantes
            DarkP(idump)=DarkPmeas
 	    DarkT(idump)=DarkTmeas
          endif
          do iang=1,nbrang
            canalP_corr(idump,iang)=canalP_brut(idump,iang)-DarkP(idump)
            canalP_corr(idump,iang)=max(0.,canalP_corr(idump,iang))
            canalT_corr(idump,iang)=canalT_brut(idump,iang)-DarkT(idump)
            canalT_corr(idump,iang)=max(0.,canalT_corr(idump,iang))
          enddo
c         Correction du premier et du dernier angle
c	  Cela augmente la valeur intégrée sur tous les angles car il y a 
c 	  beaucoup de valeurs initiales et finales nulles !
          if(canalP_corr(idump,1).eq.0.)
     .            canalP_corr(idump,1)= canalP_corr(idump,2)
          if(canalT_corr(idump,1).eq.0.)
     .          canalT_corr(idump,1)= canalT_corr(idump,2)
          do iang =2,nbrang
            if(canalP_corr(idump,iang).eq.0.)
     .          canalP_corr(idump,iang)= canalP_corr(idump,iang-1)
            if(canalT_corr(idump,iang).eq.0.)
     .          canalT_corr(idump,iang)= canalT_corr(idump,iang-1)
          enddo
        enddo
c
c ---   Ajout 2012 : les dumps mauvais sont sautés à partir de la 
c 	lecture dans polar.input
 	if (n_bad_dumpsT.ne.0)then
 	  do jdump = 1,n_bad_dumpsT
 	    idump = num_bad_dumpT(jdump)
            do iang = 1,nbrang
              canalT_corr(idump,iang)=canalT_corr(idump-1,iang)
            enddo
          enddo
 	endif
 	if (n_bad_dumpsP.ne.0)then
 	  do jdump = 1,n_bad_dumpsP
 	    idump = num_bad_dumpP(jdump)
            do iang = 1,nbrang
              canalP_corr(idump,iang)=canalP_corr(idump-1,iang)
            enddo
          enddo
 	endif
c
c ----  Moyennes glissantes AMONT en temps
c 	Si nmoy = 1 (ou 0), pas de moyenne
c	On commence par moyenner chaque angle sur un temps donné.
c 	1 mn = 15 dumps
 	if(nmoy.gt.1)then
 	  do iang = 1,nbrang
 	    do idump = 1,ndump
 	      workin(idump) = canalP_corr(idump,iang)
 	    enddo
            call moyglis(workin,ndump,nmoy,workout)
 	    do idump = 1,ndump
 	      canalP_corr(idump,iang) = workout(idump)
 	    enddo
 	  enddo
 	  do iang = 1,nbrang
 	    do idump = 1,ndump
 	      workin(idump) = canalT_corr(idump,iang)
 	    enddo
            call moyglis(workin,ndump,nmoy,workout)
 	    do idump = 1,ndump
 	      canalT_corr(idump,iang) = workout(idump)
 	    enddo
 	  enddo
 	endif
c
c       calcul du rapport T sur P cycle par cycle
        do idump=1,ndump
          do iang=1,nbrang
            rap_TsurP(idump,iang)=0.
            if(canalT_corr(idump,iang).ne.0.)
     .      rap_TsurP(idump,iang)=canalT_corr(idump,iang)/
     .                          canalP_corr(idump,iang)
          enddo
        enddo

c       Normalisation du canal P aux intensités du canal T
        do idump=1,ndump
          moyP=0.
          moyT=0.
          do iang=1,nbrang
            moyP=moyP+ canalP_corr(idump,iang)
            moyT=moyT+ canalT_corr(idump,iang)
          enddo
          moyP=moyP/float(nbrang)
          moyT=moyT/float(nbrang)
          do iang=1,nbrang
            if(moyP.ne.0)
     .        canalP_norm(idump,iang)= canalP_corr(idump,iang)*moyT/moyP
c    .        canalP_norm(idump,iang)= canalP_corr(idump,iang)
          enddo
        enddo
c
c 	Ce qui suit est un reliquat du passé quand je ne savais plus
c 	où faire les moyennes. Je garde ce transfert de nom quoique 
c 	"moy" est plutôt mal choisi...  (jl, fevrier 2014)
        do iang = 1,nbrang
          do idump=1,ndump
            canalP_moy(idump,iang)=canalP_norm(idump,iang)
          enddo
          do idump=1,ndump
            canalT_moy(idump,iang)= canalT_corr(idump,iang)
          enddo
          do idump=1,ndump
            rap_TsurP_moy(idump,iang)= rap_TsurP(idump,iang)
          enddo
        enddo
c       
c       calcul des intensités moyennes dump par dump
        do idump=1,ndump
          intbrutP(idump)=0.
          intbrutT(idump)=0.

          intcorrP(idump)=0.
          intcorrT(idump)=0.

          intnormP(idump)=0.

          intmoy(idump)=0.
          rapmoy(idump)=0.

	  nang = 1
          do iang=1,nbrang
 	    if(canalP_brut(idump,iang).ne.0.)then
              intbrutP(idump)=intbrutP(idump)+canalP_brut(idump,iang)
              intbrutT(idump)=intbrutT(idump)+canalT_brut(idump,iang)
              intcorrP(idump)=intcorrP(idump)+canalP_corr(idump,iang)
              intcorrT(idump)=intcorrT(idump)+canalT_corr(idump,iang)
              intnormP(idump)=intnormP(idump)+canalP_norm(idump,iang)
              intmoy(idump)=intmoy(idump)+canalP_moy(idump,iang)
              rapmoy(idump)=rapmoy(idump)+rap_TsurP_moy(idump,iang)
 	      nang = nang+1
 	    endif
          enddo
          intbrutP(idump)=intbrutP(idump)/float(nang)
          intbrutT(idump)=intbrutT(idump)/float(nang)
          intcorrP(idump)=intcorrP(idump)/float(nang)
          intcorrT(idump)=intcorrT(idump)/float(nang)
          intnormP(idump)=intnormP(idump)/float(nang)
          intmoy(idump)=intmoy(idump)/float(nang)
          rapmoy(idump) =rapmoy(idump)/float(nang)
        enddo
c
c       Ecriture de toutes les  mesures dans lect.output
c       Le dernier dump peut être corrompu 
c       si l'instrument s'est arrêté en cours de route
c       write(6,*)'Ecriture'
c
        open(20,file='../'//dirin(1:lchaine(dirin))//
     .       '/lect_spp_modsynch.output')
c
 	write(20,*)quelle_campagne,'    num de campagne'
        if (kelanalyse.eq.1)then        ! toute la manip
          pointeur = 0
          do idump = 1,ndump-1       ! pas le  dernier dump
            pointeur = pointeur +1
          enddo
          write(20,1015)pointeur,nmoy
          write(20,1010)ifix(an),ifix(mois),ifix(jour),dec
          write(20,1020)
          do idump = 1,ndump-1       ! pas le  dernier dump
            write(20,1030)
     .        idump,hrloc(idump),ut(idump),chideg(idump),deltaT(idump),
     .        intbrutP(idump),intbrutT(idump),
     .        intcorrP(idump),intcorrT(idump),
     .        intnormP(idump),
     .        intmoy(idump),rapmoy(idump),
     .        azym(idump),elevation(idump),Btemper(idump),Ftemper(idump)
          enddo
c         Ecriture des intensités mesurées et corrigees
          do idump=1,ndump-1
            write(20,1040) idump,chideg(idump),ut(idump),hrloc(idump),
     .            DarkP(idump)
            do iang=1,nbrang
              write(20,1050)iang,utdump(idump,iang),angdeg(iang),
     .          canalP_brut(idump,iang),canalT_brut(idump,iang),
     .          canalP_corr(idump,iang),canalT_corr(idump,iang),
     .          canalP_norm(idump,iang),
     .          canalP_moy(idump,iang),canalT_moy(idump,iang),
     .          rap_TsurP(idump,iang)
            enddo
          enddo
          close(20)
          write(6,1105)pointeur,ndump
        elseif (kelanalyse.eq.2)then
          pointeur=0
          do idump=1,ndump-1
            if(chideg(idump).ge.108)then
              pointeur = pointeur +1
            endif
          enddo
          write(20,1015)pointeur,nmoy
          write(20,1010)ifix(an),ifix(mois),ifix(jour),dec
          write(20,1020)
          do idump = 1,ndump-1       ! pas le  dernier dump
            if(chideg(idump).ge.108)then
            write(20,1030)
     .        idump,hrloc(idump),ut(idump),chideg(idump),deltaT(idump),
     .        intbrutP(idump),intbrutT(idump),
     .        intcorrP(idump),intcorrT(idump),
     .        intnormP(idump),
     .        intmoy(idump),rapmoy(idump),
     .        azym(idump),elevation(idump),Btemper(idump),Ftemper(idump)
            endif
          enddo
c         Ecriture des intensités mesurées et corrigees
          do idump=1,ndump-1
            if(chideg(idump).ge.108)then
              write(20,1040) idump,chideg(idump),ut(idump),hrloc(idump),
     .          DarkP(idump)
              do iang=1,nbrang
              write(20,1050)iang,utdump(idump,iang),angdeg(iang),
     .          canalP_brut(idump,iang),canalT_brut(idump,iang),
     .          canalP_corr(idump,iang),canalT_corr(idump,iang),
     .          canalP_norm(idump,iang),
     .          canalP_moy(idump,iang),canalT_moy(idump,iang),
     .          rap_TsurP(idump,iang)
              enddo
            endif
          enddo
          close(20)
          write(6,1100)pointeur,ndump
        elseif (kelanalyse.eq.3)then
          pointeur=0
          do idump=1,ndump-1
            if(ut(idump).ge.utdeb .and. ut(idump).le.utfin)then
              pointeur = pointeur +1
            endif
          enddo
          write(20,1015)pointeur,nmoy
          write(20,1010)ifix(an),ifix(mois),ifix(jour),dec
          write(20,1020)
          pointeur=0
          do idump = 1,ndump-1
            if(ut(idump).ge.utdeb .and. ut(idump).le.utfin)then
              pointeur = pointeur +1
            write(20,1030)
     .        idump,hrloc(idump),ut(idump),chideg(idump),deltaT(idump),
     .        intbrutP(idump),intbrutT(idump),
     .        intcorrP(idump),intcorrT(idump),
     .        intnormP(idump), intmoy(idump),rapmoy(idump),
     .        azym(idump),elevation(idump),Btemper(idump),
     .        Ftemper(idump),pointeur
            endif
          enddo
c         Ecriture des intensités mesurées et corrigees
          do idump=1,ndump-1
            if(ut(idump).ge.utdeb .and. ut(idump).le.utfin)then

              write(20,1040) idump,chideg(idump),ut(idump),hrloc(idump),
     .          DarkP(idump)
              do iang=1,nbrang
                write(20,1050)iang,utdump(idump,iang),angdeg(iang),
     .          canalP_brut(idump,iang),canalT_brut(idump,iang),
     .          canalP_corr(idump,iang),canalT_corr(idump,iang),
     .          canalP_norm(idump,iang),
     .          canalP_moy(idump,iang),canalT_moy(idump,iang),
     .          rap_TsurP(idump,iang)
              enddo
            endif
          enddo
          close(20)
          write(6,1120)pointeur,ndump,utdeb,utfin
        endif
c
c --- 	calcul du Dark éventuellement
 	mode_dark= .false.
 	if (DarkPmeas.lt.0)then
 	  mode_dark= .true.
 	  DarkPmeas = 0.
 	  pointeur = 0
          do idump = 1,ndump
            if(ut(idump).ge.utdeb .and. ut(idump).le.utfin)then
 	      pointeur = pointeur+1
              DarkPmeas = DarkPmeas+intbrutP(idump)
 	    endif
 	  enddo
          DarkPmeas = DarkPmeas/float(pointeur)
 	endif
c
 	if (DarkTmeas.lt.0)then 
 	  mode_dark= .true.
 	  DarkTmeas = 0.
 	  pointeur = 0
          do idump = 1,ndump
            if(ut(idump).ge.utdeb .and. ut(idump).le.utfin)then
 	      pointeur = pointeur+1
              DarkTmeas = DarkTmeas+intbrutT(idump) 
 	    endif
 	  enddo
          DarkTmeas = DarkTmeas/float(pointeur)
 	endif
c
 	if(mode_dark)then
 	  write(6,1130) DarkPmeas,DarkTmeas
 	endif
c
 	if (idess.ne.0)then
c 	Tracé de l'intensité
c       call GR_EXEC('dev image white')
        call GR_EXEC('SET /DEFAULT')
        call GR_EXEC('SET PLOT_PAGE PORTRAIT')
        call GR_EXEC('PENCIL /DEFAULT') 
        call GR_EXEC('TICKSPACE 0 0 0 0')
        call GR_EXEC('SET FONT DU')
c
        call GR_EXEC('SET CHARACTER .8')
        call GR_EXEC('SET VIEWPORT .20 .90 .20 .90')
        call GR_EXEC('label "UT" /X')
        call GR_EXEC('label "Intensity [counts]" /Y')
        call gr4_give('X',ndump,ut)
        call mnmx(ut,ndump,utmin,utmax,0)
 	utmin = max (utmin,utdeb)
 	utmax = min (utmax,utfin)
c
        write(entete_date,1051)dirin
        call GR_EXEC(entete_date)
1051    format('DRAW TEXT -0.0 -2.3 "',a12,'" 2 0.00  /BOX 2')
c
        call GR_EXEC('pen 0 /weight 3')

c 	Canal P
        call GR_EXEC('SET VIEWPORT .23 .90 .31 .67')
        call GR_EXEC('SET CHARACTER .7')
        call GR_EXEC('label "Polarized" /Y')
        call GR_EXEC('SET CHARACTER .5')
c
        call GR_EXEC('pen 0 /weight 3')
 	jdump = 0
 	do idump=1,ndump
 	  if(ut(idump).ge.utmin .and. ut(idump).le.utmax)then
 	    jdump = jdump + 1
 	    workin(jdump)=intbrutP(idump)
 	  endif
 	enddo
        call mnmx(workin,jdump,ymin,ymax,0)
        call gr_limi(4,utmin,utmax,ymin,ymax)
        call GR_EXEC('SET VIEWPORT .22 .90 .55 .67')
        call GR_EXEC('LABEL "Raw" /Y')
        call GR_EXEC('SET VIEWPORT .20 .90 .55 .67')
        call GR_EXEC('BOX N')
        call gr4_give('Y',ndump,intbrutP)
        call GR_EXEC('pen 2 /weight 3')
        call GR_EXEC('connect')
c
 	jdump = 0
 	do idump=1,ndump
 	  if(ut(idump).ge.utmin .and. ut(idump).le.utmax)then
 	    jdump = jdump + 1
 	    workin(jdump)=intcorrP(idump)
 	  endif
 	enddo
        call mnmx(workin,jdump,ymin,ymax,0)
        call GR_EXEC('pen 0 /weight 3')
        call gr_limi(4,utmin,utmax,ymin,ymax)
        call GR_EXEC('SET VIEWPORT .22 .90 .43 .55')
        call GR_EXEC('LABEL "Corrected" /Y')
        call GR_EXEC('SET VIEWPORT .20 .90 .43 .55')
        call GR_EXEC('BOX N')
        call gr4_give('Y',ndump,intcorrP)
        call GR_EXEC('pen 2 /weight 3')
        call GR_EXEC('connect')

        call GR_EXEC('pen 0 /weight 3')
 	jdump = 0
 	do idump=1,ndump
 	  if(ut(idump).ge.utmin .and. ut(idump).le.utmax)then
 	    jdump = jdump + 1
 	    workin(jdump)=intnormP(idump)
 	  endif
 	enddo
        call mnmx(workin,jdump,ymin,ymax,0)
        call gr_limi(4,utmin,utmax,ymin,ymax)
        call GR_EXEC('SET VIEWPORT .22 .90 .31 .43')
        call GR_EXEC('LABEL "Normalized" /Y')
        call GR_EXEC('SET VIEWPORT .245 .90 .31 .43')
        call GR_EXEC('LABEL "to Ref" /Y')
        call GR_EXEC('SET VIEWPORT .20 .90 .31 .43')
        call GR_EXEC('BOX N')
        call gr4_give('Y',ndump,intnormP)
        call GR_EXEC('pen 2 /weight 3')
        call GR_EXEC('connect')
c
c      Moyenne glissante
c
        call GR_EXEC('pen 0 /weight 3')
 	jdump = 0
 	do idump=1,ndump
 	  if(ut(idump).ge.utmin .and. ut(idump).le.utmax)then
 	    jdump = jdump + 1
 	    workin(jdump)=intmoy(idump)
 	  endif
 	enddo
        call mnmx(workin,jdump,ymin,ymax,0)
        call gr_limi(4,utmin,utmax,ymin,ymax)
        call GR_EXEC('SET VIEWPORT .22 .90 .18 .30')
        call GR_EXEC('LABEL "Final working" /Y')
        call GR_EXEC('SET VIEWPORT .245 .90 .18 .30')
        call GR_EXEC('LABEL "intensity" /Y')
        call GR_EXEC('SET VIEWPORT .20 .90 .18 .30')
        call GR_EXEC('BOX')
        call gr4_give('Y',ndump,intmoy)
        call GR_EXEC('pen 0 /weight 3')
        call GR_EXEC('connect')
        call GR_EXEC('pen 0 /weight 3')
c

c 	Canal T
        call GR_EXEC('SET VIEWPORT .23 .90 .68 .92')
        call GR_EXEC('SET CHARACTER .7')
        call GR_EXEC('label "Reference" /Y')
        call GR_EXEC('SET CHARACTER .5')
c
        call GR_EXEC('pen 0 /weight 3')
 	jdump = 0
 	do idump=1,ndump
 	  if(ut(idump).ge.utmin .and. ut(idump).le.utmax)then
 	    jdump = jdump + 1
 	    workin(jdump)=intbrutT(idump)
 	  endif
 	enddo
        call mnmx(workin,jdump,ymin,ymax,0)
        call gr_limi(4,utmin,utmax,ymin,ymax)
        call GR_EXEC('SET VIEWPORT .22 .90 .80 .92')
        call GR_EXEC('LABEL "raw" /Y')
        call GR_EXEC('SET VIEWPORT .20 .90 .80 .92')
        call GR_EXEC('BOX N')
        call gr4_give('Y',ndump,intbrutT)
        call GR_EXEC('pen 3 /weight 3 /dashed 1')
        call GR_EXEC('connect')

        call GR_EXEC('pen 0 /weight 3 /dashed 1')
 	jdump = 0
 	do idump=1,ndump
 	  if(ut(idump).ge.utmin .and. ut(idump).le.utmax)then
 	    jdump = jdump + 1
 	    workin(jdump)=intcorrT(idump)
 	  endif
 	enddo
        call mnmx(workin,jdump,ymin,ymax,0)
        call gr_limi(4,utmin,utmax,ymin,ymax)
        call GR_EXEC('SET VIEWPORT .22 .90 .68 .80')
        call GR_EXEC('LABEL "Corrected" /Y')
        call GR_EXEC('SET VIEWPORT .20 .90 .68 .80')
        call GR_EXEC('BOX N')
        call gr4_give('Y',ndump,intcorrT)
        call GR_EXEC('pen 3 /weight 3 /dashed 1')
        call GR_EXEC('connect')
        call GR_EXEC('pen 0 /weight 3')

        write(entete,1053)dirin
        call GR_EXEC(entete)
1053    format
     .   ('HARDCOPY "../',a12,'/intensity" /dev eps color /OVERWRITE')
c
        if(idess.eq.2)call gmaster_loop('')
 	endif			! sur idess
c
 	if(mode_dark)then
 	  stop
 	endif
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine asz(day,UT,alat,along,hrloc,dec,chi,chideg)
c
        real hrloc,chideg,day,chi,alat
c
c       ENTREES 
c       day = numero du jour dans l'annee
c       UT = heure universelle decimale (10h30 = 10.5)
c       alat = latitude en degres
c       along = longitude en degres
c
c       SORTIES
c       hrloc = heure locale decimale (10h30 = 10.5)
c       dec= equatorial declination of the sun
c       chi = angle solaire zenithal en radian
c       chideg = angle solaire zenithal en degres
c
        pi= 3.141593
        pio180 = 3.141593/180.0
c
        hrloc = amod(24.+UT+along/15.,24.)
c
        dec = 23.5 * sin(2.0*pi*(day-80.0)/365.0)
        clat = cos ( pio180 * alat)
        slat = sin ( pio180 * alat)
        cdec = cos ( pio180 * dec)
        sdec = sin ( pio180 * dec)
        a = slat * sdec 
        b = clat * cdec 
        cx = a - b * cos(pio180*hrloc*15.0)
        chi = acos(cx)
        chideg=chi/pio180
c
        return
        end
c       
