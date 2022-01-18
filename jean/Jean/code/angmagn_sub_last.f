c
c----------------------------------------------------------------------
c
 	subroutine angle_apparent(ifile,an,jour,Lat_instru,Long_instru,
     .      Elev_instru,Azimut_instru,eta,angle_B_los_deg,Lat_P,
     .      Long_P,d_instru_obs,PH)
c
 	implicit none
c
c --- Inputs
c 	latitude du point où est l'instrument (degres)
 	real Lat_instru
C 	Longitude du point où est l'instrument (degres)
 	real Long_instru
c 	Elevation (degres) et azimut (degres) de l'instrument. 
c 	Azimut: 0 N, 90 E, 180 S, 270 W
 	real Elev_instru, Azimut_instru
c 	date : decimal year (year+month/12.0-0.5 or
c                  year+day-of-year/365 or ../366 if leap year)c
c 	Better between 1900 and 2025 (otherwise, a warning is given)
 	real date
 	double precision date_dd,alt_dd,Lat_P_dd,Long_P_dd
c 	numéro du fichier de sortie
 	integer ifile
 	real an,jour 
c --- Outputs
c 	angle apparent sur instru du champ magnétique au point où se fait 
c 	l'emission de la raie rouge
 	real eta
c 	Distance from instru to the point of observation H
 	real d_instru_obs
c 	Angle between line of sight and B where the red line is emitted
 	real angle_B_los_rad,angle_B_los_deg
 	real cos_angle_B_los_rad
c
c --- Internal variables
c 	Azimut interne: 0 S, 90 W, 180 N, 270 E
 	real Azimut_int
c 	Les directions des axes...
c 	x = altitude up
c 	y = East
c 	z = North
 	integer x,y,z
 	data x /1/,y /2/,z /3/
 	integer Up,East,North
 	data Up /1/,East /2/,North /3/
c 	radius of the Earth (km)
 	real Rt
 	data Rt /6370./
c 	Altitude of the red line max intensity (km)
 	real PH
c	data PH /220./
c 	Vector OP between the center of the Earth and P
 	real OP(3),normOP 		! normOP doit être egal à Rt
c 	Vector AH between instru and P
 	real AH(3)
 	real normeta,etarad,sineta,coseta
c	rotation matrixes
 	real S1(3,3),S2(3,3)
 	real S3(3,3),S4(3,3)
c	latitude Lat_P et longitude Long_P au point P et leurs sin/cos
 	real Lat_P_rad,Long_P_rad
 	real Lat_P,Long_P		! en degres
 	real sinLat_P,cosLat_P,sinLong_P,cosLong_P
c 	Valeurs du champ:
C       BABS   MAGNETIC FIELD STRENGTH IN GAUSS
C       BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
C       TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
C       POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
C       AND DOWNWARD. 
 	real Bnorth_H,Beast_H,Bup_H,Babs_H,dimo
c 	B proj on instru
 	real Bnorth,Beast,Bup
        integer decl,declmn,sv_decl,incl,inclmn,sv_incl,
     .       HorizB,sv_HorizB,Bx,sv_Bx,By,sv_By,Bz,sv_Bz,Btot,sv_Btot
c
 	real pi
 	real a, b, c, delta
 	real Elev_instru_rad,Azimut_instru_rad
 	real Lat_instru_rad,Long_instru_rad
 	integer icasLat_,icasLong_
 	real facteur
c
 	date= an+jour/366.
 	
c	write(6,*)
c	write(6,1010)Lat_instru,Long_instru,Elev_instru,Azimut_instru
 	write(ifile,1010)Lat_instru,Long_instru,Elev_instru,Azimut_instru
 	pi = 4.*atan(1.)
 	Lat_instru_rad = Lat_instru * pi/180.		! Latitude
 	Long_instru_rad = Long_instru * pi/180.		! Longitude
c	on se place dans un repere avec 0 au sud
 	Azimut_int = Azimut_instru + 180.
 	if(Azimut_int.ge.360.) Azimut_int = Azimut_int-360.
 	Azimut_instru_rad = Azimut_int * pi/180.
 	Elev_instru_rad = Elev_instru * pi/180.
c 	Computes the distance d_instru_obs between instru and the point 
c 	H of observation
 	a = 1.
 	b = 2.*Rt*sin(Elev_instru_rad)
 	c = -PH*(PH+2.*Rt)
	delta = b*b - 4.*a*c
 	d_instru_obs = (-b+sqrt(delta))/2.
c
 	facteur = Rt**2 + 
     . 		d_instru_obs**2 +2.*Rt*d_instru_obs*sin(Elev_instru_rad)
 	facteur = Rt/sqrt(facteur)
 	OP(x)=facteur*(Rt+d_instru_obs*sin(Elev_instru_rad))
 	OP(y)=facteur*
     . 		(-d_instru_obs*sin(Azimut_instru_rad)*cos(Elev_instru_rad))
 	OP(z)=facteur*
     . 		(-d_instru_obs*cos(Azimut_instru_rad)*cos(Elev_instru_rad))
c 
 	normOP  = sqrt(OP(x)**2+OP(y)**2+OP(z)**2)
c	
 	sinLat_P  = OP(z)/normOP
 	cosLat_P  = sqrt(OP(x)**2+OP(y)**2)/normOP
 	sinLong_P = OP(y)/sqrt(OP(y)**2+OP(x)**2)
 	cosLong_P = OP(x)/sqrt(OP(y)**2+OP(x)**2)
c
 	Lat_P_rad = atan2(sinLat_P,cosLat_P) + Lat_instru_rad
 	Lat_P = Lat_P_rad*180./pi
 	Long_P_rad = atan2(sinLong_P,cosLong_P) + Long_instru_rad
 	Long_P = Long_P_rad*180./pi
c	write(6,*)'norme /OP/ = ',normOP,' (Rt=6370)'
c
 	alt_dd = dble(PH)
 	date_dd = dble(date)
 	Lat_P_dd=dble(Lat_P)
 	Long_P_dd=dble(Long_P)
 	call igrf_sub(date_dd,alt_dd,Lat_P_dd,Long_P_dd,
     .	     decl,declmn,sv_decl,incl,inclmn,sv_incl,HorizB,sv_HorizB,
     .       Bx,sv_Bx,By,sv_By,Bz,sv_Bz,Btot,sv_Btot)
c
c       write(6,1030)decl,declmn,sv_decl,incl,inclmn,sv_incl,
c    .       HorizB,sv_HorizB,Bx,sv_Bx,By,sv_By,Bz,sv_Bz,Btot,sv_Btot

 	Bnorth_H= float(Bx)
 	Beast_H = float(By)
 	Bup_H = -float(Bz)		! IGRF is positive downward
 	Babs_H = float(Btot)
c 
c 	Now, projection on instru:
 	Bnorth= Bnorth_H		!Initialization
 	Beast = Beast_H			!Initialization
 	Bup   = Bup_H			!Initialization

  	S1(x,x) = cos(Long_instru_rad-Long_P_rad)
 	S1(x,y)	= sin(Long_instru_rad-Long_P_rad)
 	S1(x,z) = 0.
  	S1(y,x) = -sin(Long_instru_rad-Long_P_rad)
 	S1(y,y)	= cos(Long_instru_rad-Long_P_rad)
 	S1(y,z) = 0.
  	S1(z,x) = 0.
 	S1(z,y)	= 0.
 	S1(z,z) = 1.
c
c 	rotation autour de (-sin Long_,cos Long_,0)
  	S2(x,x) = sin(Long_instru_rad)**2+
     .            cos(Long_instru_rad)**2*cos(Lat_instru_rad-Lat_P_rad)
 	S2(x,y)	= -sin(Long_instru_rad)*cos(Long_instru_rad)*
     .            (1.-cos(Lat_instru_rad-Lat_P_rad))
 	S2(x,z) = cos(Long_instru_rad)*sin(Lat_instru_rad-Lat_P_rad)
  	S2(y,x) = -sin(Long_instru_rad)*cos(Long_instru_rad)*
     .            (1.-cos(Lat_instru_rad-Lat_P_rad))
 	S2(y,y)	= cos(Long_instru_rad)**2+
     .            sin(Long_instru_rad)**2*cos(Lat_instru_rad-Lat_P_rad)
 	S2(y,z) = sin(Long_instru_rad)*sin(Lat_instru_rad-Lat_P_rad)
  	S2(z,x) = -cos(Long_instru_rad)*sin(Lat_instru_rad-Lat_P_rad)
 	S2(z,y)	= -sin(Long_instru_rad)*sin(Lat_instru_rad-Lat_P_rad)
 	S2(z,z) = cos(Lat_instru_rad-Lat_P_rad)
c
c 	Rotation de instru pour aller dans le même plan A: autour de 
c 	pi-Azimut_instru_rad
  	S3(x,x) = 1.
 	S3(x,y)	= 0.
 	S3(x,z) = 0.
  	S3(y,x) = 0.
 	S3(y,y)	= cos(pi-Azimut_instru_rad)
 	S3(y,z) = sin(pi-Azimut_instru_rad)
  	S3(z,x) = 0.
 	S3(z,y)	= -sin(pi-Azimut_instru_rad)
 	S3(z,z) = cos(pi-Azimut_instru_rad)
c
c 	rotation autour de (0, sin Azimut_instru_rad, -cos Azimut_instru_rad)
 	S4(x,x) = cos(Elev_instru_rad)
 	S4(x,y) = cos(Azimut_instru_rad)*sin(Elev_instru_rad)
 	S4(x,z)	= sin(Azimut_instru_rad)*sin(Elev_instru_rad)
  	S4(y,x) = -cos(Azimut_instru_rad)*sin(Elev_instru_rad)
  	S4(y,y) = sin(Azimut_instru_rad)**2+
     .            cos(Azimut_instru_rad)**2*cos(Elev_instru_rad)
 	S4(y,z)	= -sin(Azimut_instru_rad)*cos(Azimut_instru_rad)*
     .            (1.-cos(Elev_instru_rad))
 	S4(z,x) = -sin(Azimut_instru_rad)*sin(Elev_instru_rad)
  	S4(z,y) = -sin(Azimut_instru_rad)*cos(Azimut_instru_rad)*
     .            (1.-cos(Elev_instru_rad))
 	S4(z,z)	= cos(Azimut_instru_rad)**2+
     .            sin(Azimut_instru_rad)**2*cos(Elev_instru_rad)
c
c
 	Bnorth=S1(x,x)*Bnorth+S1(x,y)*Beast+S1(x,z)*Bup
 	Beast=S1(y,x)*Bnorth+S1(y,y)*Beast+S1(y,z)*Bup
 	Bup=S1(z,x)*Bnorth+S1(z,y)*Beast+S1(z,z)*Bup
c 
 	Bnorth=S2(x,x)*Bnorth+S2(x,y)*Beast+S2(x,z)*Bup
 	Beast=S2(y,x)*Bnorth+S2(y,y)*Beast+S2(y,z)*Bup
 	Bup=S2(z,x)*Bnorth+S2(z,y)*Beast+S2(z,z)*Bup
c 
 	Bnorth=S3(x,x)*Bnorth+S3(x,y)*Beast+S3(x,z)*Bup
 	Beast=S3(y,x)*Bnorth+S3(y,y)*Beast+S3(y,z)*Bup
 	Bup=S3(z,x)*Bnorth+S3(z,y)*Beast+S3(z,z)*Bup
c 
 	Bnorth=S4(x,x)*Bnorth+S4(x,y)*Beast+S4(x,z)*Bup
 	Beast=S4(y,x)*Bnorth+S4(y,y)*Beast+S4(y,z)*Bup
 	Bup=S4(z,x)*Bnorth+S4(z,y)*Beast+S4(z,z)*Bup
c 
 	normeta  = sqrt(Bnorth**2+Beast**2)
c	
c 	les deux lignes suivantes pour un repère avec 0 au nord
 	coseta = Beast/normeta
 	sineta = Bnorth/normeta
c 	les deux lignes suivantes pour un repère avec 0 au sud
c	sineta = Beast/normeta
c	coseta = Bnorth/normeta
 	etarad = atan2(sineta,coseta)
c	etarad = modulo(etarad,pi)
 	eta = etarad*180./pi
c
c 	Calcul de l'angle entre la ligne de visee et le champ 
c 	magnétique en P
 	AH(Up) = d_instru_obs*sin(Elev_instru_rad)
  	AH(East) = d_instru_obs*
     . 		cos(Elev_instru_rad)*sin(Azimut_instru_rad)
  	AH(North) = d_instru_obs*
     . 		cos(Elev_instru_rad)*cos(Azimut_instru_rad)
c
 	cos_angle_B_los_rad=
     .      (Bup_H*AH(Up)+Bnorth_H*AH(North)+Beast_H*AH(East))/
     .      (Babs_H*d_instru_obs)
c 	ci-dessous, B projeté dans le repère A.
  	cos_angle_B_los_rad=
     .      (Bup*AH(Up)+Bnorth*AH(North)+Beast*AH(East))/
     .      (Babs_H*d_instru_obs)
 	angle_B_los_rad=acos(cos_angle_B_los_rad)
 	angle_B_los_deg = angle_B_los_rad *180./pi
c	if(angle_B_los_deg.gt.90.) 
c    .     angle_B_los_deg=180.-angle_B_los_deg
c
 	write(ifile,1040) Lat_P,Long_P
 	write(ifile,1020)d_instru_obs
 	write(ifile,1050)eta
 	write(ifile,1060) angle_B_los_deg
c
1000 	format('Angle apparent du champ',1f10.2,' degres')
1010    format('Lat instru =',1f8.2,4x,'Long instru =',1f8.2,/,
     .         'Elevation =',1f8.2,'  Azimut =',1f8.2)
1020    format('Distance instru - obs. point',1f8.2,' km')
1030    format ( ' Declination (+ve east)',t35,
     .          'D =',I5,' deg',I4,' min',4X,'SV =',I4,' min/yr',/,
     .          ' Inclination (+ve down)',t35,
     .          'I =',I5,' deg',I4,' min',4X,'SV =',I4,' min/yr',/,
     .          ' Horizontal intensity',t35,
     .          'H =',I8,' nT     ',5X,'SV =',I4,' nT/yr',/,
     .          ' North component',t35,
     .          'X =',I8,' nT     ',5X,'SV =',I4,' nT/yr',/,
     .          ' east component',t35,
     .          'Y =',I8,' nT     ',5X,'SV =',I4,' nT/yr',/,
     .          ' Vertical component (+ve down)',t35,
     .          'Z =',I8,' nT     ',5X,'SV =',I4,' nT/yr',/,
     .          ' Total intensity',t35,
     .          'F =',I8,' nT     ',5X,'SV =',I4,' nT/yr',/,
     .          ' SV is secular variation (annual rate of change)')
1040 	format('Lat P = ',1f8.2,'   Long P = ',1f8.2)
1050    format('--> Apparent angle of B on instru = '1f8.2)
1060    format('--> Angle between B and line of sight in P = '1f8.2)
c
	return
 	end
