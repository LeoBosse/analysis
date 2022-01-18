c
c----------------------------------------------------------------------
c
 	subroutine angle_apparent(ifile,an,jour,Lat_SPP,Long_SPP,Elev_SPP,
     .       Azimut_SPP,eta,angle_B_los_deg,Bup_P,Bnorth_P,Beast_P,AH,
     .       Babs_P,d_SPP_obs,ichamp,alpha,Lat_P,Long_P)
c
 	implicit none
c
c --- Inputs
c 	latitude du point où est l'instrument (degres)
 	real Lat_SPP
C 	Longitude du point où est l'instrument (degres)
 	real Long_SPP
c 	Elevation (degres) et azimut (degres) de l'instrument. 
c 	Azimut: 0 N, 90 E, 180 S, 270 W
 	real Elev_SPP, Azimut_SPP
c 	date : decimal year (year+month/12.0-0.5 or
c                  year+day-of-year/365 or ../366 if leap year)c
c 	Better between 1900 and 2025 (otherwise, a warning is given)
 	real date
 	double precision date_dd,alt_dd,Lat_P_dd,Long_P_dd
c 	numéro du fichier de sortie
 	integer ifile
 	real an,jour 
 	integer ichamp
c --- Outputs
c 	angle apparent sur SPP du champ magnétique au point où se fait 
c 	l'emission de la raie rouge
 	real eta
c 	Distance from SPP to the point of observation H
 	real d_SPP_obs
c 	Angle between line of sight and B where the red line is emitted
 	real angle_B_los_rad,angle_B_los_deg
 	real cos_angle_B_los_rad
 	real alpha 	! angle entre B et Bup
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
 	data PH /220./
c	data PH /300./
c 	Vector OP between the center of the Earth and P
 	real OP(3),normOP 		! normOP doit être egal à Rt
c 	Vector AH between SPP and P
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
 	real Babs_ref
 	real Bnorth_P,Beast_P,Bup_P,Babs_P,dimo
c 	B proj on SPP
 	real Bnorth,Beast,Bup
        integer decl,declmn,sv_decl,incl,inclmn,sv_incl,
     .       HorizB,sv_HorizB,Bx,sv_Bx,By,sv_By,Bz,sv_Bz,Btot,sv_Btot
c
 	real pi
 	real a, b, c, delta
 	real Elev_SPP_rad,Azimut_SPP_rad,Lat_SPP_rad,Long_SPP_rad
 	integer icasLat_,icasLong_
 	real facteur
c
 	if(Elev_SPP.le.1. .or. Elev_SPP.ge.179.)then
 	  write(6,*)'elevation must be > 1°'
 	  stop
 	endif
c
 	date= an+jour/366.
 	Babs_ref = sqrt(49694.**2 + 6360.**2) 	! sqrt(Bup**2 + Bnorth**2) pour ichamp=0
 	
 	write(6,1010)Lat_SPP,Long_SPP,Elev_SPP,Azimut_SPP
 	write(ifile,1010)Lat_SPP,Long_SPP,Elev_SPP,Azimut_SPP
 	pi = 4.*atan(1.)
 	Lat_SPP_rad = Lat_SPP * pi/180.		! Latitude
 	Long_SPP_rad = Long_SPP * pi/180.		! Longitude
c	on se place dans un repere avec 0 au sud
 	Azimut_int = Azimut_SPP + 180.
 	if(Azimut_int.ge.360.) Azimut_int = Azimut_int-360.
 	Azimut_SPP_rad = Azimut_int * pi/180.
 	Elev_SPP_rad = Elev_SPP * pi/180.
c 	Computes the distance d_SPP_obs between SPP and the point 
c 	H of observation
 	a = 1.
 	b = 2.*Rt*sin(Elev_SPP_rad)
 	c = -PH*(PH+2.*Rt)
	delta = b*b - 4.*a*c
 	d_SPP_obs = (-b+sqrt(delta))/2.
c
 	facteur = Rt**2 + d_SPP_obs**2 +2.*Rt*d_SPP_obs*sin(Elev_SPP_rad)
 	facteur = Rt/sqrt(facteur)
 	OP(x)=facteur*(Rt+d_SPP_obs*sin(Elev_SPP_rad))
 	OP(y)=facteur*(-d_SPP_obs*sin(Azimut_SPP_rad)*cos(Elev_SPP_rad))
 	OP(z)=facteur*(-d_SPP_obs*cos(Azimut_SPP_rad)*cos(Elev_SPP_rad))
c 
 	normOP  = sqrt(OP(x)**2+OP(y)**2+OP(z)**2)
c	
 	sinLat_P  = OP(z)/normOP
 	cosLat_P  = sqrt(OP(x)**2+OP(y)**2)/normOP
 	sinLong_P = OP(y)/sqrt(OP(y)**2+OP(x)**2)
 	cosLong_P = OP(x)/sqrt(OP(y)**2+OP(x)**2)
c
 	Lat_P_rad = atan2(sinLat_P,cosLat_P) + Lat_SPP_rad
 	Lat_P = Lat_P_rad*180./pi
 	Long_P_rad = atan2(sinLong_P,cosLong_P) + Long_SPP_rad
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

 	Beast_P = float(By)
 	Bup_P = -float(Bz)+float(ichamp)	! IGRF is positive downward
c	Bnorth_P= float(Bx)
 	Bnorth_P= sqrt(Babs_ref**2 - Bup_P**2) 
 	Babs_P = float(Btot)
 	alpha = 180./3.14159*atan(Bup_P/Bnorth_P)
c 
c 	Now, projection on SPP:
 	Bnorth= Bnorth_P		!Initialization
 	Beast = Beast_P			!Initialization
 	Bup   = Bup_P			!Initialization

  	S1(x,x) = cos(Long_SPP_rad-Long_P_rad)
 	S1(x,y)	= sin(Long_SPP_rad-Long_P_rad)
 	S1(x,z) = 0.
  	S1(y,x) = -sin(Long_SPP_rad-Long_P_rad)
 	S1(y,y)	= cos(Long_SPP_rad-Long_P_rad)
 	S1(y,z) = 0.
  	S1(z,x) = 0.
 	S1(z,y)	= 0.
 	S1(z,z) = 1.
c
c 	rotation autour de (-sin Long_,cos Long_,0)
  	S2(x,x) = sin(Long_SPP_rad)**2+
     .            cos(Long_SPP_rad)**2*cos(Lat_SPP_rad-Lat_P_rad)
 	S2(x,y)	= -sin(Long_SPP_rad)*cos(Long_SPP_rad)*
     .            (1.-cos(Lat_SPP_rad-Lat_P_rad))
 	S2(x,z) = cos(Long_SPP_rad)*sin(Lat_SPP_rad-Lat_P_rad)
  	S2(y,x) = -sin(Long_SPP_rad)*cos(Long_SPP_rad)*
     .            (1.-cos(Lat_SPP_rad-Lat_P_rad))
 	S2(y,y)	= cos(Long_SPP_rad)**2+
     .            sin(Long_SPP_rad)**2*cos(Lat_SPP_rad-Lat_P_rad)
 	S2(y,z) = sin(Long_SPP_rad)*sin(Lat_SPP_rad-Lat_P_rad)
  	S2(z,x) = -cos(Long_SPP_rad)*sin(Lat_SPP_rad-Lat_P_rad)
 	S2(z,y)	= -sin(Long_SPP_rad)*sin(Lat_SPP_rad-Lat_P_rad)
 	S2(z,z) = cos(Lat_SPP_rad-Lat_P_rad)
c
c 	Rotation de SPP pour aller dans le même plan A: autour de 
c 	pi-Azimut_SPP_rad
  	S3(x,x) = 1.
 	S3(x,y)	= 0.
 	S3(x,z) = 0.
  	S3(y,x) = 0.
 	S3(y,y)	= cos(pi-Azimut_SPP_rad)
 	S3(y,z) = sin(pi-Azimut_SPP_rad)
  	S3(z,x) = 0.
 	S3(z,y)	= -sin(pi-Azimut_SPP_rad)
 	S3(z,z) = cos(pi-Azimut_SPP_rad)
c
c 	rotation autour de (0, sin Azimut_SPP_rad, -cos Azimut_SPP_rad)
 	S4(x,x) = cos(Elev_SPP_rad)
 	S4(x,y) = cos(Azimut_SPP_rad)*sin(Elev_SPP_rad)
 	S4(x,z)	= sin(Azimut_SPP_rad)*sin(Elev_SPP_rad)
  	S4(y,x) = -cos(Azimut_SPP_rad)*sin(Elev_SPP_rad)
  	S4(y,y) = sin(Azimut_SPP_rad)**2+
     .            cos(Azimut_SPP_rad)**2*cos(Elev_SPP_rad)
 	S4(y,z)	= -sin(Azimut_SPP_rad)*cos(Azimut_SPP_rad)*
     .            (1.-cos(Elev_SPP_rad))
 	S4(z,x) = -sin(Azimut_SPP_rad)*sin(Elev_SPP_rad)
  	S4(z,y) = -sin(Azimut_SPP_rad)*cos(Azimut_SPP_rad)*
     .            (1.-cos(Elev_SPP_rad))
 	S4(z,z)	= cos(Azimut_SPP_rad)**2+
     .            sin(Azimut_SPP_rad)**2*cos(Elev_SPP_rad)
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
 	sineta = Beast/normeta
 	coseta = Bnorth/normeta
 	etarad = atan2(sineta,coseta)
 	etarad = modulo(etarad,pi)
	eta = etarad*180./pi
c
c 	Calcul de l'angle entre la ligne de visee et le champ 
c 	magnétique en P
 	AH(Up) = d_SPP_obs*sin(Elev_SPP_rad)
 	AH(North) = d_SPP_obs*cos(Elev_SPP_rad)*sin(Azimut_SPP_rad)
 	AH(East) = d_SPP_obs*cos(Elev_SPP_rad)*cos(Azimut_SPP_rad)
 	cos_angle_B_los_rad=
     .      (Bup_P*AH(Up)+Bnorth_P*AH(North)+Beast_P*AH(East))/
     .      (Babs_P*d_SPP_obs)
 	angle_B_los_rad=acos(cos_angle_B_los_rad)
 	angle_B_los_deg = angle_B_los_rad *180./pi
 	if(angle_B_los_deg.gt.90.) 
     .     angle_B_los_deg=180.-angle_B_los_deg
c
 	write(6,1040) Lat_P,Long_P
 	write(6,1020)d_SPP_obs
 	write(6,1050)eta
 	write(6,1060) angle_B_los_deg
 	write(ifile,1040) Lat_P,Long_P
 	write(ifile,1020)d_SPP_obs
 	write(ifile,1050)eta
 	write(ifile,1060) angle_B_los_deg
c
1000 	format('Angle apparent du champ',1f10.2,' degres')
1010    format('Lat SPP =',1f8.2,4x,'Long SPP =',1f8.2,/,
     .         'Elevation SPP =',1f8.2,'  Azimut SPP =',1f8.2)
1020    format('Distance SPP - obs. point',1f8.2,' km')
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
1050    format('--> Apparent angle of B on SPP = '1f8.2)
1060    format('--> Angle between B and line of sight in P = '1f8.2)
c
	return
 	end
