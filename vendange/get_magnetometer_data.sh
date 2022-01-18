#!/bin/bash


# This script will download magnetometer data from the http://flux.phys.uit.no website and save them where you are in your terminal.

P=$1 #Tromso: tro2a / NyAlesund: nal1a

Y=$2 #Year
M=$3 #Month
D=$4 #Day of the month

PSW="ResUseNoCom" #Password used to connect to the site

DT="10sec" #delta t for the data.

wget -O $Y$M$D 'http://flux.phys.uit.no/cgi-bin/mkascii.cgi?site='$P'&year='$Y'&month='$M'&day='$D'&res='$DT'&pwd='$PSW'&format=html&comps=DHZ&getdata=+Get+Data+'
