#!/bin/bash

# This script will download all-sky camera images from the http://tid.uio.no website and save them where you are in your terminal.
#Stay in this directory.

P=$1 #Skibotn: skn4 / NyAlesund: nya6

Y=$2 #Year
M=$3 #Month 01-12 (2 digits)
D=$4 #Day 01-31 (2 digits)

l=$5 #5577 or 6300 wavelength (A)


wget -r -np "http://tid.uio.no/plasma/aurora/$1/$5/$2/$2$3$4/" > log 2>&1 &
