#!/bin/bash


folder="ptcu/20190307_Skibotn/"
./vendange.py $folder"br_a164_e45" -s False
./vendange.py $folder"mr_a164_e45"  -s False
./vendange.py $folder"rr_a164_e90"  -s False
./vendange.py $folder"vr_a164_e45"  -s False
./vendange.py $folder"br_rot_e45"  -s False
./vendange.py $folder"mr_rot_e45"  -s False

#$folder"vr_a164_e90"  -s False
#$folder"vr_a164_e90_densite"  -s False
#$folder"vr_a164_e90_densite_suite"  -s False
#$folder"vr_a164_e90_suite"  -s False

folder="ptcu/20190308_Skibotn/"
./vendange.py $folder"vr_a164_e45"  -s False

folder="ptcu/20190309_Skibotn/"
./vendange.py $folder"br_rot_e45"  -s False

folder="spp/20190303_NyAlesund/"
./vendange.py $folder"rot_lente_e30.in" -s False
