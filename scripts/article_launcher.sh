#!/bin/bash


folder="ptcu/20190307_Skibotn/"
../vendange/vendange.py.py $folder"br_a164_e45" -s False
../vendange/vendange.py.py $folder"mr_a164_e45"  -s False
../vendange/vendange.py.py $folder"rr_a164_e90"  -s False
../vendange/vendange.py.py $folder"vr_a164_e45"  -s False
../vendange/vendange.py.py $folder"br_rot_e45"  -s False
../vendange/vendange.py.py $folder"mr_rot_e45"  -s False

#$folder"vr_a164_e90"  -s False
#$folder"vr_a164_e90_densite"  -s False
#$folder"vr_a164_e90_densite_suite"  -s False
#$folder"vr_a164_e90_suite"  -s False

folder="ptcu/20190308_Skibotn/"
../vendange/vendange.py.py $folder"vr_a164_e45"  -s False

folder="ptcu/20190309_Skibotn/"
../vendange/vendange.py.py $folder"br_rot_e45"  -s False

folder="spp/20190303_NyAlesund/"
../vendange/vendange.py.py $folder"rot_lente_e30.in" -s False
