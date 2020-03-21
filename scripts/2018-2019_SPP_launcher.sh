#!/bin/bash

./vendange.py spp/20181113_NyAlesund/ -s False
./vendange.py spp/20181114_NyAlesund/ -s False
./vendange.py spp/20181115_NyAlesund/ -s False
./vendange.py spp/20181204_NyAlesund/ -s False
./vendange.py spp/20181119_NyAlesund/ -s False
./vendange.py spp/20181206_NyAlesund/ -s False
./vendange.py spp/20181208_NyAlesund/ -s False
./vendange.py spp/20181209_NyAlesund/ -s False
./vendange.py spp/20181207_NyAlesund/ -s False
./vendange.py spp/20181212_NyAlesund/ -s False
./vendange.py spp/20181214_NyAlesund/ -s False
./vendange.py spp/20181217_NyAlesund/ -s False
./vendange.py spp/20190128_NyAlesund/ -s False
./vendange.py spp/20190129_NyAlesund/ -s False
./vendange.py spp/20190131_NyAlesund/ -s False
./vendange.py spp/20190201_NyAlesund/ -s False
./vendange.py spp/20190202_NyAlesund/ -s False


folder="spp/20190301_NyAlesund/"
./vendange.py $folder -s False
./vendange.py $folder"a00_e38.in"  -s False
./vendange.py $folder"a180_e85.in"  -s False
./vendange.py $folder"rot_lente_e30.in"  -s False

folder="spp/20190303_NyAlesund/"
./vendange.py $folder -s False
./vendange.py $folder"a-90_e30.in" -s False
./vendange.py $folder"a180_e38.in" -s False
./vendange.py $folder"a180_e85.in" -s False
./vendange.py $folder"calib_angle_0deg.in" -s False
./vendange.py $folder"rot_lente_e30.in" -s False


folder="spp/20190307_Skibotn/"
./vendange.py $folder -s False
./vendange.py $folder"br_a164_e45.in"  -s False
./vendange.py $folder"br_a164_e90.in"  -s False
./vendange.py $folder"mr_a164_e90.in"  -s False
./vendange.py $folder"mr_a164_e45.in"  -s False
./vendange.py $folder"vr_a164_e45.in"  -s False
./vendange.py $folder"vr_a164_e90.in"  -s False
./vendange.py $folder"rr_a164_e90.in"  -s False
