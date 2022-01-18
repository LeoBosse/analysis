#!/bin/bash

folder="ptcu/20191021_Sob/"
../vendange/main.py $folder"mm_a225_e45" -s False
../vendange/main.py $folder"mm_rot_e45" -s False
../vendange/main.py $folder"mm_rot_a00" -s False
../vendange/main.py $folder"mm_rot_a90" -s False

folder="ptcu/20191022_Sob/"
../vendange/main.py $folder"or_a00_e45" -s False
../vendange/main.py $folder"or_a90_e45" -s False
../vendange/main.py $folder"or_a180_e45" -s False
../vendange/main.py $folder"or_a180_e90" -s False
../vendange/main.py $folder"or_rot_a00" -s False
../vendange/main.py $folder"or_rot_a90" -s False
../vendange/main.py $folder"or_rot_e45" -s False
../vendange/main.py $folder"ro_a00_e45" -s False
../vendange/main.py $folder"ro_rot_a00" -s False
../vendange/main.py $folder"ro_rot_a90" -s False
../vendange/main.py $folder"ro_rot_e45" -s False

folder="ptcu/20191023_Sob/"
../vendange/main.py $folder"bo_a225_e45" -s False
../vendange/main.py $folder"bo_rot_a00" -s False
../vendange/main.py $folder"bo_rot_a90" -s False
../vendange/main.py $folder"bo_rot_e45" -s False
../vendange/main.py $folder"vo_a225_e45" -s False
../vendange/main.py $folder"vo_rot_a00" -s False
../vendange/main.py $folder"vo_rot_a90" -s False
../vendange/main.py $folder"vo_rot_e45" -s False
