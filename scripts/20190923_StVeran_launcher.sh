#!/bin/bash

folder="ptcu/20190923_stveran/"
../vendange/main.py $folder"vm_a68_e36" -s False

folder="ptcu/20190925_stveran/"
../vendange/main.py $folder"bm_a255_e44" -s False
../vendange/main.py $folder"mm_a255_e44" -s False
../vendange/main.py $folder"vm_a255_e44" -s False
../vendange/main.py $folder"rm_a255_e44" -s False
../vendange/main.py $folder"bm_rot_e44" -s False
../vendange/main.py $folder"mm_rot_e44" -s False
../vendange/main.py $folder"vm_rot_e44" -s False
../vendange/main.py $folder"rm_rot_e44" -s False
../vendange/main.py $folder"vm_west" -s False

folder="ptcu/20190926_stveran/"
../vendange/main.py $folder"rm_a00_e45" -s False
../vendange/main.py $folder"rm_a331_e27" -s False
../vendange/main.py $folder"vm_a00_e45" -s False
../vendange/main.py $folder"vm_a331_e27" -s False
../vendange/main.py $folder"vm_a338_e22" -s False
