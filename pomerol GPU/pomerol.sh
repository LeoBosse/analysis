#!/usr/bin/sh
# -*-coding:utf-8 -*

N=$1

python3 setup.py build_ext --inplace

mpirun -np $N python3 main.py
