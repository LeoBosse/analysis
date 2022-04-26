#!/bin/bash

cd ~/iasb/analysis/data

pwd

for f in `ls -d ./ptcu/202112*_Corbel/`:
do
  echo $f"vm_a35_e45_5hz"
  cd ~/iasb/analysis/src/vendange/

  ./main.py $f"/vm_a35_e45_5hz" -s False

  cd ~/iasb/analysis/data

done
