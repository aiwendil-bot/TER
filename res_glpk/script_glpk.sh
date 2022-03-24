#!/bin/bash

for file in ../modeles/*.lp
do  
  s=${file##*/}
  glpsol --lp $file --log "${s%.lp}.log"
done
