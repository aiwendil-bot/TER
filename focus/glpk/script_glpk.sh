#!/bin/bash

for file in ../*.lp
do  
  s=${file##*/}
  glpsol --lp $file --log "${s%.lp}.log" -o "${s%.lp}.txt"
done
