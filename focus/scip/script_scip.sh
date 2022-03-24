#!/bin/bash

for file in ../*.lp
do  
  s=${file##*/}
  scip -f $file -q -l "${s%.lp}.log"
done