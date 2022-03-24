#!/bin/bash

for file in ../modeles/*.lp
do  
  s=${file##*/}
  scip -f $file -q -l "${s%.lp}.log"
done