#!/bin/bash
#This script assumes you have octave installed and main.c has been compiled

echo START OF FILE > resultsOfsimbt.txt

runs=0
while [ $runs -lt 20 ]; do
    echo Run: $runs
    octave simbt.m >> resultsOfsimbt.txt
    ./driver.x output.txt
    let runs+=1
done
