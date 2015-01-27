#!/bin/bash
#This script assumes you have octave installed and main.c has been compiled

echo START OF FILE > resultsOfsimbt.txt

SLOT=4
MASTER=2.0
SLAVE=0.3

MAXRUNS=50

runs=0
echo SLOT: $SLOT MASTER: $MASTER SLAVE: $SLAVE >> resultsOfsimbt.txt
while [ $runs -lt $MAXRUNS ]; 
do
    echo Run: $runs
    echo Run: $runs >> resultsOfsimbt.txt
    octave --silent simbt.m $SLOT $MASTER $SLAVE >> 00sim
    ./driver.x output.txt >> 00driver
    let runs+=1
done

diff -w 00sim 00driver >> resultsOfsimbt.txt

rm 00sim 00driver
