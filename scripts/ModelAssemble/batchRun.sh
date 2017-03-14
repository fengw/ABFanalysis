#!/bin/bash
opt=$1    # 1 is to calculate, and 2 is plot

wrkpth=/Users/fengw/work/Project/ABFanalysis/scripts/ModelAssemble

if [ $opt == 1 ]; then
    for T in 2.00 3.00 5.00 10.00
    do 
        for Gflag in 0 1
        do 
    	python $wrkpth/CalculateSigma.py -m 54 -r $2 $3 $4 $5 -p $T -S 1.00_1.00 -G $Gflag
        done
    done
else 
    for model in CS11 CS13.1 CS13.2 CSbbp1D CS14.2CVMH CS14S4.26 CS15.4
    do 
	python $wrkpth/StackPlots.py singleModel $model 
    done 
    python $wrkpth/StackPlots.py allModel    # for all comparision
fi 

