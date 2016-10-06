#!/bin/bash 
# Test the sensitivity of C, B, and A factors on the disaggregation
opt=$1  # pre, or ga 
if [ $opt == 'pre' ]; then
    # preparation (download CyberShake database, compute NGA flatinfo and directivity factors)
    for mflag in 54
    do 
	./ABFanalysis.py Preparation $mflag 1 1
    done 
fi 

if [ $opt == 'ga' ]; then
    # Gkxmfs analysis:
    for mflag in 54
    do 
	./ABFanalysis.py GkxmfsAnalysis $mflag 0 0 
    done 
fi 

if [ $opt == 'bmap' ]; then 
    # GMT plots
    cd ./GkxmfsAnalysis 
    for mflag in 54
    do 
	./model_Bs_SiteModel.gmt 35 8 4 5 3.00 $mflag Ref 0 1.00_1.00 OnlyBs
	./model_Bs_periods.gmt 35 6 4 7 $mflag CyberShake 0 1.00_1.00
	./model_ks.gmt 35 5 3 1 3.00 BA Cks $mflag CyberShake 1.00_1.00 0   # 0 means plot all sources k(s), where Cks can be variance map (directivity)

    done 
fi 
