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

if [ $opt == 'bcmap' ]; then 
    # GMT plots
    cd ./GkxmfsAnalysis 
    for mflag in 54
    do 
	./model_Bs_SiteModel.gmt 35 5 3 1 3.00 $mflag BA All LocalModel 0 0.60 
	./model_ks.gmt 35 5 3 1 3.00 BA Cks $mflag CyberShake 0 0.60 
    done 
fi 
