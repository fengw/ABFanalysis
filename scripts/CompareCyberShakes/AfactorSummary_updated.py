#!/usr/bin/env python 

import os, sys
import numpy as np 
import matplotlib.pyplot as plt 
import datetime 
from pyABF import *

import matplotlib.colors as colors
import matplotlib.cm as cmx

# multiple lines plot colors
def LineColorCoding(N,cmap='jet'):
    """
    Generate colormap for N lines/points from the given colormap 
    Later used in loop for each line/point by colorVal = scalarMap.to_rgba(i) i from 0
    """
    colormap_name = cmap
    cm = plt.get_cmap(colormap_name)
    cNorm  = colors.Normalize(vmin=0, vmax=N-1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    return scalarMap

modelDict = {
	'CS-LA1.0':[35,5,3,1],   
	'CS-LA13.4.a':[35,7,4,1],  
	'CS-LA13.4.b':[35,7,4,4], 
	'CS-LA14.2.b':[35,8,4,8],
	'CS-LA14.2.c':[35,6,4,7],
	'CS-LA14.2.d':[35,8,4,5],
	'CS-LA15.4':[36,8,6,5], 
	}

modelKeys = ['CS-LA1.0','CS-LA13.4.a','CS-LA13.4.b','CS-LA14.2.b','CS-LA14.2.c','CS-LA14.2.d','CS-LA15.4']
Nmodels = len(modelKeys) 
scalarMap = LineColorCoding(Nmodels)

outpth = '/Users/fengw/work/Project/ABFanalysis/scripts/CompareCyberShakes/outputs'
pltpth = '/Users/fengw/work/Project/ABFanalysis/products/Afactor' 

if not os.path.exists(pltpth): 
    os.mkdir(pltpth) 

outFile = os.path.join(outpth, 'Model_Afactor_summary.dat')

opt = sys.argv[1]
if opt == 'Afactor': 
    # organize A factor and plot A factor of all CyberShake models

    # this file will include all models (NGA and CS) at various periods
    fidout = open(outFile, 'w') 

    fig = plt.figure(1)
    fig.clf() 
    ax = fig.add_subplot(111) 

    NGA08ave = []
    NGA14ave = []
    outData = {}
    for imodel in range(len(modelKeys)): 
	colorVal = scalarMap.to_rgba(imodel) 

	model = modelKeys[imodel]
	erfID, sgtID, rupVarID, velID = modelDict[model] 
	fidout.write('#%s; %s, %s, %s, %s\n'%(model, erfID, sgtID, rupVarID, velID)) 
	
	outData[model] = []
	print 'processing %s'%model 

	if model in ['CS-LA14.2.c','CS-LA14.2.d','CS-LA15.4', ]:
	    inPath = '/Users/fengw/work/Project/ABFanalysis/scripts/map_input/Model_Rups54' 
	else: 
	    inPath = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups54' 
	prefix = 'ERF%s_SGT%s_RupVar%s_Vel%s'%(erfID, sgtID, rupVarID, velID) 
	inFullPath = os.path.join(inPath, prefix, 'Gkxmfs0/A/Sigma1.00_1.00') 

	for T in ['1.00', '2.00','3.00','5.00','10.00']: 
	    if T=='1.00' and model!='CS-LA15.4': 
		continue 

	    inFile0 = 'CyberShake.NGAs.%s.a'%T 
	    fidin = open(os.path.join(inFullPath, inFile0), 'r') 
	    line = fidin.readline() 
	    if '#' in line: 
		line = fidin.readline() 
	    spl = line.strip().split()
	    str1 = []
	    tmp1 = []; tmp2 = []
	    for i in range(len(spl)): 
		if i == 4:
		    outData[model].append(float(spl[4]))
		    str1.append(spl[4])
		else: 
		    # NGA all based on CVMS4.26
		    if model == 'CS-LA13.4.a': 
			tmp1.append(float(spl[4])-float(spl[i])) 
		    if model == 'CS-LA15.4':
			tmp2.append(float(spl[4])-float(spl[i]))
		    
		    str1.append(str(float(spl[4])-float(spl[i])))
	    fidout.write(','.join(tuple(str1))+'\n') 
	    fidin.close() 

	    if model == 'CS-LA13.4.a': 
		NGA08ave.append(sum(np.exp(tmp1))/4.) 
	    if model == 'CS-LA15.4':
		NGA14ave.append(sum(np.exp(tmp2))/4.) 

	if model != 'CS-LA15.4':
	    ax.loglog([2,3,5,10],np.exp(outData[model]),'-',color=colorVal,marker='o',lw=2, label=model) 
	else:
	    ax.loglog([1,2,3,5,10],np.exp(outData[model]),'-',color=colorVal,marker='o',lw=2, label=model) 

    fidout.close()

    # plot NGA08 
    ax.loglog([2,3,5,10], NGA08ave,'k--',lw=1.5,label='NGA08ave') 
    ax.loglog([1,2,3,5,10], NGA14ave,'k-o',lw=1.5,label='NGA14ave') 

    from matplotlib.ticker import LogLocator, FormatStrFormatter
    minorLocator = LogLocator(subs = np.linspace(2,10,8,endpoint=False))
    ax.xaxis.set_minor_locator( minorLocator)

    ax.set_xlim([0.9,11])
    ax.set_ylim([0.01,0.2])

    ax.grid(b=True,which='minor') 
    ax.set_xlabel('Periods (s)')
    ax.set_ylabel('A factor (g)')
    ax.legend(loc=0, numpoints=1)
    ax.yaxis.set_ticks_position('left') 
    
    fig.savefig(os.path.join(pltpth, 'Models_Afactor_summary.png'),dpi=600)

if opt == 'afactor':
    # plot residual a factor with NGA models (reference) 
    inFile = outFile
    fid = open(inFile, 'r')
    line = fid.readline()  
 
    # read in the summary information and do the residual plots (among models maybe)

    fid.close() 


