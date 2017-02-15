#!/usr/bin/env python 

import os 
import numpy as np 
import matplotlib.pyplot as plt 

modelDict = {
	'CS11':(35,5,3,1),
	'CS13.4a':(35,7,4,1),
	'CS13.4b':(35,7,4,4),
	'CS14.2a':(35,8,4,5), 
	'CS14.2b':(35,6,4,7)
	}

modelKeys = ['CS11','CS13.4a','CS13.4b','CS14.2a','CS14.2b'] 
clrs = ['r','k','b','g','m'] 

outpth = '/Users/fengw/work/Project/ABFanalysis/scripts/CompareCyberShakes/outputs'
outFile = os.path.join(outpth, 'afactor_summary.dat') 
fidout = open(outFile, 'w') 

fig = plt.figure(1)
fig.clf() 
ax = fig.add_subplot(111) 

if 1:
    NGA08Data = {}
    NGA14Data = {}
    outData = {}
    for imodel in range(len(modelKeys)): 
	model = modelKeys[imodel]
	erfID, sgtID, rupVarID, velID = modelDict[model] 
	fidout.write('#%s; %s, %s, %s, %s\n'%(model, erfID, sgtID, rupVarID, velID)) 
	
	outData[model] = []
	NGA08Data[model] = []
	NGA14Data[model] = []
	print 'processing %s'%model 

	if model in ['CS14.2a','CS14.2b']: 
	    inPath = '/Users/fengw/work/Project/ABFanalysis/scripts/map_input/Model_Rups54' 
	else: 
	    inPath = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/Model_Rups54' 
	prefix = 'ERF%s_SGT%s_RupVar%s_Vel%s'%(erfID, sgtID, rupVarID, velID) 
	inFullPath = os.path.join(inPath, prefix, 'Gkxmfs0/A/Sigma1.00_1.00') 

	for T in ['2.00','3.00','5.00','10.00']: 
	    inFile0 = 'CyberShake.NGAs.%s.a'%T 
	    fidin = open(os.path.join(inFullPath, inFile0), 'r') 
	    line = fidin.readline() 
	    if '#' in line: 
		line = fidin.readline() 
	    spl = line.strip().split()
	    str1 = []
	    for i in range(len(spl)): 
		if i == 4:
		    outData[model].append(float(spl[4]))
		    str1.append(spl[4])
		else: 
		    if model in ['CS11','CS13.4a','CS13.4b']: 
			NGA08Data[model].append(float(spl[4])-float(spl[i])) 
	            if model in ['CS14.2a','CS14.2b']: 
			NGA14Data[model].append(float(spl[4])-float(spl[i])) 
		    str1.append(str(float(spl[4])-float(spl[i])))
	    fidout.write(','.join(tuple(str1))+'\n') 
	    fidin.close() 
	ax.loglog([2,3,5,10],np.exp(outData[model]),clrs[imodel]+'-',marker='o',lw=1.5, label=model) 
	#ax1.loglog([2,3,5,10],np.exp(outData[model]),clrs[imodel]+'-',marker='o',label=model) 

    fidout.close()

    from matplotlib.ticker import LogLocator, FormatStrFormatter
    minorLocator = LogLocator(subs = np.linspace(2,10,8,endpoint=False))
    ax.xaxis.set_minor_locator( minorLocator)

    ax.set_xlim([1.8,11]) 
    ax.set_ylim([0.005,0.1])
    #ax1.set_ylim([0.005,0.1])

    ax.grid(b=True,which='minor') 
    ax.set_xlabel('Periods (s)')
    #ax1.set_ylabel('ln A') 
    ax.set_ylabel('A factor (g)')
    ax.legend(loc=0)
    ax.yaxis.set_ticks_position('left') 

    #subs = [-3, -4, -5]

    #majorLocator = LogLocator(subs = np.exp(subs) )
    #majorLocator = LogLocator(subs = subs )
    #majorFormatter = FormatStrFormatter('$e^{%.2f}$')
    #ax1.yaxis.set_major_locator( majorLocator )
    #ax1.yaxis.set_major_formatter(majorFormatter)
    #ax1.yaxis.set_ticks_position('left')
    fig.savefig(os.path.join(outpth, 'afactor_summary.png'),dpi=600)
    plt.show() 

if 0: 
    # plot residual a factor
    import sys 
    pathname = os.path.dirname(sys.argv[0]) 
    inFile = os.path.join(pathname, 'afactor_summary.dat')
    outpth = os.path.join(pathname, 'outputs') 

    fid = open(inFile, 'r')
    line = fid.readline() 
    outData = {}
    while line: 
        if '#' in line: 
            spl = line.strip().split(';')[0][1:]
            outData[spl] = []
            for i in range(4): 
                line = fid.readline()
                spl1 = line.strip().split(',')
                ngaAve = 0.25*(sum([float(tmp) for tmp in spl1[:4]]))
                outData[spl].append(float(spl1[4])-ngaAve)  # residual a value (ln)
        line = fid.readline() 
    fid.close()
    
    # plot 
    for imodel in range(len(modelKeys)): 
        modelName = modelKeys[imodel] 
        tmpData = outData[modelName] 
        if modelName in ['CS11','CS13.4a','CS13.4b']: 
            NGAtype = 'NGA08ave'
        else: 
            NGAtype = 'NGA14ave'
        #ax.loglog([2,3,5,10],np.exp(outData[modelName]),clrs[imodel]+'-',marker='o',lw=1.5, label=modelName+'-'+NGAtype)
        ax.plot([2,3,5,10],outData[modelName],clrs[imodel]+'-',marker='o',lw=1.5, label=modelName+'-'+NGAtype)
        
    ax.grid(b=True,which='minor') 
    ax.grid(b=True,which='major')
    ax.set_xlabel('Periods (s)')
    ax.set_xlim([1.8,11])
    ax.set_ylabel('residual a factor')
    ax.legend(loc=0)    
    fig.savefig(os.path.join(pathname,'residualAfactor.png'),dpi=600)
