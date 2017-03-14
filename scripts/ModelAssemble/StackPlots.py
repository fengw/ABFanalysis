#!/usr/bin/env python 
import os, sys 

import numpy as np 
import matplotlib.pyplot as plt

opt = sys.argv[1]   # singalModel 

mappedNameDict = {
	'CS11': 'CS-LA1.0', 
	'CS13.1': 'CS-LA13.4a',
	'CS13.2': 'CS-LA13.4b',
	'CSbbp1D': 'CS-LA14.2.b',
	'CS14.2CVMH': 'CS-LA14.2.c',
	'CS14S4.26': 'CS-LA14.2.d',
	'CS15.4': 'CS-LA15.4',
	}

# general periods
xticks = [1,2,3,4,5] 
xticks = np.array( xticks ) 
allPeriods = [1.0, 2.0, 3.0, 5.0, 10.0]
allperiods = np.array(allPeriods) 

G_colorsF = ['r','#FF8C00', '#FFFF00', '#7CFC00','#4682B4']  # stack face color
#G_colorsE = ['none','none','none','none','none']   # stack edgecolor
G_colorsE = 5*['k']
G_hatchs = ['','','','','']
G_labels = [r'${\sigma}^2_B$',r'${\bar\sigma}^2_C$',r'${\bar\sigma}^2_D$',r'${\bar\sigma}^2_M$',r'${\bar\sigma}^2_F$ or ${\sigma}^2_T$'] 

g_colors = ['r','#FF8C00', '#FFFF00', '#7CFC00','#4682B4'] 
g_labels = [r'${\sigma}^2_b$',r'${\bar\sigma}^2_c$',r'${\bar\sigma}^2_d$',r'${\bar\sigma}^2_m$',r'${\bar\sigma}^2_f$'] 
g_hatchs = ['+','x','o','/','.'] 

plt.rc('font',family='Arial')
fsize=14
plt.rcParams['xtick.labelsize']=fsize
plt.rcParams['ytick.labelsize']=fsize
plt.rcParams['axes.labelsize']=fsize

pfmt='eps' 
pfmt='png'

# Inputs 
wrk = '/Users/fengw/work/Project/ABFanalysis/scripts/ModelAssemble/' 
pltpth0 = '/Users/fengw/work/Project/ABFanalysis/products/Variance' 
if not os.path.exists(pltpth0): 
    os.mkdir(pltpth0) 

if opt == 'singleModel': 
    
    model=sys.argv[2]    # use old convention
    
    updatedModelName = mappedNameDict[model] 
    if model in ['CS14S4.26', 'CS14.2CVMH', 'CS15.4',]:
	ngaModelVersion = 'NGA14' 
    else: 
	ngaModelVersion = 'NGA08' 
    
    pltpth = os.path.join(pltpth0, updatedModelName)
    if not os.path.exists(pltpth): 
	os.mkdir(pltpth) 
    
    startIndex = 0 
    if model != 'CS15.4': 
        periods = [2.0,3.0,5.0,10.0] 
        startIndex = 1
    else: 
	periods = allPeriods
    Np = len(periods)

    # basic plots given CyberShake model (absolute and residual with NGA models for the given CyberShake model)
    inpth = wrk + 'inputs/ABFvariances_%s/%s'%(ngaModelVersion, model)  

    G_NGA = []; G_CS = [] 
    g_NGA = []; g_CS = []
    for it in xrange( len(periods) ): 
	T = '%.2f'%periods[it]
	
	if model in ['CSbbp1D','CS14S4.26','CS14.2CVMH', 'CS15.4']:
	    infile_G = inpth + '/absolute_T%s.csv'%(T) 
	    infile_g = inpth + '/residual_T%s.csv'%(T) 
	    Ginputs = np.loadtxt( infile_G, usecols=range(1,7), skiprows=1, delimiter=',') 
	    ginputs = np.loadtxt( infile_g, usecols=range(1,7), skiprows=1, delimiter=',') 
	else:
	    infile_G = inpth + '/absolute_T%s.txt'%(T) 
	    infile_g = inpth + '/residual_T%s.txt'%(T) 
	    Ginputs = np.loadtxt( infile_G ) 
	    ginputs = np.loadtxt( infile_g ) 
	
	# Absolute:
	# CS 
	X_CS = Ginputs[5,:]**2  # labels: B,C,D,M,F,G  
	#X_CS = np.sqrt(X_CS)   # to variance
	G_CS.append( X_CS ) 

	# NGA mean 
	X_NGA = np.average(Ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
	#X_NGA = np.sqrt( X_NGA )
	G_NGA.append( X_NGA )  

	# residuals 
	# NGA mean 
	X_NGA = np.average(ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
	#X_NGA = np.sqrt( X_NGA )
	g_NGA.append( X_NGA )  
    
    G_CS = np.array(G_CS ) # Nt,Nfactors
    G_NGA = np.array(G_NGA) 
    g_NGA = np.array(g_NGA) 

    # =======================
    # 1. Variance buget for G
    fig = plt.figure(1, (12,8)) 
    plt.rc(('xtick.major'),pad=10)
    ax = fig.add_subplot(111) 
    G_NGA = np.array( G_NGA )
    G_CS = np.array(G_CS ) # Nt,Nf
    Nt,Nf = G_CS.shape
    Nf = Nf - 1
    width = 0.15; dx=width/2.
    p = []
    bottom1 = np.zeros( Np )
    bottom2 = np.zeros( Np )

    #plt.grid(axis='y',color='grey',linestyle='-')
    plt.grid(axis='y')
    for ifactor1 in xrange(Nf):
	ifactor = Nf - ifactor1 - 1
	if ifactor != Nf-1: 
	    bottom1 += G_CS[:,ifactor+1]
	    bottom2 += G_NGA[:,ifactor+1]
	    p0 = ax.bar( xticks[startIndex:]+dx, G_CS[:,ifactor], width, bottom=bottom1, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	    ax.bar( xticks[startIndex:]-dx-width, G_NGA[:,ifactor], width, bottom=bottom2, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	else: 
	    p0 = ax.bar( xticks[startIndex:]+dx, G_CS[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
	    ax.bar( xticks[startIndex:]-dx-width, G_NGA[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
	p.append( p0[0] ) 

    # used for distinguish NGA and CS
    bottom1 = bottom1 + G_CS[:,ifactor] 
    bottom2 = bottom2 + G_NGA[:,ifactor] 

    p = tuple(p[::-1]) 
    #plt.legend( p, tuple(G_labels), loc=0, fontsize=12 ).draw_frame(False)
    ax.legend( p, tuple(G_labels), loc=9, fontsize=fsize, bbox_to_anchor=(1.07,0.6) ).draw_frame(False)
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_ylim( [0,1.8] )
    ax.set_xlim( [0.5,5.7] )
    #ax.set_title('Variance Buget for G')
    
    # xticks 
    #plt.xticks(xticks, ['NGAnmean CS11\n%.1f'%T for T in periods])
    plt.xticks(xticks, ['%.1f'%T for T in periods])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')
    
    yticks = list(np.arange( 0,2.0, 0.2) )
    plt.yticks(yticks, ['%.1f'%tmp for tmp in yticks] )
    #ax.set_title(updatedModelName)

    # text 
    bottom = max( [bottom1[0], bottom2[0]] )
    ax.text( xticks[startIndex+1]+dx, bottom+0.1, updatedModelName,ha='left' )
    ax.text( xticks[startIndex+1]-dx-width, bottom+0.1, '%s\nmean'%ngaModelVersion, ha='center' )
    fig.savefig( os.path.join(pltpth,'%s_Variance_Absolute.%s'%(updatedModelName,pfmt)), format=pfmt, dpi=300 )

    # =======================
    # Variance buget for g
    fig = plt.figure(2,(12,8)) 
    ax = fig.add_subplot(111) 
    g_NGA = np.array( g_NGA )
    Nt,Nf = g_NGA.shape
    Nf = Nf - 1
    width = 0.4; dx=width/2.
    p = []
    bottom2 = np.zeros( Np )
    plt.grid(axis='y')
    for ifactor1 in xrange(Nf): 
	ifactor = Nf - ifactor1 - 1
	if ifactor != Nf-1: 
	    bottom2 += g_NGA[:,ifactor+1]
	    p0 = ax.bar( xticks[startIndex:]-dx, g_NGA[:,ifactor], width, bottom=bottom2, color=g_colors[ifactor], edgecolor='none' )
	else: 
	    p0 = ax.bar( xticks[startIndex:]-dx, g_NGA[:,ifactor], width, color=g_colors[ifactor], edgecolor='none' ) 
	p.append( p0[0] ) 
    
    p0 = ax.plot( xticks[startIndex:], G_NGA[:,4], 'k--^', lw=3, ms=10 ) # sigma_T
    p = p[::-1]
    p.append( p0[0] ) 
    sigmaTlabel=r'NGA ${\sigma}^2_T$' 
    g_labels.append( sigmaTlabel ) 
    
    # used for distinguish NGA and CS
    bottom2 = bottom2 + g_NGA[:,ifactor] 

    p = tuple(p) 
    #print g_labels 
    ax.legend( p, tuple(g_labels), loc=9, fontsize=fsize, bbox_to_anchor=(1.02,0.6) ).draw_frame(False)
    ax.set_title(updatedModelName)

    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_ylim( [0, 0.7] )
    ax.set_xlim( [0.5,5.7])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')
    #ax.set_title('Variance Buget for g')

    # xticks 
    plt.xticks(xticks, ['%.1f'%T for T in periods])
    yticks = list(np.arange( 0,0.8, 0.1) )
    plt.yticks(yticks, ['%.1f'%tmp for tmp in yticks] )
    
    fig.savefig( os.path.join(pltpth,'%s_Variance_Residual.%s'%(updatedModelName,pfmt)), format=pfmt, dpi=300 )


    # =========================================
    # Variance summary for NGA08 and CS models 
    fig = plt.figure(3,(12,8)) 
    ax = fig.add_subplot(111) 
    ax.plot( periods, G_NGA[:,4], 'g--o', lw=3, ms=10, mec='none', label='NGA ${\sigma}^2_T$' )
    ax.plot( periods, g_NGA[:,5], color='#FF8C00', linestyle='--', marker='o', lw=3, ms=10, mec='none', label=r'%s-NGA ${\bar\sigma}^2_g$'%(updatedModelName))
    ax.plot( periods, G_NGA[:,5], 'b-^', lw=3, ms=10, mec='none', label=r'NGA ${\bar\sigma}^2_G$' ) 
    ax.plot( periods, G_CS[:,5], 'r-s', lw=3, ms=10, mec='none', label=r'%s ${\bar\sigma}^2_G$'%(updatedModelName))
    ax.legend(loc=0, fontsize=fsize, numpoints=1, bbox_to_anchor=(0.85,0.6)).draw_frame(False)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')

    ax.set_xlim( [1,12] )
    ax.set_ylim( [0,1.8])
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    plt.xticks( periods, ['%.1f'%T for T in periods] ) 
    yticks = list(np.arange( 0,1.8, 0.2) )
    plt.yticks(yticks, ['%.1f'%tmp for tmp in yticks] )
    plt.grid(axis='y')
    ax.set_title('Variance Summary for %s and %s'%(ngaModelVersion,updatedModelName))
    fig.savefig( os.path.join(pltpth, '%s_VarianceSummary.%s'%(updatedModelName,pfmt)), format=pfmt )


if opt == 'allModel': 
    
    pltpth = os.path.join(pltpth0, 'allModels')
    if not os.path.exists(pltpth): 
	os.mkdir(pltpth) 
    
    # plot all models (NGA08mean, NGA14mean, all CS models for absolute and residual) (bar) 
    fig = plt.figure(1, (12,8)) 
    fig.clf()
    
    plt.rc(('xtick.major'),pad=10)
    plt.grid(axis='y') 

    ax = fig.add_subplot(111) 
    width = 0.07 
    
    ibar = 2
    p = []

    G_NGA1 = []; G_NGA2 = []
    g_NGA1 = []; g_NGA2 = []
    for model in ['CS11','CS13.1','CS13.2', 'CSbbp1D', 'CS14.2CVMH', 'CS14S4.26', 'CS15.4']: 
	
	if model in ['CS14S4.26', 'CS14.2CVMH', 'CS15.4',]:
	    ngaModelVersion = 'NGA14' 
	else: 
	    ngaModelVersion = 'NGA08' 
	
	startIndex = 0 
	if model != 'CS15.4': 
	    periods = [2.0,3.0,5.0,10.0] 
	    startIndex = 1
	else: 
	    periods = allPeriods 
	
	Np = len(periods)
	updatedModelName = mappedNameDict[model] 

	# basic plots given CyberShake model (absolute and residual with NGA models for the given CyberShake model)
	inpth = wrk + 'inputs/ABFvariances_%s/%s'%(ngaModelVersion, model)  

	G_CS = []; g_CS = []
	for it in xrange( len(periods) ): 
	    T = '%.2f'%periods[it]
	    
	    if model in ['CSbbp1D','CS14S4.26','CS14.2CVMH', 'CS15.4']:
		infile_G = inpth + '/absolute_T%s.csv'%(T) 
		infile_g = inpth + '/residual_T%s.csv'%(T) 
		Ginputs = np.loadtxt( infile_G, usecols=range(1,7), skiprows=1, delimiter=',') 
		ginputs = np.loadtxt( infile_g, usecols=range(1,7), skiprows=1, delimiter=',') 
		
		if model == 'CS15.4':
		    X_NGA = np.average(Ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
		    G_NGA2.append( X_NGA )  
		    X_NGA = np.average(ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
		    g_NGA2.append( X_NGA )  
	    	
	    else:
		infile_G = inpth + '/absolute_T%s.txt'%(T) 
		infile_g = inpth + '/residual_T%s.txt'%(T) 
		Ginputs = np.loadtxt( infile_G ) 
		ginputs = np.loadtxt( infile_g ) 
		
		if model == 'CS11': 
		    X_NGA = np.average(Ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
		    G_NGA1.append( X_NGA )  
		    
		    X_NGA = np.average(ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
		    g_NGA1.append( X_NGA )  
	    
	    # CS 
	    X_CS = Ginputs[5,:]**2  # labels: B,C,D,M,F,G  
	    G_CS.append( X_CS ) 

	G_CS = np.array(G_CS ) # Nt,Nfactors
	
	# =======================
	# 1. Variance buget for G
	Nt,Nf = G_CS.shape
	Nf = Nf - 1
	
	bottom = np.zeros(Np) 
	for ifactor1 in xrange(Nf):
	    ifactor = Nf - ifactor1 - 1
	    if ifactor != Nf-1: 
		bottom += G_CS[:,ifactor+1]
		p0 = ax.bar( xticks[startIndex:]+ibar*width, G_CS[:,ifactor], width, bottom=bottom, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	    else: 
		p0 = ax.bar( xticks[startIndex:]+ibar*width, G_CS[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
	    if ibar == 2: 
		p.append( p0[0] ) 
	ibar += 1 
	#ax.text( xticks[1]+ibar*width, bottom[1]+G_CS[1,ifactor]+0.1, updatedModelName, ha='center',rotation=90 )
	ax.text( xticks[1]+ibar*width, 1.6, updatedModelName, ha='center',rotation=90, fontsize=8)
    
    # add NGA08 and NGA14
    G_NGA1 = np.array(G_NGA1)
    G_NGA2 = np.array(G_NGA2)
    g_NGA1 = np.array(g_NGA1) 
    g_NGA2 = np.array(g_NGA2) 
    
    bottom = np.zeros(4) 
    for ifactor1 in xrange(Nf):
	ifactor = Nf - ifactor1 - 1
	if ifactor != Nf-1: 
	    bottom += G_NGA1[:,ifactor+1]
	    ax.bar( xticks[1:], G_NGA1[:,ifactor], width, bottom=bottom, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	else: 
	    ax.bar( xticks[1:], G_NGA1[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
    ax.text( xticks[1], 1.6, 'NGA08mean', ha='center',rotation=90,fontsize=8 )
    
    bottom = np.zeros(5) 
    for ifactor1 in xrange(Nf):
	ifactor = Nf - ifactor1 - 1
	if ifactor != Nf-1: 
	    bottom += G_NGA2[:,ifactor+1]
	    ax.bar( xticks+width, G_NGA2[:,ifactor], width, bottom=bottom, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	else:
	    ax.bar( xticks+width, G_NGA2[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
    ax.text( xticks[1]+width, 1.6, 'NGA14mean', ha='center',fontsize=8, rotation=90 )

    # finish plots
    p = tuple(p[::-1]) 
    ax.legend( p, tuple(G_labels), loc=9, fontsize=fsize, bbox_to_anchor=(1.02,0.6) ).draw_frame(False)
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_ylim( [0,1.8] )
    ax.set_xlim( [0.5,6.5] )
    #ax.set_title('Variance Buget for G')
    
    # xticks 
    #plt.xticks(xticks, ['NGAnmean CS11\n%.1f'%T for T in periods])
    plt.xticks(xticks, ['%.1f'%T for T in allPeriods])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')
    
    yticks = list(np.arange( 0,2.0, 0.2) )
    plt.yticks(yticks, ['%.1f'%tmp for tmp in yticks] )

    fig.savefig( os.path.join(pltpth, 'AllModels_Variance_Absolute.%s'%(pfmt)), format=pfmt, dpi=300 )

