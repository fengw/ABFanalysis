#!/usr/bin/env python 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

model=sys.argv[1]   # CS11; CS11_CS13.2

model1 = model 

ngaModelVersion = 'NGA14' 

periods = [2.0,3.0,5.0,10.0] 
xticks = list( np.linspace(1,4,4) )
if 1:
    if model == 'CS14Si26':
	periods = [3.0, 5.0, 10.0] 
	xticks = [1,2,3] 
    else:
	# other CyberShake Models has simulations at T=2.0
	periods = [2.0,3.0,5.0,10.0] 
	xticks = [1,2,3,4] 

Np = len(periods) 
periods = np.array( periods ) 
xticks = np.array( xticks ) 

if 1:
    G_colorsF = ['r','#B46400', '#FFFF00', '#7CFC00','#4682B4'] 
    G_colorsF = ['r','#B46400', '#FFFF00', '#7CFC00','k']
    G_colorsF = ['r','#FF8C00', '#FFFF00', '#7CFC00','#4682B4'] 
    G_colorsE = ['none','none','none','none','none']
    G_hatchs = ['','','','','']
else:
    G_colorsE = ['k','k','k','k','k']
    G_colorsF = ['none','none','none','none','none']
    G_hatchs = ['+','x','o','/','.'] 

G_labels = [r'${\sigma}^2_B$',r'${\bar\sigma}^2_C$',r'${\bar\sigma}^2_D$',r'${\bar\sigma}^2_M$',r'${\bar\sigma}^2_F$ or ${\sigma}^2_T$'] 

g_colors = ['r','#B46400', '#FFFF00', '#7CFC00','k'] 
g_colors = ['r','#FF8C00', '#FFFF00', '#7CFC00','#4682B4'] 
g_labels = [r'${\sigma}^2_b$',r'${\bar\sigma}^2_c$',r'${\bar\sigma}^2_d$',r'${\bar\sigma}^2_m$',r'${\bar\sigma}^2_f$'] 
g_hatchs = ['+','x','o','/','.'] 

plt.rc('font',family='Arial')
fsize=14
plt.rcParams['xtick.labelsize']=fsize
plt.rcParams['ytick.labelsize']=fsize
plt.rcParams['axes.labelsize']=fsize

pfmt='png'
pfmt='eps' 

# Inputs 
wrk = '/Users/fengw/work/Project/ABFanalysis/scripts/ModelAssemble/' 

if model in ['CS11','CS13.1','CS13.2', 'CS14S4.26', 'CS14.2CVMH']: 
   
    # basic plots given CyberShake model (absolute and residual with NGA models for the given CyberShake model)
    inpth = wrk + 'inputs/ABFvariances_%s/%s'%(ngaModelVersion, model)  

    G_NGA = []; G_CS = [] 
    g_NGA = []; g_CS = []
    for it in xrange( len(periods) ): 
	T = '%.2f'%periods[it]
	
	if model in ['CS14S4.26','CS14.2CVMH']:
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
	    p0 = ax.bar( xticks+dx, G_CS[:,ifactor], width, bottom=bottom1, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	    ax.bar( xticks-dx-width, G_NGA[:,ifactor], width, bottom=bottom2, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor], hatch=G_hatchs[ifactor] )
	else: 
	    p0 = ax.bar( xticks+dx, G_CS[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
	    ax.bar( xticks-dx-width, G_NGA[:,ifactor], width, color=G_colorsF[ifactor], edgecolor=G_colorsE[ifactor],hatch=G_hatchs[ifactor] ) 
	p.append( p0[0] ) 

    # used for distinguish NGA and CS
    bottom1 = bottom1 + G_CS[:,ifactor] 
    bottom2 = bottom2 + G_NGA[:,ifactor] 

    p = tuple(p[::-1]) 
    #plt.legend( p, tuple(G_labels), loc=0, fontsize=12 ).draw_frame(False)
    ax.legend( p, tuple(G_labels), loc=9, fontsize=fsize, bbox_to_anchor=(1.02,0.6) ).draw_frame(False)
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_ylim( [0,1.8] )
    ax.set_xlim( [0.5,4.7] )
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
    ax.set_title(model) 

    # text 
    bottom = max( [bottom1[0], bottom2[0]] )
    ax.text( xticks[0]+dx, bottom+0.1, model,ha='left' )
    ax.text( xticks[0]-dx-width, bottom+0.1, '%s\nmean'%ngaModelVersion, ha='center' )
    fig.savefig( './plots/%s_Variance_Absolute.%s'%(model,pfmt), format=pfmt, dpi=300 )

    # =======================
    # Variance buget for g
    fig = plt.figure(2,(12,8)) 
    ax = fig.add_subplot(111) 
    g_NGA = np.array( g_NGA )
    print g_NGA
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
	    p0 = ax.bar( xticks-dx, g_NGA[:,ifactor], width, bottom=bottom2, color=g_colors[ifactor], edgecolor='none' )
	else: 
	    p0 = ax.bar( xticks-dx, g_NGA[:,ifactor], width, color=g_colors[ifactor], edgecolor='none' ) 
	p.append( p0[0] ) 
    
    p0 = ax.plot( xticks, G_NGA[:,4], 'k--^', lw=3, ms=10 ) # sigma_T
    p = p[::-1]
    p.append( p0[0] ) 
    sigmaTlabel=r'NGA ${\sigma}^2_T$' 
    g_labels.append( sigmaTlabel ) 
    
    # used for distinguish NGA and CS
    bottom2 = bottom2 + g_NGA[:,ifactor] 

    p = tuple(p) 
    #print g_labels 
    ax.legend( p, tuple(g_labels), loc=9, fontsize=fsize, bbox_to_anchor=(1.02,0.6) ).draw_frame(False)
    ax.set_title(model) 

    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_ylim( [0, 0.7] )
    ax.set_xlim( [0.5,4.7])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')
    #ax.set_title('Variance Buget for g')

    # xticks 
    plt.xticks(xticks, ['%.1f'%T for T in periods])
    yticks = list(np.arange( 0,0.8, 0.1) )
    plt.yticks(yticks, ['%.1f'%tmp for tmp in yticks] )
    
    fig.savefig( './plots/%s_Variance_Residual.%s'%(model,pfmt), format=pfmt, dpi=300 )


    # =========================================
    # Variance summary for NGA08 and CS models 
    fig = plt.figure(3,(12,8)) 
    ax = fig.add_subplot(111) 
    ax.plot( periods, G_NGA[:,4], 'g--o', lw=3, ms=10, mec='none', label='NGA ${\sigma}^2_T$' )
    ax.plot( periods, g_NGA[:,5], color='#FF8C00', linestyle='--', marker='o', lw=3, ms=10, mec='none', label=r'%s-NGA ${\bar\sigma}^2_g$'%(model))
    ax.plot( periods, G_NGA[:,5], 'b-^', lw=3, ms=10, mec='none', label=r'NGA ${\bar\sigma}^2_G$' ) 
    ax.plot( periods, G_CS[:,5], 'r-s', lw=3, ms=10, mec='none', label=r'%s ${\bar\sigma}^2_G$'%(model))
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
    ax.set_title('Variance Summary for %s and %s'%(ngaModelVersion,model))
    fig.savefig( './plots/%s_VarianceSummary.%s'%(model,pfmt), format=pfmt )

    plt.show() 

else: 
    # Plot to show the changes in those variance by changing CyberShake models
    NGAmodels = ['AS','BA','CB','CY']
    inga = 0 
    for imodel in xrange( len(NGAmodels) ): 
	if model == NGAmodels[imodel]: 
	    inga = imodel 
	    break 

    # for absolute: AS.S, AS.H, CS11, and CS13.H 
    #models = ['CS11','CS13.2'] 
    #models = ['CS13.1', 'CS14Si26']
    models = model1.split('_') 
    print models
    G1_CS = {}
    G1_NGA = {}
    g1_CS = {}
    g1_NGA = {}
    for model in models: 
	inpth = wrk + 'inputs/ABFvariances/%s'%model  

	G_NGA = []; G_CS = [] 
	g_NGA = []; g_CS = []
	for it in xrange( len(periods) ): 
	    T = '%.2f'%periods[it]
	    infile_G = inpth + '/absolute_T%s.txt'%(T) 
	    infile_g = inpth + '/residual_T%s.txt'%(T) 
	    
	    Ginputs = np.loadtxt( infile_G ) 
	    ginputs = np.loadtxt( infile_g ) 
	    
	    # Absolute:
	    # CS 
	    X_CS = Ginputs[5,:]**2  # labels: B,C,D,M,F,G  
	    G_CS.append( X_CS ) 

	    # AS
	    if model1 != 'NGAmean':
		X_NGA = Ginputs[inga,:]**2 
	    else: 
		X_NGA = np.average(Ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 

	    G_NGA.append( X_NGA )  

	    # residuals 
	    if model1 != 'NGAmean':
		X_NGA = ginputs[inga,:]**2 
	    else: 
		X_NGA = np.average(ginputs[:4,:]**2, axis=0, weights=np.ones(4)/(1.0*4)) 
	    g_NGA.append( X_NGA )  

        G1_CS[model] = np.array(G_CS) 
        G1_NGA[model] = np.array(G_NGA)
        g1_CS[model] = np.array(g_CS) 
        g1_NGA[model] = np.array(g_NGA)

    G_NGA = np.array( G_NGA )
    G_CS = np.array(G_CS ) # Nt,Nf
    Nt,Nf = G_CS.shape
    Nf = Nf - 1
    
    # =======================
    # Variance buget for G
    fig = plt.figure(1) 
    ax = fig.add_subplot(111) 
    
    width = 0.3; dx=0.1; offset = dx/2
    p = []
    bottoms = np.zeros( (4,Np) )

    for ifactor1 in xrange(Nf): 
	ifactor = Nf-ifactor1-1
	if ifactor != Nf-1: 
	    model = models[0]
	    bottoms[0,:] += G1_NGA[model][:,ifactor+1]
	    p0 = ax.bar( xticks-offset-dx-2*width, G1_NGA[model][:,ifactor], width, bottom=bottoms[0,:], color=G_colorsF[ifactor], edgecolor='none' )
	    bottoms[1,:] += G1_CS[model][:,ifactor+1]
	    ax.bar( xticks-offset-width, G1_CS[model][:,ifactor], width, bottom=bottoms[1,:], color=G_colorsF[ifactor], edgecolor='none' )
	    model = models[1]
	    bottoms[2,:] += G1_NGA[model][:,ifactor+1]
	    ax.bar( xticks+offset, G1_NGA[model][:,ifactor], width, bottom=bottoms[2,:], color=G_colorsF[ifactor], edgecolor='none' )
	    bottoms[3,:] += G1_CS[model][:,ifactor+1]
	    ax.bar( xticks+offset+dx+width, G1_CS[model][:,ifactor], width, bottom=bottoms[3,:], color=G_colorsF[ifactor], edgecolor='none' )
	else: 
	    model = models[0]
	    p0 = ax.bar( xticks-offset-dx-2*width, G1_NGA[model][:,ifactor], width, color=G_colorsF[ifactor], edgecolor='none' )
	    ax.bar( xticks-offset-width, G1_CS[model][:,ifactor], width, color=G_colorsF[ifactor], edgecolor='none' )
	    model = models[1]
	    ax.bar( xticks+offset, G1_NGA[model][:,ifactor], width, color=G_colorsF[ifactor], edgecolor='none' )
	    ax.bar( xticks+offset+dx+width, G1_CS[model][:,ifactor], width, color=G_colorsF[ifactor], edgecolor='none' )
	p.append( p0[0] ) 

    # used for distinguish NGA and CS
    model = models[0]
    bottoms[0,:] += G1_NGA[model][:,ifactor]
    bottoms[1,:] += G1_CS[model][:,ifactor]
    model = models[1]
    bottoms[2,:] += G1_NGA[model][:,ifactor]
    bottoms[3,:] += G1_CS[model][:,ifactor]

    p = tuple(p)[::-1] 
    ax.legend( p, tuple(G_labels), loc=0, fontsize=fsize )
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_ylim( [0,2] )
    #ax.set_title('Variance Buget for G')

    # xticks 
    plt.xticks(xticks, ['%.1f'%T for T in periods])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')
    
    # text 
    bottom0 = max( bottoms[:,0] )
    ha = 'left'
    #ax.text( xticks[1]-dx-offset-2*width+offset, bottom0+0.1, '%s.S'%model1, ha=ha, fontsize=fsize )
    #ax.text( xticks[1]-offset-width+offset, bottom0+0.1, 'CS11', ha=ha, fontsize=fsize)
    #ax.text( xticks[1]+offset+offset, bottom0+0.1, '%s.H'%model1, ha=ha,fontsize=fsize )
    #ax.text( xticks[1]+offset+dx+width+offset, bottom0+0.1, 'CS13.2', ha=ha,fontsize=fsize )
    plt.grid(axis='y')
    fig.savefig( './plots/%s_Variance_Absolute.%s'%(model1, pfmt), format=pfmt )

    # =======================
    # Variance buget for g
    fig = plt.figure(2,(12,8)) 
    ax = fig.add_subplot(111) 
    g_NGA = np.array( g_NGA )
    Nt,Nf = g_NGA.shape
    Nf = Nf - 1
    width = 0.15; dx=0.1; offset = dx/2
    p = []
    bottoms = np.zeros( (2,Np) )
    for ifactor1 in xrange(Nf): 
	ifactor = Nf - ifactor1 - 1
	if ifactor != Nf-1: 
	    model = models[0]
	    bottoms[0,:] += g1_NGA[model][:,ifactor+1]
	    p0 = ax.bar( xticks-offset-width, g1_NGA[model][:,ifactor], width, bottom=bottoms[0,:], color=g_colors[ifactor], edgecolor='none' )
	    model = models[1]
	    bottoms[1,:] += g1_NGA[model][:,ifactor+1]
	    ax.bar( xticks+offset, g1_NGA[model][:,ifactor], width, bottom=bottoms[1,:], color=g_colors[ifactor], edgecolor='none' )
	else: 
	    model = models[0]
	    p0 = ax.bar( xticks-offset-width, g1_NGA[model][:,ifactor], width, color=g_colors[ifactor], edgecolor='none' )
	    model = models[1]
	    ax.bar( xticks+offset, g1_NGA[model][:,ifactor], width, color=g_colors[ifactor], edgecolor='none' )
	p.append( p0[0] ) 

    # used for distinguish NGA and CS
    model = models[0]
    bottoms[0,:] += g1_NGA[model][:,ifactor]
    model = models[1]
    bottoms[1,:] += g1_NGA[model][:,ifactor]

    p = p[::-1]
    
    # plot NGA sigma_T
    p0 = ax.plot( xticks, G1_NGA['CS11'][:,4], 'k--^', lw=3, ms=10 )
    p.append( p0[0])
    p = tuple(p)
    #sigmaTlabel=r'%s ${\sigma}^2_T$'%(model1) 
    sigmaTlabel=r'${\sigma}^2_T$' 
    g_labels.append( sigmaTlabel ) 
    
    plt.legend(p, tuple(g_labels), fontsize=fsize, numpoints=1, bbox_to_anchor=(1.1,0.6)).draw_frame(False)
    ax.set_xlabel('Period (s)')
    ax.set_ylabel('Excitation Variance')
    ax.set_xlim( [0.55,5] )
    ax.set_ylim( [0,0.7] )
    #ax.set_title('Variance Buget for g')

    # xticks 
    plt.xticks(xticks, ['%.1f'%T for T in periods])
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False) 
    ax.yaxis.set_ticks_position('left') 
    ax.xaxis.set_ticks_position('bottom')
    
    bottom0 = max( bottoms[:,0] )
    fsize = 12; ha = 'left'
    ax.text( xticks[1]-offset-width+offset, bottoms[0,1], '%s-%s'%(models[0],ngaModelVersion), ha=ha, fontsize=fsize )
    ax.text( xticks[1]+offset+offset, bottoms[1,1], '%s-%s'%(models[1],ngaModelVersion), ha=ha, fontsize=fsize)
    plt.grid(axis='y')
    fig.savefig( './plots/%s_Variance_Residual.%s'%(model1, pfmt), format=pfmt, dpi=300 )
    plt.show() 
