#!/usr/bin/env python 
import os, sys
import optparse 

import numpy as np 
import matplotlib.pyplot as plt 
from scipy import stats 

from pyABF import *
from pynga import * 

from my_util.numer import interp 

# Basic parameters
parser = optparse.OptionParser() 
parser.add_option('-m','--model-flag', dest='mflag', type='string', help='model selection') 
parser.add_option('-r','--rupture-model-id', dest='rup_model_id', type='int', nargs=4, help='rupture model specfication') 
parser.add_option('-p','--period', dest='T', type='float', help='period') 
parser.add_option('-S','--sigma', dest='sigma', type='string', help='parameters controling the hypocenter distribution') 
parser.add_option('-G','--Gflag', dest='Gflag', type='int', help='parameters controling the absolute or residual')
#parser.add_option('-s','--source-id', dest='sid', type='int', help='SourceID' ) 
# example: CalculateSigma.py -m 54 -r 35 8 4 5 -p 3.00 -S 1.00_1.00 -G 0 

# CyberShake Dictionary
CSdict = {'(35,5,3,1)':'CS11',
          '(35,7,4,1)':'CS13.1',
	  '(35,7,4,4)':'CS13.2',
	  '(35,8,4,5)':'CS14S4.26',
	  '(35,6,4,7)':'CS14.2CVMH',
	  '(35,8,4,8)':'CSbbp1D',
	  }

# use in check inputs
parser.print_help() 

(options, args) = parser.parse_args() 

mflag = options.mflag 
rup_model_id = options.rup_model_id 
erf_id, sgt_id, rvid, vel_id = rup_model_id
T = '%3.2f'%options.T
sigma = '%s'%options.sigma 
Gflag = options.Gflag
#sid = options.sid 

# Inputs
mapin = '/Users/fengw/work/Project/ABFanalysis/scripts/map_input/'
Gpth0 = mapin + 'Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs0/'%(mflag, erf_id, sgt_id, rvid, vel_id)
Epth = Gpth0 + 'Ekxms/Sigma%s'%(sigma)

Gpth = mapin + 'Model_Rups%s/ERF%s_SGT%s_RupVar%s_Vel%s/Gkxmfs/'%(mflag, erf_id, sgt_id, rvid, vel_id)
Fpth = Gpth + 'SFks/Sigma%s'%(sigma)    # only for CyberShake (column 3)
Dpth = Gpth + 'Dkxs/Sigma%s'%(sigma)
Cpth = Gpth + 'Cks/Sigma%s'%(sigma)
Bpth = Gpth + 'Bs/Sigma%s'%(sigma) 

mapout0 = '/Users/fengw/work/Project/ABFanalysis/scripts/map_input/' 
mapout1 = mapout0 + 'Model_Rups%s'%mflag 
mapout2 = mapout1 + '/ERF%s_SGT%s_RupVar%s_Vel%s'%(erf_id,sgt_id,rvid,vel_id ) 
mapout3 = mapout2 + '/Gkxmfs' 
for f in [mapout0, mapout1, mapout2, mapout3,]:
    if not os.path.exists( f ): 
	os.mkdir(f) 

# interpolation parameters
eps = 0.001
smooth = 0.01
method = {'name':'exp','smooth':smooth}
site_locs1 = './SiteInterpMesh.txt'
sites = np.loadtxt( site_locs1 )
SiteLon1D = sites[:,0]
SiteLat1D = sites[:,1]

# read in source info 
srcfile = Gpth0 + 'SourceInfo' 
srids = np.loadtxt( srcfile ) 
sids = srids[:,0] 
Areas = srids[:,1]

# read in hypocenter info 
hypofile = Gpth0 + 'SourceRuptureHypoInfo' 
inputs = np.loadtxt( hypofile ) 
sidhypo = {} 
for i in xrange( inputs.shape[0] ): 
    sidhypo[str(int(inputs[i,0]))] = inputs[i,2]   # number of hypocenters

# source weighting function 
wrk = '/Users/fengw/work/Project/ABFanalysis'
if mflag[0] == '5': 
    # Test same disaggregation condition, but different Ntop (number of sources used) 
    IML = '0.3_G'   # consider rock sites or choose large IML value
    Str0 = 'DisaggSources_ERF%s_'%(erf_id)
    Str1 = 'DisaggIML_' + IML 
    Str2 = '_SA_3sec.txt'

    Ntops = [5,10,15,20,25,30]
    if mflag[1] in ['1','2','3','4','5','6']: 
	Ntop = Ntops[int(mflag[1])-1] 
    DisaggSourceFile = wrk + '/metadata/Disaggregation/DisaggSources/UniformSiteDistribution/' + Str0 + Str1 + Str2 
    sids, srcPDF = Disaggregation.GetSourceIDs(DisaggSourceFile,Ntop=Ntop)

# sigma values output path
CSmodel = CSdict['(%s,%s,%s,%s)'%(erf_id,sgt_id,rvid,vel_id)]
print 'CyberShake Model:', CSmodel
erf_id, sgt_id, rvid, vel_id = rup_model_id
soutpth = '/Users/fengw/work/Project/ABFanalysis/scripts/ModelAssemble/inputs/ABFvariances_NGA14/%s'%CSmodel
if not os.path.exists(soutpth):
    os.mkdir(soutpth) 
if Gflag: 
    soutfile = soutpth + '/absolute_T%s.csv'%T
    print 'Absolute sigmas'
else: 
    soutfile = soutpth + '/residual_T%s.csv'%T
    print 'Residual sigmas'

# Sigma_T after SC08 correction 
NGAmodel = ['CB','BA','CY','AS']
sigmaT = []
for nga in NGAmodel: 
    SC08 = SC08_model( nga + '08' )
    ind = (np.array(SC08.periods) == float(T)).nonzero()[0]
    tau = SC08.tau[ind] 
    sigma0 = SC08.sigma0[ind] 
    sigmaT0 = np.sqrt( tau**2 + sigma0**2 )
    #print nga, sigma0**2, tau**2, sigmaT**2
    print nga, sigma0, tau, sigmaT0 
    sigmaT.append( sigmaT0**2 )   # variance unit

# ===================
# Sigma_F 
# ===================
print 'Compute sigma_U...'
sigma_F_ks = []
for isid in xrange( len(sids) ):
    
    sid = str(int(sids[isid])) 
    AreaTmp = Areas[isid]
    Nh = int(sidhypo[sid])

    # read in Magnitude info for given sid 
    Ffile = Fpth + '/' + sid + '/CyberShake.NGAs.%s.Source%s.Sfks'%(T,sid)
    inputs = np.loadtxt( Ffile ) 
    lons0 = inputs[:,0] 
    lats0 = inputs[:,1]
    sigma_F_ks.append( inputs[:,2] )   # only CyberShake has sigma_F^2 (careful here)
    
sigma_F_ks = np.array(sigma_F_ks)  # dimension: [Nsid, Nsta,Nmodels], type: variance

# compute the average of CyberShake
pdf_s = np.repeat( 1./sigma_F_ks.shape[1], sigma_F_ks.shape[1] )
sigma_F_CS = np.average( np.average( sigma_F_ks, axis=0, weights = srcPDF ), axis=0, weights = pdf_s )
print 'sigma_F', mflag, T, rup_model_id, np.sqrt( sigma_F_CS )     # attention: printed are all sigmas and what are written in table is sigma as well

sigmaT.append(sigma_F_CS)
sigmaT.append(0)   # for reference model
sigma_F = np.array(sigmaT)
print sigma_F.shape

# ===================
# Ekxms for sigma_M
# ===================
print 'Compute sigma_M...'
mapout4 = mapout3 + '/SEks'
mapout = mapout4 + '/Sigma%s'%(sigma) 

for f in [mapout4, mapout, ]: 
    if not os.path.exists( f ): 
	os.mkdir(f) 

sigma_M_ks = []
for isid in xrange( len(sids) ):
    
    sid = str(int(sids[isid])) 
    AreaTmp = Areas[isid]
    Nh = int(sidhypo[sid])

    # Calculate the hypocenter weighting function (pdf_x)
    xh = (np.arange( Nh )+1.0) / (Nh+1.0)
    spl = sigma.strip().split('_')
    alpha, beta = float(spl[0]), float(spl[1]) 
    prob0 = stats.beta.pdf(xh,alpha,beta)     # use scipy.stats
    pdf_x = (prob0/sum(prob0)).tolist()    # to make sure sum(prob) = 1

    # read in Magnitude info for given sid 
    #Mfile = Epth + '/' + sid + '/SourceRuptureMwInfo' 
    #Mws = np.loadtxt( Mfile )
    tmp_Mws = glob.glob( Epth+'/'+sid +'/M*' )
    Mws = []
    for tmp_Mw in tmp_Mws: 
	spl = tmp_Mw.strip().split('/')[-1] 
	Mws.append( float( spl[1:] ) ) 

    mu = 3.87 + np.log10(AreaTmp) * 1.05   # somerville 2006 M-A relationship
    std = 0.2
    prob0 = 1./(np.sqrt(2.*np.pi)*std) * np.exp(-0.5*((Mws-mu)/std)**2)
    pdf_m = (prob0 / sum(prob0)).tolist()   # weighting should be summed up to 1
    
    tmp_ekxms = []
    for ihypo in xrange( Nh ): 
	tmp1 = []
	for Mw in Mws: 
	    Mw0 = str(Mw) 
	    Efile = Epth + '/'+sid+'/'+'M%.2f'%Mw +'/' +'Period%s.Hypo%s.ekxms'%(T,ihypo) 

	    inputs = np.loadtxt( Efile ) 
	    lons0 = inputs[:,0] 
	    lats0 = inputs[:,1]
	    e1 = inputs[:,2:8]   # CS-CB, CS-BA, CS-CY, CS-AS, CS, CS-Ref 
	    if Gflag:
		# absolute
		ekxms0 = np.zeros( e1.shape )
		ekxms0 = -( e1[:,:6]-e1[:,4:5] ) # NGA and Ref 
		ekxms0[:,4:5] = e1[:,4:5]   # CS 
	    else: 
		ekxms0 = e1 

	    tmp1.append( ekxms0 ) 
	tmp_ekxms.append( tmp1 ) 
    tmp_ekxms = np.array( tmp_ekxms )   # [Nhypo, NMw, Nsta, Nmodels]
    
    # standard deviation and average over hypocenters 
    tmp_sd_xs = np.sqrt( np.average( tmp_ekxms**2, axis=1,weights=pdf_m) )
    tmp_sd0 = np.average( tmp_sd_xs, axis=0, weights = pdf_x ) 
    tmp_vd0 = np.average( tmp_sd_xs**2, axis=0, weights = pdf_x ) 
    tmp_sd = []; tmp_vd = []
    for imodel in xrange( tmp_sd0.shape[1] ):
	NewValue = interp( lons0, lats0, tmp_sd0[:,imodel], SiteLon1D, SiteLat1D, eps=eps, method=method )
        tmp_sd.append( NewValue ) 
	NewValue = interp( lons0, lats0, tmp_vd0[:,imodel], SiteLon1D, SiteLat1D, eps=eps, method=method )
        tmp_vd.append( NewValue ) 
    
    sigma_M_ks.append( tmp_vd )    # variance unit

    pafile = mapout + '/CyberShake.NGAs.%s.Source%s.Seks'%(T,sid)
    fid = open( pafile, 'w' ) 
    Nrow = len(SiteLon1D)
    for irow in xrange( Nrow ): 
	fid.write( '%s %s %s %s %s %s %s %s\n'%(SiteLon1D[irow],SiteLat1D[irow], \
		   tmp_sd[0][irow], 
		   tmp_sd[1][irow], 
		   tmp_sd[2][irow], 
		   tmp_sd[3][irow], 
		   tmp_sd[4][irow], 
		   tmp_sd[5][irow] )) 
    fid.close() 

sigma_M_ks = np.array(sigma_M_ks)  # dimension: [Nsid, Nmodel, Nsta ]

# compute the average 
pdf_s = np.repeat( 1./sigma_M_ks.shape[2], sigma_M_ks.shape[2] )
sigma_M = np.average( np.average( sigma_M_ks, axis=0, weights = srcPDF ), axis=1, weights = pdf_s )
print 'sigma_M', mflag, T, rup_model_id, np.sqrt( sigma_M )

# ===================
# Dkxs for sigma_D
# ===================
print 'Compute sigma_D...'
mapout4 = mapout3 + '/SDksD'   # distinguish with SDks (residual sigma)
mapout = mapout4 + '/Sigma%s'%(sigma) 

for f in [mapout4, mapout, ]: 
    if not os.path.exists( f ): 
	os.mkdir(f) 

sigma_D_ks = []
for isid in xrange( len(sids) ):
    
    sid = str(int(sids[isid])) 
    AreaTmp = Areas[isid]
    Nh = int(sidhypo[sid])

    # Calculate the hypocenter weighting function (pdf_x)
    xh = (np.arange( Nh )+1.0) / (Nh+1.0)
    spl = sigma.strip().split('_')
    alpha, beta = float(spl[0]), float(spl[1]) 
    prob0 = stats.beta.pdf(xh,alpha,beta)     # use scipy.stats
    pdf_x = (prob0/sum(prob0)).tolist()    # to make sure sum(prob) = 1

    tmp_dkxs = []
    for ihypo in xrange( Nh ): 
	Dfile = Dpth + '/'+sid+'/'+'CyberShake.NGAs.%s.Source%s.Ih%s.dkxs'%(T,sid,ihypo) 

	inputs = np.loadtxt( Dfile ) 
	lons0 = inputs[:,0] 
	lats0 = inputs[:,1]
	d1 = inputs[:,4:10]   # CS-CB, CS-BA, CS-CY, CS-AS, CS, CS-Ref 
	if Gflag: 
	    dkxs0 = np.zeros( d1.shape )
	    dkxs0 = -( d1[:,:6]-d1[:,4:5] ) # NGA and Ref 
	    dkxs0[:,4:5] = d1[:,4:5]   # CS 
	else: 
	    dkxs0 = d1

	tmp_dkxs.append( dkxs0 ) 
    tmp_dkxs = np.array( tmp_dkxs ) 

    # standard deviation and average over hypocenters 
    tmp_sd = np.sqrt( np.average( tmp_dkxs**2, axis=0,weights=pdf_x) )
    sigma_D_ks.append( tmp_sd**2 )   # variance 

    pafile = mapout + '/CyberShake.NGAs.%s.Source%s.SdksD'%(T,sid)
    fid = open( pafile, 'w' ) 
    Nrow = len(SiteLon1D)
    for irow in xrange( Nrow ): 
	fid.write( '%s %s %s %s %s %s %s %s\n'%(SiteLon1D[irow],SiteLat1D[irow], \
		   tmp_sd[irow, 0], 
		   tmp_sd[irow, 1], 
		   tmp_sd[irow, 2], 
		   tmp_sd[irow, 3], 
		   tmp_sd[irow, 4], 
		   tmp_sd[irow, 5], 
		    )) 
    fid.close() 

# compute average 
sigma_D_ks = np.array(sigma_D_ks) 
pdf_s = np.repeat( 1./sigma_D_ks.shape[1], sigma_D_ks.shape[1] )
sigma_D = np.average( np.average( sigma_D_ks, axis=0, weights = srcPDF ), axis=0, weights = pdf_s )
print 'sigma_D', mflag, T, rup_model_id, np.sqrt( sigma_D )

# ==================
# Sigma_C 
# ==================
print 'Compute Sigma_C...'
mapout4 = mapout3 + '/SCs'   # distinguish with SDks (residual sigma)
mapout = mapout4 + '/Sigma%s'%(sigma) 

for f in [mapout4, mapout, ]: 
    if not os.path.exists( f ): 
	os.mkdir(f) 

tmp_ks = []
for isid in xrange( len(sids) ):
    
    sid = str(int(sids[isid])) 
    AreaTmp = Areas[isid]
    Nh = int(sidhypo[sid])

    Cfile = Cpth + '/CyberShake.NGAs.%s.Source%s.cks'%(T,sid) 

    inputs = np.loadtxt( Cfile ) 
    lons0 = inputs[:,0] 
    lats0 = inputs[:,1]
    c1 = inputs[:,5:11]   # CS-CB, CS-BA, CS-CY, CS-AS, CS, CS-Ref 
    if Gflag:
	cks0 = np.zeros( c1.shape )
	cks0 = -( c1[:,:6]-c1[:,4:5] ) # NGA and Ref 
	cks0[:,4:5] = c1[:,4:5]   # CS 
    else: 
	cks0 = c1 

    tmp_ks.append( cks0 ) 

tmp_ks = np.array( tmp_ks )

# standard deviation and average over hypocenters 
tmp_sd = np.sqrt( np.average( tmp_ks**2, axis=0, weights=srcPDF) )
sigma_C_s = tmp_sd**2   # variance 

pafile = mapout + '/CyberShake.NGAs.%s.SCs'%T
fid = open( pafile, 'w' ) 
Nrow = len(SiteLon1D)
for irow in xrange( Nrow ): 
    fid.write( '%s %s %s %s %s %s %s %s\n'%(SiteLon1D[irow],SiteLat1D[irow], \
	       tmp_sd[irow, 0], 
	       tmp_sd[irow, 1], 
	       tmp_sd[irow, 2], 
	       tmp_sd[irow, 3], 
	       tmp_sd[irow, 4], 
	       tmp_sd[irow, 5], 
		)) 
fid.close() 

# compute average 
pdf_s = np.repeat( 1./sigma_C_s.shape[0], sigma_C_s.shape[0] )
sigma_C = np.average( sigma_C_s, axis=0, weights = pdf_s )
print 'sigma_C', mflag, T, rup_model_id, np.sqrt( sigma_C )


# ==================
# Sigma_B 
# ==================
print 'Compute Sigma_B...'
mapout4 = mapout3 + '/SB'   # distinguish with SDks (residual sigma)
mapout = mapout4 + '/Sigma%s'%(sigma) 

for f in [mapout4, mapout, ]: 
    if not os.path.exists( f ): 
	os.mkdir(f) 

Bfile = Bpth + '/'+'CyberShake.NGAs.%s.bs'%(T) 

inputs = np.loadtxt( Bfile ) 
lons0 = inputs[:,0] 
lats0 = inputs[:,1]
b1 = inputs[:,5:11]   # CS-CB, CS-BA, CS-CY, CS-AS, CS, CS-Ref 
if Gflag:
    bs0 = np.zeros( b1.shape )
    bs0 = -( b1[:,:6]-b1[:,4:5] ) # NGA and Ref 
    bs0[:,4:5] = b1[:,4:5]   # CS 
else: 
    bs0 = b1  

# standard deviation and average over hypocenters 
pdf_s = list(np.repeat(1./bs0.shape[0],bs0.shape[0])) 
tmp_sd = np.sqrt( np.average( bs0**2, axis=0, weights=pdf_s ) )
sigma_B = tmp_sd**2 

pafile = mapout + '/CyberShake.NGAs.%s.SbsB'%(T)
fid = open( pafile, 'w' ) 
fid.write( '%s %s %s %s %s %s\n'%( \
	   tmp_sd[0], 
	   tmp_sd[1], 
	   tmp_sd[2], 
	   tmp_sd[3], 
	   tmp_sd[4], 
	   tmp_sd[5], 
	    )) 
fid.close() 

print 'sigma_B', mflag, T, rup_model_id, np.sqrt( sigma_B )

# write sigmas into file (not varirance, but above are all variance)
sfid = open(soutfile,'w')
if Gflag:
    print 'sigma_T(abf)', mflag, T, rup_model_id, np.sqrt( sigma_B+sigma_C+sigma_D+sigma_M+np.array([sigmaT[0], sigmaT[1], sigmaT[2],sigmaT[3],sigmaT[4],sigmaT[5]]))
    sigma_G = sigma_B+sigma_C+sigma_D+sigma_M+np.array([sigmaT[0], sigmaT[1], sigmaT[2],sigmaT[3],sigmaT[4],sigmaT[5]])
    sfid.write('Model, sigma_B, sigma_C, sigma_D, sigma_M, sigma_F, sigma_G\n')
    Var_ABF = np.array([sigma_B, sigma_C, sigma_D, sigma_M, sigma_F, sigma_G])    # each has size (6, but we want 0:5)
else: 
    print 'sigma_T(abf)', mflag, T, rup_model_id, np.sqrt( sigma_B+sigma_C+sigma_D+sigma_M+np.array([sigma_F_CS,sigma_F_CS,sigma_F_CS,sigma_F_CS, 0, sigma_F_CS]) )
    sigma_f = np.array([sigma_F_CS,sigma_F_CS,sigma_F_CS,sigma_F_CS, 0, sigma_F_CS])   #  NGA has no contribution to residual f variance!!!
    sigma_G = sigma_B+sigma_C+sigma_D+sigma_M+sigma_f
    sfid.write('Model, sigma_b, sigma_c, sigma_d, sigma_m, sigma_f, sigma_g\n') 
    Var_ABF = np.array([sigma_B, sigma_C, sigma_D, sigma_M, sigma_f, sigma_G])    # each has size (6, but we want 0:5)

nrow, ncol = Var_ABF.shape

# compute NGArms [3,1,2,0]
NGA_rms = []
for irow in range(nrow):
    sigma_tmp = Var_ABF[irow,:4]   # NGAs
    NGA_rms.append(np.sqrt(sum(sigma_tmp)/len(sigma_tmp)))    # root mean square

sigma_ABF = np.sqrt(Var_ABF)   
for icol in [3,1,0,2]:
    sfid.write('%s,%s,%s,%s,%s,%s,%s\n'%(NGAmodel[icol],sigma_ABF[0,icol],sigma_ABF[1,icol],sigma_ABF[2,icol], sigma_ABF[3,icol], sigma_ABF[4,icol], sigma_ABF[5,icol]))

# write out the NGA-rms (for Gflag = 0, this is CS-NGArms)
sfid.write('NGA-rms, %s, %s, %s, %s, %s, %s\n'%(NGA_rms[0],NGA_rms[1],NGA_rms[2],NGA_rms[3],NGA_rms[4],NGA_rms[5]))

# write out CyberShake
icol = 4
sfid.write('%s, %s,%s,%s,%s,%s,%s\n'%(CSmodel, sigma_ABF[0,icol],sigma_ABF[1,icol],sigma_ABF[2,icol], sigma_ABF[3,icol], sigma_ABF[4,icol], sigma_ABF[5,icol]))
sfid.close() 
