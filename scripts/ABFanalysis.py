#!/usr/bin/env python
"""
Averaging-Based Factorization of CyberShake Runs
Compare with NGA2 models
"""
from pyABF import * 

# ========================================
# User defined
opt = sys.argv[1]     # Disaggregation, Preparation, Calculation, Examination, Regression
wrk = '/Users/fengw/work/Project/ABFanalysis'   # work path

# =============================================================================
# CyberShake Model Parameters (specify the CyberShake study, one flag, 1 digit)
CyberShakeModelParameters = """
erf_id = 35  # UCERF2
erf_id = 36  # UCERF2 (finer fault grid)

sgt_id = 5   # GP v3
sgt_id = 6   # AWP
sgt_id = 7   # GP v3.0.3
sgt_id = 8   # AWP-GPU?  (CyberShake study 14.2)

rup_scenario_id = 3  # GP 2007
rup_scenario_id = 4  # GP 2010 (available only for CyberShake 13.4)
rup_scenario_id = 6  # GP updated model 2015?

vel_id = 1    # CVM-SCEC
vel_id = 2    # CVM-H v11.2
vel_id = 4    # CVM-H v11.9.1  (updated CVM-H) 
vel_id = 5    # CVM-S4.26   (CVM-SCEC model updates from Po's inversion)
vel_id = 8    # BBP 1D velocity
"""

# user change fill the password of CyberShake database for security purpose
cybershk_database_password = '****'   # this should goes to github 

# CyberShake study/model indicator (currently focused versions)
#erf_id, sgt_id, rup_scenario_id, vel_id = 35, 5, 3, 1    # original CyberShake Models 
#erf_id, sgt_id, rup_scenario_id, vel_id = 35, 7, 4, 4    # CS13.2 (CVMH)
#erf_id, sgt_id, rup_scenario_id, vel_id = 35, 8, 4, 5    # CS14S4.26 (CVM-S4.26)
erf_id, sgt_id, rup_scenario_id, vel_id = 36, 8, 6, 5    # updated CyberShake Models

rup_model_ids = ( erf_id, sgt_id, rup_scenario_id, vel_id )
print 'CyberShake Study: %s'%str(rup_model_ids)

# broadband simulations (for different CyberShake phases)
# bb could be NULL (None), '0.5','10'   in Hz
bb = None    # just low frequency (early definition) 
bb = 0       # include all frequency band (later you can specify the frequencies you interested in the analysis step)

# Target NGA models (two flags, 4 digits)
NGAmodels=['CB','BA','CY','AS'] 
modelVersion = '2008' 
modelVersion = '2014'

if modelVersion == '2014': 
    NGAs={'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1,1,1)},\
	  'BA':{'NewCoefs':None,'terms':(1,1,1)},\
	  'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
	  'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)}, \
	  'SC':{'NewCoefs':{'CB':None,'BA':None,'CY':None,'AS':None},\
	  'predictors':{'cuteps':{'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[200,250],'m_cut':[5.6,6.0]},}} \
	  }
    
    # Reference models (two flags, 4 digits)
    Reference={'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1,1,1)},\
	  'BA':{'NewCoefs':None,'terms':(1,1,1)},\
	  'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
	  'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)}, \
	  'SC':{'NewCoefs':{'CB':None,'BA':None,'CY':None,'AS':None},\
	  'predictors':{'cuteps':{'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[200,250],'m_cut':[5.6,6.0]},}} \
	  }
else: 
    NGAs = { \
	    'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	    'BA':{'NewCoefs':None,'terms':(1,1,1)},\
	    'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	    'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
	    'SC':{'NewCoefs':{'CB':None,'BA':None,'CY':None,'AS':None},\
	    'predictors':{'cuteps':{'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[200,250],'m_cut':[5.6,6.0]},}} \
	    } 

    # Reference models (two flags, 4 digits)
    Reference = { \
		    'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
		    'BA':{'NewCoefs':None,'terms':(1,1,1)},\
		    'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
		    'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
		    'SC':{'NewCoefs':{'CB':None,'BA':None,'CY':None,'AS':None},\
		    'predictors':{'cuteps':{'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[200,250],'m_cut':[5.6,6.0]},}} \
		} 

# Default parameters for hypocenter distribution
hypoPDF = {'CyberShake':{'Beta':[(1.0,1.0),]},\
	   'CB':1,'BA':1,'CY':1,'AS':1}    # uniform 
Ref = 'BA00'; Reference['BA']['terms'] = (1,1,int(Ref[2]))    
srcPDF = None   # default source distribution is uniform

mflag = sys.argv[2]  # default is '0' for workflow test
print 'Model Flag: ', mflag
if mflag == '0':
    # used for debug of the workflow and some special cases
    sids = [255,10]
    srcPDF = [0.5,0.5 ]  # test srcPDF
    #sids = [93,88] 
    #srcPDF = [1./2.]*2 

    #hypoPDF = {'CyberShake':{'Beta':[(1.0,1.0),]}, 'CB':1,'BA':1,'CY':1,'AS':1}      # gaussian sigma == 0.0: uniform distribution for hypocenters
    hypoPDF = {'CyberShake':{'Beta':[(1.0,1.0),(2.0,2.0),]}, 'CB':0,'BA':0,'CY':0,'AS':0}   

if mflag[0] == '5': 
    # this is the basis (54) for ABF analysis: Use top 20 sources from disaggregation 
    # disaggragation files are existing in metadata folder ../metadata
    IML = '0.3_G'   # consider rock sites or choose large IML value
    Str0 = 'DisaggSources_ERF%s_'%(35)
    Str1 = 'DisaggIML_' + IML 
    Str2 = '_SA_3sec.txt'

    Ntops = [5,10,15,20,25,30]
    if mflag[1] in ['1','2','3','4','5','6']: 
	Ntop = Ntops[int(mflag[1])-1] 
    try: 
	if mflag[2] == '1':
	    print 'test'
	    hypoPDF = {'CyberShake':{'Beta':[(1.0,1.0),]}, 'CB':0,'BA':0,'CY':0,'AS':0}   
    except:
	hypoPDF = {'CyberShake':{'Beta':[(1.0,1.0),(2.0,2.0),(0.5,0.5),(2.0,3.5),(1.5,0.5),(0.5,1.5),(3.5,2.0)]},\
		   'CB':1,'BA':1,'CY':1,'AS':1}   

    DisaggSourceFile = wrk + '/metadata/Disaggregation/DisaggSources/UniformSiteDistribution/' + Str0 + Str1 + Str2 
    sids, srcPDF = Disaggregation.GetSourceIDs(DisaggSourceFile,Ntop=Ntop)
    
    tmpPDF = {}; sum1 = 0
    for isid in xrange( len(sids) ):
	tmpPDF[str(sids[isid])] = srcPDF[isid] 
	sum1 = sum1 + srcPDF[isid] 

    print 'Conditional hypocenter distribution:'
    print hypoPDF['CyberShake'] 
    print [[nga, hypoPDF[nga]] for nga in NGAmodels ]

print 'Selected Source IDs: ', sids
sids0, rids = RupSelect( sids, cybershk_database_password, erf_id=erf_id )   # get ruptures given source ID
if len(rids) != 0:
    print 'Source/Rupture sets:'
    for irup in xrange( len(sids) ):
	print 'Source ID: ', sids[irup], ', Rupture (Same Area): ', rids[irup] 
else:
    print 'No ruptures satisfy your criteria of selection'
    raise ValueError

yes_no = 'y'
if yes_no == 'y':
    pass
else:
    sys.exit('No execution!\nChange your rupture-selection criterion to get good ruptures !')

# =====================================
# Initial setting and metafile loading
# =====================================
if opt == 'Preparation':
    flagO = int(sys.argv[3])   # prepare for OpenSHA flatfiles or not (prepare NGA input parameters): 0, 1
    flagD = int(sys.argv[4])   # prepare for Directivity or not (prepare directivity input parameters): 0, 1

    print '='*50
    print 'MetaData Preparation for ABF analysis'
    print '='*50
    
    CN = CyberShakeNGA.cybershk_nga(wrk,cybershk_database_password, rup_model_ids=rup_model_ids, sids=sids,periods=[3.0,],mflag=mflag,bb=bb,NGAs=NGAs,Reference=Reference,ngaModelVersion=modelVersion)
    
    # Extract database and ruptures (controls the first digit in mflag)
    print '='*30
    BlockName = 'Extract Database and Ruptures for CyberShake study (%s,%s,%s,%s)...'%(erf_id,sgt_id,rup_scenario_id,vel_id)
    start_time = HourMinSecToSec(BlockName=BlockName)
    for irup in xrange( len(sids) ):
	sid = sids[irup]
	
	print 'Working on extraction for Source %s...'%sid
	metafile1 = CN.RupsMeta + '%s/'%sid + \
		                    'meta_rup_%s_%s_%s.py'%(CN.erf_id,CN.rup_scenario_id,sid)
	
	for rid in rids[irup]:
	    metafile = CN.CSmeta + '%s/'%sid + \
				    'meta_%s_%s_%s_%s_%s_%s.py'%(CN.erf_id,CN.sgt_id,CN.rup_scenario_id,CN.vel_id,sid,rid)

	    print 'Rupture %s'%rid
	    if not os.path.exists( metafile ) or not os.path.exists( metafile1 ):
		# time consuming part
		meta_rup, meta = CN.extract_database( sid, rid )
		CN.cybershk_sites_rups( sid, meta_rup, meta )   # for GMT plot
		print '# of available sites: ', len( meta.sites_info.keys() )
	    else: 
		meta_rup, meta = CN.extract_database( sid, rid )
		print '# of available sites: ', len( meta.sites_info.keys() )

    end_time = HourMinSecToSec(BlockName='CyberShake Extracting Finished!')
    hour,min,sec = SecToHourMinSec(end_time-start_time,BlockName='CyberShake Extraction')
    print '='*30+'\n'

    # NGA input parameters (M, distance, etc.) from OpenSHA run outputs
    if flagO:
	print '='*30
	BlockName='Extract flatinfo from OpenSHA for rupture sets...'
	start_time = HourMinSecToSec(BlockName=BlockName)

	sids0 = []; rids0 = []
	for irup in xrange( len(sids) ):
	    sid = sids[irup]
	    try: 
		OpenSHA_output = OpenSHA_nga_files(CN.NGAmeta, sid, rids[irup][-1], 3.0, SiteName=None,erf_id=35)
	    except:
		rids0.append( rids[irup] )
		sids0.append( sid )
	    
	if len(sids0) == 0:
	    print 'NGA Flatinfo for all ruptures exists...'
	else:
	    print 'NGA flatinfo extraction for missing ruptures...'
	    CN = CyberShakeNGA.cybershk_nga(wrk,cybershk_database_password,rup_model_ids=rup_model_ids, sids=sids0,periods=[3.0,],mflag=mflag,NGAs=NGAs,Reference=Reference,ngaModelVersion=modelVersion)
	    CN.OpenSHA_nga_cpt([3.0,],rids0,SiteName=None)    # periods dependent (compute at four periods)
	
	end_time = HourMinSecToSec(BlockName='OpenSHA NGA MetaData Extracting Finished')
	hour,min,sec = SecToHourMinSec(end_time-start_time,BlockName='OpenSHA NGA MetaData Extracting')
	print '='*30 + '\n'
    
    if flagD:
	# =======================================================================
	# Get OpenSHA and SC08 (<=> get flatfiles for ngaP and ngaD calculation)
	# =======================================================================
	# Directivity flatinfo using Matlab isochrone theory
	print '='*30
	BlockName = 'Extract flatinfo from SC08 Matlab Once for rupture sets...'
	start_time = HourMinSecToSec(BlockName=BlockName)

	sids0 = []; rids0 = []
	for irup in xrange( len(sids) ):
	    sid = sids[irup]
	    rid = rids[irup][-1]
	    ftmpname = glob.glob(CN.fD_output + 'ERF%s_RupVar%s/SourceID%s_RuptureID%s/hypo*_%s6.txt'%(erf_id,rup_scenario_id,sid,rid,'BA'))
	    if len(ftmpname) == 0:
		sids0.append( sid )
		rids0.append( rids[irup] )
	    else:
		pass
	
	if len(sids0) == 0:
	    print 'Directivity Flatinfo for all ruptures exists...'
	else:
	    print 'Directivity info extraction for missing ruptures...'
	    CN = CyberShakeNGA.cybershk_nga(wrk,cybershk_database_password,rup_model_ids=rup_model_ids, sids=sids0, periods=[3.0,],mflag=mflag,NGAs=NGAs,Reference=Reference,ngaModelVersion=modelVersion)
	    CN.SC08_cpt(rids0)         # just use one period
       
	end_time = HourMinSecToSec(BlockName='SC08 MetaData Extracting Finished')
	hour,min,sec = SecToHourMinSec(end_time-start_time,BlockName='SC08 MetaData Extracting')
	print '='*30+'\n'

# generate fault surface projection (dash lines) (to modify the plots that have fault trace, you can add the dashed downdip lines)
if opt == 'DowndipSurface': 
    RupsMeta = wrk + '/metadata/Ruptures/ERF%s_RupVar%s/'%(erf_id,rup_scenario_id)
    outpth = wrk + '/scripts/map_input/rups_map'
    prefix = 'cybershk.%s.%s'%(erf_id, rup_scenario_id)
    suffix = 'downdip'
    Nsrc = len(sids)
    #Nsrc = 1      # test
    for irup in xrange( Nsrc ):
	sid = sids[irup]
	print 'Working on extraction for Source %s...'%sid
	
	metafile1 = RupsMeta + '%s/'%sid + \
		                    'meta_rup_%s_%s_%s.py'%(erf_id,rup_scenario_id,sid)
        meta_rups = load(metafile1).rups_info 
	faults = meta_rups['fault'] 
	MR = meta_rups['MR'] 
	Ndd, Nas = MR[2:4]    # points defines the fault surface
	faults = load(metafile1).rups_info['fault']
	
	# extract
	faults = np.array(faults).reshape((3,Ndd,Nas))
	bottoms = faults[:,-1,:]
	side1 = faults[:,:,0]
	side2 = faults[:,:,-1]

	# write into files
	fid = open(outpth + '/%s.%s.%s'%(prefix, sid, suffix),'w')
	for ip in xrange(bottoms.shape[1]):
	    fid.write('%s %s\n'%(bottoms[0,ip],bottoms[1,ip]))
	fid.write('>\n') 
	for ip in xrange(side1.shape[1]):
	    fid.write('%s %s\n'%(side1[0,ip],side1[1,ip]))
	fid.write('>\n') 
	for ip in xrange(side2.shape[1]):
	    fid.write('%s %s\n'%(side2[0,ip],side2[1,ip])) 
	fid.close() 


if opt == 'GkxmfsAnalysis':
    # main operation of ABF analysis
    # simple Statistical analysis to connect CyberShake with NGA (uncertainties) 
    NGAcpt = int( sys.argv[3] )     # for multiple periods, just set as 0 (will be recomputed internally)
    
    periods = [3.0,]
    periods = [5.0,10.0]       # For the new rupture generator source there is no calculation for 2.0s SA
    periods = [3.0, 5.0, 10.0]
    periods = [2.0,3.0,5.0,10.0]  # CS11, CS13.2, CS14S4.26 have those 4 periods
    #periods = [3.0, ]   #ERF=36 has 1.0 sec
    print 'Do calculation for periods:', periods
    CN = CyberShakeNGA.cybershk_nga(wrk,cybershk_database_password,rup_model_ids=rup_model_ids, sids=sids,periods = periods, mflag=mflag,NGAs=NGAs,Reference=Reference,ngaModelVersion=modelVersion)
    
    print '='*30
    start_time0 = HourMinSecToSec(BlockName='CyberShake (%s,%s,%s,%s) Statistical Analysis starts'%(erf_id,sgt_id,rup_scenario_id,vel_id))
    print '='*30
    
    print 'Compute CyberShake/NGA ...'
    if NGAcpt: 
	NGAcpt = True
    else: 
	NGAcpt = False 

    MBF_D = False 
    if NGAcpt:
	#CyberShakeRvar, ngaP, ngaD, sites_info, rups_info, Sources = CN.IM_cpt(NGAcpt=NGAcpt)
	MBF_D = True
	CyberShakeRvar, ngaP, ngaD, sites_info, rups_info, Sources = CN.IM_cpt(NGAcpt=NGAcpt, MBF_D=MBF_D)
    else: 
	CyberShakeRvar, sites_info, rups_info, Sources = CN.IM_cpt(NGAcpt=NGAcpt)
	ngaP = None; ngaD = None
    print 'Finish computing CyberShake/NGA ...'

    if MBF_D: 
	sys.exit('Finished')

    Debug = False
    CN.GkxmfsAnalysis( CyberShakeRvar, sites_info, rups_info, Sources, ngaP0k=ngaP, ngaD0k=ngaD, srcPDF=srcPDF, hypoPDF=hypoPDF, Ref=Ref, Debug=Debug )

    print '='*30
    end_time0 = HourMinSecToSec(BlockName='CyberShake Statistical Analysis Ends')
    hour,min,sec = SecToHourMinSec(end_time0-start_time0,BlockName='CyberShake Statistical Analysis')
    print '='*30 + '\n'

if opt == 'PlotAkE0':
    refmodel=sys.argv[3]   # CB, BA, CY, AS, CyberShake, or NoRef
    sigma=sys.argv[4]

    # simple Statistical analysis to connect CyberShake with NGA (uncertainties) 
    periods = [3.0,]
    periods = [3.0,5.0,10.0]
    periods = [2.0,3.0,5.0,10.0]
    CN = CyberShakeNGA.cybershk_nga(wrk,cybershk_database_password,rup_model_ids=rup_model_ids, sids=sids,periods=periods,mflag=mflag,NGAs=NGAs,Reference=Reference,ngaModelVersion=modelVersion)
    CN.PlotAkE0(sigma, Ref=refmodel)

sys.exit('Finished!')

