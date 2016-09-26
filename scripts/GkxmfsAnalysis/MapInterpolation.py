#!/usr/bin/env python 
""" 
Input uninterpolated maps and do the interpolation based on the SiteInterpMesh 
""" 
import os, sys 
import numpy as np 
from my_util.numer import interp 

# interpolation (used in GMT)
infile = sys.argv[1]    # file full path with name 
outfile = sys.argv[2] 

inputs = np.loadtxt( infile ) 
SiteLon = inputs[:,0] 
SiteLat = inputs[:,1] 
Value = inputs[:,2] 

site_locs1 = './SiteInterpMesh.txt' 
sites = np.loadtxt( site_locs1 ) 
SiteLon1D = sites[:,0] 
SiteLat1D = sites[:,1] 

eps = 0.001 
smooth = 0.01
method = {'name':'exp','smooth':smooth}
newValue = interp( SiteLon, SiteLat, Value, SiteLon1D, SiteLat1D, eps=eps, method=method ) 

fid = open( outfile, 'w' ) 
for i in xrange( len(SiteLon1D) ): 
    fid.write( '%s %s %s\n'%(SiteLon1D[i], SiteLat1D[i], newValue[i]) )
fid.close() 

