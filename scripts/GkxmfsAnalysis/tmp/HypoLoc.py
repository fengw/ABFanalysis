#!/usr/bin/env python 
import os, sys
import numpy as np 
from pynga.utils import * 

def To2PI(angle): 
    return (angle+2*np.pi)%(2*np.pi) 

if 1:
    hypoFile = sys.argv[1] 
    rupFile = sys.argv[2] 
    x0 = float(sys.argv[3])
    y0 = float(sys.argv[4])
    reg = sys.argv[5] 
    scale1 = sys.argv[6]
    scale2 = sys.argv[7]
    
    # outputfile 
    Afile = sys.argv[8]
    Ofile = sys.argv[9] 
else: 
    mapin = '/Users/fengw/work/Project/CyberShake_analysis/scripts/map_input/'
    hypoFile = mapin + 'rups_map/cybershk.35.3.90.hypo'
    rupFile = mapin + 'rups_map/cybershk.35.3.90' 
    x0 = 1.0 
    y0 = 1.0 
    scale2 = '2i'
    scale1 = '6i'
    reg = '-R-121/-115/32.5/36.5'
    Afile = mapin + 'arrow_file'
    Ofile = mapin + 'offset_file' 

scale = float( scale1[:-1] ) 
unit = scale1[-1]   # inch


hypoloc = np.loadtxt( hypoFile, usecols=(0,1) ) 
ruploc = np.loadtxt( rupFile ) 

spl = reg.strip().split('/') 
lat0 = float(spl[2]) 
lon0 = float(spl[0][2:]) 
lat1 = float(spl[3]) 
lon1 = float(spl[1] )
dlon0 = lon1-lon0 
dlat0 = lat1-lat0 
if dlon0 > dlat0: 
    lscale = dlon0 
else: 
    lscale = dlat0 

# compute average Strike 
Nh = len(hypoloc)
Nseg = len(ruploc)-1
AveSt = 0.0
for iseg in xrange( Nseg ): 
    loc1 = ruploc[iseg] 
    loc2 = ruploc[iseg+1] 
    r,h,v,az = LonLatToAngleDistance( loc1, loc2, CalcRadius=False, CalcDist=False, Azimuth0to2PI=True) 
    AveSt += az 
AveSt /= Nseg 
AveSt = AveSt * 180. / np.pi 

if 1:
    drot = 90./Nh
    if 0< AveSt < 90: 
	rot = np.arange(Nh) * drot + 180. 
    if 90< AveSt < 180.: 
	rot = np.arange(Nh) * drot
    if 180. < AveSt < 270: 
	rot = np.arange( Nh )[::-1] * drot + 270 
    if 270. < AveSt < 360: 
	rot = np.arange(Nh)[::-1] * drot
if 0: 
    drot = 180/Nh
    if 0< AveSt < 180: 
	rot = np.arange(Nh) * drot + AveSt + 180
    if 180. < AveSt < 360: 
	rot = np.arange( Nh ) * drot + AveSt - 180
if AveSt == 0 or AveSt == 180: 
    rot = np.repeat(90, Nh)
if AveSt == 90 or AveSt == 270: 
    rot = np.repeat(0, Nh)

if 0: 
    if 90<AveSt<180:
	AveSt += 180. 
    rot = np.repeat( AveSt + 90, Nh )

rot = rot*np.pi/180.
rot = To2PI( rot ) 

# write to output files
fid1 = open( Afile, 'w' )
fid2 = open( Ofile, 'w' ) 
fmt1 = '%s %s %s %s\n' 
fmt2 = '-Xa%s%s -Ya%s%s\n' 
in2cm = 2.6  # inch to cm
for ih in xrange( Nh ): 
    hloc = [hypoloc[ih,0],hypoloc[ih,1],0.0]
    az = rot[ih]
    hD = 50   # km or 0.5 degree
    vD = 0
    vector = [az, hD, vD] 
    hloc2 = EndLocation( hloc, vector )   # end of arrow 
   
    # compute arrow direction (relative to horizontal and counter-clockwise in degree) 
    r,h,v,az1 = LonLatToAngleDistance( hloc2, hloc, Azimuth0to2PI=True) 
    AveSt = (2*np.pi-az1+np.pi/2.+2*np.pi)%(2*np.pi) 
    hlon2, hlat2 = hloc2[:2]
    arrowL = r * 180./np.pi * scale / lscale * in2cm
    fid1.write( fmt1%( hlon2, hlat2, AveSt*180./np.pi, arrowL) )

fid1.close() 

# compute offset for each hypocenter 
x = np.zeros(Nh)
y = np.zeros(Nh)
for ih in xrange( Nh ): 

    hD = arrowL = 50   # km or 0.5 degree
    vD = 0
    
    hloc = [hypoloc[ih,0],hypoloc[ih,1],0.0]
    az = rot[ih]
    vector = [az, hD, vD] 
    hloc2 = EndLocation( hloc, vector )   # end of arrow 

    hlon2, hlat2 = hloc2[:2]

    x[ih] = (hlon2-lon0)*scale/ dlon0 #- float(scale2[:-1])/4.
    y[ih] = (hlat2-lat0)*scale/ dlat0 #- float(scale2[:-1])/4.

    fid2.write( fmt2%(x[ih],unit,y[ih],unit) )
fid2.close() 


