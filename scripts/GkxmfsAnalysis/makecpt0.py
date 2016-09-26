#!/usr/bin/env python
"""
Make cpt for GMT grdimages (simple one)
You can extend its functionality
"""
import numpy as np
from my_util.indpy import bisearch


# Colormap Library
colormaps = { 
	# define color knots
	#     red, green, blue value for each knots (from small value to large value)
	'byr': [ (0,0,1),
	         (0,0,2),
		 (0,2,2),
		 (2,2,2),
		 (2,2,0),
		 (2,0,0),
		 (1,0,0),
	       ], 

	'b2r': [ 
	          (0,0,0.5),   
		  (0,0.25,0.9),
		  (1,1,1),  # white (doesn't change)
		  (0.9,0,0),
		  (0.5,0,0),   
		  ], 
	
	'bwr0': [ # reference (this one is better)
	          # you should know if adding 
	          (0,0,0.5),   
		  (0,0.25,0.75),
		  (1,1,1),  # white (doesn't change)
		  (0.9,0,0),
		  (0.5,0,0),   
		  ], 
	
	'bwr': [ 
	          (0,0,0.5),   
		  (0,0.25,0.9),
		  (1,1,1),  # white (doesn't change)
		  (1,0,0),
		  (0.5,0,0),   
		  ], 

	}


def MakeGMTcpt( minv, maxv, Np, cmapname, cv, sigma, cptfile=None, invert=0, logscale=0):

    if minv >= maxv:
	print 'min value must be less then max value...'
	raise ValueError

    # make GMT cpt files from cmap library 
    cmap = colormaps[cmapname] 
    Nc0 = len(cmap) 
    if sigma == 0: 
	cmap = np.array(cmap,'f').T
	if invert: 
	    cmap = cmap[:,::-1]

	# normalize
	cmap /= max( 1.0, cmap.max() )
	
	r,g,b = cmap 

	Nc0 = cmap.shape[1]   # number of known colors
	
	# no central confidences
	v0 = np.linspace( minv, maxv, Nc0 )   # known stuff
	x = np.linspace( 0.0, 1.0, Nc0 )
	xi = np.linspace( 0.0, 1.0, Np ) 
	
	# interpolate colormap
	r = np.interp( xi, x, r )
	g = np.interp( xi, x, g )
	b = np.interp( xi, x, b )

	v = np.interp( xi, x, v0 )   # linear interpolation 
        Nc = Np-1 

    else: 
	# normalize sigma (in 0,1) for grid 
	sigma0 = (sigma-minv)/(maxv-minv) - 0.5 
	sigma0 = abs(sigma0)

	Nc0 = len(cmap) + 2
        
	# correct
	v01 = np.linspace( minv, cv-sigma, (Nc0-1)/2 )
	v02 = np.linspace( cv+sigma, maxv, (Nc0-1)/2 )
	x01 = np.linspace( 0.0, 0.5-sigma0, (Nc0-1)/2 )
	x02 = np.linspace( 0.5+sigma0, 1.0, (Nc0-1)/2 )
	x = np.r_[x01, 0.5, x02] 
	v0 = np.r_[v01, cv, v02]
	
	xi01 = np.linspace( 0.0, 0.5-sigma0, Np )
	xi02 = np.linspace( 0.5+sigma0, 1.0, Np )
        xi = np.r_[xi01, 0.5, xi02] 
	v = np.interp( xi, x, v0 )   # linear interpolation 
        
	Nc0 = len(cmap)
	index = (Nc0-1)/2
	cmap.insert( index, cmap[index] )
	cmap = np.array(cmap,'f').T
	if invert: 
	    cmap = cmap[:,::-1]

	# normalize
	cmap /= max( 1.0, cmap.max() )
	r,g,b = cmap
	Nc0 = cmap.shape[1] 
	x01 = np.linspace( 0.0, 0.5-sigma0, Nc0/2 )
	x02 = np.linspace( 0.5, 1.0, Nc0/2 )
	x = np.r_[x01, x02] 
	r = np.interp( xi, x, r )
	g = np.interp( xi, x, g )
	b = np.interp( xi, x, b )

	Nc = 2*Np 

    if logscale: 
	v = np.exp(v)  # actural ratio

    fmt = '%-5.3f%6.0f%6.0f%6.0f    '
    fmt = '%-5.3f%6.0f%6.0f%6.0f    '
    fmt = 2*fmt
    fmt = fmt+'\n'
    #fmt = '%-10f %3.0f %3.0f %3.0f     %-10f %f %3.0f %3.0f\n'
    if cptfile == None:
	# test formats
	cpt = ''
	for i in xrange( Nc ):
	    cpt += fmt % (
		    v[i], 255*r[i], 255*g[i], 255*b[i],
		    #v[i+1], 255*r[i+1], 255*g[i+1], 255*b[i+1],  # continuous colrmap
		    v[i+1], 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
		    )
	print cpt
    else:
	# write into file
	fid = open( cptfile, 'w' )
	head = '# COLOR_MODEL = RGB\n'
	fid.write(head)
	for i in xrange( Nc ):
	    fid.write( fmt % (
		    v[i], 255*r[i], 255*g[i], 255*b[i],
		    v[i+1], 255*r[i], 255*g[i], 255*b[i],
		    ))
	fid.write('B  0  0  0\nF  255  255  255\nN 128 128 128 \n')
	fid.close()



if __name__ == '__main__':
    import os, sys 
    maxv = float(sys.argv[1]) # maximum value
    minv = float(sys.argv[2]) # minimum value
    cv = float( sys.argv[3] ) # central value
    Np = int(sys.argv[4] )  # points of number (odd or even depends on minv and maxv)
    sigma = float( sys.argv[5] )
    cmapname = sys.argv[6]    # bwr, etc name from colormap library
    cptfile = sys.argv[7]
    logscale = int(sys.argv[8]) 

    if cptfile == 'None': 
	cptfile = None
    
    MakeGMTcpt(minv,maxv,Np,cmapname,cv, sigma, cptfile=cptfile, logscale=logscale)


