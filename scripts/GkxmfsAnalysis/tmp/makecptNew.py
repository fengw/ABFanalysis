#!/usr/bin/env python
"""
Make cpt for GMT grdimages (simple one)
You can extend its functionality
"""
import numpy as np


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

	'bwwr': [ (0,1,1,1), 
	          (0,1,1,0),
		  (1,1,1,0), 
	       ], 
	
	
	'bwr0': [ # reference (this one is better)
	          # you should know if adding 
	          (0,0,0.5),   
		  (0,0.25,0.75),
		  (1,1,1),  # white (doesn't change)
		  (1,0,0),
		  (0.5,0,0),   
		  ], 
	
	'bwr': [ 
	          (0,0,0.5),   
		  (0,0,0.75),
		  (0,0,1),
		  (1,1,1),  # white (doesn't change)
		  (1,0,0),
		  (0.75,0,0),
		  (0.5,0,0),   
		  ], 

	}



def MakeGMTcpt(minv,maxv,Np,cmapname,cv=0,sigma=0.01,cptfile=None,logscale=0,invert=0):
    
    if minv >= maxv:
	print 'min value must be less then max value...'
	raise ValueError

    # make GMT cpt files from cmap library 
    cmap = colormaps[cmapname] 
    cmap = np.array(cmap,'f').T
    if invert: 
	cmap = cmap[:,::-1]

    # normalize
    cmap /= max( 1.0, cmap.max() )
    print cmap 

    # assign to colors
    r,g,b = cmap 
    Nc = cmap.shape[1]   # number of known colors

    # no central confidences
    v0 = np.linspace( minv, maxv, Nc )   # known stuff
    x = np.linspace( 0.0, 1.0, Nc )
    xi = np.linspace( 0.0, 1.0, Np ) 
    
    # interpolate colormap
    r = np.interp( xi, x, r )
    g = np.interp( xi, x, g )
    b = np.interp( xi, x, b )

    sigma = 0.01    # in loge space, the significance of the effects
    if not logscale: 
	v = np.interp( xi, x, v0 )   # linear interpolation 
    else: 
	e = np.exp(1)
	if minv > 0:
	    v = np.logspace( np.log(minv), np.log(maxv), Np, base = e ) 
	elif minv == 0: 
	    # 0 in [minv,maxv]extend the Np
	    v = np.r_[minv, np.logspace( np.log(sigma), np.log(maxv), Np-1, base=e)]
	elif minv < 0 and maxv > 0:
	    dv = (maxv-minv)/(Np-1)
	    N1 = (-sigma-minv)/dv + 1
	    N2 = (maxv-sigma)/dv + 1
	    v = np.r_[ \
		    -np.logspace(np.log(sigma),np.log(abs(minv)),N1, base=e)[::-1], 0, np.logspace(np.log(sigma),np.log(maxv),N2,base=e) ]
	    
	    # N1 + N2 + 1 = Np; N1, N2 depends on points in each side
	    # ...
	elif maxv == 0: 
	    v = np.r_[ -np.logspace(np.log(sigma),np.log(abs(minv)),Np-1,base=e)[::-1], 0]
	elif maxv < 0: 
	    v = -np.logspace( np.log(abs(maxv)), np.log(abs(minv)), base=e )
    
    cmap = r,g,b 

    fmt = '%-5.3f%6.0f%6.0f%6.0f    '
    fmt = 2*fmt
    fmt = fmt+'\n'

    if cptfile == None:
	# test GMT formats without writing into cpt file
	cpt = ''
	for i in xrange( Np-1 ):
	    cpt += fmt % (
		    v[i], 255*r[i], 255*g[i], 255*b[i],
		    v[i+1], 255*r[i+1], 255*g[i], 255*b[i],         # discontinous case
		    #v[i+1], 255*r[i+1], 255*g[i+], 255*b[i+1],     # continuous case
		    )
        print '# COLOR_MODEL = RGB'	
        print cpt
	print 'B  0  0  0\nF  255  255  255\nN  128  128  128\n'

    else:
	# write into file
	fid = open( cptfile, 'w' )
	head = '# COLOR_MODEL = RGB\n'
	fid.write(head)
	for i in xrange( Np-1 ):
	    fid.write( fmt % (
		    v[i], 255*r[i], 255*g[i], 255*b[i],
		    v[i+1], 255*r[i+1], 255*g[i], 255*b[i],
		    ))
	fid.write('B  0  0  0\nF  255  255  255\nN  128  128  128\n')
	fid.close()


if __name__ == '__main__':
    import os, sys 

    maxv = float(sys.argv[1]) # maximum value
    minv = float(sys.argv[2]) # minimum value
    Np = int(sys.argv[3] )  # points of number (odd or even depends on minv and maxv)
    cmapname = sys.argv[4]    # bwr, etc name from colormap library
    cptfile = sys.argv[5]
    invert = int(sys.argv[6])    # reverse the colormap 
    logscale = int(sys.argv[7]) 

    if cptfile == 'None': 
	cptfile = None
    
    MakeGMTcpt(minv,maxv,Np,cmapname,cptfile=cptfile,logscale=logscale,invert=invert)

