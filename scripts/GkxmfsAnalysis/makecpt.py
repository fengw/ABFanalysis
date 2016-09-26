#!/usr/bin/env python
"""
Make cpt for GMT grdimages (simple one)
You can extend its functionality
"""
import numpy as np

# Colormap Library
colormaps = { 
	# define color knots
	#     red, green, blue value for each knots in the colorbar (from small value to large value)
	
	# For symmetric plots
	'bwr': [ 
	          (0,0,0.5),   
		  #(0,0.25,0.75),
		  (0, 0.25, 0.9), 
		  (1,1,1),  # white (doesn't change)
		  (0.9,0,0),
		  #(1,0,0),
		  (0.5,0,0),   
		  ], 
        # retrieve the simplified colorscale from GMT globe cpt 
	'wrb': [
	    (1,1,1),
	    (1,0,0),
	    (0,0.25,0.75),
	    (0,0,0.5),
	    ]

	}


def MakeGMTcpt( minv, maxv, Np, cmapname, cv, sigma, cptfile=None, invert=0, logscale=0):

    if minv >= maxv:
	print 'min value must be less then max value...'
	raise ValueError

    # make GMT cpt files from cmap library 
    cmap = colormaps[cmapname] 
    Nc0 = len(cmap) 
    
    # GMT cpt file format
    fmt = '%-5.3f%6.0f%6.0f%6.0f    '
    fmt = 2*fmt
    fmt = fmt+'\n'
    
    if cv != 'None':
	# case1  (negative to positive, symmetry about cv) (sigma == 0 or a very small value)
	if abs(maxv-cv) != abs(minv-cv): 
	    print 'central value %s should be in the middle between the interval [%s,%s]'%(cv,minv,maxv)
	    raise ValueError 

	if not np.mod( Np, 2 ): 
	    print '# points need to be odd number'
	    raise ValueError 
	
	#if sigma > (maxv-minv)/500: 
	#    sigma = (maxv-minv)/500 

	cmap = np.array(cmap,'f').T
	if invert: 
	    cmap = cmap[:,::-1]

	# normalize
	cmap /= max( 1.0, cmap.max() )
	
	r,g,b = cmap 

	Nc0 = cmap.shape[1]   # number of known colors
	
	v0 = np.linspace( minv, maxv, Nc0 )   # known stuff
	x = np.linspace( 0.0, 1.0, Nc0 )
	xi = np.linspace( 0.0, 1.0, Np ) 

	# interpolate colormap
	r = np.interp( xi, x, r )
	g = np.interp( xi, x, g )
	b = np.interp( xi, x, b )

	v = np.interp( xi, x, v0 )   # linear interpolation 
	Nc = Np
	
	if logscale: 
	    v = np.exp(v)  # actural ratio

	if cptfile == None:
	    # test formats and print
	    cpt = ''
	    for i in xrange( Nc ):
		if i < Nc/2: 
		    cpt += fmt % (
				v[i], 255*r[i], 255*g[i], 255*b[i],
				v[i+1], 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
				)
		elif i > Nc/2:
		    cpt += fmt % (
			    v[i-1], 255*r[i], 255*g[i], 255*b[i],
			    v[i], 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
			    )
		else:
		    # use the centeral color for the central value 
		    cpt += fmt % (
			    v[i], 255*r[i], 255*g[i], 255*b[i],
			    v[i]+sigma, 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
			    )
	    print cpt
	else:
	    # write into file
	    fid = open( cptfile, 'w' )
	    head = '# COLOR_MODEL = RGB\n'
	    fid.write(head)
	    for i in xrange( Nc ):
		if i < Nc/2: 
		    if i == Nc/2-1:
			fid.write( fmt % (
				    v[i], 255*r[i], 255*g[i], 255*b[i],
				    v[i+1]-sigma, 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
				    ))
			v[i+1] -= sigma
		    else: 
			fid.write( fmt % (
				    v[i], 255*r[i], 255*g[i], 255*b[i],
				    v[i+1], 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
				    ))

		elif i > Nc/2:
		    fid.write( fmt % (
			    v[i-1], 255*r[i], 255*g[i], 255*b[i],
			    v[i], 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
			    ) )
		else:
		    fid.write( fmt % (
			    v[i], 255*r[i], 255*g[i], 255*b[i],
			    v[i]+sigma, 255*r[i], 255*g[i], 255*b[i],      #  discretized colrmap
			    ) )
		    v[i] += sigma 

	    fid.write('B  0  0  0\nF  255  255  255\nN  128  128  128\n')
	    fid.close()
    
    else: 
	# case2: regular colormap (no symmetry, for example [0,1] for the D-map sigma)
	cmap = np.array(cmap,'f').T
	if invert: 
	    cmap = cmap[:,::-1]

	# normalize
	cmap /= max( 1.0, cmap.max() )
	
	r,g,b = cmap 

	Nc0 = cmap.shape[1]   # number of known colors
	
	v0 = np.linspace( minv, maxv, Nc0 )   # known stuff
	x = np.linspace( 0.0, 1.0, Nc0 )
	xi = np.linspace( 0.0, 1.0, Np ) 

	# interpolate colormap
	r = np.interp( xi, x, r )
	g = np.interp( xi, x, g )
	b = np.interp( xi, x, b )

	v = np.interp( xi, x, v0 )   # linear interpolation 
	Nc = Np-1
	
	if logscale: 
	    v = np.exp(v)  # actural ratio

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
	    fid.write('B  0  0  0\nF  255  255  255\nN  128  128  128\n')
	    fid.close()



if __name__ == '__main__':
    import os, sys 
    maxv = float(sys.argv[1]) # maximum value
    minv = float(sys.argv[2]) # minimum value
    cv = sys.argv[3]  # central value 
    if cv != 'None':
	cv = float(cv) 

    Np = int(sys.argv[4] )  # points of number (odd or even depends on minv and maxv)
    sigma = float( sys.argv[5] )
    cmapname = sys.argv[6]    # bwr, etc name from colormap library
    cptfile = sys.argv[7]
    logscale = int(sys.argv[8]) 

    if cptfile == 'None': 
	cptfile = None
    
    MakeGMTcpt(minv,maxv,Np,cmapname,cv, sigma, cptfile=cptfile, logscale=logscale)


