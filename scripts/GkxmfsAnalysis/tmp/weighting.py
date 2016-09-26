#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

Nh = 30
xh = np.arange( Nh )
mu = Nh/2. 
mu0 = mu
std = np.std(xh)
s = Nh/2.

prob = 1./(2.*np.pi)**0.5/std * np.exp( -0.5*((xh-mu)/std)**2)
prob0 = 0.5/s *(1+np.cos((xh-mu0)/s*np.pi))

print sum(prob),sum(prob0)
fig = plt.figure(1)
ax = fig.add_subplot(111)

if 0:
    line1 = ax.plot( xh, prob, 'r' )
    line2 = ax.plot( xh, prob0, 'b' )
    plt.legend( [line1,line2],['Gaussian','Cosine'],loc=0)
else:
    ax.plot( xh, prob0 )
    ax.set_xlabel( 'x' )
    ax.set_ylabel( 'w(x)' )
    ax.set_title( 'Cosine weighting function')
plt.show()

