#!/usr/bin/python

#plot "./nImp4/U_4/EC-TDA_GFResponse_4.00" u 1:3 w l, "./nImp2/NoReopt/EC-TDA_GFResponse_4.00" u 1:3 w l, "./CDMFT/ldos_u4.00_ep0.05.dat" u 1:2 w l

import os
from numpy import *
import sys
from pylab import *

#for U in arange(0.00,6.01,0.5):
#
#    print 'U: ',U
#    print 'Filename: EC-TDA_DDResponse_%g' % str(U)


#params = {'legend.fontsize': 8}
#rcParams.update(params)

ax1 = subplot(1,1,1)
data=mlab.load('nImp4/Hybrid',usecols=[0,2],unpack=True)
data2=mlab.load('nImp2/EC-TDA_GFResponse_4.00',usecols=[0,2],unpack=True)

#ax1.plot(data2[0],data2[1],linewidth=2,label='2-site',color='b')
ax1.plot(data[0],data[1],linewidth=2,label='4-site',color='b')
ax1.legend(loc=0)
xlim(-4.6,4.6)
xlabel(r'$\omega$ / t')
ylabel('Density of States')
#ylim(0,2)
#show()
savefig('Honeycomb.eps')
