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

ax1 = subplot(1,1,1,title='4-impurity Honeycomb DoS')
data2=mlab.load('EC-TDA_GFResponse_2.00',usecols=[0,2],unpack=True)
data4=mlab.load('EC-TDA_GFResponse_4.00',usecols=[0,2],unpack=True)
data6=mlab.load('EC-TDA_GFResponse_6.00',usecols=[0,2],unpack=True)
data8=mlab.load('EC-TDA_GFResponse_8.00',usecols=[0,2],unpack=True)
data10=mlab.load('EC-TDA_GFResponse_10.00',usecols=[0,2],unpack=True)

#ax1.plot(data2[0],data2[1],linewidth=2,label='2-site',color='b')
ax1.plot(data2[0],data2[1],linewidth=2,label='U=2')
ax1.plot(data4[0],data4[1],linewidth=2,label='U=4')
ax1.plot(data6[0],data6[1],linewidth=2,label='U=6')
ax1.plot(data8[0],data8[1],linewidth=2,label='U=8')
ax1.plot(data10[0],data10[1],linewidth=2,label='U=10')
ax1.legend(loc=0)
xlim(-3.8,3.8)
xlabel(r'$\omega$ / t')
ylabel('Local Density of States')
#ylim(0,2)
#show()
savefig('Honeycomb.eps')
