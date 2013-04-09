#!/usr/bin/python

import os
from numpy import *
import sys
from pylab import *

params = {'legend.fontsize': 14}
rcParams.update(params)

nImp4_Zero=0.012
nImp2_Opt_Zero = 0.0122719 

data=mlab.load('GAPS',usecols=[0,1,2,3,4],unpack=True)
data_BA=mlab.load('analyticgap_g100000_L100000.dat',usecols=[0,1],unpack=True)
clf()
ax1 = subplot(111)
#suptitle('1D hubbard spectral gap (Lattice = 624 sites)')
ylabel('Spectral gap / t',fontsize='14')
xlabel('U / t',fontsize='14')
ax1.plot(data_BA[0],data_BA[1],linewidth=2,label='Bethe Ansatze',color='b')
#ax1.plot(data[0],data[1],linewidth=1,label='2-site DMET (Reopt GS)',color='g')
ax1.plot(data[0],data[1],linewidth=1,marker='+',label='2-site DMET',color='g')
ax1.plot(data[0],data[3]-nImp4_Zero,linewidth=1,marker='x',label='4-site DMET',color='r')
ax1.plot(data[0],data[4]-nImp2_Opt_Zero,linewidth=1,marker='x',label='2-site DMET: 2-particle spectrum',color='k')
xlim(0,10)
ylim(0,7.0)
legend(loc=0)
show()
#savefig('Hubbard_Gap.eps')
