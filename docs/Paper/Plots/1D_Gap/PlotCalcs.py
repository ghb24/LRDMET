#!/usr/bin/python

import os
from numpy import *
import sys
from pylab import *

params = {'legend.fontsize': 11}
rcParams.update(params)

nImp4_Zero=0.012
nImp2_DDZero = 0.0122719 
nImp4_DDZero = 0.01198
nImp6_Zero = 0.01232

data=mlab.load('GAPS',usecols=[0,1,2,3,4,5,6,7],unpack=True)
data_BA=mlab.load('analyticgap_g100000_L100000.dat',usecols=[0,1],unpack=True)
data_CDMFT_26=mlab.load('CDMFT/data_gap/gap_26.dat',usecols=[0,1,2],unpack=True)
data_CDMFT_28=mlab.load('CDMFT/data_gap/gap_28.dat',usecols=[0,1,2],unpack=True)
data_CDMFT_46=mlab.load('CDMFT/data_gap/gap_46.dat',usecols=[0,1,2],unpack=True)
data_CDMFT_48=mlab.load('CDMFT/data_gap/gap_48.dat',usecols=[0,1,2],unpack=True)
data_CDMFT_66=mlab.load('CDMFT/data_gap/gap_66.dat',usecols=[0,1,2],unpack=True)
clf()
ax1 = subplot(111)
#suptitle('1D hubbard spectral gap (Lattice = 624 sites)')
ylabel('Error in Spectral gap / t',fontsize='14')
xlabel('U / t',fontsize='14')
#ax1.plot(data_BA[0],data_BA[1],linewidth=2,label='Bethe Ansatze',color='b')
##ax1.plot(data[0],data[1],linewidth=1,label='2-site DMET (Reopt GS)',color='g')
#ax1.plot(data[0],data[4]-nImp2_Opt_Zero,linewidth=1,marker='o',label='2-site DMET: Neutral excitation gap',color='k')
#ax1.plot(data[0],data[1],linewidth=1,marker='+',label='2-site DMET: Quasiparticle gap',color='g')
#ax1.plot(data[0],data[3]-nImp4_Zero,linewidth=1,marker='x',label='4-site DMET: Quasiparticle gap',color='r')
#ax1.plot(data_CDMFT_26[0],data_CDMFT_26[1],linewidth=1,label='CDMFT(2,6)',marker='s')
##ax1.plot(data_CDMFT_28[0],data_CDMFT_28[1],linewidth=1,label='CDMFT(2,8)',marker='d')
#ax1.plot(data_CDMFT_46[0],data_CDMFT_46[1],linewidth=1,label='CDMFT(4,6)',marker='d')
##ax1.plot(data_CDMFT_48[0],data_CDMFT_48[1],linewidth=1,label='CDMFT(4,8)')
##ax1.plot(data_CDMFT_66[0],data_CDMFT_66[1],linewidth=1,label='CDMFT(6,6)')

#ax1.plot(data_BA[0],data_BA[1],linewidth=2,label='Bethe Ansatze',color='b')
#ax1.plot(data[0],data[1],linewidth=1,label='2-site DMET (Reopt GS)',color='g')
ax1.plot(data[0],data[1]-data[7],linewidth=2,marker='s',label='2-site DMET: Quasiparticle gap')#,color='g')
ax1.plot(data[0],(data[4]-nImp2_DDZero)-data[7],linewidth=2,marker='o',label='2-site DMET: Neutral gap')#,color='k')
ax1.plot(data[0],(data[3]-nImp4_Zero)-data[7],linewidth=2,marker='s',label='4-site DMET: Quasiparticle gap')#,color='r')
ax1.plot(data[0],(data[5]-nImp4_DDZero)-data[7],linewidth=2,marker='s',label='4-site DMET: Neutral gap')#,color='r')
ax1.plot(data[0],(data[6]-nImp6_Zero)-data[7],linewidth=2,marker='s',label='6-site DMET: Quasiparticle gap')#,color='r')
ax1.plot(data_CDMFT_26[0],data_CDMFT_26[1]-data_CDMFT_26[2],linewidth=1,label='2-site, 6-bath CDMFT',marker='x')
#ax1.plot(data_CDMFT_28[0],data_CDMFT_28[1],linewidth=1,label='CDMFT(2,8)',marker='d')
ax1.plot(data_CDMFT_46[0],data_CDMFT_46[1]-data_CDMFT_46[2],linewidth=1,label='4-site, 6-bath CDMFT',marker='+')
ax1.plot(data_CDMFT_66[0],data_CDMFT_66[1]-data_CDMFT_66[2],linewidth=1,marker='d',mfc='none',color='y',mec='y',label='6-site, 6-bath CDMFT')
#ax1.plot(data_CDMFT_48[0],data_CDMFT_48[1]-data_CDMFT_48[2],linewidth=1,marker='d',mfc='none',mec='r',label='4-site, 8-bath CDMFT')
ax1.axhline(0.0,color='k')
#ax1.plot(data_CDMFT_48[0],data_CDMFT_48[1],linewidth=1,label='CDMFT(4,8)')
#ax1.plot(data_CDMFT_66[0],data_CDMFT_66[1],linewidth=1,label='CDMFT(6,6)')
xlim(0,10)
ylim(-0.5,0.3)
legend(loc=0,ncol=2)
show()
#savefig('Hubbard_Gap.eps')
