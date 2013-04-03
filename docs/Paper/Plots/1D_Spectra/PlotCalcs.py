#!/usr/bin/python

import os
from numpy import *
import sys
from pylab import *

params = {'legend.fontsize': 8}
rcParams.update(params)

fac2 = 1.8 
fac4 = 1.3 
fac6 = 1.3 
fac8 = 1.3 

ax1 = subplot(2,2,1) #,title='1D 624 site Hubbard (APBC) DoS')
data_DMFT=mlab.load('CDMFT/ldos_u2.00_ep0.05.dat',usecols=[0,1],unpack=True)
data=mlab.load('nImp4/EC-TDA_GFResponse_2.00',usecols=[0,2],unpack=True)
ax1.plot(data_DMFT[0],data_DMFT[1]/fac2,linewidth=1,label='CDMFT (6,6)',color='r')
#ax1.plot(data_NoReopt2[0],data_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='r')
ax1.plot(data[0],data[1],linewidth=2,label='4-site DMET',color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax1.annotate('U = 2',xy=(-2,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax1.get_yticklabels(),visible=False)
xlim(-3,3)

ax2 = subplot(2,2,2) #,title='1D 624 site Hubbard (APBC) DoS')
data_DMFT=mlab.load('CDMFT/ldos_u4.00_ep0.05.dat',usecols=[0,1],unpack=True)
data=mlab.load('nImp4/EC-TDA_GFResponse_4.00',usecols=[0,2],unpack=True)
ax2.plot(data_DMFT[0],data_DMFT[1]/fac4,linewidth=1,color='r')
#ax1.plot(data_NoReopt2[0],data_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='r')
ax2.plot(data[0],data[1],linewidth=2,color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax2.annotate('U = 4',xy=(-2,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax2.get_yticklabels(),visible=False)
xlim(-2.8,2.8)


ax3 = subplot(2,2,3) #,title='1D 624 site Hubbard (APBC) DoS')
data_DMFT=mlab.load('CDMFT/ldos_u6.00_ep0.05.dat',usecols=[0,1],unpack=True)
data=mlab.load('nImp4/EC-TDA_GFResponse_6.00',usecols=[0,2],unpack=True)
ax3.plot(data_DMFT[0],data_DMFT[1]/fac6,linewidth=1,color='r')
#ax1.plot(data_NoReopt2[0],data_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='r')
ax3.plot(data[0],data[1],linewidth=2,color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax3.annotate('U = 6',xy=(-2,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax3.get_yticklabels(),visible=False)
xlim(-3.52,3.52)
xlabel(r'$\omega$ / t')


ax4 = subplot(2,2,4) #,title='1D 624 site Hubbard (APBC) DoS')
data_DMFT=mlab.load('CDMFT/ldos_u8.00_ep0.05.dat',usecols=[0,1],unpack=True)
data=mlab.load('nImp4/EC-TDA_GFResponse_8.00',usecols=[0,2],unpack=True)
ax4.plot(data_DMFT[0],data_DMFT[1]/fac8,linewidth=1,color='r')
#ax1.plot(data_NoReopt2[0],data_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='r')
ax4.plot(data[0],data[1],linewidth=2,color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax4.annotate('U = 8',xy=(-2,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax4.get_yticklabels(),visible=False)
xlim(-4.36,4.36)
xlabel(r'$\omega$ / t')

subplots_adjust(hspace=0.17,wspace=0.015)

#show()
savefig('1D_Hub_Spectra.eps')

#clf()

