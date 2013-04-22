#!/usr/bin/python

import os
from numpy import *
import sys
from pylab import *

params = {'legend.fontsize': 14}
rcParams.update(params)

#fac2 = 1.8 
#fac4 = 1.3 
#fac6 = 1.3 
#fac8 = 1.3 

ax1 = subplot(2,2,1) #,title='1D 624 site Hubbard (APBC) DoS')
#data_DMFT=mlab.load('CDMFT/ldos_u2.00_ep0.05.dat',usecols=[0,1],unpack=True)
#data=mlab.load('nImp2/Reopt/EC-TDA_DDResponse_2.00',usecols=[0,2],unpack=True)
data=mlab.load('nImp2/EC-TDA_DDResponse_2.00',usecols=[0,2],unpack=True)
ax1.plot(data[0],data[1],linewidth=2,label='2-site DMET',color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax1.annotate('U = 2',xy=(-15,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax1.get_yticklabels(),visible=False)
#ax1.axhline(0)
#ax1.yaxis.tick_left()
xticks( (0, 1, 2, 3, 4), ('0','1','2','3','4') )
#yticks( (0, 1, 2), ('0','1','2') )
xlim(0,3.38)
ylim(0,1.0)

ax2 = subplot(2,2,2) #,title='1D 624 site Hubbard (APBC) DoS')
#data=mlab.load('nImp2/Reopt/EC-TDA_DDResponse_4.00',usecols=[0,2],unpack=True)
data=mlab.load('nImp2/EC-TDA_DDResponse_4.00',usecols=[0,2],unpack=True)
ax2.plot(data[0],data[1],linewidth=2,color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax2.annotate('U = 4',xy=(-25,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax2.get_yticklabels(),visible=False)
xticks( (0, 1, 2, 3, 4), ('0','1','2','3','4') )
#ax2.yaxis.tick_right()
#yticks( (0, 1, 2, 3), ('0','1','2','3') )

xlim(0,4.5)
ylim(0,1.0)
#ax2.axhline(0)


ax3 = subplot(2,2,3) #,title='1D 624 site Hubbard (APBC) DoS')
#data=mlab.load('nImp2/Reopt/EC-TDA_DDResponse_6.00',usecols=[0,2],unpack=True)
data=mlab.load('nImp2/EC-TDA_DDResponse_6.00',usecols=[0,2],unpack=True)
ax3.plot(data[0],data[1],linewidth=2,color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax3.annotate('U = 6',xy=(-25,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax3.get_yticklabels(),visible=False)
#ax3.yaxis.tick_left()
#yticks( (0, 1, 2), ('0','1','2') )
xlim(0,6.1)
ylim(0,0.5)
#ax3.axhline(0)
xlabel(r'$\omega$ / t')


ax4 = subplot(2,2,4) #,title='1D 624 site Hubbard (APBC) DoS')
#data=mlab.load('nImp2/Reopt/EC-TDA_DDResponse_8.00',usecols=[0,2],unpack=True)
data=mlab.load('nImp2/EC-TDA_DDResponse_8.00',usecols=[0,2],unpack=True)
ax4.plot(data[0],data[1],linewidth=2,color='b')
#ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
ax4.annotate('U = 8',xy=(-25,-9),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='15')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
setp(ax4.get_yticklabels(),visible=False)
#ax4.axhline(0)
#ax4.yaxis.tick_right()
#yticks( (0, 1), ('0','1') )
xlim(0,7.85)
ylim(0,0.3)
xlabel(r'$\omega$ / t')

subplots_adjust(hspace=0.17,wspace=0.04)

#show()
savefig('1D_Hub_DD.eps')

#clf()

