#!/usr/bin/python

import os
#import matplotlib
from numpy import *
import sys
from pylab import *

#for U in arange(0.00,6.01,0.5):
#
#    print 'U: ',U
#    print 'Filename: EC-TDA_DDResponse_%g' % str(U)


params = {'legend.fontsize': 8}
rcParams.update(params)

fac = 1.91

#ax1 = subplot(5,1,1,title='1D 624 site Hubbard (APBC) DoS')
#data_DMFT=mlab.load('CDMFT/ldos_u0.00_ep0.05.dat',usecols=[0,1],unpack=True)
#data_NoReopt2=mlab.load('nImp2/NoReopt/EC-TDA_GFResponse_0.00',usecols=[0,2],unpack=True)
#data_Reopt2=mlab.load('nImp2/Reopt/EC-TDA_GFResponse_0.00',usecols=[0,2],unpack=True)
##data_NoReopt3=mlab.load('nImp3/NoReopt/Refresh/EC-TDA_GFResponse_0.00',usecols=[0,2],unpack=True)
##data_Reopt3=mlab.load('nImp3/Reopt/EC-TDA_GFResponse_0.00',usecols=[0,2],unpack=True)
#ax1.plot(data_DMFT[0],data_DMFT[1]/fac,linewidth=2,label='CDMFT (6,6)',color='b')
##ax1.plot(data_NoReopt2[0],data_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='r')
#ax1.plot(data_Reopt2[0],data_Reopt2[1],linewidth=1,label='2-site',color='r')
##ax1.plot(data_NoReopt3[0],data_NoReopt3[1],linewidth=1,linestyle='dashed',label='3-site',color='g')
##ax1.plot(data_Reopt3[0],data_Reopt3[1],linewidth=1,label='3-site ReoptGS',color='g')
#ax1.annotate('U = 0',xy=(-5,26),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='16')
#ax1.legend(loc=2)
#setp(ax1.get_xticklabels(),visible=False)
#setp(ax1.get_yticklabels(),visible=False)
#xlim(-7,7)

ax2 = subplot(2,2,1)
data=mlab.load('nImp4/EC-TDA_GFResponse_4.00',usecols=[0,2],unpack=True)
ax2.plot(data[0],data[1],linewidth=2,color='b')
ax2.annotate('U = 4t',xy=(-5,-10),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='16')
xlim(-3.587,3.587)
ylabel(r'$A(\omega)$')
#setp(ax2.get_xticklabels(),visible=False)
setp(ax2.get_yticklabels(),visible=False)
#xlim(-7,7)
ax3 = subplot(2,2,2)
data=mlab.load('nImp4/EC-TDA_GFResponse_6.00',usecols=[0,2],unpack=True)
ax3.plot(data[0],data[1],linewidth=2,color='b')
ax3.annotate('U = 6t',xy=(-5,-10),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='16')
xlim(-3.457,3.457)
#setp(ax2.get_xticklabels(),visible=False)
setp(ax3.get_yticklabels(),visible=False)
#xlim(-7,7)
ax4 = subplot(2,2,3)
data=mlab.load('nImp4/EC-TDA_GFResponse_8.00',usecols=[0,2],unpack=True)
ax4.plot(data[0],data[1],linewidth=2,color='b')
ax4.annotate('U = 8t',xy=(-5,-10),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='16')
xlim(-6,6)
xlabel(r'$\omega$ / t')
ylabel(r'$A(\omega)$')
#setp(ax2.get_xticklabels(),visible=False)
setp(ax4.get_yticklabels(),visible=False)
#xlim(-7,7)
ax5 = subplot(2,2,4)
data=mlab.load('nImp4/EC-TDA_GFResponse_10.00',usecols=[0,2],unpack=True)
ax5.plot(data[0],data[1],linewidth=2,color='b')
ax5.annotate('U = 10t',xy=(-5,-10),xycoords='axes points',xytext=None, textcoords='offset points',horizontalalignment='right', verticalalignment='bottom',fontsize='16')
#setp(ax2.get_xticklabels(),visible=False)
setp(ax5.get_yticklabels(),visible=False)
xlim(-6.853,6.853)
xlabel(r'$\omega$ / t')
subplots_adjust(hspace=0.1,wspace=0.075)

#show()
savefig('2DHub_Spectra.eps')

#clf()
#
#for U in [1.0,2.0,4.0,6.0,8.0]:
#    data_DMFT=mlab.load('CDMFT/ldos_u%s0_ep0.05.dat' % str(U),usecols=[0,1],unpack=True)
#    data_NoReopt2=mlab.load('nImp2/NoReopt/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#    data_Reopt2=mlab.load('nImp2/Reopt/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#    data_NoReopt3_Ref=mlab.load('nImp3/NoReopt/Refresh/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#    data_NoReopt3_Down=mlab.load('nImp3/NoReopt/RampDown/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#    data_Reopt3_Ref=mlab.load('nImp3/Reopt/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#    data_Reopt3_Down=mlab.load('nImp3/Reopt/RampDown/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#    data_NoReopt3_AntiTile=mlab.load('nImp3/NoReopt/Refresh/TileU/EC-TDA_GFResponse_%s0' % str(U),usecols=[0,2],unpack=True)
#
#    clf()
#    ax1 = subplot(111)
#    suptitle('Hubbard(APBC) 624 sites DoS: U = %s' % str(U))
#    ylabel(r'-$\Im [A_{1,1}]$')
#    xlabel(r'$\omega$ / t')
#    ax1.plot(data_DMFT[0],data_DMFT[1]/fac,linewidth=2,label='CDMFT(6,6)',color='b')
#    ax1.plot(data_NoReopt2[0],data_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='g')
#    ax1.plot(data_Reopt2[0],data_Reopt2[1],linewidth=1,label='2-site (Reopt GS)',color='g')
#    ax1.plot(data_NoReopt3_Ref[0],data_NoReopt3_Ref[1],linewidth=1,label='3-site',linestyle='dashed',color='r')
#    ax1.plot(data_Reopt3_Ref[0],data_Reopt3_Ref[1],linewidth=1,label='3-site (Reopt GS)',color='r')
#    ax1.plot(data_NoReopt3_Down[0],data_NoReopt3_Down[1],linewidth=1,label='3-site (Ramp Down)',linestyle='dashed',color='c')
#    ax1.plot(data_Reopt3_Down[0],data_Reopt3_Down[1],linewidth=1,label='3-site (Ramp Down: Reopt GS)',color='c')
#    ax1.plot(data_NoReopt3_AntiTile[0],data_NoReopt3_AntiTile[1],linewidth=1,label='3-site, Anti-tile',color='r')
#    legend(loc=0)
##    ylim(-5,30)
##    savefig('Plot_%s.eps' % str(U))
#    show()
#    #savefig('GF_Hub_OBC_8Site_%s.png' % str(U))
#
##    dataDD_Exact=mlab.load('Exact/Exact_DDResponse_%s0' % str(U),usecols=[0,2],unpack=True)
##    dataDD_Reopt1=mlab.load('Reopt/nImp1/EC-TDA_DDResponse_%s0' % str(U),usecols=[0,2],unpack=True)
##    dataDD_Reopt2=mlab.load('Reopt/nImp2/EC-TDA_DDResponse_%s0' % str(U),usecols=[0,2],unpack=True)
##    dataDD_NoReopt1=mlab.load('NoReopt/nImp1/EC-TDA_DDResponse_%s0' % str(U),usecols=[0,2],unpack=True)
##    dataDD_NoReopt2=mlab.load('NoReopt/nImp2/EC-TDA_DDResponse_%s0' % str(U),usecols=[0,2],unpack=True)
##
##    clf()
##    ax1 = subplot(111)
##    suptitle('Hubbard(OBC) 8 sites DD response: U = %s' % str(U))
##    ylabel(r'-$\Im [A_{1,1}]$')
##    xlabel(r'$\omega$ / t')
##    ax1.plot(dataDD_Exact[0],dataDD_Exact[1],linewidth=2,label='Exact',color='b')
##    ax1.plot(dataDD_NoReopt1[0],dataDD_NoReopt1[1],linewidth=1,label='1-site',linestyle='dashed',color='g')
##    ax1.plot(dataDD_NoReopt2[0],dataDD_NoReopt2[1],linewidth=1,label='2-site',linestyle='dashed',color='r')
##    ax1.plot(dataDD_Reopt1[0],dataDD_Reopt1[1],linewidth=1,label='1-site ReoptGS',color='g')
##    ax1.plot(dataDD_Reopt2[0],dataDD_Reopt2[1],linewidth=1,label='2-site ReoptGS',color='r')
##    legend(loc=0)
##    ylim(-5,30)
##    savefig('Plot_%s.eps' % str(U))
#    #show()
#    #savefig('DD_Hub_OBC_8Site_%s.png' % str(U))
