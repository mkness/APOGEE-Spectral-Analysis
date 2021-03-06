#!/usr/bin/python
import numpy
from numpy import savetxt
import matplotlib
from matplotlib import pyplot
import scipy
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
s.set_family('serif')
s.set_size(14)
from matplotlib import rc
rc('text', usetex=False)
rc('font', family='serif')
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
s.set_family('serif')
rcParams["xtick.labelsize"] = 14
rcParams["ytick.labelsize"] = 14
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
s = matplotlib.font_manager.FontProperties()
majorLocator   = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(5)
yminorLocator   = MultipleLocator(10)
yminorLocator2   = MultipleLocator(25)
xminorLocator   = MultipleLocator(5)
yminorLocator   = MultipleLocator(5)
ymajorLocator   = MultipleLocator(50)
xmajorLocator   = MultipleLocator(10)
rcParams['figure.figsize'] = 15.0, 10.0
fig, temp = pyplot.subplots(5,1, sharex=True, sharey=False)
ax1 = temp[0]
ax2 = temp[1]
ax3 = temp[2]
ax4 = temp[3]
ax5 = temp[4]
#plot(dataall[:, 0, 0], 1. * coeffs[:, 0]) # mean spectra  
#plot(dataall[:, 0, 0], 1. * coeffs[:, 1]) # teff 
#plot(dataall[:, 0, 0], 1. * coeffs[:, 2]) # log g 
#plot(dataall[:, 0, 0], 1. * coeffs[:, 3]) # feh 

ax1.plot(dataall[:, 0, 0], 1. * coeffs[:, 0],color = 'k' ,linewidth = 2) # median
ax2.plot(dataall[:, 0, 0], 1. * coeffs[:, 1], color = 'green' ,linewidth = 2) # teff 
ax3.plot(dataall[:, 0, 0], 1. * coeffs[:, 2], color = 'blue',linewidth = 2 ) #  g 
ax4.plot(dataall[:, 0, 0], 1. * coeffs[:, 3],color = 'red' ,linewidth = 2) # feh  
ax5.plot(dataall[:, 0, 0], 1. * coeffs[:, 3], color = 'red',linewidth = 2) # feh  
ax5.plot(dataall[:, 0, 0], 1. * coeffs[:, 2], color = 'blue',linewidth = 2) #  g 
ax5.plot(dataall[:, 0, 0], 1000. * coeffs[:, 1], color = 'green',linewidth = 2) # teff 
l1a  = 15390
l2a = 15394
l1b = 15693
l2b = 15700
l1c = 15956
l2c = 15960
l1d = 16202
l2d = 16210
l1e = 16116
l2e = 16122
l1f = 1666.5
l2f = 16172.5

ax1.vlines(15697, -1,2, linestyle = 'dashed', linewidth = 2) 
ax2.vlines(15697, -1,2, linestyle = 'dashed', linewidth = 2) 
ax3.vlines(15697, -1,2, linestyle = 'dashed', linewidth = 2) 
ax4.vlines(15697, -1,2, linestyle = 'dashed', linewidth = 2) 
ax5.vlines(15697, -1,2, linestyle = 'dashed', linewidth = 2) 

ax2.set_xlim(l1b - 20, l2b + 20 ) 
ax2.set_ylim(-0.001 ,0.001) 
ax1.set_ylim(0.6 ,1.2) 
ax3.set_ylim(-0.5 ,0.5) 
ax4.set_ylim(-0.5 ,0.5) 
ax5.set_ylim(-0.5 ,0.5) 
ax1.set_xlim(l1b - 20, l2b + 20 ) 
ax2.set_xlim(l1b - 20, l2b + 20 ) 
ax3.set_xlim(l1b - 20, l2b + 20 ) 
ax4.set_xlim(l1b - 20, l2b + 20 ) 
ax5.set_xlim(l1b - 20, l2b + 20 ) 
ax1.text(l1b-19, 1.1, "mean spectra" , fontsize = 12) 
ax2.text(l1b-19, 0.00010, "Teff coeff" , fontsize = 12) 
ax3.text(l1b-19, 0.3, "log g coeff" , fontsize = 12) 
ax4.text(l1b-19, 0.3, "[Fe/H]  coeff" , fontsize = 12) 
ax5.text(l1b-19, 0.3, "[Fe/H]  coeff, log g coeff, Teff coeff*1000" , fontsize = 12) 
ax1.set_title("REGION 2 USED FOR [Fe/H] INDEX") 

ax1.axvspan(l1b, l2b, facecolor='c', alpha=0.1)
ax2.axvspan(l1b, l2b, facecolor='c', alpha=0.1)
ax3.axvspan(l1b, l2b, facecolor='c', alpha=0.1)
ax4.axvspan(l1b, l2b, facecolor='c', alpha=0.1)
ax5.axvspan(l1b, l2b, facecolor='c', alpha=0.1)
ax5.set_xlabel("Wavelength $\AA$", fontsize = 20) 
ax1.set_ylabel("coeff a0", fontsize = 20) 
ax2.set_ylabel("coeff a1", fontsize = 20) 
ax3.set_ylabel("coeff a2", fontsize = 20) 
ax4.set_ylabel("coeff a3", fontsize = 20) 
ax5.set_ylabel("coeff a1,a2,a3", fontsize = 20) 

fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
