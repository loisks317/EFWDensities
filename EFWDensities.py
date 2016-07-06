# EFWDensities.py
#
# look at stats behind EFW densities
#
# LKS, December 2015. Less than 1 month to South Africa!
#
#
# imports
import numpy as np
from spacepy import pycdf
import glob
import os
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
import matplotlib.dates as dates
from scipy.interpolate import interp1d
import pickle 
import h5py
import spacepy.datamodel as dm
import pandas as pd
from matplotlib import pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from numpy import ma
from mpl_toolkits.axisartist.grid_finder import FixedLocator, \
     MaxNLocator, DictFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
from numpy import ma
#
#
# parameters and starting conditions
date1='20130101'
date2='20150401'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
nLbins=21
nMLTbins=48
LBins=np.linspace(1.5, 6.5, nLbins)#+0.125
MLTbins=np.linspace(0, 24, nMLTbins)
dir=['A', 'B']

#MLAT_bins=[-20, -15, -10, -5, 0, 5, 10, 15, 20]
name=['rbsp-a', 'rbsp-b']
name2=['rbspa', 'rbspb']
# pickle function
#
#
# begin loop
Edens=[[[] for i in range(nLbins)] for j in range(nMLTbins)]
mEdens=[[[] for i in range(nLbins)] for j in range(nMLTbins)]
for idir in range(len(dir)):
    dt1=dt0    
    while dt1 < dt2:
        try:
                date=datetime.datetime.strftime(dt1, '%Y%m%d')
                os.chdir('/Users/loisks/Documents/ResearchProjects/ChargingProject/EFW_L3_'+dir[idir])
                efw_file='*'+date+'*.cdf'
                gemf=glob.glob(efw_file)
                pyf=pycdf.CDF(gemf[0])
                eDens=pyf['density'][...]
                MLT=np.swapaxes(np.array(pyf['mlt_lshell_mlat'][...]),1,0)[0]
                L=np.swapaxes(np.array(pyf['mlt_lshell_mlat'][...]), 1,0)[1]
                #
                # now sort
                for imlt in range(nMLTbins):
                    for iL in range(nLbins):
                        L_indexes=np.where((L >= 1.25+.125+0.25*iL) & (L < 1.25+0.125+0.25*(iL+1)))[0]
                        MLT_indexes=np.where((MLT >= 0.5*imlt) & (MLT < 0.5*(imlt+1)))[0]
                        matches=list(set(L_indexes) & set(MLT_indexes))
                        Edens[imlt][iL]+=list(eDens[matches])
        except(IndexError):
                print dt1
        dt1=dt1+datetime.timedelta(days=1)

# now plot
lenEdens=[[[] for i in range(nLbins)] for j in range(nMLTbins)]
stdEdens=[[[] for i in range(nLbins)] for j in range(nMLTbins)]
os.chdir('/Users/loisks/Documents/ResearchProjects/EFWDensities')
# take the median
for imlt in range(nMLTbins):
    for iL in range(nLbins):
        try:
          temp=np.array(Edens[imlt][iL])
          temp[temp < 0]=np.nan
          mEdens[imlt][iL]=np.nanmedian(temp)
          lenEdens[imlt][iL]=len(temp[~np.isnan(temp)])
          stdEdens[imlt][iL]=np.std(temp[~np.isnan(temp)])/mEdens[imlt][iL]
        except:
            mEdens[imlt][iL]=np.nan
            lenEdens[imlt][iL]=np.nan
            stdEdens[imlt][iL]=np.nan




def dual_half_circle(center, radius, angle=90, ax=None, colors=('w','k'),
                     **kwargs):
    from matplotlib.patches import Wedge
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]
#
#
# intitialize the figure
MLTbins=MLTbins*(15.5*np.pi/180.)
fig=plt.figure()
ax=fig.add_subplot(111, polar=True)
datah_m=ma.masked_invalid(np.array(mEdens).transpose())
X,Y=np.meshgrid(MLTbins, LBins)
ax.set_ylim(0, 6.5)
dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
vmin=1
vmax=4
col=ax.pcolormesh( X, Y, np.log10(datah_m), cmap='jet', vmin=vmin, vmax=vmax )
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
xT2=['', '2', '', '4', '', '6']
#xT2=['', '1', '', '2', '', '3']
xL2=['00','', '06','', '12','', '18']
cb.set_label('cm$^{-3}$', fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=25)    
cb.ax.tick_params(labelsize=35) 
ax.set_yticklabels(xT2, fontsize=30)
ax.set_xticklabels(xL2, fontsize=30)
ax.grid(True)
plt.draw()
subdir_name='EFW_L_MLT_Densities'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('EFWDensities_L_MLT_AllActivity.png')
plt.close(fig)
os.chdir('..')

#
# at L < 4
fig=plt.figure()
ax=fig.add_subplot(111, polar=True)
datah_m=ma.masked_invalid(np.array(mEdens).transpose())
X,Y=np.meshgrid(MLTbins, LBins)
ax.set_ylim(0, 4)
dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
vmin=500
vmax=3000
col=ax.pcolormesh( X, Y, datah_m, cmap='jet', vmin=vmin, vmax=vmax )
#cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
cb = plt.colorbar(col, cax = cbaxes)
plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
#xT2=['', '2', '', '4', '', '6']
xT2=['', '1', '', '2', '', '3']
xL2=['00','', '06','', '12','', '18']
cb.set_label('cm$^{-3}$', fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=25)    
cb.ax.tick_params(labelsize=35) 
ax.set_yticklabels(xT2, fontsize=30)
ax.set_xticklabels(xL2, fontsize=30)
ax.grid(True)
plt.draw()
subdir_name='EFW_L_MLT_Densities'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('EFWDensities_Lless4_MLT_AllActivity.png')
plt.close(fig)
os.chdir('..')
#
#
# number of points
fig=plt.figure()
ax=fig.add_subplot(111, polar=True)
datah_m=ma.masked_invalid(np.array(lenEdens).transpose())
X,Y=np.meshgrid(MLTbins, LBins)
ax.set_ylim(0, 4)
dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
vmin=2
vmax=4
col=ax.pcolormesh( X, Y, np.log10(datah_m), cmap='BrBG', vmin=vmin, vmax=vmax )
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
#cb = plt.colorbar(col, cax = cbaxes)
plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
#xT2=['', '2', '', '4', '', '6']
xT2=['', '1', '', '2', '', '3']
xL2=['00','', '06','', '12','', '18']
cb.set_label('Number of Points', fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=25)    
cb.ax.tick_params(labelsize=35) 
ax.set_yticklabels(xT2, fontsize=30)
ax.set_xticklabels(xL2, fontsize=30)
ax.grid(True)
plt.draw()
subdir_name='EFW_L_MLT_Densities'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('NoPoints_Lless4_MLT_AllActivity.png')
plt.close(fig)
os.chdir('..')
#


#
#
# standard Deviation
fig=plt.figure()
ax=fig.add_subplot(111, polar=True)
datah_m=ma.masked_invalid(np.array(stdEdens).transpose())
X,Y=np.meshgrid(MLTbins, LBins)
ax.set_ylim(0, 4)
dual_half_circle((0,0), 1.0, angle=90, ax=None, colors=('w','k'))
cbaxes = fig.add_axes([0.75, 0.15, 0.03, 0.75])
vmin=.1
vmax=1
col=ax.pcolormesh( X, Y, datah_m, cmap='BrBG_r', vmin=vmin, vmax=vmax )
#cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1))
cb = plt.colorbar(col, cax = cbaxes)
plt.subplots_adjust(right=0.73, top=0.8, bottom=0.15)
#xT2=['', '2', '', '4', '', '6']
xT2=['', '1', '', '2', '', '3']
xL2=['00','', '06','', '12','', '18']
cb.set_label('$\sigma_{rel}$', fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=25)    
cb.ax.tick_params(labelsize=35) 
ax.set_yticklabels(xT2, fontsize=30)
ax.set_xticklabels(xL2, fontsize=30)
ax.grid(True)
plt.draw()
subdir_name='EFW_L_MLT_Densities'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name, 0777) 
os.chdir(subdir_name)#
fig.set_size_inches(13,9)
plt.savefig('stddev_Lless4_MLT_AllActivity.png')
plt.close(fig)
os.chdir('..')
