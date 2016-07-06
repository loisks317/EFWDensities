# EFWKp_load.py
#
# load int he processed EFW data sorted by
# kp, MLT, L, and MLAT
#
# LKS, December 2015
#
# imports
import numpy as np
from spacepy import pycdf
import glob
from numpy import ma
import os
import pandas as pd
import datetime
import spacepy.pybats.kyoto as spk
import itertools as itert
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
import matplotlib.dates as dates
import pickle 
from dateutil.relativedelta import relativedelta 
import h5py
import spacepy.datamodel as dm
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
import matplotlib.pyplot as plt
os.chdir('/Users/loisks/Desktop/Functions/')
import pickling as pickling
os.chdir('/Users/loisks/Documents/ResearchProjects/EFWDensities/')

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
nLbins=21
nMLATbins=31
nMLTbins=48
LBins=np.linspace(1.5, 6.5, nLbins)#+0.125
MLTbins=np.linspace(0, 2*np.pi, nMLTbins)
MLATbins=np.linspace(-20, 10, nMLATbins)+0.5
for iKp in range(8):
    kpData=[[[[] for k in range(nMLATbins)] for i in range(nLbins)] for j in range(nMLTbins)]
    mKpData=[[[] for i in range(nLbins)] for j in range(nMLTbins)]
    os.chdir('/Users/loisks/Documents/ResearchProjects/EFWDensities/Sorted_Kp/')
    #EFiles=glob.glob('Kp='+str(iKp)+'_'+str(iKp+1)+'*')
    data1=pickling.hdf5_data_open('Kp='+str(iKp)+'_'+str(iKp+1)+'_A',nLbins, nMLTbins)
    data2=pickling.hdf5_data_open('Kp='+str(iKp)+'_'+str(iKp+1)+'_B',nLbins, nMLTbins)
    for imlt in range(nMLTbins):
        for iL in range(nLbins):
            #for imlat in range(nMLATbins):
                mKpData[imlt][iL]=np.nanmedian(np.array(list(data1[imlt][iL])+list(data2[imlt][iL])))

    # plot the data
    fig=plt.figure()
    ax=fig.add_subplot(111, polar=True)
    datah_m=ma.masked_invalid(np.array(mKpData).transpose())
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
    xL2=['00','', '06','', '12','', '18']
    cb.set_label('N$_e$ cm$^{-3}$', fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=25)    
    cb.ax.tick_params(labelsize=35) 
    ax.set_yticklabels(xT2, fontsize=30)
    ax.set_xticklabels(xL2, fontsize=30)
    ax.grid(True)
    plt.draw()
    os.chdir('..')
    subdir_name='Kp_Figures'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name, 0777) 
    os.chdir(subdir_name)#
    fig.set_size_inches(13,9)
    plt.savefig('EFWDensities_Kp_'+str(iKp)+'_'+str(iKp+1)+'.png')
    plt.close(fig)
    os.chdir('..')
