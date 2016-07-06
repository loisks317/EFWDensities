# EFWKp.py
#
# looks at EFW densities with Kp
# try to establish MLT and Kp radial distance of plasmapause
#
# LKS December 2015
#
#
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
os.chdir('/Users/loisks/Documents/Functions')
import pickling as pickling
os.chdir('/Users/loisks/Desktop/')
#
#
# parameters and starting conditions
date1='20130101'
#date2='20130201'
date2='20150501'
dt0=datetime.datetime.strptime(date1, '%Y%m%d')
dt2=datetime.datetime.strptime(date2, '%Y%m%d')
nLbins=21
nMLTbins=48
nMLATbins=31
LBins=np.linspace(1.5, 6.5, nLbins)#+0.125
MLTbins=np.linspace(0, 2*pi, nMLTbins)
MLATbins=np.linspace(-20, 10, nMLATbins)+0.5
dir=['A', 'B']

name=['rbsp-a', 'rbsp-b']
name2=['rbspa', 'rbspb']
# pickle function
#
#
# begin loop

for iKp in range(8):
  Edens=[[[] for i in range(nLbins)] for j in range(nMLTbins)]
  for idir in range(len(dir)):
    dt1=dt0    
    while dt1 < dt2:
        try:
                date=datetime.datetime.strftime(dt1, '%Y%m%d')
                os.chdir('/Users/loisks/Desktop/liemohn10/loisks/EFW_L3_'+dir[idir])
                efw_file='*'+date+'*.cdf'
                gemf=glob.glob(efw_file)
                pyf=pycdf.CDF(gemf[0])
                epoch=pd.DatetimeIndex(pyf['epoch'][...])
                MLML=pyf['mlt_lshell_mlat'][...]
                dMLT=pd.DataFrame({'MLT':np.swapaxes(MLML, 1,0)[0]}, index=epoch)
                dL=pd.DataFrame({'L':np.swapaxes(MLML, 1,0)[1]}, index=epoch)
                dMLAT=pd.DataFrame({'MLAT':np.swapaxes(MLML,1,0)[2]},index=epoch)
                deDens=pd.DataFrame({'eDens':pyf['density'][...]}, index=epoch)
                rt = pd.period_range(date,periods=1441, freq='T').to_timestamp()
                eDens=np.array(deDens['eDens'].resample('1min', how='median').reindex(index=rt,fill_value=np.nan))
                L=np.array(dL['L'].resample('1min', how='median').reindex(index=rt,fill_value=np.nan))
                MLT=np.array(dMLT['MLT'].resample('1min', how='median').reindex(index=rt,fill_value=np.nan))
                MLAT=np.array(dMLAT['MLAT'].resample('1min', how='median').reindex(index=rt,fill_value=np.nan))
                # get Kp
                kyoto=spk.fetch('kp', dt1, dt1)
                day=datetime.datetime.strftime(dt1,'%d') # need to get day of month
                kp_arr=kyoto['kp'][(int(day)-1)*8: int(day)*8]
                kpC=0
                kpM=0
                nKP=np.zeros(1440)
                for iKP in range(1440):
                    nKP[iKP]=kp_arr[kpM]
                    kpC+=1
                    if kpC == 180:
                        kpC=0
                        kpM+=1

                Kp=np.array(nKP)
                #
                # great, now everything is on same time scale
                # now sort
                for imlt in range(nMLTbins):
                   for iL in range(nLbins):
                      #for imlat in range(nMLATbins):
                        L_indexes=np.where((L >= 1.25+.125+0.25*iL) & (L < 1.25+0.125+0.25*(iL+1)))[0]
                        MLT_indexes=np.where((MLT >= 0.5*imlt) & (MLT < 0.5*(imlt+1)))[0]
                        Kp_indexes=np.where((Kp >= iKp) & (Kp < iKp+1))[0]
                        #MLAT_indexes=np.where((MLAT >= -20-0.5+imlat) & (MLAT < -20-0.5+(imlat+1)))[0]
                        
                        matches=list(set(L_indexes) & set(MLT_indexes) & set(Kp_indexes)) #& set(MLAT_indexes))
                        Edens[imlt][iL]+=list(eDens[matches])
                        #Edens[imlt][iL][imlat]+=list(eDens[matches])
        except:
                print dt1
        dt1=dt1+datetime.timedelta(days=1)
    print 'saving first file'
    os.chdir('/Users/loisks/Desktop/ResearchProjects/EFWDensities')
    subdir_name='Sorted_Kp'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name, 0777) 
    pickling.hdf5_data_save(Edens,'Kp='+ str(iKp)+'_'+str(iKp+1)+'_'+ dir[idir], subdir_name,nLbins,nMLTbins)
    os.chdir('..')




        


