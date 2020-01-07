# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:00:11 2019

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import math
from scipy.spatial import distance

#import os to handle operating system
import os
#============================================================================= 
#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Images\\live_SDC_images\\Kinesin-13mNG\\data\\linescan_data\\"

df = pd.DataFrame()

#all of the means and SEM for the data are condensed into a single excel file
#with distance from the falgellar tip indicated

df = pd.read_excel(datadir + '190422_kinesin13mNG_taxol_dmso_intensity_profile_exp_fits.xlsx')

#setup a range of lengths between 0-1.2 to plot the function from Jane
x = np.linspace(0,1.3,1e3)

#for each flagellar pair, write a specific function with the fit parameters Jane
#dervied in Mathematica

afd_y = 1023.12 * np.cosh(2.66569 * (-1.2 + x))
aft_y = 735.56 * np.cosh(2.50137 * (-1.2 + x))
cfd_y = 1489.83 * np.cosh(2.81018 * (-1.2 + x))
cft_y = 692.984 * np.cosh(3.37861 * (-1.2 + x))
pfd_y = 1071.2 * np.cosh(2.77688 * (-1.2 + x))
pft_y = 1232.22 * np.cosh(2.10299 * (-1.2 + x))

#==============================================================================

#plot each of the functions with the measured data and errors

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,afd_y)
    plt.errorbar(df['distance'], df['AF_dmso_mean'], df['AF_dmso_sem'], \
                 linestyle='None', marker='o', label='AF_DMSO', capsize=5)
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=24)
    plt.ylabel('Intensity (au)', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,1.3])
    plt.ylim([0,1.2e4])    
    plt.legend(loc='best', fontsize=16)
    plt.tight_layout()
#    plt.savefig('190909_kinesin13mNG_decay_AF_DMSO.svg')

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,aft_y)
    plt.errorbar(df['distance'], df['AF_taxol_mean'], df['AF_taxol_sem'], \
                 linestyle='None', marker='o', label='AF_Taxol', capsize=5)
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=24)
    plt.ylabel('Intensity (au)', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,1.3])
    plt.ylim([0,1.2e4])    
    plt.legend(loc='best', fontsize=16)
    plt.tight_layout()
#    plt.savefig('190909_kinesin13mNG_decay_AF_taxol.svg')

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,cfd_y)
    plt.errorbar(df['distance'], df['CF_dmso_mean'], df['CF_dmso_sem'],\
                 linestyle='None', marker='o', label='CF_DMSO', capsize=5)
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=24)
    plt.ylabel('Intensity (au)', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,1.3])
    plt.ylim([0,1.2e4])    
    plt.legend(loc='best', fontsize=16)
    plt.tight_layout()
#    plt.savefig('190909_kinesin13mNG_decay_CF_DMSO.svg')

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,cft_y)
    plt.errorbar(df['distance'], df['CF_taxol_mean'], df['CF_taxol_sem'], \
                 linestyle='None', marker='o', label='CF_Taxol', capsize=5)
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=24)
    plt.ylabel('Intensity (au)', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,1.3])
    plt.ylim([0,1.2e4])    
    plt.legend(loc='best', fontsize=16)
    plt.tight_layout()
#    plt.savefig('190909_kinesin13mNG_decay_CF_taxol.svg')

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,pfd_y)
    plt.errorbar(df['distance'], df['PF_dmso_mean'], df['PF_dmso_sem'],\
                 linestyle='None', marker='o', label='PF_DMSO', capsize=5)
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=24)
    plt.ylabel('Intensity (au)', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,1.3])
    plt.ylim([0,1.2e4])    
    plt.legend(loc='best', fontsize=16)
    plt.tight_layout()
#    plt.savefig('190909_kinesin13mNG_decay_PF_DMSO.svg')

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,pft_y)
    plt.errorbar(df['distance'], df['PF_taxol_mean'], df['PF_taxol_sem'], \
                 linestyle='None', marker='o', label='PF_Taxol', capsize=5)
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=24)
    plt.ylabel('Intensity (au)', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,1.3])
    plt.ylim([0,1.2e4])    
    plt.legend(loc='best', fontsize=16)
    plt.tight_layout()
#    plt.savefig('190909_kinesin13mNG_decay_PF_taxol.svg')

#=============================================================================
#now import the summary data for the lamba fits (mean, sem) and plot those on 
#a bar chart

df2 = pd.DataFrame()

df2 = pd.read_excel( datadir + "190910_kinesin13mNG_lambda_fits.xlsx")

x_pos = np.arange(len(df2['fp']))
la = df2['lambda']
la_err = df2['lambda_err']*1.96


#first plot everything on the same plot
with sns.axes_style('white'):
    fig, ax = plt.subplots(figsize=(5,5))
    plt.bar(x_pos[:], la[:], align='center', \
            edgecolor='black', linewidth=3)
    plt.errorbar(x_pos[:], la[:], yerr=la_err[:], \
                 linestyle='None', marker='.', label='AF_Taxol', capsize=7, color='k')    
    ax.set_ylabel(u'lambda (${\mu}m$)', fontsize=30)
    ax.set_xticks(x_pos[:])
    ax.set_xticklabels(df2['fp'][:], fontsize=30)
#    plt.ylim([0,0.1])
    plt.tight_layout()
#    plt.savefig('190910_AF_PF_CF_lambda_mean_err.svg')  


with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.plot(x,pft_y, label='PF_taxol', color='g')
    plt.plot(x,pfd_y, label='PF_DMSO', color='g')
    plt.plot(x,aft_y, label='AF_taxol', color='b')
    plt.plot(x,afd_y, label='AF_DMSO', color='b')
    plt.plot(x,cft_y, label='CF_taxol', color='y')
    plt.plot(x,cfd_y, label='CF_DMSO', color='y')
    plt.legend(loc='best', fontsize=16)
    
    

