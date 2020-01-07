# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 11:40:14 2019

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
from scipy import stats
import glob
import statsmodels.stats.api as sms
from sklearn import metrics
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
#import os to handle operating system
import os
#============================================================================= 
#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\Taxol-IFT Expts\\190130_taxol_tubulin_saturation\\"

data = pd.read_excel(datadir + "\\190205_taxol_tubulin_saturation_analysis_na.xlsx")

data = data.dropna()

#first find the mean for each time point for each flagellum

AF_mean = data[data['fp'] == 'AF'].groupby('tx').length.mean()

CF_mean = data[data['fp'] == 'CF'].groupby('tx').length.mean()

PF_mean = data[data['fp'] == 'PF'].groupby('tx').length.mean()

VF_mean = data[data['fp'] == 'VF'].groupby('tx').length.mean()



with sns.axes_style('white'):
    plt.figure(figsize=(6,6))
    sns.lineplot(data[data['fp'] == 'AF']['tx'],\
                     (data[data['fp'] == 'AF']['length']-AF_mean[0]), label = 'AF')
    sns.lineplot(data[data['fp'] == 'CF']['tx'],\
                     (data[data['fp'] == 'CF']['length']-CF_mean[0]), label = 'CF')    
    sns.lineplot(data[data['fp'] == 'PF']['tx'],\
                     (data[data['fp'] == 'PF']['length']-PF_mean[0]), label = 'PF')    
    sns.lineplot(data[data['fp'] == 'VF']['tx'],\
                     (data[data['fp'] == 'VF']['length']-VF_mean[0]), label = 'VF')        
    plt.ylabel(u'Change in flagellar length(${\mu}m$)', fontsize=20)
    plt.xlabel('Taxol treatment time (hrs)', fontsize=20)
    plt.xticks(np.arange(0,6,1.0))    
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.xlim([0,5])
    plt.ylim([0,6])    
    plt.legend(loc='best', fontsize=10)
    plt.tight_layout()
#    plt.savefig('190808_taxol_saturation_length_changes.svg')
   




