# -*- coding: utf-8 -*-
"""
Created on Thu May 17 14:29:50 2018

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
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
#This script is used to import the full axoneme flourescence profile of AF and
#PF flagella and to plot this mean flourescence profile with 95% CI shaded. It
#also measures the integrated intensity for 2um of the membrane-bound and 
#cytoplasmic axoneme regions and compares the ratios of flourescence intensity
#between these two regions. These values are then used in Eq3 to compare the 
#relative differences in flourescence intensity in these two regions.
#=============================================================================
#first establish the data direction to import the data
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Flagellar_pore_intensity_measurements\\data\\"

#initalize data frame to append all data 
data = pd.DataFrame()

#import the data from the data directory
data = pd.read_excel(datadir + '180521_AF_PF_full_axo_summary_stats.xlsx')

#Plot the means and 95%CI for each flagellum; zero denotes the BB
with sns.axes_style('white'):
    plt.figure(figsize=(7,5))
    plt.plot(data['AF_length'], data['AF_mean']/data['AF_mean'].max(), \
             label='Anterior', color='b')
    plt.fill_between(data['AF_length'], \
                     data['AF_mean']/data['AF_mean'].max()-data['AF_ci'], \
                     data['AF_mean']/data['AF_mean'].max()+data['AF_ci'], \
                     alpha=0.5, edgecolor='#3F7F4C', facecolor='#90CAF9',\
                     linewidth=0)
    plt.plot(data['PF_length'], data['PF_mean']/data['PF_mean'].max(), \
             label='Posteriolateral', color='r')
    plt.fill_between(data['PF_length'], \
                     data['PF_mean']/data['PF_mean'].max()-data['PF_ci'], \
                     data['PF_mean']/data['PF_mean'].max()+data['PF_ci'], \
                     alpha=0.5, edgecolor='#3F7F4C', facecolor='#E57373', \
                     linewidth=0)
    plt.ylabel('Relative Intensity', fontsize=32)
    plt.xlabel(u'Length (${\mu}m$)', fontsize=32)
    plt.title('Full axoneme flourescence profile of AF and PF')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(loc='best')
    plt.xlim([0, 18])

#=============================================================================
# Compare the integrated intensity between CA and membrane bound regions
# of the axoneme

# first define the CA region for integartion, use 2-4um

PF_cyto = metrics.auc(data['PF_length'][15:24], data['PF_mean'][15:24])
AF_cyto = metrics.auc(data['AF_length'][15:31], data['AF_mean'][15:31])

PF_cyto = metrics.auc(data['PF_length'][15:24], data['PF_mean'][15:24])


# now do the same for the membrane region, use 12-14um
PF_mem = metrics.auc(data['PF_length'][90:99], data['PF_mean'][90:99])
AF_mem = metrics.auc(data['AF_length'][90:106], data['AF_mean'][90:106])

#print the ratio of cyto to mem

print(PF_cyto/PF_mem)
print(AF_cyto/AF_mem)


  
    
    
    
    
    