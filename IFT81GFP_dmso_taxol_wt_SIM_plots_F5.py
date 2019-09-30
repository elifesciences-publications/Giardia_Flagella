# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 10:17:11 2017

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
#The purpose of this script is to generate final plots of the SIM data for cells
#untreated, or treated with DMSO or Taxol.
#============================================================================= 

#now reimport the data with annotations (taxol, dmso, or untreated cells)
data = pd.DataFrame()

data = pd.read_csv('E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Nikon_SIM\\SIM_data\\170514_IFT81GFP_taxol_SIM_data.csv')

cont_data = pd.DataFrame()
cont_data = pd.read_csv("E:\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Nikon_SIM\\SIM_data\\170405_SIM_predictions_v2.csv")

#============================================================================= 
#calculate the average length changes between control and taxol tx'd cells
cont = data[data['tx']=='control']
tx = data[data['tx']=='taxol']

cont2 = cont_data[cont_data['regime']== 1]


AF_taxol = tx[tx['flagellar_pair']=='AF']
PF_taxol = tx[tx['flagellar_pair']=='PF']
CF_taxol = tx[tx['flagellar_pair']=='CF']

AF_cont = cont[cont['flagellar_pair']=='AF']
PF_cont = cont[cont['flagellar_pair']=='PF']
CF_cont = cont[cont['flagellar_pair']=='CF']


#============================================================================= 
#plot length changes for each treatment group
with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    sns.barplot(x=data['flagellar_pair'], y=data['flagellar_length'],
               hue=data['tx'], capsize=.2, linewidth=2, palette='Set1',\
                       edgecolor=".2")
    plt.ylabel(r'Flagellar Length ($\mu m$)', fontsize=20)
    plt.xlabel('Flagellar Pair', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.legend(loc='best', title='Treatment')
    plt.ylim([0, 18])
#    plt.title("Total IFT intensity increases with flagellar length", fontsize=16)
    plt.tight_layout()
#    plt.savefig('170707_taxol_length_changes.svg') 

#=============================================================================

#use linear regression to determine r2

#CONTROL
#strategy one:
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = \
scipy.stats.linregress(cont2['flagellar_length'], cont2['auc_int'])
print("r-squared:", r_value**2)

#strategy two:
from scipy import stats
def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2
sns.jointplot(cont2['flagellar_length'], cont2['auc_int'], kind="reg", \
              stat_func=r2)
sns.regplot(cont2['flagellar_length'], cont2['auc_int'])
#different plot; replaces scatter plots with density estimates
sns.jointplot(cont2['flagellar_length'], cont2['auc_int'], kind="kde", space=0, 
              color="b", stat_func=r2)

#DMSO CONTROL
#strategy one:
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = \
scipy.stats.linregress(cont['flagellar_length'], cont['auc_int'])
print("r-squared:", r_value**2)

#strategy two:
from scipy import stats
def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2
sns.jointplot(cont['flagellar_length'], cont['auc_int'], kind="reg", \
              stat_func=r2)
sns.regplot(cont['flagellar_length'], cont['auc_int'])
#different plot; replaces scatter plots with density estimates
sns.jointplot(cont['flagellar_length'], cont['auc_int'], kind="kde", space=0, 
              color="b", stat_func=r2) 

#TAXOL
#strategy one:
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = \
scipy.stats.linregress(tx['flagellar_length'], tx['auc_int'])
print("r-squared:", r_value**2)

#strategy two:
from scipy import stats
def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2
sns.jointplot(tx['flagellar_length'], tx['auc_int'], kind="reg", stat_func=r2)
sns.regplot(tx['flagellar_length'], tx['auc_int'])
#different plot; replaces scatter plots with density estimates
sns.jointplot(tx['flagellar_length'], tx['auc_int'], kind="kde", space=0, 
              color="b", stat_func=r2) 


#=============================================================================
#plot the untreated control cells
with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.xlim([0,20])
    plt.ylim([0,30000])
    sns.regplot(cont_data['flagellar_length'], cont_data['auc_int'], \
                scatter_kws={"s": 50}, label='Control')
    plt.xlabel(r'Flagellar Length ($\mu m$)', fontsize=20)
    plt.ylabel('Integrated Intensity (au)', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.legend(loc='best')
#    plt.title("Total IFT intensity increases with flagellar length", fontsize=16)
    plt.tight_layout() 
#    plt.savefig('180804_total_integrated_intensity_equilibirum_lengths.svg')

#Plot the taxol and dmso cells
with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    plt.xlim([0,20])
    plt.ylim([0,30000])
    sns.regplot(cont['flagellar_length'], cont['auc_int'], \
                scatter_kws={"s": 50}, label='DMSO, R2=0.93')
    sns.regplot(tx['flagellar_length'], tx['auc_int'], scatter_kws={"s": 80}, \
                marker='o', label='Taxol, R2=0.92')
    plt.xlabel(r'Flagellar Length ($\mu m$)', fontsize=20)
    plt.ylabel('Integrated Intensity (au)', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.legend(loc='best')
#    plt.title("Total IFT intensity increases with flagellar length", fontsize=16)
    plt.tight_layout() 
#    plt.savefig('180804_total_integrated_intensity_Taxol_lengths.svg')







