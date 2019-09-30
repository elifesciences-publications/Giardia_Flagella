# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:19:22 2016

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
#The purpose of this script is to analyze line scan data from SIM images of
#IFT81GFP.
#==============================================================================
#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Nikon_SIM\\SIM_data\\161012\\"
files = np.array(os.listdir(datadir))
files = files[np.array(['csv' in f for f in files])]
#============================================================================= 
#define the patterns that glob will screen for
flagellar_pair = np.array(['AF', 'PF', 'CF'])
cell_number = np.array(['1','2','3','4','5','6','7','8','9','10'])
#============================================================================= 
#initalize data frame to append all data 
df = pd.DataFrame()
#use glob to find the files based on the patterns set above
for j, cn in enumerate(cell_number):
    for i, fp in enumerate(flagellar_pair):
        intensity = glob.glob(datadir + '*' + 'IFT81GFP_SIM_' + '*' + fp + '*' + cn + '.csv')[0]
        int_data = pd.read_csv(intensity)
        #measure the auc for each flagellum
        int_auc = metrics.auc(int_data['X'], int_data['Y'])
        #determine the maximum length of each flagellum
        flagellar_length = max(int_data['X'])
        #calculate the total auc/length
        int_per_leng = int_auc/max(int_data['X'])
        #append all data to a single dataframe for plotting
        df = df.append([[cn, fp, int_auc, int_per_leng, flagellar_length]],\
                       ignore_index=True)
#now name the columns
df.columns = np.array(['cell_number', 'flagellar_pair', 'auc_int',\
                       'int_per_length', 'flagellar_length'])

#now plot using a scatter plot
plt.scatter(df['flagellar_length'], df['auc_int'])
plt.scatter(df['flagellar_length'], df['int_per_length'])

#make a scatterplot using seaborn
sns.lmplot('flagellar_length', 'auc_int', data=df)

#use linear regression to determine r2

#strategy one:
from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = \
scipy.stats.linregress(df['flagellar_length'], df['auc_int'])
print("r-squared:", r_value**2)

#strategy two:
from scipy import stats
def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2
sns.jointplot(df['flagellar_length'], df['auc_int'], kind="reg", stat_func=r2)
sns.regplot(df['flagellar_length'], df['auc_int'])
#different plot; replaces scatter plots with density estimates
sns.jointplot(df['flagellar_length'], df['auc_int'], kind="kde", space=0, 
              color="b", stat_func=r2)


#==============================================================================

#prepare and save final figures for this data:
with sns.axes_style('white'):  
    plt.figure()
    sns.jointplot(df['flagellar_length'], df['auc_int'], kind="reg", stat_func=r2)
    plt.xlabel('Flagellar Length (um)', fontsize=20)
    plt.ylabel('Total Intensity (AUC)', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.ylim([0, 30000])
    plt.title("Total IFT intensity increases with flagellar length", fontsize=16)
    plt.tight_layout()
#    plt.savefig('161019_IFT81GFP_SIM_AUC_intensity_v1.png')

   
