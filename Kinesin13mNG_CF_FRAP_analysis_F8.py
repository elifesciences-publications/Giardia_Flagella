# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:41:39 2017

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from math import pi
import scipy
import scipy.stats
import statsmodels.stats.api as sms
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
#try to do the fit using scipy.optimize_curvefit
from scipy.optimize import curve_fit
#import os to handle operating system
import os
#=============================================================================
#Notes: this script is used to load in the individual csv files containing
#FRAP recovery data derived from FIJI ROI analysis of bleached regions through
#time with accompanying background ROI data. Once loaded the background intensity
#is subtracted from each frame of the timeseries. Finally, the data is normalized
#to the pre-bleach intensity. The output is a dataframe with each representing
#a single cell and each column representing a single time point during the 
#acqusition. Transposing the dataframe allows for the columns to become labelled
#with cell_number and to calculate the average recovery for all cells. Finally 
# we can calcualte the 95% confidence interval and plot the recovery data.
#=============================================================================
#initialize a new dataframe
CF_df = pd.DataFrame()

CF_df2 = pd.DataFrame()

#set patterns for glob to search for
cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell5', 'cell6','cell7', \
                        'cell8', 'cell9', 'cell10', 'cell11', 'cell12','cell13','cell14',\
                        'cell15', 'cell16', 'cell17', 'cell18','cell19'])

#setup an array of time values to interate the for loop through
time = np.linspace(0,179,180)
#set the directory for the data
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Images\\live_SDC_images\\Kinesin-13mNG\\data\\FRAP_data\\"
files = np.array(os.listdir(datadir))
#use glob to find the correct files then read those files to subtract background
#fluoro, then for each time point calculate the percent fluoro relative to the
#inital intensity; finally append the data into a final dataframe; this format
#gives us rows of fluoro data for each time and lets us take an average to plot
#instead of individual traces
for j, cn in enumerate(cell_number):
#    AF_int_data = glob.glob(datadir + 'AF_FRAP\\' + '*' + '_Kinesin13mNG_FRAP_' + '*' + cn + '.xlsx')[0]
#    PF_int_data = glob.glob(datadir + 'PF_FRAP\\' + '*' + '_Kinesin13mNG_FRAP_' + '*' + cn + '.xlsx')[0]
    CF_int_data = glob.glob(datadir + 'CF_FRAP\\' + '*' + '_Kinesin13mNG_FRAP_' + '*' + cn + '.xlsx')[0]
    
#    AF_int_data_2 = pd.read_excel(AF_int_data)
#    PF_int_data_2 = pd.read_excel(PF_int_data)
    CF_int_data_2 = pd.read_excel(CF_int_data)
    
#    AF_delta_data = AF_int_data_2['intensity'].sub(AF_int_data_2['bckgrd'], \
#                                 axis=0).div(AF_int_data_2['bleach']/ \
#                                       AF_int_data_2['bleach'][0])
#    PF_delta_data = PF_int_data_2['intensity'].sub(PF_int_data_2['bckgrd'], \
#                                 axis=0).div(PF_int_data_2['bleach']/ \
#                                       PF_int_data_2['bleach'][0])
    CF_delta_data = CF_int_data_2['intensity'].sub(CF_int_data_2['bckgrd'], \
                                 axis=0).div(CF_int_data_2['bleach']/ \
                                       CF_int_data_2['bleach'][0])    
    for i in time:
#        AF_per_diff2 = 1 - (AF_delta_data[0] - AF_delta_data) / AF_delta_data[0]
#        PF_per_diff2 = 1 - (PF_delta_data[0] - PF_delta_data) / PF_delta_data[0]
        CF_per_diff2 = 1 - (CF_delta_data[0] - CF_delta_data) / CF_delta_data[0]

#    AF_df = AF_df.append([AF_per_diff2], ignore_index=True)
#    PF_df = PF_df.append([PF_per_diff2], ignore_index=True)
    CF_df = CF_df.append([CF_per_diff2], ignore_index=True)
    
#AF_df2 = AF_df.T
#PF_df2 = PF_df.T
CF_df2 = CF_df.T

#output the data 
CF_df2.to_excel(datadir + '\\Kinesin13mNG_FRAP_CF_data.xlsx', index=False)

#==============================================================================
#add summary statistics to the dataframe

CF_df2['int_avg'] = CF_df2.mean(numeric_only=True, axis=1)
CF_df2['int_SEM'] = CF_df2.sem(numeric_only=True, axis=1)

#now make new dataframe for plotting data, only go to t=180
CF_plot = pd.DataFrame()
CF_plot = CF_df2[:180]

CF_plot['time'] = np.linspace(-5,174,180)

#==============================================================================

#write a function to fit the data and extract the diffusion constant


def diffusion_fit (t, I, D):
    '''
    This function is going to use the Ellenberg method for determining the
    diffusion constant from our experimental data.
    
    '''
    diffusion = I * (1 - (0.5**2 * (0.5**2 + 4 * pi * D * t)**-1 ) **0.5)
    return diffusion 

#=============================================================================  
CF_fit = CF_plot[5:]

CF_fit['avg_int_norm'] = CF_plot['int_avg']-CF_plot['int_avg'].min()

#=============================================================================    


with sns.axes_style('white'):
    plt.figure(figsize=(7,6))
    plt.plot(CF_fit['time'], diffusion_fit(CF_fit['time'],cf_I, cf_D), \
             'k', label='Fit')
    plt.scatter(CF_fit['time'], CF_fit['avg_int_norm'], s=100, color='y', edgecolor='k', label='Avg Intensity')
    plt.fill_between(CF_fit['time'],(CF_fit['avg_int_norm'] - 1.96*CF_fit['int_SEM']), \
                     (CF_fit['avg_int_norm'] + 1.96*CF_fit['int_SEM']), color='gray', alpha=0.15)    
    plt.xlabel('Time (seconds)', fontsize=24)
    plt.ylabel('Recovery', fontsize=24)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.legend(loc='lower right', fontsize=16)
    plt.xlim([0,180])
    plt.ylim([0,0.5])
    plt.tight_layout()
#    plt.savefig('190306_kinesin13mNG_CF_only_FRAP.svg')
    


