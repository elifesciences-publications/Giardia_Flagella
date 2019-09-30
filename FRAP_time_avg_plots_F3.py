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
from scipy.stats import sem
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
# This script picks up where 'FRAP_single_cell_fits_analysis_plot.py" left off.
#It calculates the time averaged recovery from all cells for each regions 
#and then plots these with the measured values from the script above.
#Finally, it performs linear regression on the inital phase of recovery to 
#compare the slopes of the FP and CA regions.
#=============================================================================

#Use the diffusion_fit equation to plot the fit to the time averaged data
#It is simpler to have one for each the FP and CA regions to plot both datasets

#Equantion for FP region (w and D fixed based on measurements)
def diffusion_fit_fp (t, C, b):
    '''
    This function is the same as 'diffusion_fit' in the previous script, but
    now the measured values are plugged in for D and w.
    
    '''
    diffusion_fp = b + (C*1.5 - b) * (1 - (1.5**2 * (1.5**2 + 4 * pi * 0.055 * t)**-1 ) **0.5)
    return diffusion_fp 

#Equantion for CA region (w and D fixed based on measurements)
def diffusion_fit_ca (t, C, b):
    '''
    This function is the same as 'diffusion_fit' in the previous script, but
    now the measured values are plugged in for D and w
    
    '''
    diffusion_ca = b + (C*1 - b) * (1 - (1**2 * (1**2 + 4 * pi * 0.018 * t)**-1 ) **0.5)
    return diffusion_ca 


#============================================================================= 
#set the directory for the PF_CA data
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\IFT FRAP\\data\\AFFP_PFCA_data\\"
files = np.array(os.listdir(datadir))
#initialize the dataframe for the data
df_ca = pd.DataFrame()
PFCA_df = pd.DataFrame()
        
#============================================================================= 

#first use the PF CA data
pf_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell5', 'cell6',\
                        'cell8', 'cell9', 'cell10', 'cell11', 'cell12',\
                        'cell13', 'cell16', 'cell17', 'cell18', 'cell19',\
                        'cell20', 'cell21', 'cell22', 'cell23', 'cell24', \
                        'cell25', 'cell26', 'cell27', 'cell28', 'cell29', \
                        'cell30', 'cell31', 'cell32', 'cell33', 'cell34', \
                        'cell35'])
    
for j, pcn in enumerate(pfca_cell_number):
    #use glob to read the data from xlsx files into the dataframe    
    pfca_data = glob.glob(datadir + '\\PF_data\\' + '*' + '_IFT81NG_FRAP_' + \
    '*' + pcn + '.xlsx')[0]
    pfca_data_2 = pd.read_excel(pfca_data)

    #correct for photbleaching and conduct background subtraction    
    pfca_corr_data = (pfca_data_2['int'] - pfca_data_2['bckgrd']) / \
    ((pfca_data_2['bleach']-pfca_data_2['bckgrd'])/\
      (pfca_data_2['bleach'][0]-pfca_data_2['bckgrd']))

    #setup time interval for fitting    
    pfca_time = np.linspace(0,299,300)
    
    #append the data for time averaged plotting
    df_ca = df_ca.append([pfca_corr_data], ignore_index=True)
           

#now transpose the data and change the column names to cell number
PFCA_df = df_ca.T
PFCA_df.columns = pfca_cell_number
#add a column of averages (do this before adding the time column!)
PFCA_df['average_int'] = df_ca.mean()

#add a column of time
PFCA_df['time'] = np.linspace(0,299,300)

#now plot the data
plt.plot(PFCA_df['time'], PFCA_df['average_int'])

    
#============================================================================ 

#now do the PF FP data
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\IFT FRAP\\data\\PFFP_AFFP_data\\"
files = np.array(os.listdir(datadir))
p0 = [5000, 1000]#inital guesses for curvefit

f = 300 #final frame to fit

#initialize dataframes for data
df_fp = pd.DataFrame()
PFFP_df = pd.DataFrame()
#============================================================================= 

#now try with the PF CA data
pf_fp_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell6', \
                        'cell7', 'cell8', 'cell9', 'cell10', 'cell11', \
                        'cell12', 'cell13', 'cell14', 'cell15', 'cell16',\
                        'cell17', 'cell18', 'cell19','cell20', 'cell21', \
                        'cell22','cell23', 'cell24', 'cell25', 'cell26'])
    
for j, pfcn in enumerate(pf_fp_cell_number):
    #use glob to read the data from xlsx files into the dataframe        
    pffp_data = glob.glob(datadir + '\\PF_data\\' + '*' + '_IFT81NG_FRAP_' + \
    '*' + pfcn + '.xlsx')[0]
    pffp_data_2 = pd.read_excel(pffp_data)
    
    #correct for photbleaching and conduct background subtraction    
    pffp_corr_data = (pffp_data_2['int'] - pffp_data_2['bckgrd']) / \
    ((pffp_data_2['bleach']-pffp_data_2['bckgrd'])/\
      (pffp_data_2['bleach'][0]-pffp_data_2['bckgrd']))

    #setup time interval for fitting            
    pffp_time = np.linspace(0,599,600)
    
    #append the data for time averaged plotting    
    df_fp = df_fp.append([pffp_corr_data], ignore_index=True)
    
#now transpose the data and change the column names to cell number
PFFP_df = df_fp.T
PFFP_df.columns = pf_fp_cell_number
#add a column of averages (do this before adding the time column!)
PFFP_df['average_int'] = df_fp.mean()

#add a column of time
PFFP_df['time'] = np.linspace(0,599,600)


#=============================================================================
#now make linear fits for initial recovery phase    

#first fit the FP data
fp_linear = 8 #final frame to fit
  
from scipy.stats import linregress
slope_fp, intercept_fp, r_value_fp, p_value_fp, std_err_fp = \
scipy.stats.linregress(PFFP_df['time'][1:fp_linear],\
         (diffusion_fit_fp(PFFP_df['time'][1:fp_linear]-PFFP_df['time'][1], *pf_fp_f)-\
         diffusion_fit_fp(PFFP_df['time'][1:fp_linear]-PFFP_df['time'][1], *pf_fp_f).min())/\
          diffusion_fit_fp(PFFP_df['time'][1:fp_linear]-PFFP_df['time'][1], *pf_fp_f).max())
print("r-squared_fp:", r_value_fp**2, "slope_fp:", slope_fp)  


ca_linear = 20 #final frame to fit

slope_ca, intercept_ca, r_value_ca, p_value_ca, std_err_ca = \
scipy.stats.linregress(PFCA_df['time'][1:ca_linear],\
         (diffusion_fit_ca(PFCA_df['time'][1:ca_linear]-PFCA_df['time'][1], *pf_ca_f)-\
         diffusion_fit_ca(PFCA_df['time'][1:ca_linear]-PFCA_df['time'][1], *pf_ca_f).min())/\
          diffusion_fit_ca(PFCA_df['time'][1:ca_linear]-PFCA_df['time'][1], *pf_ca_f).max())
print("r-squared_ca:", r_value_ca**2, "slope_ca:", slope_ca)

#print the ratio of slopes between the two regions
print(slope_fp / slope_ca)

#=============================================================================
#generate final plot of the data for both regions

with sns.axes_style('white'):
    plt.figure()
    plt.plot(PFCA_df['time'][1:f]-PFCA_df['time'][1],\
         (diffusion_fit_ca(PFCA_df['time'][1:f]-PFCA_df['time'][1], *pf_ca_f)-\
         diffusion_fit_ca(PFCA_df['time'][1:f]-PFCA_df['time'][1], *pf_ca_f).min())/\
          diffusion_fit_ca(PFCA_df['time'][1:f]-PFCA_df['time'][1], *pf_ca_f).max(),\
                          color='k',label='CA')
         
    plt.plot(PFCA_df['time'][1:f],(intercept_ca + slope_ca*PFCA_df['time'][1:f]),\
             'k--', label='ca_linear_regression')
    
    plt.scatter(PFCA_df['time'][1:f]-PFCA_df['time'][1],(PFCA_df['average_int'][1:f] -\
         diffusion_fit_ca(PFCA_df['time'][1:f]-PFCA_df['time'][1], *pf_ca_f).min())/\
          diffusion_fit_ca(PFCA_df['time'][1:f]-PFCA_df['time'][1], *pf_ca_f).max(),color='m',\
                          edgecolor='k', label='time averaged intensity')
    plt.xlabel('Time (seconds)', fontsize=24)
    plt.ylabel('Recovery', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)    
    plt.legend(loc='lower right')
    plt.xlim([0,300])
    plt.ylim([0,1])
#    plt.savefig('181029_PFCA_linearfit_norm_scatter.svg')    

with sns.axes_style('white'):
    plt.figure()    
    plt.plot(PFFP_df['time'][1:f]-PFFP_df['time'][1],\
         (diffusion_fit_fp(PFFP_df['time'][1:f]-PFFP_df['time'][1], *pf_fp_f)-\
         diffusion_fit_fp(PFFP_df['time'][1:f]-PFFP_df['time'][1], *pf_fp_f).min())/\
          diffusion_fit_fp(PFFP_df['time'][1:f]-PFFP_df['time'][1], *pf_fp_f).max(),\
                          color='k', label='FP')
         
    plt.plot(PFFP_df['time'][1:f],(intercept_fp + slope_fp*PFFP_df['time'][1:f]),\
             'k--', label='fp_linear_regression')
    
    plt.scatter(PFFP_df['time'][1:f]-PFFP_df['time'][1],(PFFP_df['average_int'][1:f] -\
         diffusion_fit_fp(PFFP_df['time'][1:f]-PFFP_df['time'][1], *pf_fp_f).min())/\
          diffusion_fit_fp(PFFP_df['time'][1:f]-PFFP_df['time'][1], *pf_fp_f).max(), \
                           color='g', edgecolor='k',label='time averaged intensity')          
    
    plt.xlabel('Time (seconds)', fontsize=24)
    plt.ylabel('Recovery', fontsize=24)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)        
    plt.legend(loc='lower right')
    plt.xlim([0,300])
    plt.ylim([0,1])
#    plt.savefig('181029_PFFP_linearfit_norm_scatter.svg')   





