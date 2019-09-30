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
#This script has multiple parts and is likely written ineffeciently. First, it
#loads in the individual xlsx files containing all data on the FRAP experiments
#on individual cells (ROI intensity, background intensity, etc), computes the
#bleach and background corrected values for ROIs and fit the diffusion equation
#to the data in order to obtain D, Co, and b. The means and sem of these
#values for each cell are then exported to a xlsx file and reimported to plot
#summary stats (inefficient). 
#=============================================================================

#First, write a function to fit the data and extract the diffusion constant


def diffusion_fit (t, C, b, w, D):
    '''
    This function is going to use the Ellenberg method for determining the
    diffusion constant from our experimental data.
    t = time (minutes)
    C = inital flouresecene intensity of ROI (au)
    b = bleach depth (intensity following bleach of ROI, au)
    w = width of bleached region (um)
    D = effective diffusion coefficient (um^2/min)
    
    '''
    diffusion = b + (C*w - b) * (1 - (w**2 * (w**2 + 4 * pi * D * t)**-1 ) **0.5)
    return diffusion 

#============================================================================= 

#set the directory for glob to serach through
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\IFT FRAP\\data\\AFFP_PFCA_data\\"
files = np.array(os.listdir(datadir))
p0 = [5000, 1000, 1, 1] #inital guesses for curvefit
f = 125 #final frame to perform fit

#initialize dataframe for the loop below
PF_df = pd.DataFrame()
#============================================================================= 

#first do this with the PF CA data
pf_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell5', 'cell6',\
                        'cell8', 'cell9', 'cell10', 'cell11', 'cell12',\
                        'cell13', 'cell16', 'cell17', 'cell18', 'cell19',\
                        'cell20', 'cell21', 'cell22', 'cell23', 'cell24', \
                        'cell25', 'cell26', 'cell27', 'cell28', 'cell29', \
                        'cell30', 'cell31', 'cell32', 'cell33', 'cell34', \
                        'cell35'])
    
for j, pcn in enumerate(pf_cell_number):
    #use glob to read the data from xlsx files into the dataframe
    pf_data = glob.glob(datadir + '\\PF_data\\' + '*' + '_IFT81NG_FRAP_' + '*' +\
                        pcn + '.xlsx')[0]
    pf_data_2 = pd.read_excel(pf_data)
    
    #correct for photbleaching and conduct background subtraction
    pf_data_2['int_corr'] = (pf_data_2['int'] - pf_data_2['bckgrd']) / \
    ((pf_data_2['bleach']-pf_data_2['bckgrd'])/(pf_data_2['bleach'][0]-\
      pf_data_2['bckgrd']))
    
    #setup time interval for fitting
    pf_data_2['time'] = np.linspace(0,5,300)
    
    #fit the data with curvefit; ignore inital intensity frame during fit to only
    #fit the recovery phase    
    pf_f,pf_cov = curve_fit(diffusion_fit,\
                            pf_data_2['time'][1:f]-pf_data_2['time'][1],\
                                     pf_data_2['int_corr'][1:f], p0)
    
    #append the individual cell data to a dataframe for statistical analyses;
    #column order is cell number, C, b, w, D
    PF_df = PF_df.append([[pcn, pf_f[0], pf_f[1], pf_f[2], pf_f[3]]], \
                         ignore_index=True)
    

#============================================================================= 
#Now repeat the above for the FP data
#set the directory for the PF FP data
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\IFT FRAP\\data\\PFFP_AFFP_data\\"
files = np.array(os.listdir(datadir))
p0 = [5000, 1000, 1, 1]#inital guesses for curvefit
f = 125 #final frame to fit

PF_fp_df = pd.DataFrame()
#df_fp = pd.DataFrame()
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
    pffp_data_2['int_corr'] = (pffp_data_2['int'] - pffp_data_2['bckgrd']) / \
    ((pffp_data_2['bleach']-pffp_data_2['bckgrd'])/(pffp_data_2['bleach'][0]-\
      pffp_data_2['bckgrd']))
    
    pffp_data_2 = pffp_data_2.dropna()

    #setup time interval for fitting    
    pffp_data_2['time'] = np.linspace(0,5,300)
        
    #fit the data with curvefit; ignore inital intensity frame during fit to only
    #fit the recovery phase    
    pffp_f,pffp_cov = curve_fit(diffusion_fit, \
                                pffp_data_2['time'][1:f]-pffp_data_2['time'][1],\
                                           pffp_data_2['int_corr'][1:f], p0)
    
    #append the individual cell data to a dataframe for statistical analyses;
    #column order is cell number, C, b, w, D      
    PF_fp_df = PF_fp_df.append([[pfcn, pffp_f[0], pffp_f[1], pffp_f[2], \
                                 pffp_f[3]]], ignore_index=True)
    
#============================================================================= 
#output the fit data to a csv file for calculating summary stats
#PF_fp_df.to_csv(datadir + 'PF_FP_singlecellfits.csv', index=False)
#
#PF_df.to_csv(datadir + 'PF_CA_singlecellfits.csv', index=False)

#============================================================================= 
#Now import a single xlsx file with the summary stats for the effective 
#diffusion constant

#initalize a dataframe for the data   
param_data = pd.DataFrame()
#import the data to the above data frame
param_data = pd.read_excel('E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\IFT FRAP\\data\\181009_FRAP_fit_parameters.xlsx')

#setup aspects of the data to plot as a bar chart
x_pos = np.arange(len(param_data['region']))
diff_app = param_data['D_mean']
diff_app_error = param_data['D_95CI']


#first plot everything on the same plot
with sns.axes_style('white'):
    fig, ax = plt.subplots(figsize=(5,6))
    plt.bar(x_pos[:], diff_app[:], yerr=diff_app_error[:], align='center', \
            edgecolor='black', capsize=10, linewidth=3)
    plt.errorbar(x_pos[:], diff_app[:], yerr=diff_app_error[:], capsize=10, \
                 elinewidth=3, markeredgewidth=3, fmt='.', color='k')
    ax.set_ylabel(u'Diffusion coefficient (${\mu}m^2$/sec)', fontsize=24)
    ax.set_xticks(x_pos[:])
    ax.set_xticklabels(param_data['region'][:], fontsize=20)
#    plt.ylim([0,0.1])
    plt.tight_layout()
#    plt.savefig('181010_PF-FP_PF-CA_FRAP_Dapp_plot.svg') 

#=============================================================================



    





    


    
    
