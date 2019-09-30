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

#from __future__ import print_function
#from scipy.integrate import simps
#from numpy import trapz

#import os to handle operating system
import os
#==============================================================================
#The purpose of this script is to import the linescan data for Kinesin-13mNG
#cells and plot the data from the tip to the base and calculate the flagellar tip
#intensity for each flagellum.
#==============================================================================
#setup the data directory
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Images\\live_SDC_images\\Kinesin-13mNG\\data\\linescan_data\\181011_intensity_linescan\\"

#initialize dfs for all the PF data
pf_data = pd.DataFrame()
pf_df = pd.DataFrame()
pf_max = pd.DataFrame()

pf_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell5', 'cell6', \
                        'cell7','cell8', 'cell9', 'cell10', 'cell11', 'cell12',\
                        'cell13', 'cell14', 'cell15', 'cell16', 'cell18',\
                        'cell20','cell21','cell22', 'cell23','cell24','cell25',\
                        'cell26', 'cell27','cell28','cell29','cell30','cell31',\
                        'cell32','cell33','cell34','cell35'])
    
for j, pcn in enumerate(pf_cell_number):
    #use glob to pull from the data directory
    pf_data = glob.glob(datadir + '\\PF\\' + '*' + '_kinesin13mNG_linescan_PF_' + '*' + pcn + '.xlsx')[0]
    pf_data_2 = pd.read_excel(pf_data)
    #write the cell number with the corresponding linescan
    pf_data_2['cell_number'] = pcn
    #append the data together       
    pf_df = pf_df.append(pf_data_2)
    #calculate length of the flagellum
    pf_length = max(pf_data_2['488_length'])
    #substract the average background intensity
    pf_bcksub = pf_data_2['488_int']-1300  
    pf_int_max = max(pf_bcksub)

    pf_sum = pf_bcksub.sum()   
    #calculate the AUC
    pf_auc = metrics.auc(pf_data_2['488_length'], pf_bcksub)
    #append all of the data together
    pf_max = pf_max.append([[pcn, pf_length, pf_int_max, pf_auc, pf_sum]], ignore_index=True)
#label columns
pf_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc', 'int_sum'])
#add a column to ID flagellum specific data
pf_max['flagellar_pair'] = 'PF'

#initialize dfs for all the AF data
af_data = pd.DataFrame()
af_df = pd.DataFrame()
af_max = pd.DataFrame()
af_cell_number = np.array(['cell3', 'cell5', 'cell6', 'cell7', 'cell9', 'cell11',\
                        'cell13', 'cell16', 'cell17','cell19', 'cell21', 'cell23',\
                        'cell24','cell25','cell26','cell27','cell29',\
                        'cell30','cell31','cell32','cell33','cell34','cell35'])
    
for j, acn in enumerate(af_cell_number):
    #use glob to pull from the data directory    
    af_data = glob.glob(datadir + '\\AF\\' + '*' + '_kinesin13mNG_linescan_AF_' + '*' + acn + '.xlsx')[0]    
    af_data_2 = pd.read_excel(af_data)
    #write the cell number with the corresponding linescan   
    af_data_2['cell_number'] = acn
    #append the data together              
    af_df = af_df.append(af_data_2)
    #calculate length of the flagellum
    af_length = max(af_data_2['488_length'])
    #substract the average background intensity    
    af_bcksub = af_data_2['488_int']-1300
    af_int_max = max(af_bcksub)
    
    af_sum = af_bcksub.sum()
    #calculate the AUC        
    af_auc = metrics.auc(af_data_2['488_length'], af_bcksub)
    #append all of the data together
    af_max = af_max.append([[acn, af_length, af_int_max, af_auc, af_sum]], ignore_index=True)
#label columns
af_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc', 'int_sum'])
#add a column to ID flagellum specific data
af_max['flagellar_pair'] = 'AF'    
    
#initialize dfs for all the CF data
cf_data = pd.DataFrame()
cf_df = pd.DataFrame()
cf_check = pd.DataFrame()
cf_max = pd.DataFrame()
cf_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell5', 'cell7', \
                        'cell8', 'cell9', 'cell10', 'cell11', 'cell12','cell14',\
                        'cell15', 'cell16', 'cell17', 'cell18','cell19','cell20',\
                        'cell21','cell22', 'cell23','cell24','cell25','cell26', \
                        'cell27','cell28','cell29','cell30','cell31','cell32',\
                        'cell33','cell34','cell35'])
    
for j, ccn in enumerate(cf_cell_number):
    #use glob to pull from the data directory    
    cf_data = glob.glob(datadir + '\\CF\\' + '*' + '_kinesin13mNG_linescan_CF_' + '*' + ccn + '.xlsx')[0]
    cf_data_2 = pd.read_excel(cf_data)
    #write the cell number with the corresponding linescan       
    cf_data_2['cell_number'] = ccn
    #drop NA values   
    cf_df = cf_df.append(cf_data_2.dropna())
    #calculate length of the flagellum
    cf_length = max(cf_data_2['488_length'])
    #substract the average background intensity        
    cf_bcksub = cf_data_2['488_int']-1300 
    cf_int_max = max(cf_bcksub)
    
    cf_sum = cf_bcksub.sum()
    #calculate the AUC        
    cf_auc = metrics.auc(cf_data_2['488_length'], cf_bcksub)
    #append all of the data together       
    cf_max = cf_max.append([[ccn, cf_length, cf_int_max, cf_auc, cf_sum]], ignore_index=True)
#label columns
cf_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc', 'int_sum'])
#add a column to ID flagellum specific data
cf_max['flagellar_pair'] = 'CF' 

#combine all of the max data into a single df
      
dfs_all = [cf_max,pf_max,af_max]
    
max_all_data = pd.DataFrame()

max_all_data = pd.concat(dfs_all)      

#output the data 
#cf_check.to_csv(datadir + '\\190305_cf_check.csv', index=False)
#pf_df.to_csv(datadir + '\\181108_pf_check.csv', index=False)

#=============================================================================
#plot the intensity distribution from tip to base
with sns.axes_style('white'):
    plt.figure(figsize=(9,5))
    sns.lineplot(x=pf_df['488_length'], y=pf_df['488_int']-1300, label='PF', ci=68)    
    sns.lineplot(x=af_df['488_length'], y=af_df['488_int']-1300, label='AF', ci=68)    
    sns.lineplot(x=cf_df['488_length'], y=cf_df['488_int']-1300, label='CF', ci=68)     
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=30)
    plt.ylabel('Intensity (au)', fontsize=30)
#    plt.yscale('log')
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.xlim([0,16])
    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190415_kinesin13mNG_length_distribution_9x5.svg')
    
#plot the flagellar tip intensity for each flagellum
with sns.axes_style('white'):
    plt.figure(figsize=(5,6))
    sns.barplot(x=max_all_data['flagellar_pair'], y=max_all_data['max_int'],\
                capsize=.2, errcolor=".2", edgecolor=".2", \
                linewidth=2.5, ci=95)   
    plt.xlabel('Flagellar pair', fontsize=30)
    plt.ylabel('Flagellar tip intensity (au)', fontsize=30)
#    plt.yscale('log')
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
#    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190305_kinesin13mNG_max_int.svg')




    
    
  






