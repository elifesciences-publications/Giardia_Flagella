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
#==============================================================================
#The purpose of this script is to import the linescan data for Kinesin-13mNG
#cells and plot the data from the tip to the base and calculate the flagellar tip
#intensity for each flagellum.
#==============================================================================
#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Images\\live_SDC_images\\Kinesin-13mNG\\data\\linescan_data\\taxol_linescans\\projections\\"

#initialize dfs for all the CF data
cftaxol_data = pd.DataFrame()
cftaxol_df = pd.DataFrame()
cftaxol_max = pd.DataFrame()

cftaxol_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell6', 'cell7', \
                        'cell8', 'cell9', 'cell11', 'cell12', 'cell13', 'cell14',\
                        'cell15', 'cell16', 'cell17'])
    
for j, cnt in enumerate(cftaxol_cell_number):
    #use glob to pull from the data directory
    cftaxol_data = glob.glob(datadir + '\\taxol\\CF\\' + '*' + '_kinesin13mNG_linescan_CF_taxol_' + '*' + cnt + '.xlsx')[0]
    cftaxol_data_2 = pd.read_excel(cftaxol_data)
    #write the cell number with the corresponding linescan    
    cftaxol_data_2['cell_number'] = cnt
    #append the data together                         
    cftaxol_df = cftaxol_df.append(cftaxol_data_2)
    #calculate length of the flagellum
    cftaxol_length = max(cftaxol_data_2['488_length'])
    #subtract the bacakground intensity
    cftaxol_bcksub = cftaxol_data_2['488_int']-cftaxol_data_2['488_bck'].mean()
    cftaxol_int_max = max(cftaxol_bcksub)
    #calculate the AUC 
    cftaxol_auc = metrics.auc(cftaxol_data_2['488_length'], cftaxol_bcksub)
    #append all of the data together
    cftaxol_max = cftaxol_max.append([[cnt, cftaxol_length, cftaxol_int_max, cftaxol_auc]], ignore_index=True)

cftaxol_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc'])
cftaxol_max['flagellar_pair'] = 'CF_taxol'    
    
#compile all the CF data
cfdmso_data = pd.DataFrame()
cfdmso_df = pd.DataFrame()
cfdmso_max = pd.DataFrame()
cf_dmso_exp_dist = pd.DataFrame()

cfdmso_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell6', 'cell7', \
                        'cell8', 'cell9', 'cell10', 'cell11', 'cell12', 'cell13', 'cell14',\
                        'cell15', 'cell16'])
    
for j, cnd in enumerate(cfdmso_cell_number):
    #use glob to pull from the data directory
    cfdmso_data = glob.glob(datadir + '\\dmso\\CF\\' + '*' + '_kinesin13mNG_linescan_CF_dmso_' + '*' + cnd + '.xlsx')[0]
    cfdmso_data_2 = pd.read_excel(cfdmso_data)
    #write the cell number with the corresponding linescan        
    cfdmso_data_2['cell_number'] = cnd
    #append the data together                                
    cfdmso_df = cfdmso_df.append(cfdmso_data_2)
    #calculate length of the flagellum
    cfdmso_length = max(cfdmso_data_2['488_length'])
    #subtract the bacakground intensity    
    cfdmso_bcksub = cfdmso_data_2['488_int']-cfdmso_data_2['488_bck'].mean()  
    cfdmso_int_max = max(cfdmso_bcksub)
    #calculate the AUC      
    cfdmso_auc = metrics.auc(cfdmso_data_2['488_length'], cfdmso_bcksub)
    #append all of the data together    
    cfdmso_max = cfdmso_max.append([[cnd, cfdmso_length, cfdmso_int_max, cfdmso_auc]], ignore_index=True)

cfdmso_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc'])
cfdmso_max['flagellar_pair'] = 'CFdmso' 
          
#plot the spatial distribution of kinesin13mNG w/ DMSO vs Taxol tx
with sns.axes_style('white'):
    plt.figure(figsize=(9,5))
    sns.lineplot(x=cftaxol_df['488_length'], y=cftaxol_df['488_int']-cftaxol_df['488_bck'], label='CF_taxol', ci=68)    
    sns.lineplot(x=cfdmso_df['488_length'], y=cfdmso_df['488_int']-cfdmso_df['488_bck'], label='CF_dmso', ci=68)    
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=30)
    plt.ylabel('Intensity (au)', fontsize=30)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.xlim([0,20])
    plt.ylim([0,1.5e4])    
    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190107_kinesin13mNG_CF_distribution_dmso_taxol.svg')           


#=============================================================================

#compile all the AF data
aftaxol_data = pd.DataFrame()
aftaxol_df = pd.DataFrame()
aftaxol_max = pd.DataFrame()
af_taxol_exp_dist = pd.DataFrame()


aftaxol_cell_number = np.array(['cell4','cell5', 'cell6', 'cell7', 'cell8',\
                         'cell9', 'cell10', 'cell12', 'cell13', 'cell14',\
                        'cell15', 'cell16', 'cell17'])
    
for j, ant in enumerate(aftaxol_cell_number):
    #use glob to pull from the data directory    
    aftaxol_data = glob.glob(datadir + '\\taxol\\AF\\' + '*' + '_kinesin13mNG_linescan_AF_taxol_' + '*' + ant + '.xlsx')[0]
    aftaxol_data_2 = pd.read_excel(aftaxol_data)
    #write the cell number with the corresponding linescan            
    aftaxol_data_2['cell_number'] = ant
    #append the data together                                       
    aftaxol_df = aftaxol_df.append(aftaxol_data_2)
    #calculate length of the flagellum
    aftaxol_length = max(aftaxol_data_2['488_length'])
    #subtract the bacakground intensity        
    aftaxol_bcksub = aftaxol_data_2['488_int']-aftaxol_data_2['488_bck'].mean()    
    aftaxol_int_max = max(aftaxol_bcksub)
    #calculate the AUC 
    aftaxol_auc = metrics.auc(aftaxol_data_2['488_length'], aftaxol_bcksub)
    #append all of the data    
    aftaxol_max = aftaxol_max.append([[ant, aftaxol_length, aftaxol_int_max, aftaxol_auc]], ignore_index=True)

aftaxol_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc'])
aftaxol_max['flagellar_pair'] = 'AF_taxol'    
    
#compile all the AF data
afdmso_data = pd.DataFrame()
afdmso_df = pd.DataFrame()
afdmso_max = pd.DataFrame()
af_dmso_exp_dist = pd.DataFrame()

afdmso_cell_number = np.array(['cell1','cell2','cell3', 'cell4', 'cell5', 'cell6',\
                        'cell8', 'cell9','cell12','cell13','cell15','cell16'])
    
for j, dan in enumerate(afdmso_cell_number):
    #use glob to pull from the data directory    
    afdmso_data = glob.glob(datadir + '\\dmso\\AF\\' + '*' + '_kinesin13mNG_linescan_AF_dmso_' + '*' + dan + '.xlsx')[0]
    afdmso_data_2 = pd.read_excel(afdmso_data)
    #write the cell number with the corresponding linescan            
    afdmso_data_2['cell_number'] = dan
    #append the data together                                             
    afdmso_df = afdmso_df.append(afdmso_data_2)
    #calculate length of the flagellum
    afdmso_length = max(afdmso_data_2['488_length'])
    #subtract the bacakground intensity        
    afdmso_bcksub = afdmso_data_2['488_int']-afdmso_data_2['488_bck'].mean()    
    afdmso_int_max = max(afdmso_bcksub)
    #calculate the AUC  
    afdmso_auc = metrics.auc(afdmso_data_2['488_length'], afdmso_bcksub)
    #append all of the data    
    afdmso_max = afdmso_max.append([[dan, afdmso_length, afdmso_int_max, afdmso_auc]], ignore_index=True)

afdmso_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc'])
afdmso_max['flagellar_pair'] = 'AFdmso' 

with sns.axes_style('white'):
    plt.figure(figsize=(9,5))
    sns.lineplot(x=aftaxol_df['488_length'], y=aftaxol_df['488_int']-aftaxol_df['488_bck'], label='AF_taxol', ci=68)    
    sns.lineplot(x=afdmso_df['488_length'], y=afdmso_df['488_int']-afdmso_df['488_bck'], label='AF_dmso', ci=68)    
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=30)
    plt.ylabel('Intensity (au)', fontsize=30)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.xlim([0,20])
    plt.ylim([0,1.5e4])    
    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190107_kinesin13mNG_AF_distribution_dmso_taxol.svg')          

#=============================================================================
#compile all the AF data
pftaxol_data = pd.DataFrame()
pftaxol_df = pd.DataFrame()
pftaxol_max = pd.DataFrame()
pf_taxol_exp_dist = pd.DataFrame()

pftaxol_cell_number = np.array(['cell1','cell2', 'cell3', 'cell5', 'cell7',\
                         'cell10', 'cell11', 'cell12', 'cell13', \
                        'cell15', 'cell16'])
    
for j, pnt in enumerate(pftaxol_cell_number):
    #use glob to pull from the data directory    
    pftaxol_data = glob.glob(datadir + '\\taxol\\PF\\' + '*' + '_kinesin13mNG_linescan_PF_taxol_' + '*' + pnt + '.xlsx')[0]
    pftaxol_data_2 = pd.read_excel(pftaxol_data)
    #write the cell number with the corresponding linescan                
    pftaxol_data_2['cell_number'] = pnt
    #append the data together                                                    
    pftaxol_df = pftaxol_df.append(pftaxol_data_2)
    #calculate length
    pftaxol_length = max(pftaxol_data_2['488_length'])
    #subtract background intensity
    pftaxol_bcksub = pftaxol_data_2['488_int']-pftaxol_data_2['488_bck'].mean()   
    pftaxol_int_max = max(pftaxol_bcksub)
    #calculate AUC  
    pftaxol_auc = metrics.auc(pftaxol_data_2['488_length'], pftaxol_bcksub)
    #append all of the data
    pftaxol_max = pftaxol_max.append([[pnt, pftaxol_length, pftaxol_int_max, pftaxol_auc]], ignore_index=True)

pftaxol_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc'])
pftaxol_max['flagellar_pair'] = 'PF_taxol'    
    
#compile all the AF data
pfdmso_data = pd.DataFrame()
pfdmso_df = pd.DataFrame()
pfdmso_max = pd.DataFrame()
pf_dmso_exp_dist = pd.DataFrame()

pfdmso_cell_number = np.array(['cell2','cell5', 'cell6', 'cell7', 'cell10',\
                        'cell11','cell12','cell13', 'cell14', 'cell15'])
    
for j, pnd in enumerate(pfdmso_cell_number):
    #use glob to pull from the data directory    
    pfdmso_data = glob.glob(datadir + '\\dmso\\PF\\' + '*' + '_kinesin13mNG_linescan_PF_dmso_' + '*' + pnd + '.xlsx')[0]
    pfdmso_data_2 = pd.read_excel(pfdmso_data)
    #write the cell number with the corresponding linescan                    
    pfdmso_data_2['cell_number'] = pnd
    #append the data together                                                           
    pfdmso_df = pfdmso_df.append(pfdmso_data_2)
    #calculate length
    pfdmso_length = max(pfdmso_data_2['488_length'])
    #subtract background intensity   
    pfdmso_bcksub = pfdmso_data_2['488_int']-pfdmso_data_2['488_bck'].mean()    
    pfdmso_int_max = max(pfdmso_bcksub)
    #calculate AUC        
    pfdmso_auc = metrics.auc(pfdmso_data_2['488_length'], pfdmso_bcksub)
    #append all of the data
    pfdmso_max = pfdmso_max.append([[pnd, pfdmso_length, pfdmso_int_max, pfdmso_auc]], ignore_index=True)

pfdmso_max.columns = np.array(['cell_number', 'flagellar_length', 'max_int', 'int_auc'])
pfdmso_max['flagellar_pair'] = 'PFdmso' 

with sns.axes_style('white'):
    plt.figure(figsize=(9,5))
    sns.lineplot(x=pftaxol_df['488_length'], y=pftaxol_df['488_int']-pftaxol_df['488_bck'], label='PF_taxol', ci=68)    
    sns.lineplot(x=pfdmso_df['488_length'], y=pfdmso_df['488_int']-pfdmso_df['488_bck'], label='PF_dmso', ci=68)    
    plt.xlabel(u'Distance from flagellar tip (${\mu}m$)', fontsize=30)
    plt.ylabel('Intensity (au)', fontsize=30)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.xlim([0,20])
    plt.ylim([0,1.5e4])    
    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190107_kinesin13mNG_PF_distribution_dmso_taxol.svg') 

#=============================================================================

#combine all of the max data into a single df for plotting
      
dfs_all = [cfdmso_max,cftaxol_max,afdmso_max,aftaxol_max,pfdmso_max,pftaxol_max]
    
max_all_data = pd.DataFrame()

max_all_data = pd.concat(dfs_all)      

#plot length of flagella for DMSO and taxol cells
with sns.axes_style('white'):
    plt.figure(figsize=(6,7))
    sns.barplot(x=max_all_data['flagellar_pair'], y=max_all_data['flagellar_length'],\
                capsize=.2, errcolor=".2", edgecolor=".2", linewidth=2.5, ci=95) 
    plt.xlabel('Treatment', fontsize=30)
    plt.ylabel(u'Flagellar length (${\mu}m$)', fontsize=30)
#    plt.yscale('log')
    plt.xticks(rotation=90)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190107_kinesin13mNG_length_dmso_taxol.svg') 

#plot the max intensity for taxol and dmso
with sns.axes_style('white'):
    plt.figure(figsize=(5,6))
    sns.barplot(x=max_all_data['flagellar_pair'], y=max_all_data['max_int'],\
                capsize=.3, errcolor=".2", edgecolor=".2", linewidth=2, ci=95) 
    plt.xlabel('Treatment', fontsize=30)
    plt.ylabel('Maximum intensity', fontsize=30)
#    plt.yscale('log')
    plt.xticks(rotation=90)
    plt.rc('xtick', labelsize=24)
    plt.rc('ytick', labelsize=24)
    plt.legend(loc='best', fontsize=20)
    plt.tight_layout()
#    plt.savefig('190910_kinesin13mNG_max_int_dmso_taxol.svg') 

    


    


    


    
    
  






