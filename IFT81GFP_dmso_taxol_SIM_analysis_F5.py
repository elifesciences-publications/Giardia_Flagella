# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 16:08:56 2017

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
#IFT81GFP treated with either taxol or dmso for 1 hour. 
#============================================================================= 
#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Nikon_SIM\\SIM_data\\170413\\"
files = np.array(os.listdir(datadir))
files = files[np.array(['csv' in f for f in files])]
#============================================================================= 
#define the patterns that glob will screen for
flagellar_pair = np.array(['AF', 'PF', 'CF'])
cell_number = np.array(['1','2','3','4','5','6','7','8','9','10', '11', '12',\
                        '13', '14', '15', '16', '17', '18' , '19', '20', '21',\
                        '22', '23'])
#============================================================================= 
#initalize data frame to append all data 
df = pd.DataFrame()
#use glob to find the files based on the patterns set above
for j, cn in enumerate(cell_number):
    for i, fp in enumerate(flagellar_pair):
        intensity = glob.glob(datadir + '*' + 'IFT81GFP_taxol_SIM_' + '*' + fp + \
        '*' + cn + '.csv')[0]
        int_data = pd.read_csv(intensity)
        int_auc = metrics.auc(int_data['X'], int_data['Y'])
        flagellar_length = max(int_data['X'])
        int_per_leng = int_auc/max(int_data['X'])
#        print(int_data)
        df = df.append([[cn, fp, int_auc, int_per_leng, flagellar_length]],\
                       ignore_index=True)
#now name the columns
df.columns = np.array(['cell_number', 'flagellar_pair', 'auc_int',\
                       'int_per_length', 'flagellar_length'])
    
#now output the data for analysis and some quick annotations
#outputdir = 'D:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Nikon_SIM\\SIM_data\\170413'
#df.to_csv(outputdir + '\\170418_IFT81GFP_taxol_SIM_data.csv', index=False)

#============================================================================= 
