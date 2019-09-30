# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 14:30:52 2017

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
#import os to handle operating system
import os
#============================================================================= 
#The purpose of this script is to import all of the IFT injection events and
#to measure the time between each injection. 
#==============================================================================

#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\all\\"
files = np.array(os.listdir(datadir))
files = files[np.array(['txt' in f for f in files])]
#============================================================================= 
#define the patterns that glob will screen for
flagellar_pair = np.array(['AF', 'PF'])
#cell_number = np.linspace(1, 60, 60)
#cell_number = np.array(['cell1','cell2','cell3','cell4','cell5'])
cell_number = np.array(['1','3','4','5','6', '8', '9', '10', '11', '13',\
                        '15', '16', '17', '18', '19', '20', '22', '23','27',\
                        '29', '30', '32', '33', '34', '37', '39', '40','42',\
                        '43', '45', '46', '47', '48', '49', '50', '51', '52',\
                        '54', '55', '56', '57', '27b'])
direction = np.array(['forward', 'backward'])
#============================================================================= 
#initalize data frame to append all data 
df = pd.DataFrame()

#use glob to find the files based on the patterns set above 
for j, cn in enumerate(cell_number):
    for i, fp in enumerate(flagellar_pair):
        for k, direct in enumerate(direction):
            try: #we need to use try/except to iterate through particle_number w/o KeyError
                position = glob.glob(datadir + '*' + \
                                         'Particles_positions_' + '*' + direct + '*' + \
                                         fp + '*' + cn + '.txt')[0]
                position_data = pd.read_table(position)
                for line in position_data:
                    df = df.append([[cn, line, position_data[line].max(), fp, direct]], ignore_index=True)
            except KeyError:
                pass
                                    
#now name each column
df.columns = np.array(['cell_number', 'line', 'particle_maxT', 'flagellar_pair', 'direction'])

#==============================================================================

#now just get all of the time data into a single dataframe
df2 = pd.DataFrame()
df2 = df[1::2] #we are selecting all of the even rows, which contain the time data

AF = df2[df2['flagellar_pair']=='AF']
PF = df2[df2['flagellar_pair']=='PF']
AF_fwd = AF[AF['direction']=='forward']
AF_bck = AF[AF['direction']=='backward']
PF_fwd = PF[PF['direction']=='forward']
PF_bck = PF[PF['direction']=='backward']

#now sort on particlemaxT so that we avoid making unrelated comparisons that
#generate negative values
AF_fwd = AF_fwd.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])
AF_bck = AF_bck.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])
PF_fwd = PF_fwd.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])
PF_bck = PF_bck.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])

#calculate deltaT and add to the dataframe
AF_fwd['deltaT'] = AF_fwd['particle_maxT'].diff()
AF_bck['deltaT'] = AF_bck['particle_maxT'].diff()
PF_fwd['deltaT'] = PF_fwd['particle_maxT'].diff()
PF_bck['deltaT'] = PF_bck['particle_maxT'].diff()

#now recombine all of the data and drop the blank/NA values
all_data = [AF_fwd.dropna(), AF_bck.dropna(), PF_fwd.dropna(), PF_bck.dropna()]
all_data_joined = pd.concat(all_data)

#now instead of taking absolute values lets just remove all of the
# 'time_ms1' values b/c those are messing things up and don't make any sense anyway
all_data_joined = all_data_joined[all_data_joined.line != 'time_ms1']

#there are a few (~3-4) weird values (~-25,000) so just drop anything less than zero
all_data_joined = all_data_joined[all_data_joined.deltaT >= 0]

#plot a frequency histogram
#plt.hist plots all of the data so split on direction and flagellar pair again

AF_plot = all_data_joined[all_data_joined['flagellar_pair']=='AF']
PF_plot = all_data_joined[all_data_joined['flagellar_pair']=='PF']
AF_fwd_plot = AF_plot[AF_plot['direction']=='forward']
AF_bck_plot = AF_plot[AF_plot['direction']=='backward']
PF_fwd_plot = PF_plot[PF_plot['direction']=='forward']
PF_bck_plot = PF_plot[PF_plot['direction']=='backward']

fwd = all_data_joined[all_data_joined['direction']=='forward']
bck = all_data_joined[all_data_joined['direction']=='backward']

#output the data to check that everything is sorted and calculating the correct
#deltaT
#outputdir = 'D:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\'
#AF_fwd_plot.to_csv(outputdir + '\\170203_AF_anterograde_deltaT.csv', index=False)
#AF_bck_plot.to_csv(outputdir + '\\170203_AF_retrograde_deltaT.csv', index=False)
#PF_fwd_plot.to_csv(outputdir + '\\170203_PF_anterograde_deltaT.csv', index=False)
#PF_bck_plot.to_csv(outputdir + '\\170203_PF_retrograde_deltaT.csv', index=False)
#=============================================================================
























