# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 16:29:54 2016

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
datadir = "D:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\all\\"
files = np.array(os.listdir(datadir))
files = files[np.array(['txt' in f for f in files])]
#============================================================================= 
#define the patterns that glob will screen for
flagellar_pair = np.array(['AF', 'CF'])
#cell_number = np.linspace(1, 60, 60)
#cell_number = np.array(['cell1','cell2','cell3','cell4','cell5'])
cell_number = np.array(['10', '16', '18', '24', '26', '28', '30', '31', '33',\
                        '34', '36', '38', '40', '41', '44', '47', '48', '53',\
                        '56', '58', '28b', '56b'])
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
CF = df2[df2['flagellar_pair']=='CF']
AF_fwd = AF[AF['direction']=='forward']
AF_bck = AF[AF['direction']=='backward']
CF_fwd = CF[CF['direction']=='forward']
CF_bck = CF[CF['direction']=='backward']

#now sort on particlemaxT so that we avoid making unrelated comparisons that
#generate negative values
AF_fwd = AF_fwd.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])
AF_bck = AF_bck.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])
CF_fwd = CF_fwd.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])
CF_bck = CF_bck.sort_values(['cell_number', 'flagellar_pair', 'particle_maxT'])

#calculate deltaT and add to the dataframe
AF_fwd['deltaT'] = AF_fwd['particle_maxT'].diff()
AF_bck['deltaT'] = AF_bck['particle_maxT'].diff()
CF_fwd['deltaT'] = CF_fwd['particle_maxT'].diff()
CF_bck['deltaT'] = CF_bck['particle_maxT'].diff()

#now recombine all of the data and drop the blank/NA values
all_data = [AF_fwd.dropna(), AF_bck.dropna(), CF_fwd.dropna(), CF_bck.dropna()]
all_data_joined = pd.concat(all_data)

#now instead of taking absolute values lets just remove all of the
# 'time_ms1' values b/c those are messing things up and don't make any sense anyway
all_data_joined = all_data_joined[all_data_joined.line != 'time_ms1']

#there are a few (~3-4) weird values (~-25,000) so just drop anything less than zero
all_data_joined = all_data_joined[all_data_joined.deltaT >= 0]

#plot a frequency histogram
#plt.hist plots all of the data so split on direction and flagellar pair again

AF_plot = all_data_joined[all_data_joined['flagellar_pair']=='AF']
CF_plot = all_data_joined[all_data_joined['flagellar_pair']=='CF']
AF_fwd_plot = AF_plot[AF_plot['direction']=='forward']
AF_bck_plot = AF_plot[AF_plot['direction']=='backward']
CF_fwd_plot = CF_plot[CF_plot['direction']=='forward']
CF_bck_plot = CF_plot[CF_plot['direction']=='backward']

fwd = all_data_joined[all_data_joined['direction']=='forward']
bck = all_data_joined[all_data_joined['direction']=='backward']

#output the data to check that everything is sorted and calculating the correct
#deltaT
outputdir = 'D:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\'
AF_fwd_plot.to_csv(outputdir + '\\161108_AF_anterograde_deltaT.csv', index=False)
AF_bck_plot.to_csv(outputdir + '\\161108_AF_retrograde_deltaT.csv', index=False)
CF_fwd_plot.to_csv(outputdir + '\\161108_CF_anterograde_deltaT.csv', index=False)
CF_bck_plot.to_csv(outputdir + '\\161108_CF_retrograde_deltaT.csv', index=False)
#=============================================================================
 






















