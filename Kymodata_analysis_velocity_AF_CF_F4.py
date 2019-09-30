# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 14:49:23 2016

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
#The purpose of this set of scripts is to pull the IFT dynamics data from the
#kymograph analysis files generated by kmyograph clear/direct software.
#============================================================================= 
#import files to analyze
datadir = "E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\all\\"
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
            avg_intensity = glob.glob(datadir + '*' + \
            'Individual_particle_average_velocity_' + '*' + direct + '*' + \
            fp + '*' + cn + '.txt')[0]
            int_data = pd.read_table(avg_intensity)
#            print(int_data)
            df = df.append([[cn, fp, direct, int_data['particle'].count(),\
                             int_data['average_velocity'].mean()]],\
            ignore_index=True)

#now name each column
df.columns = np.array(['cell_number', 'flagellar_pair', 'direction',\
'frequency', 'avg_velocity'])
#now plot the data based on flagellar pair
AF_fwd = AF[AF['direction']=='forward']
AF_bck = AF[AF['direction']=='backward']
CF_fwd = CF[CF['direction']=='forward']
CF_bck = CF[CF['direction']=='backward']

fwd = df[df['direction']=='forward']
bck = df[df['direction']=='backward']


#=============================================================================
#Final plots for figures    
    
with sns.axes_style('white'):
    plt.figure(figsize=(4,5))
#    #forward
    sns.barplot(x=fwd['flagellar_pair'], y=fwd['avg_velocity'], capsize=0.2, \
                linewidth=2.5, edgecolor=".2", palette='Set1')
    plt.xlabel('FP', fontsize=20)
    plt.ylabel(u'Average Particle Velocity (${\mu}m$/sec)', fontsize=16)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.ylim([0,5])
#    plt.title("IFT particle size is greater in longer flagella", fontsize=16)
    plt.tight_layout()
#    plt.savefig('180802_IFT81NG_AF_CF_fwd_velocity_barplot.svg')

with sns.axes_style('white'):
    plt.figure(figsize=(4,5))
#backward
    sns.barplot(x=bck['flagellar_pair'], y=bck['avg_velocity']*-1, capsize=0.2, \
                linewidth=2.5, edgecolor=".2", palette='Set2')
    plt.xlabel('FP', fontsize=20)
    plt.ylabel(u'Average Particle Velocity (${\mu}m$/sec)', fontsize=16)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.ylim([0,5])
#    plt.title("IFT particle size is greater in longer flagella", fontsize=16)
    plt.tight_layout()
#    plt.savefig('180802_IFT81NG_AF_CF_bck_velocity_barplot.svg')

