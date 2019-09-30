# -*- coding: utf-8 -*-
"""

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy
import scipy.stats
from math import pi
import glob
import statsmodels.stats.api as sms
#import matplotlib for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
#import os to handle operating system
import os
#=============================================================================
#This script is used to import flagellar length measurements for WBC6 and 
#mNG-Btub trophozoites and plot the means and 95%CI for each of these datasets.
#Also calculate the means and sem for statistical analyses.  
#=============================================================================
#first import the mNG-Btub data

mNG_data = pd.DataFrame()

mNG_data = pd.read_excel('E:\\Dawson_Lab\\Projects\\Flagella\\Drugs\\190402_mNGBtub_cytoaxo&mem_lengths.xlsx')

#now plot the data as a barplot
with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    sns.barplot(x=mNG_data['FP'], y=mNG_data['length'], hue=mNG_data['region'],\
                capsize=.2, linewidth=3, palette='Set1', edgecolor=".2")
    plt.ylabel(r'Flagellar Length ($\mu m$)', fontsize=20)
    plt.xlabel('Flagellar Pair', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.legend(loc='best', title='Region')
    plt.title('Cytoplasmic and membrane-bound flagellar length of mNG-Btub')
    plt.ylim([0, 18])
    plt.tight_layout()
#    plt.savefig('190507_mNGBtub_full_flagellar_length.svg')

#now repeat for the WBC6 data
WBC6_data = pd.DataFrame()

WBC6_data = pd.read_csv('E:\\Dawson_Lab\\Projects\\Flagella\\Drugs\\WBC6_flagellar_length.csv')

WBC6_data = data.dropna()

AF = WBC6_data[WBC6_data['flagellar_pair']=='Anterior']
CF = WBC6_data[WBC6_data['flagellar_pair']=='Caudal']
PF = WBC6_data[WBC6_data['flagellar_pair']=='Posteriolateral']
VF = WBC6_data[WBC6_data['flagellar_pair']=='Ventral']

print(AF['len'].mean(),AF['len'].sem())
print(CF['len'].mean(),CF['len'].sem())
print(PF['len'].mean(),PF['len'].sem())
print(VF['len'].mean(),VF['len'].sem())

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    ax = sns.barplot(x=WBC6_data['flagellar_pair'], y=WBC6_data['len'], capsize=.2,\
                     linewidth=3, palette='Set1', edgecolor=".2")
    plt.ylabel(r'Flagellar Length ($\mu m$)', fontsize=20)
    plt.xlabel('Flagellar Pair', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.legend(loc='best', title='Treatment')
    plt.ylim([0, 18])
    plt.title("Membrane-bound flagellar length of WBC6", fontsize=16)
    plt.tight_layout()
#    plt.savefig('181025_WBC6_flagellar_length.svg')





