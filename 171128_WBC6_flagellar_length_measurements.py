# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:42:43 2017

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
#the goal is to analyze the pooled FRAP data from 12/15/16 and 1/6/17 and to 
#generate figures for this data
#=============================================================================
#first import the data into a dataframe
data = pd.DataFrame()

data = pd.read_csv('D:\\Dawson_Lab\\Projects\\Flagella\\Drugs\\WBC6_flagellar_length.csv')

data = data.dropna()

AF = data[data['flagellar_pair']=='Anterior']
CF = data[data['flagellar_pair']=='Caudal']
PF = data[data['flagellar_pair']=='Posteriolateral']
VF = data[data['flagellar_pair']=='Ventral']

print(AF['len'].mean(),AF['len'].sem())
print(CF['len'].mean(),CF['len'].sem())
print(PF['len'].mean(),PF['len'].sem())
print(VF['len'].mean(),VF['len'].sem())


from scipy.stats import sem
sem(AF['len'])
sem(CF['len'])
sem(PF['len'])

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    ax = sns.barplot(x=data['flagellar_pair'], y=data['len'], capsize=.2, linewidth=3, palette='Set1', edgecolor=".2")
    plt.ylabel(r'Flagellar Length ($\mu m$)', fontsize=20)
    plt.xlabel('Flagellar Pair', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.setp(ax.get_xticklabels(), rotation=45)
    plt.legend(loc='best', title='Treatment')
    plt.ylim([0, 18])
#    plt.title("Total IFT intensity increases with flagellar length", fontsize=16)
    plt.tight_layout()
    plt.savefig('181025_WBC6_flagellar_length.svg')
