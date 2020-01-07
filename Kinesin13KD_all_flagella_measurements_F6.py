# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 17:42:14 2019

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

data = pd.read_excel('E:\\Dawson_Lab\\Projects\\CRISPR-Cas9\\images\\190405_k13_KD\\190405_kinesin13_KD_analysis.xlsx')

data = data.dropna()

with sns.axes_style('white'):
    plt.figure(figsize=(5,5))
    sns.barplot(x=data['FP'], y=data['length'], hue=data['tx'], \
                capsize=.2, linewidth=3, palette='Set1', edgecolor=".2")
    plt.ylabel(r'Flagellar Length ($\mu m$)', fontsize=20)
#    plt.xlabel('Flagellar Pair', fontsize=20)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.legend(loc='best', title='Region')
    plt.ylim([0, 18])
#    plt.title("Total IFT intensity increases with flagellar length", fontsize=16)
    plt.tight_layout()
#    plt.savefig('190408_kinesin13gRNA60_all_flagellar_lengths.svg')




