# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 14:55:24 2017

@author: Shane
"""

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import scipy as sp
from scipy.optimize import curve_fit
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
#The purpose of this script is import and analyze the deltaT data for the AF,
#CF, and PF flagella. A previous script
#(Kymodata_analysis_allparticles_deltaT_calc_AF_CF_antero_F4.py and 
# Kymodata_analysis_allparticles_deltaT_calc_AF_PF_antero_F4)
#generated these data sets.
#==============================================================================
#import the AF data

AF_data = pd.read_csv('E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\170203_AF_anterograde_deltaT.csv')

dT_AF = AF_data['deltaT']

h_AF,b_AF = sp.histogram(dT_AF,bins=100,normed=True) # h=bin count, b=bin edge
t_AF   = (b_AF[:-1]+b_AF[1:])/2. # middle of the bins

f_AF = lambda t_AF,t0_AF: 1./t0_AF*sp.exp(-t_AF/t0_AF)

fit_AF = curve_fit(f_AF,t_AF,h_AF,p0=(1e3)) # p0 is the initial guess for a,b
t0_AF = fit_AF[0]  # fit[1] contains the error matrix
print (t0_AF)

plt.plot(t_AF/1000,h_AF*1000,'k.-')
plt.plot(t_AF/1000,f_AF(t_AF,t0_AF)*1000,'k.-')
plt.hist(dT_AF/1000, 100, normed=True, color='green')

#calculate the 95%CI for the fit parameters
sigma_af = np.sqrt(np.diagonal(fit_AF[1]))
from uncertainties import ufloat
af_error = ufloat(fit_AF[0], sigma_af[0])
print(af_error/1000) 

print(dT_AF.mean()/1000, dT_AF.sem()/1000)

af_auc = metrics.auc(t_AF/1000,f_AF(t_AF,t0_AF)*1000)
print(af_auc)
#=============================================================================

#import the CF data

CF_data = pd.read_csv('E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\161108_CF_anterograde_deltaT.csv')

dT_CF = CF_data['deltaT']

h_CF,b_CF = sp.histogram(dT_CF,bins=100,normed=True) # h=bin count, b=bin edge
t_CF   = (b_CF[:-1]+b_CF[1:])/2. # middle of the bins

f_CF = lambda t_CF,t0_CF: 1./t0_CF*sp.exp(-t_CF/t0_CF)

fit_CF = curve_fit(f_CF,t_CF,h_CF,p0=(1e3)) # p0 is the initial guess for a,b
t0_CF = fit_CF[0]  # fit[1] contains the error matrix
print (t0_CF)

plt.plot(t_CF/1000,h_CF*1000, color='orange')
plt.plot(t_CF/1000,f_CF(t_CF,t0_CF)*1000,'k.-')
plt.hist(dT_CF/1000, 100, normed=True, color='orange')

#calculate the 95%CI for the fit parameters
sigma_cf = np.sqrt(np.diagonal(fit_CF[1]))
from uncertainties import ufloat
cf_error = ufloat(fit_CF[0], sigma_cf[0])
print(cf_error/1000)

print(dT_CF.mean()/1000, dT_CF.sem()/1000)

cf_auc = metrics.auc(t_CF/1000,f_CF(t_CF,t0_CF)*1000)
print(cf_auc)
#=============================================================================

#import the PF data

PF_data = pd.read_csv('E:\\Dawson_Lab\\Projects\\Flagella\\IFT_GFP\\Live_IFT_Videos\\IFT_live_data\\data\\170203_PF_anterograde_deltaT.csv')

dT_PF = PF_data['deltaT']

h_PF,b_PF = sp.histogram(dT_PF,bins=100,normed=True) # h=bin count, b=bin edge
t_PF   = (b_PF[:-1]+b_PF[1:])/2. # middle of the bins

f_PF = lambda t_PF,t0_PF: 1./t0_PF*sp.exp(-t_PF/t0_PF)

fit_PF = curve_fit(f_PF,t_PF,h_PF,p0=(1e3)) # p0 is the initial guess for a,b
t0_PF = fit_PF[0]  # fit[1] contains the error matrix
print (t0_PF)

plt.plot(t_PF,h_PF*1000, color='orange')
plt.plot(t_PF/1000,f_PF(t_PF,t0_PF)*1000,'k.-')
plt.hist(dT_PF/1000, 100, normed=True, color='orange')

#calculate the 95%CI for the fit parameters
sigma_pf = np.sqrt(np.diagonal(fit_PF[1]))
from uncertainties import ufloat
pf_error = ufloat(fit_PF[0], sigma_pf[0])
print(pf_error/1000)

print(dT_PF.mean()/1000, dT_PF.sem()/1000)

pf_auc = metrics.auc(t_PF/1000,f_PF(t_PF,t0_PF)*1000)
print(pf_auc)

#pf_upper = f_PF(dT_PF/1000, (fit_PF[0] + sigma_pf))
#pf_lower = f_PF(dT_PF/1000, (fit_PF[0] - sigma_pf))

#=============================================================================

#FINAL PLOTS 

    
with sns.axes_style('white'):  
    plt.figure(figsize=(5,5))
    plt.hist(dT_AF/1000, 100, normed=True, color='blue', edgecolor='None',\
             label='Anterior')
    plt.plot(t_AF/1000,f_AF(t_AF,t0_AF)*1000, color='k')   
    
#    plt.hist(dT_CF/1000, 100, normed=True, color='orange', edgecolor='None',\
#             label='Caudal')
#    plt.plot(t_CF/1000,f_CF(t_CF,t0_CF)*1000, color='k')
#    
#    plt.hist(dT_PF/1000, 100, normed=True, color='green', edgecolor='None',\
#             label='Posteriolateral')
#    plt.plot(t_PF/1000,f_PF(t_PF,t0_PF)*1000,'k')
    
    
    plt.xlabel('Delta Time (sec)', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.legend(loc='upper right', title='Flagellar Pair', fontsize=20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.xlim([0,15])
    plt.ylim([0,1])
    plt.tight_layout()    
#    plt.savefig('190415_IFT81NG_AF_deltaT_distribution.svg')
#==============================================================================    
param_data = pd.DataFrame()

param_data = pd.read_excel('E:\\Dawson_Lab\\Projects\\Flagella\\Theory\\181009_IFT81mNG_freqdist_fit_stats.xlsx')

x_pos = np.arange(len(param_data['fp']))
fit_app = param_data['freq_mean']
fit_app_error = param_data['freq_CI']


#first plot everything on the same plot
with sns.axes_style('white'):
    fig, ax = plt.subplots(figsize=(5,6))
    plt.bar(x_pos[:], fit_app[:], yerr=fit_app_error[:], align='center', edgecolor='black', capsize=10, linewidth=3)
    plt.errorbar(x_pos[:], fit_app[:], yerr=fit_app_error[:], capsize=10, elinewidth=3, markeredgewidth=3, fmt='.', color='k')
    ax.set_ylabel('Time-lag between injections (sec)', fontsize=24)
    ax.set_xlabel('Flagellar Pair', fontsize=24)
    ax.set_xticks(x_pos[:])
    ax.set_xticklabels(param_data['fp'][:], fontsize=20)
    plt.ylim([0,1.6])
    plt.tight_layout()
#    plt.savefig('181029_AF_PF_CF_freqdist_fits_summary.svg')       
    
