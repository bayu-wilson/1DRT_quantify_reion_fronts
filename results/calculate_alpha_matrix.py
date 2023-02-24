import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
import os
from scipy.optimize import fsolve

#constants
h_eV = 4.135667e-15
nu_HI = 13.6/h_eV
beta=2.75
m_H = 1.672622e-24 #g
Y=0.24

#functions
def phi_ndot(nu_i):
    return alpha/nu_HI/(1-4**(-alpha))*(nu_i/nu_HI)**(-alpha-1)
def sigmaHI_over_sigmaHI0_numeric(nu):
    return (nu/nu_HI)**(-beta)
def sigma_over_sigmaHI0_analytic(alpha):
    num = alpha * (1 - 4**(-alpha-beta))
    den = (alpha + beta)*(1 - 4**(-alpha))
    return num/den
def find_root_sigma(alpha,*C):
    return sigma_over_sigmaHI0_analytic(alpha)-C


column_names = ["los pos[pkpc]","rho","ne","dene_dt","n_tot","T[K]","xHI","xHeI","xHeII","gamma_HI","q_lya"]
los_pos_bins = np.linspace(0,2500,11)
alpha_matrix = np.zeros((40+1,10))
for i in range(0,40+1):
    skewer_number = "{:04d}".format(i)
    dir_path = "../output_files/gasprops/sk{}_hardRun/".format(skewer_number)

    print(dir_path)

    n=1
    gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
    spectrum_path = dir_path+"n{}_spectrum.txt".format(n)
    alpha_list_i=[]
    length_list_i=[]
    while os.path.exists(gasprop_path):
        df_gasprops = pd.read_csv(gasprop_path,delim_whitespace=True,skiprows=1,names=column_names)
        los_pos = df_gasprops["los pos[pkpc]"]-df_gasprops["los pos[pkpc]"].iloc[0]

        nu,I_nu = np.loadtxt(spectrum_path).T
        sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
        avg_sigma_log = np.trapz(sigma_nu*I_nu/nu,x=nu)/np.trapz(I_nu/nu,x=nu)
        my_alpha = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma_log)[0]
        alpha_list_i.append(my_alpha)
        length_list_i.append(np.interp(0.01,df_gasprops["xHI"],los_pos))
        n+=1
        gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
        spectrum_path = dir_path+"n{}_spectrum.txt".format(n)
    cut_bins = pd.cut(length_list_i,bins=los_pos_bins)
    series = pd.DataFrame({"los pos[pkpc]":length_list_i,"alpha":alpha_list_i})["alpha"] 
    alpha_matrix[i] = series.groupby(cut_bins).median().values
    #alpha_mean_i = series.groupby(cut_bins).median().values()
    #print(alpha_mean_i)
    #alpha_mean_i = pd.DataFrame({"los pos[pkpc]":length_list_i,"alpha":alpha_list_i})["alpha"].groupby(cut_bins).mean().values()
    #alpha_matrix[i] = alpha_mean_i
    #print(alpha_mean_i)
#print(alpha_matrix)

np.savetxt(fname="matrices_alpha/alpha_matrix.txt",X = alpha_matrix)

