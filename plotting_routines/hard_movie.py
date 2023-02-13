import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
from scipy.optimize import fsolve

#constants
h_eV = 4.135667e-15
nu_HI = 13.6/h_eV
beta=2.75
m_H = 1.672622e-24 #g
Y=0.24

#Functions
#"""
#def sigma_HI_scaling(nu_ev): #divided by sigma0
#    beta=2.75
#    nu0_ev = 13.6 #eV
#    return (nu_ev/nu0_ev)**(-beta)
#def weighted_average(x,weights):
#    #x and weights must be same shape
#    return np.sum(x*weights,axis=0)/np.sum(weights,axis=0)
#def find_root_sigma(alpha,*C):
#    beta = 2.75
#    num = alpha*(1-4**(-alpha-beta))
#    den = (alpha+beta)*(1-4**(-alpha))
#    return num/den-C

#def sigmaHI_over_sigmaHI0(nu):
#    return (nu/nu_HI)**(-beta)
#def sigma_logspace_over_sigma0(alpha):
#    lognu = np.log10(nu_HI)
#    lognu4 = np.log10(4*nu_HI)
#    coef = nu_HI**(beta)*(1+alpha)/(1+alpha+beta)
#    num = np.exp(-lognu4*(1+alpha+beta)*np.log(10))-np.exp(-lognu*(1+alpha+beta)*np.log(10))
#    den = np.exp(-lognu4*(1+alpha)*np.log(10))-np.exp(-lognu*(1+alpha)*np.log(10))
#    return coef*num/den
#def find_root_sigma_log(alpha,*C):
#    return sigma_logspace_over_sigma0(alpha)-C
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
#column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]",
#                "xHI","xHeI","xHeII","xHeIII","q_lya"]
rows,cols=(2,2)

skewer_number = "0003"
dir_path = "../output_files/gasprops/sk{}_hardRun/".format(skewer_number)
n=1
fig_name_list = []
gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
spectrum_path = dir_path+"n{}_spectrum.txt".format(n)
alpha_list = []
length_list = []
while os.path.exists(gasprop_path):
    fig,ax = plt.subplots(rows,cols)
    fig.set_size_inches(6,6)
    df_gasprops = pd.read_csv(gasprop_path,delim_whitespace=True,skiprows=1,names=column_names)
    los_pos = df_gasprops["los pos[pkpc]"]-df_gasprops["los pos[pkpc]"].iloc[0]
    nH_arr = (1. - Y)*df_gasprops["rho"]/m_H
    nHI = df_gasprops["xHI"]*nH_arr
    locIF_calc =  np.interp(0.5,df_gasprops["xHI"],los_pos)
    backIF = np.interp(0.10,df_gasprops["xHI"],los_pos)
    frontIF = np.interp(0.75,df_gasprops["xHI"],los_pos)
    gamma_HI_tot = df_gasprops["gamma_HI"]

    nu,I_nu = np.loadtxt(spectrum_path).T
    nu_eV = nu*h_eV

    sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
    #sigma_avg = weighted_average(x=sigma_nu,weights=I_nu)
    avg_sigma_log = np.trapz(sigma_nu*I_nu/nu,x=nu)/np.trapz(I_nu/nu,x=nu)
    my_alpha = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma_log)
    alpha_list.append(my_alpha)
    length_list.append(np.interp(0.1,df_gasprops["xHI"],los_pos))
    n+=5
    gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
    spectrum_path = dir_path+"n{}_spectrum.txt".format(n)
    ax[0,0].plot(nu_eV/13.6,I_nu)
    ax[0,1].plot(los_pos,np.log10(nHI))
    ax[1,0].plot(los_pos,np.log10(gamma_HI_tot))
    ax[1,1].plot(length_list,alpha_list)
    
    ax[0,0].set_ylim(0,5e-23)
    ax[0,1].set_xlim(0,2.5e3)
    ax[1,0].set_xlim(0,2.5e3)
    ax[1,0].set_ylim(-13.5,-11.5)
    ax[1,1].set_xlim(0,2.5e3)
    ax[1,1].set_ylim(0,1.75)
    
    ax[1,0].set_xlabel("los pos [pkpc]")
    fig_name="movies/movie_snapshots/n{}.png".format(n)
    fig_name_list.append(fig_name)
    fig.savefig(fig_name,dpi=100)
    plt.close()

images = []
for filename in fig_name_list:
    images.append(imageio.imread(filename))
#imageio.mimsave("movies/movie_dotN{}_a={}.gif".format(dotN,spectral_index), images,duration=0.1)
imageio.mimsave("movies/movie_sk{}_hardRun.gif".format(skewer_number), images,duration=0.2)
