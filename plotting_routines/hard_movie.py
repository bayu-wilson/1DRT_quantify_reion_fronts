import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
from scipy.optimize import fsolve

#Functions
def sigma_HI_scaling(nu_ev): #divided by sigma0
    beta=2.75
    nu0_ev = 13.6 #eV
    return (nu_ev/nu0_ev)**(-beta)
def weighted_average(x,weights):
    #x and weights must be same shape
    return np.sum(x*weights,axis=0)/np.sum(weights,axis=0)
def find_root_sigma(alpha,*C):
    beta = 2.75
    num = alpha*(1-4**(-alpha-beta))
    den = (alpha+beta)*(1-4**(-alpha))
    return num/den-C

h_eV = 4.135667e-15
m_H = 1.672622e-24 #g
Y=0.24
column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]",
                "xHI","xHeI","xHeII","xHeIII","q_lya"]
rows,cols=(2,2)


skewer_number = "0000"
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

    nu,I_nu = np.loadtxt(spectrum_path).T
    nu_eV = nu*h_eV

    sigma_nu = sigma_HI_scaling(nu_eV)
    sigma_avg = weighted_average(x=sigma_nu,weights=I_nu)
    my_alpha = fsolve(func=find_root_sigma, x0=0.5,args=sigma_avg)
    alpha_list.append(my_alpha)
    length_list.append(np.interp(0.1,df_gasprops["xHI"],los_pos))
    n+=5
    gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
    spectrum_path = dir_path+"n{}_spectrum.txt".format(n)
    ax[0,0].plot(nu_eV/13.6,I_nu)
    ax[0,1].plot(los_pos,np.log10(nHI))
    ax[1,1].plot(length_list,alpha_list)

    ax[0,0].set_ylim(0,5e-23)
    ax[1,1].set_xlim(0,1e3)
    
    fig_name="movies/movie_snapshots/n{}.png".format(n)
    fig_name_list.append(fig_name)
    fig.savefig(fig_name,dpi=100)
    plt.close()

images = []
for filename in fig_name_list:
    images.append(imageio.imread(filename))
#imageio.mimsave("movies/movie_dotN{}_a={}.gif".format(dotN,spectral_index), images,duration=0.1)
imageio.mimsave("movies/movie_sk{}_hardRun.gif".format(skewer_number), images,duration=0.5)
