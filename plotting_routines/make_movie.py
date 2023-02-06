import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import imageio
import warnings
import sys

skewer_number = str(sys.argv[1])
spectral_index = str(sys.argv[2]) # "1.500"
warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")

delta_n = 50
z = 5.7
omega_b = 0.048
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rhoBaryons_0 = rho_crit_z*omega_b #g/cm3
m_H = 1.672622e-24 #g
Y = 0.24 #helium fraction
h = 6.626068e-27
c = 2.998e10
pi=np.pi
lambda_lya_cm = 1.21567e-5

#spectral_index = "1.500"
#dotN = "2.51189e+56"
#dir_path = "../output_files/gasprops/dotN{}_a={}/".format(dotN,spectral_index)
#dir_path = "../output_files/gasprops/sk{}_a={}/".format(skewer_number,spectral_index)
dir_path = "../output_files/gasprops/sk{}_hardRun/".format(skewer_number)
fontsize=10

panel_matrix = [["IF width [kpc]",r"F$_{inc}$",r"$x_{HI}$(r)",r"$n_e \times n_p$ (r)"],
                [r"log nHI(r) [cm$^{-3}$].", r"F$_{Ly\alpha}$ (r)",r"log nHI(r) [cm$^{-3}$]",r"$q^eff_{Ly\alpha}$"],
                [None,r"F$_{Ratio}$ (r)","T(r) [K]",r"j$_{Ly\alpha}$(r)"]]
column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI","xHeII","xHeIII","q_lya"]
fig_name_list = []

n=1
gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
otf_matrix = []
while os.path.exists(gasprop_path_otf):
    with open(gasprop_path_otf,'r') as f:
        line = f.readline()
    try:
        line_array = np.asarray(line.split(' \t'),float)
        otf_matrix.append(line_array)
        if line_array[1]<10:
            n_minTime = n
    except:
        otf_matrix.append(np.ones(13)*np.nan)
    n+=1
    gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
table_otf = np.array(otf_matrix)
#print(table_otf)
#mask = table_otf.T[1]>10
step_arr,t_arr, dt_arr,I_lya_otf,Finc_otf,clump_otf,locIF_otf,vIF_fd_otf,vIF_fm_otf,TIF_otf,rho_var_otf, width_otf, nH_avg_otf = table_otf.T#[mask]
F_lya_otf = I_lya_otf/(h*c/lambda_lya_cm)*4*pi

initial_grid = pd.read_csv(dir_path + "n{}_gasprops.txt".format(1),
                           skiprows=1,delim_whitespace=True,names = column_names)
initial_nh1 = (1. - Y)*initial_grid["rho"]/m_H * initial_grid["xHI"]

n=1
gasprop_path = dir_path + "n{}_gasprops.txt".format(n)
rows,cols = 3,4
while os.path.exists(gasprop_path):
    fig,ax = plt.subplots(rows,cols)
    fig.set_size_inches(12,10)
    df_gasprops = pd.read_csv(gasprop_path,delim_whitespace=True,skiprows=1,names=column_names)
    los_pos = df_gasprops["los pos[pkpc]"]-df_gasprops["los pos[pkpc]"].iloc[0]
    nH_arr = (1. - Y)*df_gasprops["rho"]/m_H
    nHI = df_gasprops["xHI"]*nH_arr
    locIF_calc =  np.interp(0.5,df_gasprops["xHI"],los_pos)
    backIF = np.interp(0.05,df_gasprops["xHI"],los_pos)
    frontIF = np.interp(0.75,df_gasprops["xHI"],los_pos)
    with open(gasprop_path,'r') as f:
        line = f.readline()
    step_n,t_n,dt_n,I_lya_n,Finc_n,clump_value_n,locIF_n,vIF_fd_n,vIF_fm_n,T_IF_n, rho_var_n, width_n, nH_avg_n = np.asarray(line.split(' \t'),float)
    vIF_fd = float(vIF_fd_n)/1e5
    vIF_fm = float(vIF_fm_n)/1e5
    j_lya = df_gasprops["q_lya"]*nHI*df_gasprops["ne"]*h*c/lambda_lya_cm/4/pi#*1e30
    
    
#     ax[0][0].plot(los_pos,df_gasprops["T[K]"])
    ax[0][0].scatter(width_otf[:n],F_lya_otf[:n],s=1)
    ax[1][0].plot(los_pos,np.log10(initial_nh1),color='gray')
    ax[1][0].plot(los_pos,np.log10(nHI))
    ax[0][1].scatter(np.log10(vIF_fm_otf[n_minTime:n]/1e5),Finc_otf[n_minTime:n],s=1)
    ax[1][1].scatter(np.log10(vIF_fm_otf[n_minTime:n]/1e5),F_lya_otf[n_minTime:n],s=1)
    ax[2][1].scatter(np.log10(vIF_fm_otf[n_minTime:n]),(F_lya_otf/Finc_otf)[n_minTime:n],s=1)
    ax[0][2].plot(los_pos-locIF_calc,df_gasprops["xHI"])
    ax[1][2].plot(los_pos-locIF_calc,nHI)
    ax[2][2].plot(los_pos-locIF_calc,df_gasprops["T[K]"])
    ax[2][0].set_visible(False)
    ax[0][3].plot(los_pos-locIF_calc,nHI*df_gasprops["ne"])
    ax[1][3].plot(los_pos-locIF_calc,df_gasprops["q_lya"])
    ax[2][3].plot(los_pos-locIF_calc,j_lya)
    
    ax[1][0].set_xlabel("skewer position [pkpc]")
    ax[1][0].axvline(locIF_calc,ls="dashed",color='k',alpha=0.5,lw=1)
    ax[1][0].axvspan(backIF, frontIF, alpha=0.1, color='red')
    ax[1][0].set_ylim(-5,-3.75)
    ax[2][1].set_xlabel("log10 vIF [km s^-1]")
    ax[2][2].set_xlabel("relative position from IF [pkpc]")
    ax[2][3].set_xlabel("relative position from IF [pkpc]")

    
    

    ax[2][2].text(-0.1,0.1,"n:                     {}\n".format(n)+
                      "otf time:           {:.1f} Myr\n".format(t_n)+
                      r"F$_{\gamma,inc}$(x$_{HI}=0.1$)"+": {:.1e} photons ".format(Finc_n)+ r"s$^{-1}$"+'\n'+
                      "clumping factor: {:.3f}\n".format(clump_value_n)+# +
                      #r"$n_H$(x$_{HI}=0.5$)"+":     {:.1e}".format(nH_value)+r" cm$^{-3}$"+'\n'+
                      "HI IF loc:            {:.1f} pkpc\n".format(locIF_calc)+
                      "HI IF vel (fd):      {:.1e}".format(vIF_fd)+r" km s$^{-1}$"+'\n'+
                      "HI IF vel (fm):     {:.1e}".format(vIF_fm)+r" km s$^{-1}$"+'\n'+
                      r"I$_{Ly\alpha}$"+":                    {:.1e}".format(I_lya_n)+" erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$"
                      ,transform=ax[2][0].transAxes)
    for j in range(rows):
    # ax[j][0].set_xlim(-10,310)
    # ax[j][1].set_xlim(0,320)
        ax[j][2].set_xlim(-25,25)
        ax[j][3].set_xlim(-25,25)
        ax[j][2].axvline(0,ls="dashed",color='k',alpha=0.5,lw=1)
        ax[j][3].axvline(0,ls="dashed",color='k',alpha=0.5,lw=1)
        for k in range(cols):
            ax[j][k].tick_params(bottom=True, top=True, left=True, right=False, direction="in")
            ax[j][k].text(0.05,0.55,panel_matrix[j][k],transform=ax[j][k].transAxes,fontsize=fontsize)
    
    fig_name="movies/movie_snapshots/n{}.png".format(n)
    n+=50
    gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
#     plt.tight_layout()
    fig_name="movies/movie_snapshots/n{}.png".format(n)
    fig_name_list.append(fig_name)
    fig.savefig(fig_name,dpi=100)
    plt.close()

images = []
for filename in fig_name_list:
    images.append(imageio.imread(filename))
#imageio.mimsave("movies/movie_dotN{}_a={}.gif".format(dotN,spectral_index), images,duration=0.1)
imageio.mimsave("movies/movie_sk{}_a={}.gif".format(skewer_number,spectral_index), images,duration=0.1)    
