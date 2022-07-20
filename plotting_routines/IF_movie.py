import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import imageio
import warnings

warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")

Y = 0.24
c = 2.998e10
m_H = 1.672622e-24
h = 6.626068e-27
h_eV = 4.135667e-15
lambda_lya_cm = 1.21567e-5
pi=np.pi
kpc_to_cm = 3.086e21
rows = 3
cols = 3
panel_matrix = [["T[K]","log xHI",r"$10^{23}$I_nu"],["log nHI","j_lya","log timestep [kyr]"],[None,"q_eff_lya","I_lya/F_inc"]]
column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI","xHeII","xHeIII","q_lya"]
fig_name_list = []

initial_grid = pd.read_csv("../output_files/gas_test_hydro_000myr.txt",skiprows=1,delim_whitespace=True)
initial_nh1 = initial_grid["nH1"]
spectral_index = 2.5
gasprop_path = "../output_files/gasprops/a={:.1f}/".format(spectral_index)
spec_path = "../output_files/incident_spectra/a={:.1f}/".format(spectral_index)

#### FOR LOOP OVER n
for n in range(10,240,5):
# for n in range(220,246,2):
    fig,ax = plt.subplots(rows,cols)
    fig.set_size_inches(10,8)
    df_gasprops = pd.read_csv(gasprop_path+"n{}_gasprops.txt".format(n),delim_whitespace=True,skiprows=1,names=column_names)
    los_pos_initial = df_gasprops["los pos[pkpc]"][0]
    df_gasprops["los pos[pkpc]"] = df_gasprops["los pos[pkpc]"]-los_pos_initial
    df_spec = pd.read_csv(spec_path+"n{}_spectrum.txt".format(n),delim_whitespace=True,names=["nu","I_nu"])#,skiprows=0,names=column_names)
    with open(gasprop_path+"n{}_gasprops.txt".format(n),'r') as f:
        line = f.readline()
    step,t,dt,I_lya,Finc,nH_value,locIF,vIF_fd,vIF_fm,T_IF = np.asarray(line.split(' \t'),float)
    dt_arr = np.zeros(n+1)
    t_arr = np.zeros(n+1)
    locIF_otf = np.zeros(n+1)
    Finc_otf = np.zeros(n+1)
    I_lya_otf = np.zeros(n+1)
    nH_otf = np.zeros(n+1)
    vIF_otf = np.zeros(n+1)
    TIF_otf  = np.zeros(n+1)
    for i in range(1,n+1):
        with open(gasprop_path+"n{}_gasprops.txt".format(i),'r') as f:
            line = f.readline()
        otf_data = np.asarray(line.split(' \t'),float)
        t_arr[i]= otf_data[1]
        dt_arr[i] = otf_data[2]*1e3 #kilo-years
        I_lya_otf[i]=otf_data[3]
        Finc_otf[i] = otf_data[4]
        nH_otf[i] = otf_data[5]
        locIF_otf[i] = otf_data[6]/kpc_to_cm -los_pos_initial
        vIF_otf[i] = otf_data[8]/1e5 #km/s
        TIF_otf[i] = otf_data[9] #K
    F_lya_otf = I_lya_otf/(h*c/lambda_lya_cm)*4*pi

    mask =(t_arr>10)&(locIF_otf<280)

    los_pos = df_gasprops["los pos[pkpc]"]
    nH_arr = (1. - Y)*df_gasprops["rho"]/m_H
    nHI = df_gasprops["xHI"]*nH_arr
    frequency,intensity = df_spec["nu"],df_spec["I_nu"]
    # locIF = locIF/kpc_to_cm-los_pos_initial
    locIF_calc =  np.interp(0.5,df_gasprops["xHI"],los_pos)
    backIF = np.interp(0.01,df_gasprops["xHI"],los_pos)
    frontIF = np.interp(0.9,df_gasprops["xHI"],los_pos)
    # temperature_IF = np.interp(0.5,df_gasprops["xHI"],df_gasprops["T[K]"],)
    # frontIF,backIF = locIF+wIF/2,locIF-wIF/2


    vIF_fd = float(vIF_fd)/1e5
    vIF_fm = float(vIF_fm)/1e5
    j_lya = df_gasprops["q_lya"]*nHI*df_gasprops["ne"]*h*c/lambda_lya_cm/4/pi*1e30

    # ax[0][0].plot(los_pos,df_gasprops["T[K]"])
    ax[0][0].scatter(np.log10(vIF_otf[mask]),TIF_otf[mask],s=5)
    ax[0][1].plot(los_pos-locIF_calc,np.log10(df_gasprops["xHI"]))
    ax[0][2].plot(frequency*h_eV/13.6,intensity*1e23)
    ax[1][0].plot(los_pos,np.log10(initial_nh1),color='grey')
    ax[1][0].plot(los_pos,np.log10(nHI))
    ax[1][1].plot(los_pos-locIF_calc,j_lya)
    ax[1][2].scatter(t_arr,np.log10(dt_arr),s=5)
    ax[2][0].set_visible(False)
    # ax[2][1].set_visible(False)
    ax[2][1].plot(los_pos-locIF_calc,df_gasprops["q_lya"])
    # ax[2][2].scatter(np.log10(vIF_otf[mask]),(I_lya_otf/nH_otf)[mask],s=5)
    ax[2][2].scatter(np.log10(vIF_otf[mask]),(F_lya_otf/Finc_otf)[mask],s=5)

    # ax[2][2].scatter(np.log10(nH_otf[mask]),np.log10(I_lya_otf[mask]),s=5)

    # ax[0][0].axvline(locIF_calc,ls="dashed",color='k',alpha=0.5,lw=1)
    # # ax[0][0].set_ylim(1e4,30000)
    # ax[0][0].set_ylim(15000,30000)
    # ax[0][0].set_xlim(2.75,4.25)
    # # ax[0][0].set_ylim(-100,13000)
    ## ax[0][0].axvspan(backIF, frontIF, alpha=0.1, color='red')
    ## ax[0][0].axvline(0,ls="dashed",color='k',alpha=0.5,lw=1)
    ## ax[0][0].axvspan(backIF-locIF_calc, frontIF-locIF_calc, alpha=0.1, color='red')

    ax[0][1].axvline(0,ls="dashed",color='k',alpha=0.5,lw=1)
    ax[0][1].axvspan(backIF-locIF_calc, frontIF-locIF_calc, alpha=0.1, color='red')
    ax[0][1].set_ylim(-4,0.1)

    ## ax[0][2].set_ylim(-1.7,-.5)
    # ax[0][2].set_ylim(0.5,5)

    ax[1][0].axvline(locIF_calc,ls="dashed",color='k',alpha=0.5,lw=1)
    ax[1][0].axvspan(backIF, frontIF, alpha=0.1, color='red')
    ax[1][0].set_ylim(-5,-3)

    ax[1][1].axvline(0,ls="dashed",color='k',alpha=0.5,lw=1)
    ax[1][1].axvspan(backIF-locIF_calc, frontIF-locIF_calc, alpha=0.1, color='red')
    ax[1][1].set_ylim(-0.1,4)

    ax[1][2].set_ylim(-1,1.5)
    ax[1][2].set_xlim(-10,105)

    ax[2][1].axvline(0,ls="dashed",color='k',alpha=0.5,lw=1)
    ax[2][1].axvspan(backIF-locIF_calc, frontIF-locIF_calc, alpha=0.1, color='red')
    ax[2][1].set_ylim(0,5e-10)

    # ax[2][2].set_ylim(0,2.75e-4)#(0,4e-4) #for Ilya/nH vs log10 vIF
    # ax[2][2].set_xlim(2.9,4.25)#1.1) #for Ilya/nH vs log10 vIF

    #### # ## ax[2][2].set_ylim(1e-10,1e-8)
    # ax[2][2].set_xlim(2.75,4.25)
    # ax[2][2].set_ylim(1.5e-13,4.5e-13)
    ## ax[2][2].set_xlim(2.75,4.25)
    ## ax[2][2].set_ylim(1e-11,2e-12)
    #### # # ax[2][2].set_xscale('log')

    ax[1][0].set_xlim(-10,310)
    for j in range(rows):
        # ax[j][0].set_xlim(-10,310)
        # ax[j][1].set_xlim(0,320)
        ax[j][1].set_xlim(-20,20)
        for k in range(cols):
            ax[j][k].tick_params(bottom=True, top=True, left=True, right=False, direction="in")
            ax[j][k].text(0.05,0.55,panel_matrix[j][k],transform=ax[j][k].transAxes)

    ax[0][0].set_xlabel("log10 vIF [km s^-1]")
    ax[0][2].set_xlabel("frequency [Ryd]")
    ax[1][0].set_xlabel("proper radial position [pkpc]")
    ax[1][1].set_xlabel("relative position from IF [pkpc]")
    ax[1][2].set_xlabel("simulation time [Myr]")
    ax[2][2].set_xlabel("log10 vIF [km s^-1]")
    # ax[2][2].set_xlabel("log10 nH [cm^-3]")
    ax[2][2].text(-0.1,0.1,"n:                     {}\n".format(n)+
                          "otf time:           {:.1f} Myr\n".format(t)+
                          r"F$_{\gamma,inc}$(x$_{HI}=0.1$)"+": {:.1e} photons ".format(Finc)+ r"s$^{-1}$"+'\n'+
                          r"$n_H$(x$_{HI}=0.5$)"+":     {:.1e}".format(nH_value)+r" cm$^{-3}$"+'\n'+
                          "HI IF loc:            {:.1f} pkpc\n".format(locIF_calc)+
                          "HI IF vel (fd):      {:.1e}".format(vIF_fd)+r" km s$^{-1}$"+'\n'+
                          "HI IF vel (fm):     {:.1e}".format(vIF_fm)+r" km s$^{-1}$"+'\n'+
                          r"I$_{Ly\alpha}$"+":                    {:.1e}".format(I_lya)+" erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$"
                          ,transform=ax[2][0].transAxes)
    # plt.subplots_adjust(wspace=0, hspace=0)
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)
    fig_name="movies/movie_snapshots/n{}.png".format(n)
    fig_name_list.append(fig_name)
    fig.savefig(fig_name,dpi=100)
    plt.close()

images = []
for filename in fig_name_list:
    images.append(imageio.imread(filename))
imageio.mimsave("movies/movie_a={}.gif".format(spectral_index), images,duration=0.5)
