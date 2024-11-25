import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
import imageio
import warnings
from functions import get_inc_LyC_flux

#user input
ud_or_fd = "fd_pp"
delta_n = 2

#ignoring annoying warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in log10")
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")

#constants
rho_crit_0 = 8.688e-30 #critical density at z=0 [g cm^-3]
Omega_b = 0.048
z=5.7
m_H = 1.672622e-24 #g
Y = 0.24
h = 6.626068e-27
c = 2.998e10
pi=np.pi
lambda_lya_cm = 1.21567e-5
kpc_to_cm = 3.086e21
z = 5.7
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
omega_b = 0.048
rhoBaryons_0 = rho_crit_z*omega_b #g/cm3

#plotting parameters
fontsize=15
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
# colors = ['red', 'orange','gold','green','blue','purple','cyan','pink']
plt.rcParams['font.family'] = "serif"
#plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
# plt.rcParams['axes.prop_cycle'] = cycler(ls=['-', '--'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size

### making on-th-fly data table
if ud_or_fd == "fd":
    dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/sk0001_a=1.500_Lum=4.0e+46/"
    R0 = 1e4
    R = 1000
    jmax = 1e-29
    zoom_bounds = 60
elif ud_or_fd == "fd_pp": #fd_pp_a\=0.0_sk0007_L\=3.08e+45
    #dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=0.0_sk0007_L=3.08e+45/"
    #dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=0.0_sk0017_L=4.13e+44/"
    ###dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=0.0_sk0019_L=4.52e+44/"
    #dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=0.0_sk0001_L=6.30e+45/"
    #dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=0.0_sk0025_L=3.09e+44/"
    #dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=1.5_sk0044_L=7.07e+44_v2/"
    dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/fd_pp_a=1.5_sk0043_L=7.26e+46_v2/"
    #"sk0056_L=5.76e+45_v2/"
    ###drwxrwsr-x 2 bwils033 uot171  57344 Sep 25 12:18 'fd_pp_a=0.0_sk0002_L=1.63e+46'
    R0 = 1e4
    R = 1500
    jmax = 1e-29
    zoom_bounds = 60

elif ud_or_fd == "ud":
    dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/ud_a=1.500_Lum=8.42e+44/"
    R0 = 1e3
    R = 1667
    jmax = 1e-29
    zoom_bounds = 10
elif ud_or_fd == 'ud_pp':
    #dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/ud_pp_a=1.5_L=5.50e+46/"
    dir_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/output_files/gasprops/ud_pp_a=1.5_L=1.51e+47/"
    R0 = 1e4
    R = 940 
    #R = 760
    #1.5 5.50e+46 759.8 1519 1.0e+04 565.7
    jmax = 1e-29
    zoom_bounds = 50
else:
    print("Something is wrong")

n=1
n_minTime=0
gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
otf_matrix = []
while os.path.exists(gasprop_path_otf):
    with open(gasprop_path_otf,'r') as f:
        line = f.readline()
    try:
        line_array = np.asarray(line.split(' \t'),float)
        otf_matrix.append(line_array)
        #if line_array[1]<10:
        #    n_minTime = n
    except:
        otf_matrix.append(np.ones(13)*np.nan)
    n+=1
    gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
table_otf = np.array(otf_matrix)
#mask = (df_fd["locIF"]>100)&(df_fd["locIF"]<1000-100)&(df_fd["F_inc"]>0)
#print(table_otf)
#mask = table_otf.T[1]>10
step_arr,t_arr, dt_arr,F_lya_otf,Finc_otf,locIF_otf,vIF_fd_otf,vIF_Flex_otf,TIF_otf, Tc_otf, width_otf, nH_center_otf, trash_otf = table_otf.T#[mask]
#print(locIF_otf/kpc_to_cm)
locs = locIF_otf/kpc_to_cm-R0
mask = (locs>100)&(locs<R-100)&(Finc_otf>0)
step_arr,t_arr, dt_arr,F_lya_otf,Finc_otf,locIF_otf,vIF_fd_otf,vIF_Flex_otf,TIF_otf, Tc_otf, width_otf, nH_center_otf, trash_otf = table_otf[mask].T
# plotting the properties along the whole line of sight

column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]",
                "xHI","xHeI","xHeII","xHeIII","q_lya"]

initial_grid = pd.read_csv(dir_path + "n{}_gasprops.txt".format(1),
                           skiprows=1,delim_whitespace=True,names = column_names)
initial_nh1 = (1. - Y)*initial_grid["rho"]/m_H * initial_grid["xHI"]

#file_name = "../output_files/gasprops/sk0001_a=1.500_Lum=1.0e+45/n100_gasprops.txt"
#n=100

fig_name_list = []
n=1
gasprop_path = dir_path + "n{}_gasprops.txt".format(n)
while os.path.exists(gasprop_path):
    gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
    #print(gasprop_path_otf)
    df = pd.read_csv(gasprop_path_otf,delim_whitespace=True,skiprows=1,names=column_names)
    r_los = df["los pos[pkpc]"]-R0
    xHI = df["xHI"]
    rho = df["rho"]
    nHI = (1 - Y)*rho/m_H*xHI
    IF_loc = np.interp(x=0.5,xp=xHI,fp=r_los)
    frontIF = np.interp(0.99,df["xHI"],r_los)
    backIF = np.interp(1e-3,df["xHI"],r_los)
    T_skew = df["T[K]"]
    jlya = df["ne"]*nHI*df["q_lya"]*h*c/lambda_lya_cm /4/pi

    fig,axs = plt.subplots(nrows=2, ncols=3, figsize=(12, 6), constrained_layout=True,
                            gridspec_kw={'width_ratios': [1, 1, 2],'height_ratios': [1, 1]})

    axs[0,0].plot(r_los,rho/rhoBaryons_0,color="k")
    axs[1,0].plot(r_los-IF_loc,xHI,color="k")
    #axs[1,0].semilogy(r_los,nHI,color="k")
    axs[0,1].semilogy(r_los-IF_loc,T_skew,color="k")
    axs[1,1].semilogy(r_los-IF_loc,jlya,color="k")

    gs = axs[0, 0].get_gridspec()
    ax_combined = fig.add_subplot(gs[:, -1])

    for ax in axs[:, -1]:
        ax.remove()
    
    if ud_or_fd == "fd_pp":
        axs[0,0].set_ylim(-15,48)
    axs[0,0].axvline(IF_loc,color='k',ls="--",alpha=0.5)
    axs[0,0].axvspan(backIF,frontIF,color="red",alpha=0.3)
    axs[0,0].set_xlabel(r"R$_{\mathrm{sim}}$[pkpc]")
    #axs[0,0].set_yscale("log")

    #axs[1,0].set_ylim(1e-4,1e+1)
    axs[1,0].set_ylim(-0.05,1.05)
    #axs[1,0].set_xlabel(r"R$_{\mathrm{sim}}$[pkpc]")

    axs[0,1].set_ylim(2e1,4e4)

    #axs[1,1].set_ylim(1e-8,1e-2)
    #axs[1,1].set_xlabel(r"R$_{\mathrm{sim}}$[pkpc]")

    axs[1,1].set_ylim(1e-35,jmax)

    ax_combined.set_xlabel(r"IF speed, log$_{10}$ v$_\mathrm{IF}$ [km/s]")
    # ax_combined.set_ylabel(r"Photon Flux ratio, $\zeta(\alpha, v_{\mathrm{IF}})$")

    if ud_or_fd == "fd_pp":
        incLyC_Flux = get_inc_LyC_flux(Lum=7.26e+46,#df_i["Lum"].values,
                               spectral_index=1.5,
                               R0_cm=R0*kpc_to_cm,
                               t_Myr=None,t0_Myr=None,eta=0.0)
        ax_combined.scatter(vIF_Flex_otf[n_minTime:n]/1e5,(F_lya_otf/incLyC_Flux)[n_minTime:n],s=1,color='k')
    else:
        ax_combined.scatter(vIF_Flex_otf[n_minTime:n]/1e5,(F_lya_otf/Finc_otf)[n_minTime:n],s=1,color='k')
    ax_combined.set_xscale("log")
    ax_combined.set_xlim(2e+2,5e+4)
    #ax_combined.set_xlim(1e+0,3e+5)
    ax_combined.set_ylim(-0.05,1)
    
    labels = [[r"Gas overdensity, $\Delta$","Temperature [K]"],[r"Neutral fraction, x$_\mathrm{HI}$", r"Ly$\alpha$ emissivity, $j_{\mathrm{Ly}\alpha}$"]]
    for i in range(2):
        for j in range(2):
            ax = axs[i,j]
            if not (i==0 and j==0):
                ax.set_xlim(-zoom_bounds,zoom_bounds)
                ax.axvline(0,color='k',ls="--",alpha=0.5)
                ax.axvspan(backIF-IF_loc,frontIF-IF_loc,color="red",alpha=0.3)
                ax.set_xlabel(r"R$_{\mathrm{sim}} - R_\mathrm{IF}$[pkpc]")
            ax.text(0.05, 0.05, f'{labels[i][j]}', transform=ax.transAxes,
                verticalalignment='bottom', bbox=dict(facecolor='grey', alpha=0.7))        
    ax_combined.text(0.05, 0.95, r"Ly$\alpha$ photon efficiency, $\zeta(\alpha, v_{\mathrm{IF}})$", transform=ax_combined.transAxes,
                    verticalalignment='top', bbox=dict(facecolor='grey', alpha=0.7))

    n+=delta_n
    gasprop_path = dir_path+"n{}_gasprops.txt".format(n)
    fig_name="movies/movie_snapshots/n{}.png".format(n)
    fig_name_list.append(fig_name)
    fig.savefig(fig_name,dpi=100,bbox_inches="tight")
    plt.close()

images = []
for filename in fig_name_list:
    images.append(imageio.imread(filename))
imageio.mimsave(f"movies/qual_movie_{ud_or_fd}_test.gif", images, duration=0.1)
#imageio.mimsave("movies/movie_dotN{}_a={}.gif".format(dotN,spectral_index), images,duration=0.1)
#imageio.mimsave("movies/movie_sk{}_a={}_Lum={}.gif".format(skewer_number,spectral_index,Lum), images,duration=0.1)

#fig.savefig("figures/test_movie.png", bbox_inches="tight")
#plt.close()
