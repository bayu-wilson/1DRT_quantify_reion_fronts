import matplotlib.pyplot as plt
import pandas as pd
from user_input import *
from cycler import cycler
import matplotlib.ticker as ticker
from scipy.interpolate import griddata
from scipy import stats

#initial figure parameters
fontsize=17
tick_dir='in'
scaler = 1.75
lw=1.5
major_tick_size = 5*scaler
minor_tick_size = 2.5*scaler
tick_width = 1.5
colors = ['red', 'orange','gold','green','blue','purple']
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
# plt.rcParams['axes.prop_cycle'] = cycler(ls=['-', '--'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['xtick.major.width'] = tick_width*1.2
plt.rcParams['xtick.minor.width'] = tick_width
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.width'] = tick_width*1.2
plt.rcParams['ytick.minor.width'] = tick_width
plt.rcParams['xtick.labelsize'] = fontsize*1.25
plt.rcParams['ytick.labelsize'] = fontsize*1.25


uniform_density = True

#input data -
#df = pd.read_csv(
#    "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/221202/otf.csv")
#        #    "/Volumes/Extreme SSD/anson/test_expanse/master_1d_rt/plotting_routines/paper/final_data/otf2.csv")
#mask = (df["t"]>10)&(df["vIF_fm"]>1e1)&(df["vIF_fm"]<1e5)
#df = df[mask]
if uniform_density:
    df = pd.read_csv("cleaned_df_ud.csv")
    #df = pd.read_csv("/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/final_results/fd_otf.csv")
    #df["locIF"]-=1e4
    #mask = ((df["locIF"])>100)&(df["locIF"]<1000-100)&(df["F_inc"]>0)&(np.isfinite(df["F_lya"]))
    min_logvIF,max_logvIF = np.log10(2.5e2),np.log10(6e4)
    #levels = np.arange(-0.3,1.2,0.1)
    levels = np.insert(np.arange(-0.1,1.0,0.1),2,[0.01,0.03])
    label_levels =  levels #[::2]
    figname = "figures/pspace_lya_eff_ud.png"
    save_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/interp_tables/ud_parameter_space_interp_table.txt"

    alpha_list = np.unique(df["alpha_i"])
    lya_eff = df["lya_eff"].values
    v_IF = np.log10(df["vIF_Flex"].values)
    alpha_i = df["alpha_i"]
    bb = np.column_stack((v_IF,alpha_i))
    title = "Uniform-density" 
else:
    #print("start")
    df = pd.read_csv("cleaned_df_fd.csv")
    #print("loaded")
    #df = pd.read_csv("/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/final_results/ud_pp_otf_Sept20.csv") 
    #pd.read_csv("/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/final_results/ud_otf.csv")
    #df["locIF"]-=1e3
    #mask = ((df["locIF"])>100)&(df["locIF"]<1667-100)&(df["F_inc"]>0)&(np.isfinite(df["F_lya"]))
    min_logvIF,max_logvIF = np.log10(2.5e2),np.log10(6e4)
    levels = np.insert(np.arange(-0.1,1.0,0.1),2,[0.01,0.03])
    label_levels = levels
    figname = "figures/pspace_lya_eff_fd.png"
    save_path = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/interp_tables/fd_parameter_space_interp_table.txt"

    title = "Fluctuating-density"

    alpha_list = np.unique(df["alpha_i"])
    bins = 40
    total_lya_eff = [] #np.array([])

    for i in range(len(alpha_list)):    
        mask = (df["vIF_Flex"]>1e2)&(df["alpha_i"]==alpha_list[i])
        mean_lya_eff,bin_edges,binnumber = stats.binned_statistic(x=np.log10(df["vIF_Flex"][mask].values),#binning lya_eff
                           values=df["lya_eff"][mask].values,#the data on which the statistic is computed
                           statistic='mean',
                           bins=bins)
        bin_widths = np.diff(bin_edges)
        logvIF_centers = bin_edges[:-1] + bin_widths / 2
        #std_lya_eff,_,_ = stats.binned_statistic(x=np.log10(total_df_fd["vIF_Flex"][mask].values),#binning lya_eff
        #                   values=total_df_fd["lya_eff"][mask].values,#the data on which the statistic is computed
        #                   statistic='std',
        #                  bins=bins)
        total_lya_eff.append([mean_lya_eff,np.ones_like(mean_lya_eff)*alpha_list[i], logvIF_centers]) #mean lya, alpha, vi
    col_stack = np.column_stack((total_lya_eff))
    #print("stacked")
    bb = np.column_stack((col_stack[2],col_stack[1]))
    lya_eff = col_stack[0]

    #alpha_list = np.unique(df["alpha_i"])
    #lya_eff = df["lya_eff"].values
    #v_IF = np.log10(df["vIF_Flex"].values)
    #alpha_i = df["alpha_i"]
    #bb = np.column_stack((v_IF,alpha_i))

#grid_V, grid_A = np.mgrid[min_logvIF:max_logvIF:8j, 0:max_alpha:8j] #min_alpha:max_alpha:10j]
grid_V, grid_A = np.mgrid[min_logvIF:max_logvIF:25j, 0:max_alpha:25j] #min_alpha:max_alpha:10j]
#grid_V, grid_A = np.mgrid[min_logvIF:np.log10(5e+4):10j, 0:max_alpha:10j] #min_alpha:max_alpha:10j]

#bb = np.column_stack((v_IF,alpha_i))
grid_lya_eff = griddata(bb,lya_eff,(grid_V,grid_A),method='linear')
extent=(10**(min_logvIF),10**(max_logvIF),
        min_alpha,max_alpha)
#print(fratgrid)
#levels = np.arange(-0.3,1.2,0.1)
fig, ax = plt.subplots(1,figsize=(8,8))
ax.set_xscale("log")
contour = ax.contourf(10.**grid_V, grid_A, grid_lya_eff, levels=levels,extent=extent,
    cmap="Oranges",alpha=0.9)
contour2 = ax.contour(10.**grid_V, grid_A, grid_lya_eff, levels=levels,colors='k')
ax.clabel(contour2, label_levels, inline=True, fontsize=16, fmt='%.2f',colors='k')
ax.set_xlabel(r'I-front Velocity, v$_{\mathrm{IF}}$ [km/s]')
ax.set_ylabel(r'Incident Spectral Index, $\alpha_{\mathrm{IF}}$')
ax.set_title(title)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(ticker.MultipleLocator(delta_bin_alpha))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.yaxis.set_ticks_position('both')


fig.savefig(figname, bbox_inches="tight",dpi=200)
plt.close()

#### Saving interpolation table
X = np.column_stack((bb.T[0],bb.T[1],lya_eff))
#np.savetxt(fname="../results/interp_tables/test1_Dec5_2022.txt",X=X)
#np.savetxt(fname="../results/final_results/fd_parameter_space_interp_table.txt",X=X)
np.savetxt(fname=save_path,X=X)

