import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.patheffects as mpe
from scipy.interpolate import interp1d
from functions import get_inc_LyC_flux
from constants import kpc_to_cm
import user_input as ui
from matplotlib.lines import Line2D
from scipy import stats

#initial figure parameters
#outline=mpe.withStroke(linewidth=4, foreground='black')
#s = 5
#fontsize=15
color_list = ["red","orange","gold","green","blue","purple","pink","black","brown","cyan"]
fontsize=15
tick_dir='in'
scaler = 1.75
lw=1.5

labels = [r"$\alpha$={:.1f}".format(float(i)) for i in ui.bincenters_alpha]
                # Line2D([0], [0], lw=0)]
#labels.insert(-3, '')

major_tick_size = 3*scaler #5*scaler
minor_tick_size = 1.5*scaler #2.5*scaler
tick_width = 1 #1.5
#colors = ['red', 'orange','gold','green','blue','purple']
plt.rcParams['font.family'] = "serif"
#plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
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
plt.rcParams['xtick.labelsize'] = fontsize#*1.25
plt.rcParams['ytick.labelsize'] = fontsize#*1.25
custom_lines = [Line2D([0], [0], color=color_list[0], lw=9, marker=None),
                Line2D([0], [0], color=color_list[1], lw=9, marker=None),
                Line2D([0], [0], color=color_list[2], lw=9, marker=None),
                Line2D([0], [0], color=color_list[3], lw=9, marker=None),
                #Line2D([0], [0],lw=0,marker=None), 

                Line2D([0], [0], color=color_list[4], lw=9, marker=None),
                Line2D([0], [0], color=color_list[5], lw=9, marker=None),
                Line2D([0], [0], color=color_list[6], lw=9, marker=None)]
                #Line2D([0], [0], color='k', lw=3, ls='--')]



### paths
dir_home = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/"
dir_input_params = dir_home+"parameters_input/input_params/"
dir_results = dir_home+"results/final_results/"

### load in uniform density data
R0_pkpc=1e4 #pkpc
df_ud_pp = pd.read_csv(dir_results+"ud_pp_otf_Sept20.csv") #results/final_results/ud_pp_otf.csv 
params_ud_pp = np.loadtxt(dir_input_params+"planeParallel_params.txt")
df_ud_pp["locIF"]-=R0_pkpc
mask=(df_ud_pp["locIF"]>0)&(df_ud_pp["F_inc"]>0)&(df_ud_pp["vIF_Flex"]>0)#&(df_pp["locIF"]>100)
# (df_pp["Lum"]==L_list[4])&(df_pp["locIF"]>100)
df_ud_pp = df_ud_pp[mask]

### load in fluctuating density data
df_fd_pp = pd.read_csv(dir_results+"fd_pp_otf_v3.csv")#"fd_pp_otf.csv")
#params_fd_pp = np.loadtxt(dir_input_params+"fd_planeParallel_params.txt")
params_fd_pp = np.loadtxt(dir_input_params+"fd_April_params.txt")
alpha_list = df_fd_pp["alpha_i"].unique()
df_fd_pp["locIF"]-=R0_pkpc
df_fd_pp = df_fd_pp[((df_fd_pp["locIF"])>100)&(df_fd_pp["F_inc"]>0)]

fig,ax = plt.subplots(1,2,figsize=(10,5))

newcols = list(df_fd_pp.columns)#.append("Lya Eff")
newcols.append("lya_eff")
total_df_fd = pd.DataFrame(columns=newcols)
### plotting fluctuating density
for i in range(0,len(params_fd_pp),3):
    color_i = color_list[i//91] #color_list[i] #color_list[i//78]
    alpha_i = params_fd_pp[i,0] #alpha_list[i]#fd_pp_params[i,0]
    skewer_i = int(params_fd_pp[i,1]) #np.asarray(fd_pp_params.T[1],int)
    Lum_i = params_fd_pp[i,2]
   
    # the dataframe for the i'th skewer at a certain alpha only. Ignore 94 and 95 because they have weird signatures.
    df_i = df_fd_pp[(df_fd_pp["alpha_i"]==alpha_i)&(df_fd_pp["skewer"]==skewer_i)]#&(skewer_i not in [94,95])] 
                 
    max_locIF = np.min([np.max(df_i["locIF"]),1500]) #the max lenght of the skewer is 1.5 Mpc. But sometimes the I-front doesn't get that far. 
    #We just want a little buffer in case anything weird is going on
    df_i = df_i[df_i["locIF"]<max_locIF-100] #100 because it is a larger than the width of a wide I-front 
    
    ax[0].scatter(df_i["vIF_Flex"],df_i["T_avg"],s=1,color=color_i,label="alpha={:.1f}".format(alpha_i),alpha=0.3)
    incLyC_Flux = get_inc_LyC_flux(Lum=df_i["Lum"].values,
                               spectral_index=alpha_i,
                               R0_cm=R0_pkpc*kpc_to_cm,
                               t_Myr=None,t0_Myr=None,eta=0.0)

    ax[1].scatter(df_i["vIF_Flex"],df_i["F_lya"]/incLyC_Flux,s=1,color=color_i, alpha = 0.3)

    lya_eff_fd = df_i["F_lya"]/incLyC_Flux
    df_i["lya_eff"] = lya_eff_fd
    total_df_fd = pd.concat((total_df_fd,df_i))

#print("Lum fd :", df_i["Lum"].values)

#print("alpha_i:", alpha_i)
#print("F_lya:", df_i["F_lya"].values)
#print("F_lyC:",incLyC_Flux)
#print("Lya flux:", df_i["F_lya"])
#print("LyC Flux:", incLyC_Flux)

### cleaning uniform density data and adding it to total_df_ud
newcols = list(df_ud_pp.columns)#.append("Lya Eff")
newcols.append("lya_eff")
total_df_ud = pd.DataFrame(columns=newcols)
for i in range(len(params_ud_pp)):
    # color_i = color_list[i//7]
    alpha_i = params_ud_pp[i,0]
    #color_i = color_list[np.where(alpha_list==alpha_i)[0][0]]
    Lum_i = np.char.mod("%.2e",(params_ud_pp[i,1]))
    R_sim_i = params_ud_pp[i,2]
    t0_i =  params_ud_pp[i,6]
    mask_i = (
        df_ud_pp["alpha_i"] == alpha_i) & (
        np.char.mod("%.2e",df_ud_pp["Lum"]) == Lum_i)&(
        # asdf["locIF"]<R_sim_i*0.95) & (
        df_ud_pp["locIF"]<R_sim_i*0.90) & (
        df_ud_pp["locIF"]>R_sim_i*0.20)  #& (
    #print(R_sim_i, Lum_i, alpha_i, np.unique(df_ud_pp["Lum"]))    
    df_i = df_ud_pp[mask_i].copy()
    
    # ax[0].scatter(qwer["vIF_Flex"],qwer["T_avg"],color=color_i,s=5)
    # ax[0].plot(qwer["vIF_Flex"],qwer["T_avg"],color=color_i,ls='dashed',lw=2)
    incLyC_Flux_ud = get_inc_LyC_flux(Lum=df_i["Lum"].values,
                               spectral_index=alpha_i,
                               R0_cm=R0_pkpc*kpc_to_cm,
                               t_Myr=df_i["t"],
                               t0_Myr=t0_i,eta=2)
    lya_eff_ud = df_i["F_lya"]/incLyC_Flux_ud
    df_i["lya_eff"] = lya_eff_ud
    total_df_ud = pd.concat((total_df_ud,df_i))
#print(total_df_ud)
#print(df_i["Lum"].values)
#print(df_i["t"])
#print(t0_i)

#print(df_ud_pp)
#print(total_df_ud)

#plotting uniform density
for i in range(len(alpha_list)):
    color_i = color_list[i]
    df_i = total_df_ud[total_df_ud["alpha_i"] == alpha_list[i]]
    vIF_interpolate = np.logspace(np.log10(min(df_i["vIF_Flex"])),np.log10(max(df_i["vIF_Flex"])),100)
    f_interp_T_avg = interp1d(df_i["vIF_Flex"], df_i["T_avg"], kind='linear', fill_value='extrapolate')
    f_interp_lya_eff = interp1d(df_i["vIF_Flex"], df_i["lya_eff"], kind='linear', fill_value='extrapolate')
    T_avg_interpolated = f_interp_T_avg(vIF_interpolate)
    lya_eff_interpolated = f_interp_lya_eff(vIF_interpolate)
    ax[0].plot(vIF_interpolate,T_avg_interpolated,color="k",ls='dotted',lw=lw)#,path_effects=[outline])
    ax[1].plot(vIF_interpolate,lya_eff_interpolated,color="k",ls='dotted',lw=lw)#,path_effects=[outline])

custom_lines.append(Line2D([0], [0], color='k', lw=2,ls='dotted', marker=None))
labels.append(r"Uni. $\rho$")

bins=20
flip_flop=1
for i in range(len(alpha_list)):
    cleaned_df_fd = total_df_fd[total_df_fd["alpha_i"] == ui.bincenters_alpha[i]]
    mean_lya_eff,bin_edges,binnumber = stats.binned_statistic(x=np.log10(cleaned_df_fd["vIF_Flex"].values),#binning lya_eff
                       values=cleaned_df_fd["lya_eff"].values,#the data on which the statistic is computed
                       statistic='mean',
                       bins=bins)
    bin_widths = np.diff(bin_edges)
    logvIF_centers = bin_edges[:-1] + bin_widths / 2
    std_lya_eff,_,_ = stats.binned_statistic(x=np.log10(cleaned_df_fd["vIF_Flex"].values),#binning lya_eff
                       values=cleaned_df_fd["lya_eff"].values,#the data on which the statistic is computed
                       statistic='std',
                      bins=bins)
    #offset = (1+flip_flop*0.005)
    offset =  (1+flip_flop*0)
    ax[1].errorbar(x=10**logvIF_centers*offset,y=mean_lya_eff,yerr=std_lya_eff, color=color_list[i],
        ls='none', marker='o', capsize=5, capthick=1, ecolor='black',
        markeredgecolor='black',markersize=6)
    flip_flop*=-1
    if alpha_list[i] == ui.bincenters_alpha[3]:
        #*m, b = np.polyfit(np.log10(cleaned_df_fd["vIF_Flex"]), cleaned_df_fd["lya_eff"], 4)
        *m, b = np.polyfit(logvIF_centers, mean_lya_eff, 6)
        y_fit = (m[0]*logvIF_centers**6
                + m[1]*logvIF_centers**5
                + m[2]*logvIF_centers**4
                + m[3]*logvIF_centers**3
                + m[4]*logvIF_centers**2
                + m[5]*logvIF_centers
                + b)
        print("first number is highest degree polynomial (5)")
        print([f"{j:.4e}" for j in m])
        print(f"{b:.4e}",b)
        print(logvIF_centers)
        #print("Coefficients:{} {}".format(m,b))
        #ax[1].plot(10**logvIF_centers,y_fit,color="k",lw=lw,ls='dashed')
        percent_accuracy = (mean_lya_eff-y_fit)/mean_lya_eff*100
        print(percent_accuracy)
        print(max(np.abs(percent_accuracy)))
#custom_lines.append(Line2D([0], [0], color="k", lw=lw,ls='dashed', marker=None))
#labels.append(r"Best fit")

ax[0].legend(custom_lines,labels,loc="upper left",ncol=2,frameon=False,fontsize=fontsize-2)


ax[0].set_xscale('log')
# ax[0].set_yscale('log')
#ax[0].set_xlim(1e2,1e5)
ax[0].set_xlim(2e2,9e4)
ax[0].set_ylim(1e4,4.5e4)
ax[0].set_xlabel(r"I-front speed, $v_\mathrm{IF}$ [km/s]",fontsize=fontsize)
ax[0].set_ylabel(r"I-front averaged temperature, $\left<T\right>_\mathrm{IF}$ [K]",fontsize=fontsize)


ax[1].set_xscale('log')
ax[1].set_xlim(2e2,9e4)
#ax[1].set_xlim(1e2,1e5)
ax[1].set_ylim(-0.05,1.0)
ax[1].set_xlabel(r"I-front speed, $v_\mathrm{IF}$ [km/s]",fontsize=fontsize)
ax[1].set_ylabel(r"Ly$\alpha$ production efficiency, $\zeta(\alpha,v_\mathrm{IF})$",fontsize=fontsize)
plt.tight_layout()

total_df_fd.to_csv("cleaned_df_fd.csv")
total_df_ud.to_csv("cleaned_df_ud.csv")
figname = "figures/temperature_lyaEff_plot_v3.png"
plt.savefig(figname,dpi=400)



