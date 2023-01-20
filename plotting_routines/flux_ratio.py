import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler
from constants import *
from user_input import *
from matplotlib.lines import Line2D
# from functions import coll_ex_rate_H1_acc,alphaA_H2
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
#warnings.filterwarnings("ignore", message="RuntimeWarning: invalid value encountered in true_divide")
#warnings.filterwarnings("ignore", message="RuntimeWarning: Degrees of freedom <= 0 for slice")
#warnings.filterwarnings("ignore", message="RuntimeWarning: Mean of empty slice.")


#initial figure parameters
fontsize=15
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
colors = ['red', 'orange','gold','green','blue','purple','cyan','pink']
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
# plt.rcParams['axes.prop_cycle'] = cycler(ls=['-', '--'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                Line2D([0], [0],lw=0,marker=None),
                Line2D([0], [0], color=colors[3], lw=9, marker=None),
                Line2D([0], [0], color=colors[4], lw=9, marker=None),
                Line2D([0], [0], color=colors[5], lw=9, marker=None),
                Line2D([0], [0], color='k', lw=3, ls='--')]
                #Line2D([0], [0], color=colors[6], lw=9, marker=None),
                #Line2D([0], [0], color=colors[7], lw=9, marker=None)]
# Line2D([0], [0], color=colors[5], lw=9, marker=None)]
labels = [r"$\alpha$={:.1f}".format(float(i)) for i in bincenters_alpha]
                # Line2D([0], [0], lw=0)]
labels.insert(-3, '')

#input data -
df = pd.read_csv("/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/221202/otf.csv")
#    "/Volumes/Extreme SSD/anson/test_expanse/master_1d_rt/plo/final_data/otf2.csv")
mask = (df["t"]>10)&(df["vIF_fm"]>1e1)&(df["vIF_fm"]<1e5)
df = df[mask]
F_inc = df["F_inc"].values
F_lya = df["F_lya"].values
T_center = df["T_center"].values
v_IF = np.log10(df["vIF_fm"].values)
R_IF = df["width_IF"].values
nH_avg = df["nH_avg"].values
nH_center = df["nH_center"].values
C_em = df["C_em"].values
flux_ratio = F_lya/F_inc

#alpha_array = np.asarray(df["alpha"],str)
#alpha_array = np.array(["{:.3f}".format(i) for i in alpha_array])
alpha_array = df["alpha"]
print(set(alpha_array))
print(bincenters_alpha)


fig, ax = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(w=6,h=6)

#car = ['red', 'orange','gold','green','blue','purple']
flip_flop = 1
for a in range(2,nbins_alpha):
    mask = (alpha_array == bincenters_alpha[a])
    print(bincenters_alpha[a],np.sum(mask))
    flux_ratio_means = np.zeros(nbins_logvIF)
    flux_ratio_std = np.zeros(nbins_logvIF)
    for i in range(nbins_logvIF):
        mask_bin = (v_IF>binedges_logvIF[i])&(v_IF<binedges_logvIF[i+1])&(alpha_array == bincenters_alpha[a])
        flux_ratio_means[i] = np.mean(flux_ratio[mask_bin])
        flux_ratio_std[i] = np.std(flux_ratio[mask_bin])

    if bincenters_alpha[a]==1.5: #"1.500":
        m, b = np.polyfit(np.log10(v_IF[mask]), flux_ratio[mask], 1)
        y_fit = m*np.log10(bincenters_logvIF)+b
        print(m,b)
        ax.plot(bincenters_logvIF,y_fit,color='black',ls='dashed')

    offset = (1+flip_flop*0.005)
    ax.scatter(v_IF[mask],flux_ratio[mask],s=1,label="a={}".format(bincenters_alpha[a]),alpha=0.25)
    ax.errorbar(x=bincenters_logvIF*offset,y=flux_ratio_means,yerr=flux_ratio_std,
        ls='none', marker='o', capsize=5, capthick=1, ecolor='black',
        markeredgecolor='black',markersize=8)
    flip_flop*=-1

    # ax[2].scatter(np.log10(nH_avg[mask]),R_IF[mask]/kpc_to_cm,s=1)
    # ax[1][1].scatter(F_inc[mask],v_IF[mask],s=1)
    # ax11 = ax[1][1].twiny()
    # ax11.scatter(np.log10(nH_center[mask]),v_IF[mask],s=1,color=car[a])

# ax.annotate(text=r'math',xycoords='figure fraction',xy=(0.4,0.8))

ax.set_xlabel(r"IF speed, log$_{10}$ v$_{\mathrm{IF}}$ [km/s]")
ax.set_ylabel(r"Flux ratio, $\zeta(\alpha, v_{\mathrm{IF}})$")
ax.set_ylim(-0.03,1.25)
ax.set_xlim(2.5,4.65)

custom_lines.append(Line2D([0], [0], color='k', lw=3,ls='dashed', marker=None))
labels.append(r"Best fit")
#
# ax[1].set_xlabel(r"IF speed, log$_{10}$ v$_{\mathrm{IF}}$ [km/s]")
# ax[1].set_ylabel(r"Emissivity Clumping Factor, C$_e$")

ax.legend(custom_lines,labels[2:],ncol=2,frameon=False)

# plt.tight_layout()
fig.savefig("figures/flux_ratio.pdf", bbox_inches="tight")
plt.close()
# plt.show()
