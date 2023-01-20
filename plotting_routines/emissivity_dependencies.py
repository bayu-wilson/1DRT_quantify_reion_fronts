import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler
from constants import *
from user_input import *
from matplotlib.lines import Line2D
# from functions import coll_ex_rate_H1_acc,alphaA_H2

#initial figure parameters
fontsize=15
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
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
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                Line2D([0], [0], color=colors[3], lw=9, marker=None),
                Line2D([0], [0], color=colors[4], lw=9, marker=None),
                Line2D([0], [0], color=colors[5], lw=9, marker=None)]
labels = [r"$\alpha$={:.1f}".format(float(i)) for i in bincenters_alpha]
                # Line2D([0], [0], lw=0)]

#input data -
#df = pd.read_csv(
#    "/Volumes/Extreme SSD/anson/test_expanse/master_1d_rt/plotting_routines/paper/final_data/otf2.csv")
df=pd.read_csv("/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/221202/otf.csv")
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

alpha_array = np.asarray(df["alpha"],str)
alpha_array = np.array(["{:<05}".format(i) for i in alpha_array])
bincenters_alpha = np.array(["{:<05}".format(i) for i in bincenters_alpha])

### binning



fig, ax = plt.subplots(nrows=1,ncols=2)
fig.set_size_inches(w=10,h=4.5)

car = ['red', 'orange','gold','green','blue','purple']
flip_flop = 1
for a in range(2,nbins_alpha):
    mask = (alpha_array == bincenters_alpha[a])
    # print(alpha_array)
    C_em_means = np.zeros(nbins_logvIF)
    C_em_std = np.zeros(nbins_logvIF)

    for i in range(nbins_logvIF):
        mask_bin = (v_IF>binedges_logvIF[i])&(v_IF<binedges_logvIF[i+1])&(alpha_array == bincenters_alpha[a])
        C_em_means[i] = np.mean(C_em[mask_bin])
        C_em_std[i] = np.std(C_em[mask_bin])

    offset = (1+flip_flop*0.005)
    ax[0].scatter(v_IF[mask],T_center[mask],s=1)#,label="a={:.1f}".format(float(bincenters_alpha[a])))
    ax[1].scatter(v_IF[mask],C_em[mask],s=1,alpha=0.15)
    ax[1].errorbar(x=bincenters_logvIF*offset,y=C_em_means,yerr=C_em_std,
        ls='none', marker='o', capsize=5, capthick=1, ecolor='black',
        markeredgecolor='black',markersize=8)
    flip_flop*=-1

    # ax[2].scatter(np.log10(nH_avg[mask]),R_IF[mask]/kpc_to_cm,s=1)
    # ax[1][1].scatter(F_inc[mask],v_IF[mask],s=1)
    # ax11 = ax[1][1].twiny()
    # ax11.scatter(np.log10(nH_center[mask]),v_IF[mask],s=1,color=car[a])

ax[0].set_xlabel(r"IF speed, log$_{10}$ v$_{\mathrm{IF}}$ [km/s]")
ax[0].set_ylabel(r"IF temperature at center, T$_{\mathrm{c}}$[K]")

ax[1].set_xlabel(r"IF speed, log$_{10}$ v$_{\mathrm{IF}}$ [km/s]")
ax[1].set_ylabel(r"Emissivity Clumping Factor, C$_e$")

ax[0].legend(custom_lines,labels[2:],ncol=1,frameon=False)

plt.tight_layout()
fig.savefig("figures/emissivity_dependencies.pdf", bbox_inches="tight")
plt.close()
# plt.show()
