import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from functions import coll_ex_rate_H1_acc,alphaA_H2
import warnings

fontsize=15
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.prop_cycle'] = cycler(color=['k', 'k'])
plt.rcParams['axes.prop_cycle'] = cycler(ls=['-', '--'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size

temperatures = np.logspace(3,5.1,100) #K
q_eff_lya = np.array([coll_ex_rate_H1_acc(i) for i in temperatures])
a_eff_lya = np.array([alphaA_H2(i) for i in temperatures])

fig, ax = plt.subplots(nrows=1,ncols=1)
fig.set_size_inches(w=6,h=5.5)
ax.set_xscale("log")
ax.set_yscale("log")
# ax.axvline(1.5e4)
# ax.axvline(1.5e4)
ax.set_xlim(1e3,1e5)
ax.set_ylim(1e-14,3e-8)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.plot(temperatures,q_eff_lya)
ax.plot(temperatures,a_eff_lya)
ax.annotate(text=r'Coll. Exc. (q$_{\mathrm{Ly}\alpha}^{\mathrm{eff}}$)',xycoords='figure fraction',xy=(0.64,0.44))
ax.annotate(text=r'Rec. ($\alpha_{\mathrm{Ly}\alpha}^{\mathrm{eff}}$)',xycoords='figure fraction',xy=(0.22,0.45))
ax.set_xlabel("T(K)")
ax.set_ylabel(r"cm$^3$s$^{-1}$")

plt.tight_layout()
fig.savefig("figures/rates.pdf")
plt.close()
plt.show()
