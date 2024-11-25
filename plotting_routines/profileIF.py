import matplotlib.pyplot as plt
from constants import *
import pandas as pd
from cycler import cycler
import matplotlib.gridspec as gridspec


#initial figure parameters
fontsize=15
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.prop_cycle'] = cycler(color=['k','r','b'])
plt.rcParams['axes.prop_cycle'] = cycler(ls=['-'])
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size

#input data - skewer
#gasprop_path = "../output_files/gasprops/sk0001_a=1.500_Lum=1.0e+45/n260_gasprops.txt"
gasprop_path = "../output_files/gasprops/fd_pp_a=1.5_sk0002_L=1.63e+46/n175_gasprops.txt"
#sk0004_a=1.500/n70_gasprops.txt" 
#"final_data/n80_gasprops.txt"
column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI","xHeII","GammaHI","q_lya"]
df_gasprops = pd.read_csv(gasprop_path,delim_whitespace=True,skiprows=1,names=column_names)
Delta = df_gasprops["rho"]/rhoBaryons_0
nH_skewer = (1. - Y)*df_gasprops["rho"]/m_H
nHI_skewer = df_gasprops["xHI"]*nH_skewer
r_skewer = df_gasprops["los pos[pkpc]"]-1e4
T_skewer = df_gasprops["T[K]"]
xHI_skewer = df_gasprops["xHI"]
xHeI_skewer = df_gasprops["xHeI"]
jlya_skewer = df_gasprops["ne"]*nHI_skewer*df_gasprops["q_lya"]*h*c/lambda_lya_cm/4/pi
#frontIF = np.interp(0.99,df_gasprops["xHI"],r_skewer)
#backIF = np.interp(1e-3,df_gasprops["xHI"],r_skewer)

#input data - OTF
with open(gasprop_path,'r') as f:
    line = f.readline()
    line_array = np.asarray(line.split(' \t'),float)
    locIF=line_array[5]/kpc_to_cm-1e4
    #widthIF=line_array[10]/kpc_to_cm
    #backIF = frontIF-widthIF
    # print(widthIF)
    # print(line_array[6]/kpc_to_cm-1e4)


fig,ax = plt.subplots(nrows=4,ncols=1,figsize=(5,8))
ax[0].plot(r_skewer-locIF,Delta)
ax[0].axvline(0,ls='dashed')
ax[0].set_xlim(-45,45)
ax[0].set_ylabel(r"x$_j$")
ax[0].set_ylabel(r"$\Delta$")
ax[0].tick_params(labelbottom=False,top=True)
ax[0].set_yscale("log")
ax[0].annotate("", xy=(15, 10), xytext=(0, 10),
            arrowprops=dict(arrowstyle="->"))
#ax[0].set_yticks(np.arange(0.5,2.5,0.5),minor=False)

ax[1].plot(r_skewer-locIF,xHI_skewer)
ax[1].plot(r_skewer-locIF,xHeI_skewer,ls='dotted')
ax[1].axvline(0,ls='dashed')
ax[1].set_xlim(-45,45)
ax[1].set_ylim(-0.1,1.1)
ax[1].set_ylabel(r"x$_j$")
ax[1].tick_params(labelbottom=False,top=True)

ax[2].plot(r_skewer-locIF,T_skewer)
ax[2].axvline(0,ls='dashed')
ax[2].set_yticks(np.arange(0,35000,10000),minor=False)
ax[2].set_xlim(-45,45)
ax[2].set_ylabel("T[K]")
ax[2].tick_params(labelbottom=False,top=True)

ax[3].plot(r_skewer-locIF,np.log10(jlya_skewer))
ax[3].axvline(0,ls='dashed')
ax[3].set_yticks(np.arange(-33,-31+1,1),minor=False)
ax[3].set_ylim(-33.5,-30.5)
ax[3].set_xlim(-45,45)
ax[3].set_ylabel(r"log$_{10}$(j$_{\mathrm{Ly}\alpha}$)")
ax[3].set_xlabel(r"R$-$R$_{\mathrm{IF}}$ [pkpc]")
ax[3].tick_params(labelbottom=True,top=True)

fig.savefig("figures/profileIF_april24.png", bbox_inches="tight",dpi=250)
plt.close()

"""
# fig, ax = plt.subplots(nrows=2,ncols=1)
fig = plt.figure()
fig.set_size_inches(w=5,h=8)
plt.subplots_adjust(hspace=0.9)
# ax = plt.subplot2grid((4, 2), (0, 0))
# ax0 = plt.subplot2grid((4, 2), (0, 0), colspan=2)
# ax1 = plt.subplot2grid((4, 2), (1, 0), colspan=2)
# ax2 = plt.subplot2grid((4, 2), (2, 0), colspan=2)
# ax3 =  plt.subplot2grid((4, 2), (3, 0), colspan=2)

gs = gridspec.GridSpec(4,2)
# gs.update(left=0.05, right=0.48, wspace=0.05)
# gs.update(hspace=0.00)
ax0 = plt.subplot(gs[0,:])
ax1 = plt.subplot(gs[1,:])
ax2 = plt.subplot(gs[2,:])
ax3 = plt.subplot(gs[3,:])

# ax.set_xscale("log")
# ax.set_yscale("log")
ax0.set_xlim(0,1.5e3)
# ax.set_ylim(1e-14,3e-8)

#ax0.set_yticks(np.arange(0,40+20,20),minor=False)
#ax0.set_yticks(np.arange(10,50+20,20),minor=True)
#ax0.set_yticks(np.arange(-1,2,20),minor=False)

#ax0.set_yticks(np.arange(0,10+5,5),minor=False)
#ax0.set_yticks(np.arange(2.5,12.5+5,5),minor=True)
# ax0.xaxis.set_label_position('top')

ax0.xaxis.set_ticks_position('both')
ax0.yaxis.set_ticks_position('both')
ax0.set_yticks(np.arange(0,3,1),minor=False)

#ax0.plot(r_skewer,np.log10(Delta))
ax0.plot(r_skewer,np.log10(Delta))
ax0.axvline(locIF,ls='dashed')
ax0.axvspan(locIF-45,locIF+45,alpha=0.3)
#ax0.axvspan(backIF,frontIF,alpha=0.3)
# ax.plot(temperatures,a_eff_lya)
# ax.annotate(text=r'Coll. Exc. (q$_{\mathrm{Ly}\alpha}^{\mathrm{eff}}$)',xycoords='figure fraction',xy=(0.64,0.44))
# ax.annotate(text=r'Rec. ($\alpha_{\mathrm{Ly}\alpha}^{\mathrm{eff}}$)',xycoords='figure fraction',xy=(0.22,0.45))
#ax0.set_xlabel(r"R$_{\mathrm{sim}}-$R$_0$ [pkpc]")
ax0.set_xlabel(r"R [pkpc]")
ax0.set_ylabel(r"log$_\mathrm{10} (\Delta)$")
# ax0.xaxis.tick_top()

#ax0.annotate("", xy=(locIF+150, 30), xytext=(locIF, 30),
#            arrowprops=dict(arrowstyle="->"))
ax0.annotate("", xy=(locIF+150, 1), xytext=(locIF, 1),
            arrowprops=dict(arrowstyle="->"))


ax1.plot(r_skewer-locIF,xHI_skewer)
ax1.plot(r_skewer-locIF,xHeI_skewer,ls='dotted')
ax1.axvline(0,ls='dashed')
#ax1.set_yticks(np.logspace(-3,-1,3),minor=True)
#ax1.axvspan(backIF-locIF,frontIF-locIF,alpha=0.2)
ax1.set_xlim(-45,45)
ax1.set_ylim(-0.1,1.1)
ax1.set_ylabel(r"x$_j$")

#ax1.set_yscale("log")
#ax1.set_ylim(9e-5,1.5)
#ax1.set_ylabel(r"log$_{10}$(x$_i$)")

# ax2 = ax1.twinx()
ax2.plot(r_skewer-locIF,T_skewer)
#ax2.axvspan(backIF-locIF,frontIF-locIF,alpha=0.2)
ax2.axvline(0,ls='dashed')
ax2.set_yticks(np.arange(0,35000,10000),minor=False)
#ax2.set_yticks(np.arange(5,45+20,20),minor=True)
ax2.set_xlim(-45,45)
ax2.set_ylabel("T[K]")
# ax2.subplots_adjust(left=0.15)
# ax2.yaxis.labelpad = -20
# ax2.tick_params(axis='y', labelcolor="red")
# ax2.yaxis.labelpad = 20
# ax2.set_ylim(3,4.5)

ax3.plot(r_skewer-locIF,np.log10(jlya_skewer))
#ax3.axvspan(backIF-locIF,frontIF-locIF,alpha=0.2)
ax3.axvline(0,ls='dashed')
ax3.set_yticks(np.arange(-33,-31+1,1),minor=False)
ax3.set_ylim(-33.5,-30.5)
ax3.set_xlim(-45,45)
ax3.set_ylabel(r"log$_{10}$(j$_{\mathrm{Ly}\alpha}$)")
ax3.set_xlabel(r"R$-$R$_{\mathrm{IF}}$ [pkpc]")
#ax3.set_xlabel(r"(R$_{\mathrm{sim}}$-R$_0$)-R$_{\mathrm{IF}}$ [pkpc]")

pos1 = ax1.get_position()
pos2 = ax2.get_position()
pos3 = ax3.get_position()
points1 = pos1.get_points()
points2 = pos2.get_points()
points3 = pos3.get_points()
#00  left  01 bottom
#01  right 11 top
total_height = points1[1][1]-points3[0][1]
# print(total_height,points3[0][1])
loc12 = points3[0][1]+total_height*2/3
loc23 = points3[0][1]+total_height/3
# loc12 = (points1[0][1]+points2[1][1])/2
# loc23 = (points2[0][1]+points3[1][1])/2
#
points1[0][1] = loc12
points2[1][1] = loc12
points2[0][1] = loc23
points3[1][1] = loc23

#
# points1[0][1]=0.5#points2[1][1]
# points2[0][0]=0.5
# points2[1][1]=0.25
# points3[1][1]=0.25#points2[0][1]
pos1.set_points(points1)
pos2.set_points(points2)
pos3.set_points(points3)
ax1.set_position(pos1)
ax2.set_position(pos2)
ax3.set_position(pos3)

# plt.tight_layout()

# plt.show()
fig.savefig("figures/profileIF_Oct.png", bbox_inches="tight",dpi=250)
# plt.show()
plt.close()
# plt.show()


# print(line[6])
# print(line[6])

# #step, t [Myr],dt [Myr], I_lya , Finc, clump, IF loc, vIF fd, vIF fm,
# #treion, temp eff, IF width, I test, nH_avg, nH2 avg
# # print(line)]


# plt.plot(r_skewer,Delta)
# print(df_gasprops)
# plt.show()
"""
