import matplotlib.pyplot as plt
import pandas as pd
from user_input import *
from cycler import cycler
import matplotlib.ticker as ticker
from scipy.interpolate import griddata


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



#input data -
df = pd.read_csv(
    "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/results/221202/otf.csv")
        #    "/Volumes/Extreme SSD/anson/test_expanse/master_1d_rt/plotting_routines/paper/final_data/otf2.csv")
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
#alpha_array = np.array(["{:<05}".format(i) for i in alpha_array])

#get the matrix
# Ndata = len(logvIF)
# frat_matrix = np.zeros((nbins_alpha,nbins_logvIF))
# for i in range(Ndata)
# for a in range(nbins_alpha):
#     for v in range(nbins_logvIF):
#         mask = (alpha_array==bincenters_alpha[a])&(v_IF>binedges_logvIF[v])&(v_IF<binedges_logvIF[v+1])
#         frat_matrix[a][v] = np.mean(flux_ratio[mask])

grid_V, grid_A = np.mgrid[min_logvIF:max_logvIF:10j, 0:max_alpha:10j] #min_alpha:max_alpha:10j]

bb = np.column_stack((v_IF,df["alpha"].values))
fratgrid = griddata(bb,flux_ratio,(grid_V,grid_A),method='linear')
extent=(10**(min_logvIF),10**(max_logvIF),
        min_alpha,max_alpha)
#print(fratgrid)
levels = np.arange(-0.1,1.1,0.1)
fig, ax = plt.subplots(1,figsize=(8,8))
ax.set_xscale("log")
contour = ax.contourf(10.**grid_V, grid_A, fratgrid, levels=levels,extent=extent,
    cmap="inferno",alpha=0.8)
contour2 = ax.contour(10.**grid_V, grid_A, fratgrid, levels=levels,colors='k')
ax.clabel(contour2, levels[::2], inline=True, fontsize=16, fmt='%.2f',colors='k')
ax.set_xlabel(r'I-front Velocity, v$_{\mathrm{IF}}$ [km/s]')
ax.set_ylabel(r'Incident Spectral Index, $\alpha_{\mathrm{IF}}$')

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(ticker.MultipleLocator(delta_bin_alpha))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.yaxis.set_ticks_position('both')

fig.savefig("figures/Frat_parameter_space.pdf", bbox_inches="tight")
plt.close()

#### Saving interpolation table
X = np.column_stack((bb.T[0],bb.T[1],flux_ratio))
np.savetxt(fname="../results/interp_tables/test1_Dec5_2022.txt",X=X)





