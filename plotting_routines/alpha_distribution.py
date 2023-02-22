import numpy as np
import matplotlib.pyplot as plt

#initial figure parameters
fontsize=14
tick_dir='in'
lw=1.5
major_tick_size = 5
minor_tick_size = 3
colors = ['red', 'orange','gold','green','blue','purple','cyan','pink']
plt.rcParams['font.family'] = "serif"
plt.rcParams['lines.linewidth'] = lw
plt.rcParams['font.size'] = fontsize
plt.rcParams['xtick.direction'] = tick_dir
plt.rcParams['ytick.direction'] = tick_dir
plt.rcParams['xtick.major.size'] = major_tick_size
plt.rcParams['xtick.minor.size'] = minor_tick_size
plt.rcParams['ytick.major.size'] = major_tick_size
plt.rcParams['ytick.minor.size'] = minor_tick_size

#Load alpha matrix 
X = np.loadtxt("../results/matrices_alpha/alpha_matrix.txt")

#Define x-axis data (distance traveled through ionized IGM) 
los_pos_bin_edges = np.linspace(0,2500,11)/1e3 #pMpc
delta_r = np.median(los_pos_bin_edges[1:]-los_pos_bin_edges[:-1])
los_pos_bin_centers = los_pos_bin_edges[:-1]+delta_r/2

#Define what skewers to cycle through
skewer_array = np.arange(0,41,1)
skewer_array = np.delete(skewer_array,[0,2,3,4,16,18,20,21,22,23,27,28,30,36]) #these have high density spikes that produce outliers

#Plotting each skewer's alpha distribution
fig,ax = plt.subplots(1,figsize=(6,5))
for i in skewer_array:
    ax.plot(los_pos_bin_centers,X[i],alpha=0.75,color='gray')

mean_per_bin = np.nanmean(X[skewer_array],axis=0)
std_per_bin = np.nanstd(X[skewer_array],axis=0)
ax.plot(los_pos_bin_centers,mean_per_bin,lw=5,color='red')
ax.fill_between(los_pos_bin_centers, y1=mean_per_bin-std_per_bin, y2=mean_per_bin+std_per_bin,
                 alpha=0.4,color='red')

ax.set_ylabel(r"incident spectral index, $\alpha$")
ax.set_xlabel("Distance traveled through \ninhomogeneous, ionized medium [pMpc]")
ax.set_xlim(0,2.25)
ax.set_ylim(0.95,1.5)

plt.savefig("figures/alpha_distribution.pdf",bbox_inches='tight')
plt.close()

