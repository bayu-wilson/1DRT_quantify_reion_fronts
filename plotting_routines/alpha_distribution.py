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
#working_date = "Dec2023"

working_date = "march24"
X = np.loadtxt(f"../results/matrices_alpha/alpha_matrix_{working_date}.txt")
Nbins = 20

#Define x-axis data (distance traveled through ionized IGM) 
los_pos_bin_edges = np.linspace(0,2200,Nbins+1)/1e3 #pMpc
delta_r = np.median(los_pos_bin_edges[1:]-los_pos_bin_edges[:-1])
los_pos_bin_centers = los_pos_bin_edges[:-1]+delta_r/2
los_pos_bin_centers = los_pos_bin_centers*(1+5.7)*0.68


#Define what skewers to cycle through
skewer_array = np.loadtxt("/expanse/lustre/projects/uot171/bwils033/master_1d_rt/input_files/skewer_list_no_SS.txt")
skewer_array = np.array([i for i in skewer_array if i not in [94.,95.] ])
# skewer_i not in [94,95])
#skewer_array = np.arange(0,41,1)
#skewer_array = np.delete(skewer_array,[0,2,3,4,16,18,20,21,22,23,27,28,30,36]) #these have high density spikes that produce outliers

#print(los_pos_bin_centers)
#print(X[0])
#Plotting each skewer's alpha distribution
fig,ax = plt.subplots(1,figsize=(7,6))
for i,sk_number in enumerate(skewer_array):
    ax.plot(los_pos_bin_centers,X[i],alpha=0.75,color='gray')

mean_per_bin = np.nanmean(X,axis=0) #removed skewer_array when it was used for indices
std_per_bin = np.nanstd(X,axis=0)
ax.plot(los_pos_bin_centers,mean_per_bin,lw=5,color='red')
ax.fill_between(los_pos_bin_centers, y1=mean_per_bin-std_per_bin/2, y2=mean_per_bin+std_per_bin/2,
                 alpha=0.4,color='red')

ax.set_ylabel(r"incident spectral index, $\alpha$")
#ax.set_xlabel("Distance traveled through \ninhomogeneous, ionized medium [pMpc]")
ax.set_xlabel("Distance traveled through \ninhomogeneous, ionized medium [cMpc/h]")
ax.set_xlim(0.75,9.1)#8.5)
#ax.set_xlim(0,2.25)
ax.set_ylim(0.4,1.5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

yticks = np.arange(0.5,1.5+0.2,0.2)
ax.set_yticks(ticks=yticks)
xticks = np.arange(1,10,1)
ax.set_xticks(ticks=xticks)
#        ticks = np.linspace(0,bins_list[i]-1,4)
#        labels = np.asarray(np.linspace(0,24,4),int)
#        ax[i,j].scatter(x=[bins_list[i]/2],y=[bins_list[i]/2],s=25,marker="+",color="red",alpha=0.5)
#        ax[i,j].set_xticks(ticks=ticks)
#        ax[i,j].set_xticklabels(labels)

plt.tight_layout()
plt.savefig(f"figures/alpha_distribution_{working_date}.png",dpi=200) #bbox_inches='tight',dpi=200)
plt.close()

