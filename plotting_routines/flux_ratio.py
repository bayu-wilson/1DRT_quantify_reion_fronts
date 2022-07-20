import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# gasprop_names = ["los pos[kpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI","xHeII","xHeIII","emLya"]
# alpha="1.5"
# sigma="0" #kpc
# los = "11"

otf_names = ["step","t","dt","I_lya","Finc","nH_value","locIF","vIF_fd","vIF_fm","T_IF"]
N_output = 250
R_sim = 300
width_IF = 20 #an approximation
kpc_to_cm = 3.086e21
los_pos_initial = 5317
h = 6.626068e-27
c = 2.998e10
lambda_lya_cm = 1.21567e-5
pi=np.pi
fontsize = 15

output_matrix = np.zeros((10,N_output))
n=0
while n<245:
    with open("../output_files/gasprops/n{}_gasprops.txt".format(n+1),'r') as f:
        line = f.readline()
    otf_data = np.asarray(line.split(' \t'),float)
    for index in range(len(otf_data)):
        output_matrix[index][n]=otf_data[index]
    n+=1
otf_df = pd.DataFrame(output_matrix.T,columns=otf_names)
locIF = (otf_df["locIF"]/kpc_to_cm-1e4)
mask = (otf_df["t"]>10)&(locIF<(R_sim-width_IF))
F_lya_otf = otf_df["I_lya"]/(h*c/lambda_lya_cm)*4*pi
F_inc = otf_df["Finc"]
flux_ratio = (F_lya_otf/F_inc)[mask]
vIF = otf_df["vIF_fm"][mask]/1e5

#FITTING
a, b = np.polyfit(np.log10(vIF), flux_ratio, 1)
y_fit = a*np.log10(vIF)+b

# print(a,b)

fig,ax = plt.subplots(1,1)
fig.set_size_inches(6,6)
ax.scatter(np.log10(vIF),flux_ratio,s=10)
ax.plot(np.log10(vIF),y_fit,color="gray",ls="dotted",lw=3,alpha=0.8,label="{:.1f}log10(vIF)+{:.1f}".format(a,b))
ax.set_xlabel(r"log$_{10}$v$_{\mathrm{IF}}$ [km s$^{-1}$]",fontsize=fontsize)
ax.set_ylabel(r"$F^{\mathrm{emit}}_{\mathrm{Ly}\alpha}/F^{\mathrm{inc}}_{\mathrm{ion}}$",fontsize=fontsize)
ax.legend()
ax.set_xlim(2.6,3.8)
ax.set_ylim(0.2,0.47)
# plt.savefig("figures/flux_ratio_alpha1.5.pdf",bbox_inches='tight')

plt.show()
# plt.xlim(2.6,4.2)
# plt.ylim(0.15,0.5)
# otf_df["t"]
