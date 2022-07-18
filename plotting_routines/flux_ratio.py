import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# gasprop_names = ["los pos[kpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI","xHeII","xHeIII","emLya"]
alpha="1.5"
sigma="0" #kpc
los = "11"

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
plt.scatter(np.log10(vIF),flux_ratio,s=10)
plt.show()
# plt.xlim(2.6,4.2)
# plt.ylim(0.15,0.5)
# otf_df["t"]
