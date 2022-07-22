import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

matrix = np.zeros((11,11))
# matrix2= np.zeros((11,11))
alpha_list = np.linspace(0.5,2.5,11)
vIF_bins = np.logspace(2.7,3.6,11)
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
for i in range(11):
    spectral_index = alpha_list[i]
    gasprop_path = "../output_files/gasprops/a={:.1f}/".format(spectral_index)
#     spec_path = "output_files/incident_spectra/a={:.1f}/".format(spectral_index)
    output_matrix = np.zeros((10,N_output))
    n=1
    while n<249:
        with open(gasprop_path+"n{}_gasprops.txt".format(n),'r') as f:
            line = f.readline()
        otf_data = np.asarray(line.split(' \t'),float)
        for index in range(len(otf_data)):
            output_matrix[index][n]=otf_data[index]
        n+=1
    otf_df = pd.DataFrame(output_matrix.T,columns=otf_names)
    locIF = (otf_df["locIF"]/kpc_to_cm-1e4)
    vIF_fm = (otf_df["vIF_fm"]/1e5)
    mask = ((otf_df["t"]>10)&(locIF<(R_sim-width_IF))&(vIF_fm>=10**2.0)&(vIF_fm<=10**4.0))
    F_lya_otf = otf_df["I_lya"][mask]/(h*c/lambda_lya_cm)*4*pi
    F_inc = otf_df["Finc"][mask]
    nH = otf_df["nH_value"][mask]
    flux_ratio = (F_lya_otf/F_inc)
    vIF_fm = vIF_fm[mask]

    a, b = np.polyfit(np.log10(vIF_fm), flux_ratio[mask], 1)
    for j in range(11):
        matrix[i][j] = a*np.log10(vIF_bins[j])+b


fontsize=20
levels = np.arange(0,1,0.1)
# levels = np.arange(2000,12000,1000)
fig, ax = plt.subplots(1)
fig.set_size_inches(8,8)
pos = ax.imshow(matrix,origin='lower')

ax.locator_params(axis="x", nbins=11)
ax.locator_params(axis="y", nbins=11)
xlabels = [None]+["{:.1f}".format(i) for i in np.log10(vIF_bins)]
ylabels = [None]+["{:.1f}".format(i) for i in alpha_list]
cbar = fig.colorbar(pos, ax=ax, orientation="horizontal", pad=0.15)
cbar.ax.get_yaxis().labelpad = 30
cbar.ax.set_xlabel(r"$F^{\mathrm{emit}}_{\mathrm{Ly}\alpha}/F^{\mathrm{inc}}_{\mathrm{ion}}$", rotation=0,fontsize=fontsize)
ax.set_xticklabels(xlabels)
ax.set_yticklabels(ylabels)
ax.set_xlabel(r"log$_{10}$v$_{\mathrm{IF}}$ [km s$^{-1}$]",fontsize=fontsize)
ax.set_ylabel(r"Incident spectral index, $\alpha_{\mathrm{IF}}$",fontsize=fontsize)

countour = ax.contour(matrix, levels=levels,origin='lower', colors='red')

fmt = {}
levels_str = ["{:.1f}".format(i) for i in countour.levels]
for l, s in zip(countour.levels, levels_str):
    fmt[l] = s
ax.clabel(countour, countour.levels[::1], inline=True, fmt=fmt, fontsize=15)

# plt.savefig("figures/parameter_space.pdf",bbox_inches='tight')
plt.show()
