import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import options as opt

df = pd.read_csv(opt.output_dir+"otf_bayu_100myr.txt",skiprows=0,delim_whitespace=True)[4:]

fig, ax = plt.subplots(4,1)
fig.set_size_inches(6,12)
########################################
############### FIGURE 0 ###############
########################################
N_df = len(df)
H1_IF_x = df["H1_IF_x"]/3.068e21-1e4
time = df["time"]
ax[0].plot(time,time/max(time),label="time")
ax[0].plot(time,df["timestep"]/max(df["timestep"]),label="timestep")
ax[0].plot(time,df["inc_photon_flux"]/max(df["inc_photon_flux"]),label="F incident")
# plt.plot(time,df[3:N_test]["HI_IF_v_flxMthd"]/max(df[3:N_test]["HI_IF_v_flxMthd"]),label="vIF_fm")
ax[0].plot(time,H1_IF_x/max(H1_IF_x),label="location IF")
ax[0].plot(time,df["H1_IF_v"]/max(df["H1_IF_v"]),label="vIF finite diff")
ax[0].legend()
ax[0].set_xlabel("Time (Myr)")
ax[0].set_ylabel("normalized quantity")
ax[0].set_ylim(-0.1,1.1)
# plt.savefig("/Users/bayuwilson/Desktop/normalized5.pdf",dpi=200)

########################################
############### FIGURE 1 ###############
########################################

mask1 = (df["time"]>5)&(df["time"]<15)
mask2 = (df["time"]>15)&(df["time"]<25)
mask3 = (df["time"]>25)&(df["time"]<35)
mask4 = (df["time"]>35)&(df["time"]<45)
mask5 = (df["time"]>45)&(df["time"]<55)
mask6 = (df["time"]>55)&(df["time"]<65)
mask7 = (df["time"]>65)&(df["time"]<75)
mask8 = (df["time"]>75)&(df["time"]<100)
mask_list = [mask1, mask2, mask3, mask4, mask5, mask6, mask7, mask8]
label_list = ["5-15 Myr","15-25","25-35","35-45","45-55","55-65","65-75","75-100 "]
colors = ["red","orange","gold","green","blue","indigo","violet","black"]
# ax.set_xscale('log')
# ax.set_ylim(4.9e-9,5.6e-9)

for i,mask in enumerate(mask_list):
    df_masked = df[mask]
    ax[1].scatter(df_masked["nH_boundary"],df_masked["I_lya"],s=10,color=colors[i],label=label_list[i]) #starts from 0
ax[1].legend(reversed(ax[1].legend().legendHandles),reversed(label_list));
ax[1].set_xlabel("nH, hydrogen number density")
ax[1].set_ylabel("Intensity")

########################################
############### FIGURE 2 ###############
########################################

for i,mask in enumerate(mask_list):
    df_masked = df[mask]
    bin_height,bin_centers = np.histogram(df_masked["I_lya"],bins=6)
    ax[2].plot(bin_centers[:-1],bin_height,color=colors[i],label=label_list[i])
# plt.xlim(4.9e-9,5.6e-9)
ax[2].set_ylim(-1,35)
ax[2].legend(reversed(ax[2].legend().legendHandles),reversed(label_list),loc="upper left");
ax[2].set_xlabel("Intensity")
ax[2].set_ylabel("Number per bin")

########################################
############### FIGURE 3 ###############
########################################
ax[3].scatter(df["H1_IF_x"]/3.068e21-df["H1_IF_x"].values[0]/3.068e21,df["width_IF"]/3.068e21,s=5)
ax[3].set_xlabel("IF loc along LOS [kpc]")
ax[3].set_ylabel("IF width kpc")

plt.tight_layout()
fig.savefig("/Users/bayuwilson/Desktop/figures_IF/sample_skewers_constT/diagnostic.pdf")
plt.close()
# plt.show()
