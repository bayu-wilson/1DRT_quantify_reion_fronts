import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

"""
 The purpose of this program is to add density fluctuations to the 1D RT sightlines using
 sherwood simulations given to us by Matt. Matt's simulations are organized differently and
 are at a different redshift so we need to scale them correctly. First, I rescale matt's HI
 and position into units of proper distance and then choose a chunk with similar length and
 resolution as 1dRT. Then I interpolate Matt's nHI densities to the 1dRT distance coordinate
 grid.

 In order to do this I created an initial data file "gas_test_hydro_000myr.txt" which is just
 the initial grid. From the initial grid, I added matt's density fluctuations and then saved
 it to "rho_fluc_line0000.txt"
"""

make_plot = True
col_names = ["velocity [km/s]", "tau_HILya", "tau_HeI584", "tau_HeIILya", "nHI [cm^-3]",
             "nHeII [cm^-3]", "Delta_b","T [K]", "los pos [kpc/h]", "vpec [km/sec]"]
line_number = "0001"
section = "1"
path_to_dir = "/Volumes/Extreme SSD/anson/master_1d_rt/input_files/"
path_to_dir2 = "/Volumes/Extreme\ SSD/anson/master_1d_rt/input_files/"
df_matt_sim = pd.read_csv(path_to_dir+"spec_xHeII1_003_ls_line{0}.dat".format(
                line_number),skiprows=4,delim_whitespace=True,names=col_names)
z = 7.1
# temperature = 5e2
h_const = 0.68
factor = (1+z)*h_const # dividing by factor converts comoving to proper.
R = 300
m_H = 1.672622e-24 # Mass of hydrogen atom (g)
Y = 0.24
pi=np.pi
kpc_to_cm = 3.08568e21
G_grav = 6.67259e-8
H0 = h_const*100*1e5/(1000*kpc_to_cm)
rho_crit_z0 = (3*H0**2/(8*pi*G_grav))
rho_crit_z = rho_crit_z0*(1+z)**3
Omega_b = 0.048
pos_col_ckpc_overh = df_matt_sim["los pos [kpc/h]"] #original position column
pos_col_pkpc = pos_col_ckpc_overh/factor #converting to pkpc
# nHI_col_ckpc_cubed = df_matt_sim["nHI [cm^-3]"] #HI number Density comoving
rho_pkpc = df_matt_sim["Delta_b"]*rho_crit_z*Omega_b
T_pkpc = df_matt_sim["T [K]"]
if section == "1":
    initial_loc_pkpc = int(max(pos_col_pkpc)*.515) #avoiding a significant over density but essentially middle of the sightline
else:
    initial_loc_pkpc = int(max(pos_col_pkpc)*.715)
mask_chunk = (pos_col_pkpc>initial_loc_pkpc) & (
              pos_col_pkpc<initial_loc_pkpc+R) # choosing a R size pkpc chunk


pos_col_pkpc_chunk = pos_col_pkpc[mask_chunk]-initial_loc_pkpc+5317 # so position coordinate goes from 0-300 pkpc
# pos_col_pkpc_chunk = pos_col_pkpc[mask_chunk]-initial_loc_pkpc+1e4 # so position coordinate goes from 0-300 pkpc
rho_pkpc_chunk = rho_pkpc[mask_chunk]
T_pkpc_chunk = T_pkpc[mask_chunk]
# nHI_col = df_matt_sim["nHI [cm^-3]"][mask_chunk]

# df_bayu_1dRT = pd.read_csv("plotting_routines/iLya_vs_vIF/rho_flux_line000_UNIFORM_DENSITY.txt",
#                       skiprows=0, delim_whitespace=True)
df_bayu_1dRT = pd.read_csv(path_to_dir+"gas_test_hydro_000myr.txt",
                      skiprows=2, delim_whitespace=True)
rho_matt_interpolated = np.interp(x=df_bayu_1dRT["radius"],xp=pos_col_pkpc_chunk,fp=rho_pkpc_chunk)
T_matt_interpolated = np.interp(x=df_bayu_1dRT["radius"],xp=pos_col_pkpc_chunk,fp=T_pkpc_chunk)
# print(rho_matt_interpolated)

### SMOOTHING
x_vals = df_bayu_1dRT["radius"]
y_vals = rho_matt_interpolated
sigma = 3
smoothed_vals = np.zeros(y_vals.shape)
for i,x_position in enumerate(x_vals):
    kernel = np.exp(-(x_vals - x_position) ** 2 / (2 * sigma ** 2))
    # kernel = kernel / sum(kernel)
    # print( -(x_vals.values - x_position) ** 2 / (2 * sigma ** 2))
    # break;

    smoothed_vals[i] = sum(y_vals * kernel)/sum(kernel)
rho_matt_interpolated = smoothed_vals
# print(x_vals[:10])
# print(y_vals[:10])
# print(rho_matt_interpolated[:10])
### testing
kernel = np.exp(-(x_vals - x_vals[50]) ** 2 / (2 * sigma ** 2))
print(kernel)
###


if make_plot:
    delta_radius = df_bayu_1dRT["radius"]-df_bayu_1dRT["radius"][0]
    # plt.plot(delta_radius,df_bayu_1dRT["rho"],label="uniform density")
    plt.plot(delta_radius,rho_matt_interpolated,label="fluctuating density")
    plt.xlabel("Position [pkpc]")
    plt.ylabel(r"rho [g cm$^{-3}$]")
    plt.legend()
    # plt.xlim(0,20)
    # fig.savefig("/Users/bayuwilson/Desktop/figures_IF/sample_skewers/0000_section1.pdf")
    plt.savefig("/Users/bayuwilson/Desktop/figures_IF/sample_skewers_constT/skewer{}_section{}.pdf".format(
        line_number,section),dpi=250)
    # plt.show()


df_bayu_1dRT["rho"] = rho_matt_interpolated
df_bayu_1dRT["nH1"] = df_bayu_1dRT["rho"]*(1-Y)/m_H
df_bayu_1dRT["T"] = T_matt_interpolated
# df_bayu_1dRT["P"] = 1.38065e-16*df_bayu_1dRT["T"]*nHI_matt_interpolated
# df_bayu_1dRT["rho"] = nHI_matt_interpolated*m_H/(1-Y)
df_bayu_1dRT.to_csv(path_to_dir+"rho_fluc_line{0}.txt".format(line_number),sep="\t",index=False)
# !echo 'task goes here\nhi' | cat - gas_test_hydro_000myr.txt > temp && mv temp gas_test_hydro_000myr.txt
os.system("echo 'Gas Data\nTime: 0.000000e+00' | cat - "+path_to_dir2+"rho_fluc_line{0}.txt > temp".format(line_number)
+" && mv temp "+path_to_dir2+"rho_fluc_line{0}.txt".format(line_number))
