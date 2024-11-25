import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

columns = ["v [km/s]", "tau_HILya", "tau_HeI584", "tau_HeIILya", "nHI [cm-3]", "nHeII [cm-3]", "Delta_b",
           "T [K]","los pos [ckpc/h]","vpec [km/s]"]
h=6.626068e-27 #erg s
h_eV=4.135667e-15 #eV s
nu0=13.6/h_eV
alpha=1.5
pi=np.pi
chi = 0.08

Y = 0.24
m_H = 1.672622e-24
rho_crit_0 = 8.688e-30
redshift_bayu = 5.7
Omega_b = 0.048
rho_0 = rho_crit_0*Omega_b*(1+redshift_bayu)**3
nH_0 = (1 - Y)*rho_0/m_H

cm_per_kpc = 3.086e+21
cm_per_km  = 1e+5
sec_per_year = 3.154e+7
year_per_Myr = 1e+6

def get_ionizing_flux(L,r,R0): #photons/s/cm2
    return 1/h/nu0 * (alpha-1)/alpha * (4**-alpha-1)/(4**(1-alpha)-1) * L/4/pi/(R0+r)**2

def get_IF_speed(flux,n_H): #cms
    return flux/n_H/(1-chi)

def get_tsim(speed_cm_s,delta_r_cms): #input cms units, output Myr
    calibration_factor = 10/7
    return np.sum(1/speed_cm_s * delta_r_cms)/sec_per_year/year_per_Myr * calibration_factor

L_test_list = [1e+45,5e+45,4e+46]
R0_test = 1e4*cm_per_kpc
R_start_test = 50
R_test = 1000
h_const = 0.68

skewers = ["0001","0002","0003","0004","0005"]

tsim_list = []
skew_list = []
L_list = []
for sk in skewers:
    path_to_skewers = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/input_files/hydro_skewers/"
    file_path = path_to_skewers+"spec_xHeII1_007_ls_line{}.dat".format(sk)
    fileheader = open(file_path, 'r')
    redshift_mmc  = float(fileheader.readline()[12:20])
    df = pd.read_csv(file_path,skiprows=4,delim_whitespace=True,names=columns)
    R_temp = df["los pos [ckpc/h]"].values/(1+redshift_mmc)/h_const
    nH = df["Delta_b"]*nH_0
    mask_odrt = (R_temp>R_start_test)&(R_temp<R_test+R_start_test)
    r_cm = R_temp[mask_odrt]*cm_per_kpc
    delta_R_temp = np.mean(R_temp[1:]-R_temp[:-1])*cm_per_kpc

    for i,L_i in enumerate(L_test_list): 
        flux_test = get_ionizing_flux(L=L_i,r=r_cm,R0=R0_test)
        speed_test = get_IF_speed(flux = flux_test, n_H=nH[mask_odrt].values) #cm per sec
        #all_speeds.append(speed_test)
        skew_list.append(sk)
        L_list.append(L_i)
        tsim_list.append(get_tsim(speed_test,delta_R_temp))
        #print("Skewer: {}\nLuminosity: {:.1e}\ntsim: {:.1f}\n".format(skewer_number,L_i,tsim))

stack = np.column_stack((np.asarray(skew_list,int),np.array(L_list),np.array(tsim_list)))
print(stack)
#np.savetxt(fname="input_params/varySkew.txt",X=stack,fmt=('%04d','%.1e','%.2f'))
print(np.loadtxt(fname="input_params/varySkew.txt",dtype=str))







