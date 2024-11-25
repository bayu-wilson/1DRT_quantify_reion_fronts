import numpy as np
import pandas as pd

#Constants
h=6.626068e-27 #erg s
h_eV=4.135667e-15 #eV s
nu0=13.6/h_eV
pi=np.pi

h_const =  0.68
Omega_b = 0.048
z_neutral_island = 5.7
rho_crit_z = 8.688e-30*(1+z_neutral_island)**3 #g/cm3
rho_crit_z_b = rho_crit_z*Omega_b
m_H =  1.672622e-24 #g
Y = 0.24
nH_mean_density = (1-Y)*rho_crit_z_b/m_H
chi = 0.08

cm_per_kpc = 3.086e+21
cm_per_km  = 1e+5
sec_per_yr = 3.154e7
yr_per_Myr = 1e6

sigma0 = 6.34e-18
R0_pcm = 1e+4*cm_per_kpc

eV_cgs = 1.60218e-12
Eth_H1 = 1.360e+1
k_B = 1.38065e-16
def lambda_func(T,Eth):
    K2eV = k_B/eV_cgs
    eV2K = 1.0/K2eV
    xHI  = Eth*eV2K/T
    x    = 2*xHI
    return x

#Case A recombination coefficient of HII from Hui & Gnedin 1996
def alphaA_H2(T):
    x = lambda_func(T,Eth_H1)
    return 1.269e-13*np.power(x,1.503)/np.power((1.0+np.power((x/0.522),0.470)),1.923) #cm3/s

def luminosity_for_equlibibrium(nH,alpha_spectral):
    T_gas = 1e+4
    beta_sigmaHI = 2.75
    A = nH*(1+chi)/2 * alphaA_H2(T_gas) * h*nu0/sigma0
    B = alpha_spectral*(beta_sigmaHI+alpha_spectral-1)*(np.power(4,1-alpha_spectral)-2)**2 * 4*pi*R0_pcm**2
    C = (alpha_spectral-1)**2 * (4**(-alpha_spectral)-1) * (4**(1-alpha_spectral-beta_sigmaHI)-1)
    return A * B/C

def alpha_dependence(alpha):
    #for use in the following luminosity equation, use inverse in the flux and vIF
    if 1-alpha == 0:
        return 8*np.log(2)/3
    elif alpha == 0:
        return 3/2/np.log(2)
    else:
        return alpha/(alpha-1) * (4**(1-alpha)-1)/(4**-alpha-1)
    
def get_luminosity(avg_speedIF,R_skewer,delta_r,R0,nH,alpha):
    neutral_gas_per_area = np.sum(nH)*(1+chi)*delta_r
    return avg_speedIF/R_skewer * neutral_gas_per_area * alpha_dependence(alpha) *h*nu0 * 4*pi*R0**2 



R_start_pkpc = 219
R_sim_pkpc = 1500
R_start_hckpc = 1000
R_sim_hckpc = R_sim_pkpc*(1+z_neutral_island)*h_const
N_r = R_sim_pkpc*4

path_to_skewers = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/input_files/hydro_skewers/"
column_names = ["v_kms","tau_HILya","tau_HeI584","tau_HeIILya", "nHI", "nHeII", "Delta_b", "T", "los_pos", "vpec"]

max_nH_list = []
skewer_number = np.arange(0,100,1)
for i in range(len(skewer_number)):
    skewer_filename = f"spec_xHeII1_007_ls_line{i:04d}.dat"
    df = pd.read_csv(path_to_skewers+skewer_filename,skiprows=4,delim_whitespace=True,names=column_names)
    R_ckpc_over_h = df["los_pos"]
    mask = (R_ckpc_over_h>R_start_hckpc) & (R_ckpc_over_h<R_sim_hckpc+R_start_hckpc)
    df = df[mask]
    R_ckpc_over_h = R_ckpc_over_h[mask]
    Delta_cmv = df["Delta_b"]
    R_pkpc = R_ckpc_over_h/(1+z_neutral_island)/h_const #coverting ckpc/h at z-7.1 to pkpc
    nH = (1. - Y)/m_H * Delta_cmv*rho_crit_z_b
    max_nH_list.append(np.max(nH))
R_pkpc -=R_start_pkpc
log_nH_list = np.log10(max_nH_list)
mask_SS = log_nH_list>-2

sorted_pairs = sorted(zip(log_nH_list[~mask_SS], skewer_number[~mask_SS]))
log_nH_sorted,skewers_sorted = np.array(sorted_pairs).T ###only every 3
nH_sorted = 10**log_nH_sorted


max_speedIF_range = [1e+7,1e+10]
N_samples = len(skewers_sorted)
speed_sampling = np.logspace(np.log10(max_speedIF_range[0]),
                             np.log10(max_speedIF_range[1]),N_samples)

R_skew_pcm = R_sim_pkpc*cm_per_kpc
delta_r_pcm = (R_pkpc.iloc[1]-R_pkpc.iloc[0])*cm_per_kpc

tsim_array = np.logspace(np.log10(4180),np.log10(30),N_samples)

#alpha bins
alpha_min,alpha_max,alpha_bins = 0,3,7
alpha_bincenters = np.linspace(alpha_min,alpha_max,alpha_bins)

stack = []
R0_array = np.ones(N_samples)*1.e+4
R_array = np.ones(N_samples)*R_sim_pkpc
N_r_array = np.ones(N_samples)*N_r
for alpha in alpha_bincenters:
    L_test = get_luminosity(avg_speedIF=speed_sampling,
                           R_skewer=R_skew_pcm,
                           delta_r=delta_r_pcm,
                           R0=R0_pcm,
                           nH=nH_sorted,
                           alpha=1.5)
    L_max = np.max([luminosity_for_equlibibrium(nH_sorted,1.5+1e-3)*3,L_test],axis=0)

    alpha_array = np.ones(N_samples)*alpha
    stack.append(np.column_stack((alpha_array,np.asarray(skewers_sorted,int),L_max,R_array, N_r_array ,R0_array,tsim_array)))

stack = np.concatenate(stack)
print(stack.shape)
np.savetxt(fname="input_params/fd_April_params.txt",X=stack,fmt=('%.1f','%04d','%.2e','%.1f','%d','%.1f','%.1f'))
print("alpha, skewer number, L[erg/s], r_skew[pkpc], N_r, R0[pkpc], t_sim[Myr]")
print(np.loadtxt(fname="input_params/fd_April_params.txt",dtype=str))

