import numpy as np
import pandas as pd

#Constants
h=6.626068e-27 #erg s
h_eV=4.135667e-15 #eV s
nu0=13.6/h_eV
pi=np.pi

Omega_b = 0.048
z = 5.7
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rho_crit_z_b = rho_crit_z*Omega_b
m_H =  1.672622e-24 #g
Y = 0.24
nH_mean_density = (1-Y)*rho_crit_z_b/m_H
chi = 0.08

cm_per_kpc = 3.086e+21
cm_per_km  = 1e+5
sec_per_yr = 3.154e7
yr_per_Myr = 1e6

R0 = 1e4*cm_per_kpc

little_h = 0.68

#binning for vIF
#logvIF_min,logvIF_max,nbins = 1.5,5,4
#logvIF_bincenters = np.linspace(logvIF_min,logvIF_max,nbins)
#delta_logvIF = (logvIF_max-logvIF_min)/(nbins-1)
#logvIF_binedges = np.linspace(logvIF_min-delta_logvIF/2,logvIF_max+delta_logvIF/2,nbins+1)
#vIF_binedges_cms = 10**logvIF_binedges * cm_per_km #km/s *cm/km = cm/s
#vIF_bincenters_cms = 10**logvIF_bincenters * cm_per_km

#alpha bins
alpha_min,alpha_max,alpha_bins = 0,3,7
alpha_bincenters = np.linspace(alpha_min,alpha_max,alpha_bins)


#functions
def alpha_dependence(alpha):
    #for use in the following luminosity equation, use inverse in the flux and vIF
    if 1-alpha == 0:
        return 8*np.log(2)/3
    elif alpha == 0:
        return 3/2/np.log(2)
    else:
        return alpha/(alpha-1) * (4**(1-alpha)-1)/(4**-alpha-1)

def get_luminosity(max_speedIF,R0,alpha):
    #luminosity necessary to get the max speed at mean density
    return max_speedIF * h*nu0 * alpha_dependence(alpha) * 4*pi*R0**2 * nH_mean_density*(1+chi)

def get_ionizing_flux(L,R0,alpha): #photons/s/cm2
    return 1/h/nu0 *1/alpha_dependence(alpha)* L/4/pi/R0**2

def get_IF_speed(flux,nH): #cms
    #here nH is the density in each cell in a simulation
    return flux/nH/(1-chi)

def get_tsim(delta_r,vIF):
    #t =  integral of dr/vIF
    return np.sum(delta_r/vIF/sec_per_yr/yr_per_Myr) #Myr

cleaned_skewer_list=['0044', '0071', '0076', '0070', '0034', '0096', '0089', '0079',
       '0025', '0005', '0042', '0063', '0099', '0057', '0087', '0016',
       '0066', '0008', '0084', '0029', '0090', '0060', '0088', '0021',
       '0072', '0009', '0010', '0050', '0097', '0067', '0045', '0083',
       '0054', '0091', '0058', '0014', '0092', '0052', '0075', '0053',
       '0017', '0007', '0056', '0035', '0031', '0033', '0040', '0064',
       '0012', '0061', '0073', '0082', '0032', '0013', '0036', '0047',
       '0081', '0018', '0085', '0093', '0038', '0026', '0039', '0002',
       '0055', '0041', '0037', '0095', '0098', '0048', '0011', '0077',
       '0019', '0074', '0015', '0043', '0068', '0024', '0006', '0080',
       '0078', '0086', '0003', '0001', '0065', '0049', '0094', '0059',
       '0027', '0023', '0051'] #length 91, ordered by increasing max(nH) over-density

Lum_list = ['7.08e+44', '8.53e+44', '1.14e+45', '1.33e+45', '1.38e+45',
       '1.40e+45', '1.65e+45', '1.67e+45', '1.80e+45', '1.96e+45',
       '2.01e+45', '2.07e+45', '2.16e+45', '2.20e+45', '2.28e+45',
       '2.28e+45', '2.28e+45', '2.29e+45', '2.32e+45', '2.38e+45',
       '2.40e+45', '2.40e+45', '2.44e+45', '2.46e+45', '2.47e+45',
       '2.55e+45', '2.71e+45', '2.80e+45', '2.81e+45', '2.84e+45',
       '2.98e+45', '3.02e+45', '3.05e+45', '3.16e+45', '3.17e+45',
       '3.18e+45', '3.29e+45', '3.37e+45', '3.45e+45', '3.45e+45',
       '3.55e+45', '3.63e+45', '3.65e+45', '3.65e+45', '3.96e+45',
       '3.97e+45', '4.35e+45', '4.39e+45', '4.40e+45', '4.47e+45',
       '4.48e+45', '4.74e+45', '4.83e+45', '5.01e+45', '5.07e+45',
       '5.19e+45', '5.34e+45', '5.36e+45', '5.45e+45', '5.48e+45',
       '5.59e+45', '5.60e+45', '5.61e+45', '5.77e+45', '6.07e+45',
       '6.39e+45', '6.72e+45', '7.07e+45', '7.45e+45', '7.84e+45',
       '8.25e+45', '8.68e+45', '9.14e+45', '9.62e+45', '1.01e+46',
       '1.07e+46', '1.12e+46', '1.18e+46', '1.24e+46', '1.31e+46',
       '1.38e+46', '1.45e+46', '1.52e+46', '1.60e+46', '1.69e+46',
       '1.78e+46', '1.87e+46', '1.97e+46', '2.07e+46', '2.18e+46',
       '2.29e+46']

path_to_skewers = "/expanse/lustre/projects/uot171/bwils033/master_1d_rt/input_files/hydro_skewers/"
column_names = ["v_kms","tau_HILya","tau_HeI584","tau_HeIILya", "nHI", "nHeII", "Delta_b", "T", "los_pos", "vpec"]
#for i in range(len(cleaned_skewer_list)):
#max_speedIF_range = [1e+7,3e+9] #cm/s
#max_speedIF_range = [8e+8,5e+9] #cm/s
max_speedIF_range = [1e+7,5e+9] #cm/s
np.random.seed(0)
len_skewer_list = len(cleaned_skewer_list)
multiplier = 1
N_samples = len_skewer_list*multiplier
#N_samples = 10
#random_speed_sampling = 10**np.random.uniform(np.log10(max_speedIF_range[0]),np.log10(max_speedIF_range[1]),N_samples)
speed_sampling = np.logspace(np.log10(max_speedIF_range[0]),
                             np.log10(max_speedIF_range[1]),N_samples)


stack = []
R_skew = 1500 #2500 #pkpc #assume R_start = 219
N_r = R_skew*4
for alpha in alpha_bincenters:
    for i in range(N_samples):
        skewer_filename = f"spec_xHeII1_007_ls_line{cleaned_skewer_list[i]}.dat"
        df = pd.read_csv(path_to_skewers+skewer_filename,skiprows=4,delim_whitespace=True,names=column_names)
        R_ckpc_over_h = df["los_pos"]
        R_pkpc = R_ckpc_over_h/(1+5.7)/little_h #coverting ckpc/h at z-7.1 to pkpc
        Delta_baryon = df["Delta_b"]
        nH = (1. - Y)/m_H * Delta_baryon*rho_crit_z_b 
        #Lum = get_luminosity(max_speedIF=speed_sampling[i],R0=R0,alpha=alpha)
        Lum = float(Lum_list[i])
        F_lyC = get_ionizing_flux(Lum,R0,alpha=alpha)
        vIF_array = get_IF_speed(flux=F_lyC,nH=nH)
        delta_r = (R_pkpc[2]-R_pkpc[1])*cm_per_kpc
        t_sim = get_tsim(delta_r,vIF_array)
        
        stack.append([alpha,int(cleaned_skewer_list[i]),Lum,R_skew, N_r ,R0/cm_per_kpc,t_sim])
        #print(f"t_sim={t_sim:.3f} and vIF={np.median(vIF_array)/1e5:.3e}")

np.savetxt(fname="input_params/fd_planeParallel_params.txt",X=stack,fmt=('%.1f','%04d','%.2e','%.1f','%d','%.1f','%.1f'))
print("alpha, skewer number, L[erg/s], r_skew[pkpc], N_r, R0[pkpc], t_sim[Myr]")
print(np.loadtxt(fname="input_params/fd_planeParallel_params.txt",dtype=str))



