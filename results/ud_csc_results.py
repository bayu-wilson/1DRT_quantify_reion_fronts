import numpy as np
import os.path
import pandas as pd
#import matplotlib.pyplot as plt
from scipy.optimize import fsolve

N_output = 500

kpc_to_cm = 3.086e21
Y = 0.24
m_H = 1.672622e-24
h = 6.626068e-27
c = 2.998e10
lambda_lya_cm = 1.21567e-5
pi=np.pi
h_eV = 4.135667e-15
nu_HI = 13.6/h_eV
beta=2.75

#functions
def sigmaHI_over_sigmaHI0_numeric(nu):
    return (nu/nu_HI)**(-beta)
def sigma_over_sigmaHI0_analytic(alpha):
    num = alpha * (1 - 4**(-alpha-beta))
    den = (alpha + beta)*(1 - 4**(-alpha))
    return num/den
def find_root_sigma(alpha,*C):
    return sigma_over_sigmaHI0_analytic(alpha)-C

#amin = -1.000
#amax = 2.500 #2.505

np.seterr(divide='ignore', invalid='ignore')

#skewers = ["sk0000","sk0001","sk0002","sk0003"]
#Lum_list = ["4.74e+44","1.50e+45","4.74e+45","1.50e+46"]
#skewers = ["0001"]#, "0002", "0003", "0005"]
Lum_list = ["2.06e+46","1.26e+47","2.52e+47"]
Rsim_list=["340.9","945.5","1355.2"]
#skewer = "0002"

#csc_ud_L1.51e+47_N_r9456
#csc_ud_Rsim1355.2_L2.52e+47_N_r271

alpha_list = ["1.500"]

#N_r_list = [94560,  9456,  1891,  1261,   946,   630,   473,   189,    95]
#N_r_matrix = np.array([[ 34090,   3409,    682,    455,    341,    227,    170,     68],
#       [ 94560,   9456,   1891,   1261,    946,    630,    473,    189],
#       [135520,  13552,   2710,   1807,   1355,    903,    678,    271]])
N_r_matrix = np.array([[ 3409,   682,   455,   341,   227,   170],
       [ 9456,  1891,  1261,   946,   630,   473],
       [13552,  2710,  1807,  1355,   903,   678]])
num_cell_sizes = N_r_matrix.shape[1]
#[20000, 13333,  8000,  6666,  4000,  2666,  2000,  1000,   400, 200]

#otf_names = ["step","t","dt","I_lya","Finc","C_factor","locIF","vIF_fd","vIF_fm","T_IF","density_var","IF_width","nH_avg"]
otf_names = ["step", "t", "dt", "F_lya", "F_inc", "locIF","vIF_fd","vIF_Flex","T_avg","T_center","width_IF","nH_center","gamma_loc"]
total_names = ["Lum", "alpha_i", "N_r"]+otf_names

#column_names = ["los pos[pkpc]","rho","ne","dene_dt","n_tot","T[K]","xHI","xHeI","xHeII","gamma_HI","q_lya"]

print("begin for loop")
otf_total_df = pd.DataFrame(columns=total_names)
len_otf = len(otf_names)

for i in range(num_cell_sizes):
    for j in range(3): #slow, med, fast I-fronts
        N_r = N_r_matrix[j][i]
        Rsim = Rsim_list[j]
        Lum = Lum_list[j]
        spectral_index = alpha_list[0]
        n=1
        dir_path = f"../output_files/gasprops/csc_ud_alpha{spectral_index}_Rsim{Rsim}_L{Lum}_N_r{N_r}/"
        print(dir_path,flush=True)
        gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
        spectrum_path_otf = dir_path+"n{}_spectrum.txt".format(n)
        otf_matrix = np.zeros((N_output,len(otf_names)))
        while (os.path.exists(gasprop_path_otf))&(n<=N_output):
            #nu,I_nu = np.loadtxt(spectrum_path_otf).T
            #sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
            #avg_sigma_log = np.trapz(sigma_nu*I_nu/nu,x=nu)/np.trapz(I_nu/nu,x=nu)
            #my_alpha = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma_log)[0]

            with open(gasprop_path_otf,'r') as f:
                line = f.readline()
            try:
                line_array = np.asarray(line.split(' \t'),float)
                #line_array = np.insert(line_array,0,my_alpha)
                otf_matrix[n-1]=line_array #np.array([sk[2:],spectral_index]+line_list)
            except:
                otf_matrix[n-1]=np.ones(len(otf_names))*np.nan
            n+=1
            gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
            spectrum_path_otf = dir_path+"n{}_spectrum.txt".format(n)
        otf_df = pd.DataFrame(otf_matrix,columns=otf_names)[2:]
        #L = Lum_list[0]
        otf_df["Lum"] = Lum
        otf_df["N_r"] = N_r
        otf_df["alpha_i"] = spectral_index
        otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm#-1e4
        otf_df["vIF_fd"]/=1e5
        otf_df["vIF_Flex"]/=1e5
        mask_otf = np.asarray(np.ones(len(otf_df)),bool)
        otf_total_df = pd.concat([otf_total_df,otf_df],ignore_index=True)

save_dir = "final_results/"
otf_total_df.to_csv(save_dir+f"csc_ud_otf_alpha{spectral_index}.csv",index=False)
