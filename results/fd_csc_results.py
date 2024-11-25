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

#csc_fd_sk0013_L2.52e+46_N_r8333
#skewers = ["sk0000","sk0001","sk0002","sk0003"]
#Lum_list = ["4.74e+44","1.50e+45","4.74e+45","1.50e+46"]
#skewers = ["0013"]#, "0002", "0003", "0005"]
#skewers = ["0027","0049","0067","0078", "0092"]
skewers = ["0044","0060","0056","0002","0065"]
Lum_list = ["7.07e+44","2.39e+45","5.76e+45","2.89e+46","1.45e+47"]
#Lum_list = ["2.33e+46","3.06e+44","2.19e+45","1.11e+45","7.22e+46"]
#Lum_list = ["6.01e+46","7.23e+46","2.37e+46","4.2e+46","5.48e+46"]
#["2.52e+46"]
#skewer = "0002"

alpha_list = ["1.500"]
#N_r_list = [20000, 13333,  8000,  6666,  4000,  2666,  2000,  1000,   400, 200]
#N_r_list = [8333,7895,7500,6667,6000,5455,5000,3000,1500,750,300,150]
#N_r_list = [8333,7500,6000,5000,3000,1500,300]
#N_r_list = [10000,5000,4000,3333,2000,1000,500]
N_r_list = [8333,7500,6000,5000,3000,1500,300]
#N_r_list=('10000' '5000' '4000' '3333' '2000' '1000', '500')
#len_N_r="${#N_r_list[@]}"
#skewer_list=('0027' '0049' '0067' '0078' '0092')
#t_max_list=('114' '75' '229' '158' '93')
#Lum_list=('6.01e+46' '7.23e+46' '2.37e+46' '4.2e+46' '5.48e+46')


#otf_names = ["step","t","dt","I_lya","Finc","C_factor","locIF","vIF_fd","vIF_fm","T_IF","density_var","IF_width","nH_avg"]
otf_names = ["step", "t", "dt", "F_lya", "F_inc", "locIF","vIF_fd","vIF_Flex","T_avg","T_center","width_IF","nH_center","gamma_loc"]
total_names = ["skewer","Lum", "alpha_i", "N_r"]+otf_names

#column_names = ["los pos[pkpc]","rho","ne","dene_dt","n_tot","T[K]","xHI","xHeI","xHeII","gamma_HI","q_lya"]

print("begin for loop")
otf_total_df = pd.DataFrame(columns=total_names)
len_otf = len(otf_names)

for i,N_r in enumerate(N_r_list):
    spectral_index = alpha_list[0]
    for j in range(len(skewers)):
        skewer = skewers[j]
        Lum = Lum_list[j]
        n=1
        dir_path = f"../output_files/gasprops/csc_fd_sk{skewer}_L{Lum}_N_r{N_r}_v2/"
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
        otf_df["skewer"] = skewer
        otf_df["Lum"] = Lum
        otf_df["N_r"] = N_r
        otf_df["alpha_i"] = spectral_index
        otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm#-1e4
        otf_df["vIF_fd"]/=1e5
        otf_df["vIF_Flex"]/=1e5
        mask_otf = np.asarray(np.ones(len(otf_df)),bool)
        otf_total_df = pd.concat([otf_total_df,otf_df],ignore_index=True)

save_dir = "final_results/"
otf_total_df.to_csv(save_dir+"fd_csc_otf_alpha1.5_v2.csv",index=False)
