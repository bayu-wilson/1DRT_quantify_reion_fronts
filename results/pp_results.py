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

np.seterr(divide='ignore', invalid='ignore')

otf_names = ["step", "t", "dt", "F_lya", "F_inc", "locIF","vIF_fd","vIF_Flex","T_avg","T_center","width_IF","nH_center","trash"]
total_names = ["Lum", "alpha_i"]+otf_names
otf_total_df = pd.DataFrame(columns=total_names)
len_otf = len(otf_names)

#array([['0.0', '6.95e+46', '762.8', '1525', '1.0e+04', '596.3'],
#       ['0.0', '2.08e+47', '1109.4', '2218', '1.0e+04', '73.0'],

ud_or_fd = "fd"
if ud_or_fd == "ud":
    input_params = np.loadtxt(fname="../parameters_input/input_params/planeParallel_params.txt",dtype=str)
elif ud_or_fd == "fd":
    #input_params = np.loadtxt(fname="../parameters_input/input_params/fd_planeParallel_params.txt",dtype=str)
    input_params = np.loadtxt(fname="../parameters_input/input_params/fd_April_params.txt",dtype=str)

#L_list = ["5.50e+46","1.65e+47","2.75e+47"]
#L_list = ["8.25e+46","4.95e+47","9.17e+47"]
#L_list = ["5.50e+46"]
#alpha_list = ["1.5"]
#alpha_list = ["1.5","1.5","1.5"]

for i in range(0,len(input_params),3):#range(78):#range(len(input_params)):#len(L_list)):
    #for i in range(0,int(91*4),3): #range(8):
    spectral_index = input_params[i][0] #alpha_list[i]
    #print(spectral_index)
    n=1
    if ud_or_fd == "ud":
        L = input_params[i][1]#L_list[i]
        dir_path = f"../output_files/gasprops/ud_pp_a={spectral_index}_L={L}/"
    elif ud_or_fd == "fd":
        skewer = input_params[i][1]#L_list[i]
        L = input_params[i][2]#L_list[i]
        # print("alpha, skewer number, L[erg/s], r_skew[pkpc], N_r, R0[pkpc], t_sim[Myr]")
        dir_path = f"../output_files/gasprops/fd_pp_a={spectral_index}_sk{skewer}_L={L}_v2/"
    print(dir_path,flush=True)
    gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
    spectrum_path_otf = dir_path+"n{}_spectrum.txt".format(n)
    otf_matrix = np.zeros((N_output,len(otf_names)))
    #alpha_list_i=[]
    while (os.path.exists(gasprop_path_otf))&(n<=N_output):
        #nu,I_nu = np.loadtxt(spectrum_path_otf).T
        #sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
        #avg_sigma_log = np.trapz(sigma_nu*I_nu/nu,x=nu)/np.trapz(I_nu/nu,x=nu)
        #my_alpha = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma_log)[0]
        #alpha_list_i.append(my_alpha)

        with open(gasprop_path_otf,'r') as f:
            line = f.readline()
        try:
            line_array = np.asarray(line.split(' \t'),float)
            #line_list = line.split(' \t')
            #line_array = np.insert(line_array,0,my_alpha)
            otf_matrix[n-1]=line_array #np.array([sk[2:],spectral_index]+line_list)
            #if line_array[2]<10: #skip first 10 Myr
            #    n_minTime = n
        except:
            otf_matrix[n-1]=np.ones(len(otf_names))*np.nan
        n+=1
        gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
        spectrum_path_otf = dir_path+"n{}_spectrum.txt".format(n)
    if n<50:
        print("Something is wrong with the path", dir_path)
    otf_df = pd.DataFrame(otf_matrix,columns=otf_names)[2:]
    otf_df["Lum"] = L
    if ud_or_fd == "fd":
        otf_df["skewer"]=skewer 
    otf_df["alpha_i"] = spectral_index
    otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm#-1e4
    otf_df["vIF_fd"]/=1e5
    otf_df["vIF_Flex"]/=1e5
    mask_otf = np.asarray(np.ones(len(otf_df)),bool)
    otf_total_df = pd.concat([otf_total_df,otf_df],ignore_index=True)

#save_dir = "230428/"
save_dir = "final_results/"
otf_total_df.to_csv(save_dir+f"{ud_or_fd}_pp_otf_v3.csv",index=False)

print(pd.read_csv(save_dir+f"{ud_or_fd}_pp_otf_v3.csv"))
print(save_dir+f"{ud_or_fd}_pp_otf_v3.csv")
