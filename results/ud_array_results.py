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

#alpha_list = np.arange(amin,amax+0.125,0.125)
#alpha_list = np.arange(amin,amax+0.125,0.5)#0.125)
#print(alpha_list,flush=True)
#alpha_list = np.loadtxt("../parameters_input/input_params/spectral_indices.txt")
#Nv = len(alpha_list)
#logvIF_bincenters =  np.linspace(vlogmin,vlogmax,Nv)
#grid_V, grid_A = np.mgrid[vlogmin:vlogmax:complex(N_grid), amin:amax:complex(N_grid)]
#delta_logvIF = logvIF_bincenters[1]-logvIF_bincenters[0]
#logvIF_binedges = np.linspace(vlogmin-delta_logvIF/2,
#                            vlogmax+delta_logvIF/2,Nv+1)
#skewers = ["sk0000","sk0001","sk0002","sk0003","sk0004","sk0005","sk0006","sk0007","sk0008","sk0009"]
#skewers = ["sk0000","sk0001","sk0002","sk0003"]
#dotN = np.loadtxt("../input_files/dotN_Rsim.txt",dtype=str)[:,0]
#dotN = np.loadtxt("../parameters_input/input_params/flucRho.txt",dtype=str)[:,1]
#dotN_list = ["5.e+54","1.e+55","2.e+55", "4.e+55"]#,"5.e+54"]
#dotN_list = ['1.e+57','1.06e+57','1.24e+57','1.2e+57']#["1.e+57","1.2e+57","1.6e+57","2.e+57"]
#Lum_list = ["2.27e+44", "1.05e+45", "4.89e+45", "2.27e+46"]
#Lum_list = ["8.27e+43","3.84e+44","1.78e+45","8.27e+45"]
#Lum_list = ["4.54e+43","2.11e+44","9.77e+44","4.54e+45"]
#Lum_list = ["9.07e+43","4.21e+44","1.95e+45","9.07e+45"]
Lum_list = ["2.96e+44","9.36e+44","2.96e+45","9.36e+45"]

alpha_list = ["0.500", "1.000", "1.500", "2.000", "2.500"]

#otf_names = ["step","t","dt","I_lya","Finc","C_factor","locIF","vIF_fd","vIF_fm","T_IF","density_var","IF_width","nH_avg"]
otf_names = ["alpha","step", "t", "dt", "F_lya", "F_inc", "locIF","vIF_fd","vIF_Flex","T_avg","T_center","width_IF","nH_center","C_em"]
total_names = ["Lum", "alpha_i"]+otf_names

#column_names = ["los pos[pkpc]","rho","ne","dene_dt","n_tot","T[K]","xHI","xHeI","xHeII","gamma_HI","q_lya"]

print("begin for loop")
otf_total_df = pd.DataFrame(columns=total_names)
len_otf = len(otf_names)
for L in Lum_list: #looping through the skewers
    for i in range(len(alpha_list)): #looping through spectral index
        #print(sk,i)
        spectral_index = alpha_list[i]
        n=1
        dir_path = "../output_files/gasprops/ud_a={:.3f}_Lum={}/".format(float(spectral_index),L) #I shouldn't need the .3f here
        print(dir_path,flush=True)
        gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
        spectrum_path_otf = dir_path+"n{}_spectrum.txt".format(n)
        otf_matrix = np.zeros((N_output,len(otf_names)))
        #alpha_list_i=[]
        while (os.path.exists(gasprop_path_otf))&(n<=N_output):
            nu,I_nu = np.loadtxt(spectrum_path_otf).T
            sigma_nu = sigmaHI_over_sigmaHI0_numeric(nu)
            avg_sigma_log = np.trapz(sigma_nu*I_nu/nu,x=nu)/np.trapz(I_nu/nu,x=nu)
            my_alpha = fsolve(func=find_root_sigma, x0=1.0,args=avg_sigma_log)[0]
            #alpha_list_i.append(my_alpha)

            with open(gasprop_path_otf,'r') as f:
                line = f.readline()
            try:
                line_array = np.asarray(line.split(' \t'),float)
                #line_list = line.split(' \t')
                line_array = np.insert(line_array,0,my_alpha)
                otf_matrix[n-1]=line_array #np.array([sk[2:],spectral_index]+line_list)
                #if line_array[2]<10: #skip first 10 Myr
                #    n_minTime = n
            except:
                otf_matrix[n-1]=np.ones(len(otf_names))*np.nan
            n+=1
            gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
            spectrum_path_otf = dir_path+"n{}_spectrum.txt".format(n)
        otf_df = pd.DataFrame(otf_matrix,columns=otf_names)[2:]
        otf_df["Lum"] = L
        otf_df["alpha_i"] = spectral_index
        #otf_df["alpha"] = spectral_index
        #print(len(alpha_list_i[:-1]))
        #print(np.shape(otf_matrix))
        #otf_df["alpha"] = alpha_list_i[:-1]
        otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm#-1e4
        otf_df["vIF_fd"]/=1e5
        otf_df["vIF_Flex"]/=1e5
        mask_otf = np.asarray(np.ones(len(otf_df)),bool)
            #mask_otf = np.ones_like(otf_df["I_lya"]) #((otf_df["t"]>10)&(locIF[n]<(R_sim-width_IF)))#&
            #            (log_vIF_fm>=2.0)&(log_vIF_fm<=5.0))
        
        otf_total_df = pd.concat([otf_total_df,otf_df],ignore_index=True)
        
save_dir = "230428/"
otf_total_df.to_csv(save_dir+"otf40.csv",index=False)

