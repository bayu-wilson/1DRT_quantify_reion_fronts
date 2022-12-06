import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#print("hi")

N_output = 250
#R_sim = 1000
#width_IF = 20 #an approximation
kpc_to_cm = 3.086e21
# los_pos_initial = 5317
Y = 0.24
m_H = 1.672622e-24
h = 6.626068e-27
c = 2.998e10
lambda_lya_cm = 1.21567e-5
pi=np.pi

#vlogmin = 2.5
#vlogmax = 4.5
#amin = 0.500
#amax = 2.500 #2.505
#N_grid = 17
spectral_index= 1.500

np.seterr(divide='ignore', invalid='ignore')

#alpha_list = np.arange(amin,amax+0.125,0.125)
#alpha_list = np.arange(amin,amax+0.125,0.5)#0.125)
#Nv = len(alpha_list)
#logvIF_bincenters =  np.linspace(vlogmin,vlogmax,Nv)
#grid_V, grid_A = np.mgrid[vlogmin:vlogmax:complex(N_grid), amin:amax:complex(N_grid)]
#delta_logvIF = logvIF_bincenters[1]-logvIF_bincenters[0]
#logvIF_binedges = np.linspace(vlogmin-delta_logvIF/2,
#                            vlogmax+delta_logvIF/2,Nv+1)
#skewers = ["sk0000","sk0001","sk0002","sk0003","sk0004"]
#dotN = np.loadtxt("../input_files/dotN_Rsim.txt",dtype=str)[:,0]
sk_array = np.loadtxt("../parameters_input/input_params/flucRho.txt",dtype=str)[:,0]
otf_names = ["step","t","dt","I_lya","Finc","C_em","locIF","vIF_fd","vIF_fm","T_reion","temp_eff","three_avg","width_IF","nH2_avg"]
column_names = ["los pos[pkpc]","rho","ne","dne_dt","n_tot","T[K]","xHI","xHeI","xHeII","gamma_HI","q_lya"]

vIF_total = np.array([])
vIF_fd_total = np.array([])
F_lya_total = np.array([])
F_inc_total = np.array([])
alpha_total = np.array([])
c_em_total = np.array([])
#nH2_total =  np.array([])
#three_avg_total= np.array([])
width_total= np.array([])
#T_reion_total=np.array([])
#temp_eff_total = np.array([])

#density_matrix = np.zeros((Nv,Nv))
#density2_matrix = np.zeros((Nv,Nv))
#num_matrix = np.zeros((Nv,Nv))
print("running python code")
#for sk in skewers: #looping through the skewers
for sk in sk_array:    
    #for i in range(len(alpha_list)): #looping through spectral index
    #spectral_index = alpha_list[i]
    gasprop_path = "gasprops/sk{}_a={:<05}/".format(sk,spectral_index)
    output_matrix = np.zeros((len(otf_names),N_output))
    for n in range(1,N_output,1): #looping through each output file per skewer/alpha
        try:
            df_gasprops = pd.read_csv(gasprop_path+"n{}_gasprops.txt".format(n),delim_whitespace=True,skiprows=1,
                                   names=column_names)
        except FileNotFoundError:
            print("File Not Found: "+gasprop_path+"n{}_gasprops.txt".format(n))
            break
        with open(gasprop_path+"n{}_gasprops.txt".format(n),'r') as f:
            line = f.readline()
        otf_data = np.asarray(line.split(' \t'),float)
        for index in range(len(otf_data)):
            output_matrix[index][n]=otf_data[index]
        los_pos_initial = df_gasprops["los pos[pkpc]"][0]
        df_gasprops["los pos[pkpc]"] = df_gasprops["los pos[pkpc]"]-los_pos_initial
        los_pos = df_gasprops["los pos[pkpc]"]
        nH_arr = (1. - Y)*df_gasprops["rho"]/m_H
        #backIF = np.interp(0.05,df_gasprops["xHI"],los_pos)
        #frontIF = np.interp(0.75,df_gasprops["xHI"],los_pos)
        #mask_gas = (los_pos>backIF)&(los_pos<frontIF)

        otf_df = pd.DataFrame(output_matrix.T,columns=otf_names)
        locIF = (otf_df["locIF"]/kpc_to_cm-1e4)
        log_vIF_fd = np.log10(otf_df["vIF_fd"]/1e5)
        log_vIF_fm = np.log10(otf_df["vIF_fm"]/1e5)
        mask_otf = np.asarray(np.ones(len(otf_df)),bool)
        #mask_otf = np.ones_like(otf_df["I_lya"]) #((otf_df["t"]>10)&(locIF[n]<(R_sim-width_IF)))#&
        #            (log_vIF_fm>=2.0)&(log_vIF_fm<=5.0))
        F_lya_otf = otf_df["I_lya"][mask_otf]/(h*c/lambda_lya_cm)*4*pi
        F_inc = otf_df["Finc"][mask_otf]
        c_em_factor = otf_df["C_em"][mask_otf]
        #nH2_avg = otf_df["nH2_avg"][mask_otf]
        widthIF = otf_df["width_IF"][mask_otf]
        #three_avg = otf_df["three_avg"][mask_otf]
        #temp_eff = otf_df["temp_eff"][mask_otf]
        #T_reion = otf_df["T_reion"][mask_otf]
        # flux_ratio = (F_lya_otf/F_inc)
        alpha_total = np.append(alpha_total,np.ones_like(log_vIF_fm[mask_otf])*float(spectral_index))
        vIF_total = np.append(vIF_total,log_vIF_fm[mask_otf])
        vIF_fd_total = np.append(vIF_fd_total,log_vIF_fd[mask_otf])
        F_lya_total = np.append(F_lya_total,F_lya_otf)
        F_inc_total = np.append(F_inc_total,F_inc)
        c_em_total = np.append(c_em_total,c_em_factor)
        #nH2_total = np.append(nH2_total,nH2_avg)
        #three_avg_total = np.append(three_avg_total,three_avg)
        #T_reion_total = np.append(T_reion_total,T_reion)
        width_total = np.append(width_total,widthIF)
        #temp_eff_total = np.append(temp_eff_total,temp_eff)
        # frat_total = np.append(frat_total,flux_ratio)
        #for m in range(Nv):
        #    if (log_vIF_fm[n]>logvIF_binedges[m])&((log_vIF_fm[n]<logvIF_binedges[m+1])):
        #       density_matrix[i][m] += np.sum(nH_arr[mask_gas])
        #       density2_matrix[i][m] += np.sum(nH_arr[mask_gas]**2)
        #       num_matrix[i][m] += len(nH_arr[mask_gas])

#clumping_matrix = density2_matrix/density_matrix**2*num_matrix
#np.save("clumpMatrix.npy",clumping_matrix)
save_dir = "results/220919/"
np.save(save_dir+"alpha.npy", alpha_total)
np.save(save_dir+"vIF.npy", vIF_total)
np.save(save_dir+"vIF_fd.npy", vIF_fd_total)
np.save(save_dir+"FLya.npy", F_lya_total)
np.save(save_dir+"Finc.npy", F_inc_total)
np.save(save_dir+"C_em.npy", c_em_total)
#np.save(save_dir+"temp_eff.npy", temp_eff_total)
np.save(save_dir+"width.npy", width_total)
#np.save(save_dir+"T_reion.npy", T_reion_total)
#np.save(save_dir+"three_avg.npy",three_avg_total)
#np.save(save_dir+"nH2.npy",nH2_total)

print("python code complete")
