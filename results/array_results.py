import numpy as np
import os.path
import pandas as pd
#import matplotlib.pyplot as plt

N_output = 250
kpc_to_cm = 3.086e21
Y = 0.24
m_H = 1.672622e-24
h = 6.626068e-27
c = 2.998e10
lambda_lya_cm = 1.21567e-5
pi=np.pi

#amin = -1.000
#amax = 2.500 #2.505

np.seterr(divide='ignore', invalid='ignore')

#alpha_list = np.arange(amin,amax+0.125,0.125)
#alpha_list = np.arange(amin,amax+0.125,0.5)#0.125)
#print(alpha_list,flush=True)
alpha_list = np.loadtxt("../parameters_input/input_params/spectral_indices.txt")
#Nv = len(alpha_list)
#logvIF_bincenters =  np.linspace(vlogmin,vlogmax,Nv)
#grid_V, grid_A = np.mgrid[vlogmin:vlogmax:complex(N_grid), amin:amax:complex(N_grid)]
#delta_logvIF = logvIF_bincenters[1]-logvIF_bincenters[0]
#logvIF_binedges = np.linspace(vlogmin-delta_logvIF/2,
#                            vlogmax+delta_logvIF/2,Nv+1)
#skewers = ["sk0000","sk0001","sk0002","sk0003","sk0004","sk0005","sk0006","sk0007","sk0008","sk0009"]
skewers = ["sk0000","sk0001","sk0002","sk0003"]
#dotN = np.loadtxt("../input_files/dotN_Rsim.txt",dtype=str)[:,0]
#dotN = np.loadtxt("../parameters_input/input_params/flucRho.txt",dtype=str)[:,1]

#otf_names = ["step","t","dt","I_lya","Finc","C_factor","locIF","vIF_fd","vIF_fm","T_IF","density_var","IF_width","nH_avg"]
otf_names = ["step", "t", "dt", "F_lya", "F_inc", "locIF","vIF_fd","vIF_Flex","T_center","width_IF","nH_avg","nH_center","C_em"]
total_names = ["skewer","alpha"]+otf_names

#column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI","xHeII","xHeIII","q_lya"]

print("begin for loop")
otf_total_df = pd.DataFrame(columns=total_names)
len_otf = len(otf_names)
for sk in skewers: #looping through the skewers
    for i in range(len(alpha_list)): #looping through spectral index
        #print(sk,i)
        spectral_index = alpha_list[i]
        n=1
        dir_path = "../output_files/gasprops/{}_a={:.3f}/".format(sk,spectral_index) #I shouldn't need the .3f here
        #print(dir_path,flush=True)
        gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
        otf_matrix = np.zeros((N_output,len(otf_names)))
        while (os.path.exists(gasprop_path_otf))&(n<=N_output):
            with open(gasprop_path_otf,'r') as f:
                line = f.readline()
            try:
                line_array = np.asarray(line.split(' \t'),float)
                #line_list = line.split(' \t')
                otf_matrix[n-1]=line_array #np.array([sk[2:],spectral_index]+line_list)
                if line_array[3]<10: #skip first 10 Myr
                    n_minTime = n
            except:
                otf_matrix[n-1]=np.ones(len(otf_names))*np.nan
            n+=1
            gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
        otf_df = pd.DataFrame(otf_matrix,columns=otf_names)[1:]
        otf_df["skewer"] = sk[2:]
        otf_df["alpha"] = spectral_index
        otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm-1e4
        otf_df["vIF_fd"]/=1e5
        otf_df["vIF_Flex"]/=1e5
        mask_otf = np.asarray(np.ones(len(otf_df)),bool)
            #mask_otf = np.ones_like(otf_df["I_lya"]) #((otf_df["t"]>10)&(locIF[n]<(R_sim-width_IF)))#&
            #            (log_vIF_fm>=2.0)&(log_vIF_fm<=5.0))
        
        otf_total_df = pd.concat([otf_total_df,otf_df],ignore_index=True)
        
save_dir = "230413/"
otf_total_df.to_csv(save_dir+"otf.csv",index=False)

"""
vIF_total = np.array([])
vIF_fd_total = np.array([])
F_lya_total = np.array([])
F_inc_total = np.array([])
alpha_total = np.array([])
clump_total = np.array([])
nH_total =  np.array([])
rho_var_total= np.array([])
width_total= np.array([])

#density_matrix = np.zeros((Nv,Nv))
#density2_matrix = np.zeros((Nv,Nv))
#num_matrix = np.zeros((Nv,Nv))
print("running python code")
for sk in skewers: #looping through the skewers
    for i in range(len(alpha_list)): #looping through spectral index
        spectral_index = alpha_list[i]
        gasprop_path = "gasprops/dotN{}_a={:<05}/".format(dn,spectral_index)
        #print(gasprop_path) #gasprops/sk0000_a\=2.505/
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
            backIF = np.interp(0.05,df_gasprops["xHI"],los_pos)
            frontIF = np.interp(0.75,df_gasprops["xHI"],los_pos)
            mask_gas = (los_pos>backIF)&(los_pos<frontIF)

            otf_df = pd.DataFrame(output_matrix.T,columns=otf_names)
            locIF = (otf_df["locIF"]/kpc_to_cm-1e4)
            log_vIF_fd = np.log10(otf_df["vIF_fd"]/1e5)
            log_vIF_fm = np.log10(otf_df["vIF_fm"]/1e5)
            mask_otf = np.asarray(np.ones(len(otf_df)),bool)
            #mask_otf = np.ones_like(otf_df["I_lya"]) #((otf_df["t"]>10)&(locIF[n]<(R_sim-width_IF)))#&
            #            (log_vIF_fm>=2.0)&(log_vIF_fm<=5.0))
            F_lya_otf = otf_df["I_lya"][mask_otf]/(h*c/lambda_lya_cm)*4*pi
            F_inc = otf_df["Finc"][mask_otf]
            C_factor = otf_df["C_factor"][mask_otf]
            nH_avg = otf_df["nH_avg"][mask_otf]
            widthIF = otf_df["IF_width"][mask_otf]
            rho_var = otf_df["density_var"][mask_otf]
            # flux_ratio = (F_lya_otf/F_inc)
            alpha_total = np.append(alpha_total,np.ones_like(log_vIF_fm[mask_otf])*float(spectral_index))
            vIF_total = np.append(vIF_total,log_vIF_fm[mask_otf])
            vIF_fd_total = np.append(vIF_fd_total,log_vIF_fd[mask_otf])
            F_lya_total = np.append(F_lya_total,F_lya_otf)
            F_inc_total = np.append(F_inc_total,F_inc)
            clump_total = np.append(clump_total,C_factor)
            nH_total = np.append(nH_total,nH_avg)
            rho_var_total = np.append(rho_var_total,rho_var)
            width_total = np.append(width_total,widthIF)
            # frat_total = np.append(frat_total,flux_ratio)
            #for m in range(Nv):
            #    if (log_vIF_fm[n]>logvIF_binedges[m])&((log_vIF_fm[n]<logvIF_binedges[m+1])):
            #       density_matrix[i][m] += np.sum(nH_arr[mask_gas])
            #       density2_matrix[i][m] += np.sum(nH_arr[mask_gas]**2)
            #       num_matrix[i][m] += len(nH_arr[mask_gas])

#clumping_matrix = density2_matrix/density_matrix**2*num_matrix
#np.save("clumpMatrix.npy",clumping_matrix)
save_dir = "results/221027/"
np.save(save_dir+"alpha.npy", alpha_total)
np.save(save_dir+"vIF.npy", vIF_total)
np.save(save_dir+"vIF_fd.npy", vIF_fd_total)
np.save(save_dir+"FLya.npy", F_lya_total)
np.save(save_dir+"Finc.npy", F_inc_total)
np.save(save_dir+"C.npy", clump_total)
np.save(save_dir+"width.npy", width_total)
np.save(save_dir+"rho_var.npy",rho_var_total)
np.save(save_dir+"nH.npy",nH_total)

"""
#print("python code complete")
