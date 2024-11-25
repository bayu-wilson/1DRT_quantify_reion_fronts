import pandas as pd
import numpy as np
import os.path

N_output = 500
kpc_to_cm = 3.086e21
Y = 0.24
m_H = 1.672622e-24
h = 6.626068e-27
c = 2.998e10
lambda_lya_cm = 1.21567e-5
pi=np.pi


otf_names = ["step", "t", "dt", "F_lya", "F_inc", "locIF","vIF_fd","vIF_FlexRT","T_avg","T_center","width_IF","nH_center","vIF_test"]

dir_path = "../output_files/gasprops/vIF_test_sk0052/"
n=1
print(dir_path,flush=True)
gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
otf_matrix = np.zeros((N_output,len(otf_names)))
while (os.path.exists(gasprop_path_otf))&(n<=N_output):
    with open(gasprop_path_otf,'r') as f:
        line = f.readline()
    try:
        line_array = np.asarray(line.split(' \t'),float)
        otf_matrix[n-1]=line_array
        if line_array[3]<5:
            n_minTime = n
        
    except:
        otf_matrix[n-1]=np.ones(len(otf_names))*np.nan
    print(n)
    n+=1
    gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)

print(n)
#print(otf_matrix.shape)

otf_df = pd.DataFrame(otf_matrix,columns=otf_names)[1:]
otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm-1e4
otf_df["vIF_fd"]/=1e5
otf_df["vIF_FlexRT"]/=1e5
otf_df["vIF_test"]/=1e5



save_dir = "241104/"
otf_df.to_csv(save_dir+"otf.csv",index=False)

"""
otf_total_df = pd.DataFrame(columns=total_names)
len_otf = len(otf_names)
for sk in range(len(skewers)): #looping through the skewers
    for i in [0]:#range(len(alpha_list)): #looping through spectral index
        #print(sk,i)
        spectral_index = alpha_list[0]
        n=1
        dir_path = f"../output_files/gasprops/fd_pp_a=0.0_{skewers[sk]}_L={Lums[sk]}_v2/"
        #dir_path = "../output_files/gasprops/{}_a={}_test/".format(sk,spectral_index) #I shouldn't need the .3f here
        print(dir_path,flush=True)
        gasprop_path_otf = dir_path + "n{}_gasprops.txt".format(n)
        otf_matrix = np.zeros((N_output,len(otf_names)))
        while (os.path.exists(gasprop_path_otf))&(n<=N_output):
            with open(gasprop_path_otf,'r') as f:
                line = f.readline()
            try:
                line_array = np.asarray(line.split(' \t'),float)
                #line_list = line.split(' \t')
                otf_matrix[n-1]=line_array #np.array([sk[2:],spectral_index]+line_list)
                if line_array[3]<5: #skip first 10 Myr
                    n_minTime = n
            except:
                otf_matrix[n-1]=np.ones(len(otf_names))*np.nan
            n+=1
            gasprop_path_otf = dir_path+"n{}_gasprops.txt".format(n)
        otf_df = pd.DataFrame(otf_matrix,columns=otf_names)[1:]
        otf_df["skewer"] = skewers[sk][2:]
        otf_df["alpha"] = float(spectral_index)
        otf_df["locIF"]=otf_df["locIF"]/kpc_to_cm-1e4
        otf_df["vIF_fd"]/=1e5
        otf_df["vIF_FlexRT"]/=1e5
        otf_df["Lum"] = Lums[sk]
        mask_otf = np.asarray(np.ones(len(otf_df)),bool)
            #mask_otf = np.ones_like(otf_df["I_lya"]) #((otf_df["t"]>10)&(locIF[n]<(R_sim-width_IF)))#&
            #            (log_vIF_fm>=2.0)&(log_vIF_fm<=5.0))

        otf_total_df = pd.concat([otf_total_df,otf_df],ignore_index=True)

save_dir = "241104/"
otf_total_df.to_csv(save_dir+"otf.csv",index=False)
"""
