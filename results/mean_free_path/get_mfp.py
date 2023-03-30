import numpy as np
import pandas as pd
import glob

# Constants
kpc_to_cm =  3.086e21
m_H = 1.672622e-24
Y = 0.24
sigma0_LyC = 6.34e-18 #cm2
delta_r = 1 #kpc
delta_r*=kpc_to_cm #cm

# Functions
def get_tau(delta_r,nHI,idx_locIF):
    return delta_r*sigma0_LyC*np.sum(nHI[:idx_locIF])

def mean_free_path(r,delta_r,nHI,idx_locIF): #cgs
    constants = sigma0_LyC*delta_r
    tau = np.array([get_tau(delta_r,nHI,i) for i in range(idx_locIF)])
    r_nHI_subset = (r*nHI)[0:idx_locIF]
    mfp = constants*np.sum((r_nHI_subset*np.exp(-tau)))
    return mfp

# Load in sightlines
column_names = ["los pos[pkpc]","rho","ne","ne_prev","n_tot","T[K]","xHI","xHeI",
                "xHeII","xHeIII","q_lya"]
skewer_number = ["{:04d}".format(i) for i in range(0,41)]
mfp_list = []
for sk in range(len(skewer_number)):
    dir_path = "../../output_files/gasprops/sk{0}_hardRun/".format(skewer_number[sk])
    files = glob.glob(dir_path+"n*_gasprops.txt")
    final_step = max([int(f.split("/")[-1][1:-13]) for f in files])
    final_path = dir_path + "n{}_gasprops.txt".format(final_step)
    df = pd.read_csv(final_path,delim_whitespace=True,skiprows=1,names=column_names)
    los_pos = (df["los pos[pkpc]"]-df["los pos[pkpc]"][0])*kpc_to_cm
    N = len(los_pos)
    nH = (1. - Y)*df['rho']/m_H
    nHI = nH*df["xHI"]
    try:
        idx_locIF = int(np.interp(x=0.5, fp=np.arange(0,N,1), xp=df["xHI"])+0.5)
        mfp = mean_free_path(r=los_pos,delta_r=delta_r,nHI=nHI,idx_locIF=idx_locIF)/kpc_to_cm
        mfp_list.append(mfp)
        print(mfp)
    except:
        idx_locIF = N-1
        mfp = mean_free_path(r=los_pos,delta_r=delta_r,nHI=nHI,idx_locIF=idx_locIF)/kpc_to_cm
        mfp_list.append(mfp)
        print("larger then skewer: sk:{}, mfp:{}".format(sk,mfp))

print(np.mean(mfp_list))
