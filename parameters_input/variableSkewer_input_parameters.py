import numpy as np
import pandas as pd

#user-input
skewer_min = 0
skewer_max = 9
len_skewers = skewer_max-skewer_min+1
desired_vIF = 3*10**4*1e5 #[km/s][cm/km]=cm/s
R_start = 1000
R_sim = 1000

#constants
z = 5.7
omega_b = 0.048
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rhoBaryons_0 = rho_crit_z*omega_b #g/cm3
m_h = 1.672622e-24 #g
pi = np.pi
Y = 0.24 #helium fraction
chi = 0.08
nH_const = (1-Y)*rhoBaryons_0/m_h #hydrogen number density
c_cms = 3e+10 #km/s
km_to_cm =1e+5 #cm/km
kpc_to_km = 3.086e+16 #km/kpc
kpc_to_cm = 3.086e+21 #cm/kpc
myr_to_s = 1.e+6*3.154e+7 #s/myr
mpc_to_km = 3.086e19 #km/mpc
R_naught = 1e+4*kpc_to_cm #cm #to minimize geometric attenuation
t_sim_initial = 100 * myr_to_s #100myr in seconds. I will add a 100kpc cushion so the actual t_sim will be larger
R_cushion = 100 #kpc
#min_log_Ndot = 56 #Ndot = 10^56 photons/s
#max_log_Ndot = 57
bins_log_Ndot = 6
mpc_to_kpc = 1e+3 #kpc/mpc
cell_size = 0.735 #kpc #N_r*cell_size = sim_size

#functions
def get_vIF(n_dot):
    #IF speed given dot N, through density nH, at R0 distance away from source. including correction for relativistic speeds 
    v_kms_s = c_cms*n_dot/(n_dot + 4*pi*R_naught**2*c_cms*nH_const*(1+chi)) / km_to_cm #[cm/s][km/cm]=[km/s]
    return v_kms_s
def get_Rsim(n_dot):
    #size of simulation for the IF to pass through it and sample it well
    Rsim_func = get_vIF(n_dot)/mpc_to_km*t_sim_initial*mpc_to_kpc #[km/s][mpc/km][s][kpc/mpc]=[kpc]
    return Rsim_func + R_cushion
def get_tsim(n_dot):
    t_sim_cushion = (R_cushion/get_vIF(n_dot))*kpc_to_km/myr_to_s #[kpc][s/km][km/kpc][myr/s]=[Myr]
    return t_sim_cushion+t_sim_initial/myr_to_s #time to pass through cushion as well as the rest of the skewer



#parameter calculations
skewer_temp = np.arange(skewer_min,skewer_max+1,1)
skewer_arr = [str(i).zfill(4) for i in skewer_temp]
matt_col=["velocity [km/s]", "tau_HILya", "tau_HeI584", "tau_HeIILya", "nHI [cm^-3]",
    "nHeII [cm^-3]", "Delta_b", "T [K]", "los pos [kpc/h]", "vpec [km/sec]"]
dotN_approx_arr = np.zeros([len_skewers])
t_sim_arr = np.zeros([len_skewers])
for i in range(len_skewers):
    df_matt = pd.read_csv("../input_files/hydro_skewers/spec_xHeII1_007_ls_line{}.dat".format(skewer_arr[i]),delim_whitespace=True,
                       skiprows=4,names=matt_col)
    rho_b = df_matt["Delta_b"].values*rhoBaryons_0
    nH_matt = (1-Y)*rho_b/m_h
    z_mmc=7.19
    h_const=0.68
    los_pos_adjusted = df_matt["los pos [kpc/h]"]/(1+z_mmc)/h_const
    mask = (los_pos_adjusted>R_start)&(los_pos_adjusted<R_start+R_sim)
    nH_mean=np.mean(nH_matt[mask])

    dotN_approx = desired_vIF*4*pi*R_naught**2*nH_mean
    dotN_approx_arr[i]=dotN_approx

    t_sim_arr[i] = R_sim*kpc_to_cm/desired_vIF/myr_to_s
    #dotN_approx_mean+=dotN_approx
    #print(dotN_approx)
#dotN_approx_mean /= len(skewer_arr)
stack = np.column_stack((skewer_temp,dotN_approx_arr,t_sim_arr)) 
#print(stack)
np.savetxt(fname="input_params/flucRho.txt",X=stack,fmt=('%04d','%.5e','%.1f'))
print(np.loadtxt(fname="input_params/flucRho.txt",dtype=str))

#print(dotN_approx_mean)
#v_cms_s = c_cms*n_dot/(n_dot + 4*pi*R_test**2*c_cms*x[index]*(1+ 0.08)) #km/s
#dotN_arr = np.logspace(min_log_Ndot,max_log_Ndot,bins_log_Ndot) #ionizing photon rate
#R_sim_arr = get_Rsim(dotN_arr) #size of sim plus 100 pkc cushion
#N_r_arr = np.asarray(R_sim_arr/cell_size-1,int) #number of cells to maintain 0.735 pkpc cell size
#t_sim_arr = get_tsim(dotN_arr) # simulation time required for IF to pass through calculated skewer PLUS 100pkc beginning part

#stack = np.column_stack((dotN_arr,R_sim_arr,N_r_arr,t_sim_arr))
#np.savetxt(fname="input_params/constRho.txt",X=stack,fmt=('%.5e','%.5e','%d','%.1f'))
#np.loadtxt("test.txt",dtype=str)




