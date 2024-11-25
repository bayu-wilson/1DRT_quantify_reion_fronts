import numpy as np
#Constants
h=6.626068e-27 #erg s
h_eV=4.135667e-15 #eV s
nu0=13.6/h_eV
# alpha=1.5
pi=np.pi

Omega_b = 0.048
z = 5.7
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rho_crit_z_b = rho_crit_z*Omega_b
m_H =  1.672622e-24 #g
Y = 0.24
n_H = (1-Y)*rho_crit_z_b/m_H
# print(n_H)
# n_H *=1 #BAYU MAY 16, 2023. Increasing mean density by 10. This should decrease IF width by 10!
# print(n_H)
chi = 0.08

cm_per_kpc = 3.086e+21
cm_per_km  = 1e+5
sec_per_yr = 3.154e7
yr_per_Myr = 1e6

R0 = 1e4*cm_per_kpc

#binning for vIF
logvIF_min,logvIF_max,nbins = 1.5,5,4
alpha_min,alpha_max,alpha_bins = 0,3,7

logvIF_bincenters = np.linspace(logvIF_min,logvIF_max,nbins)
delta_logvIF = (logvIF_max-logvIF_min)/(nbins-1)
logvIF_binedges = np.linspace(logvIF_min-delta_logvIF/2,logvIF_max+delta_logvIF/2,nbins+1)
vIF_binedges_cms = 10**logvIF_binedges * cm_per_km #km/s *cm/km = cm/s
vIF_bincenters_cms = 10**logvIF_bincenters * cm_per_km

alpha_bincenters = np.linspace(alpha_min,alpha_max,alpha_bins)

#functions
def alpha_dependence(alpha):
    if 1-alpha == 0:
        return 8*np.log(2)/3
    elif alpha == 0:
        return 3/2/np.log(2)
    else:
        return alpha/(alpha-1) * (4**(1-alpha)-1)/(4**-alpha-1)
def get_luminosity(speedIF,R0,n_H,alpha):
    #return speedIF * h*nu0 * alpha/(alpha-1) * (4**(1-alpha)-1)/(4**-alpha-1) * 4*pi*(R0+r)**2 * n_H*(1+chi)
    return speedIF * h*nu0 * alpha_dependence(alpha) * 4*pi*(R0)**2 * n_H*(1+chi)

def time_decay_factor(t,t0,eta):
    # t,t0 must be same units
    if power_law_index<0:
        raise ValueError('Power-law index must be non-negative to decay')
        exit(1)
    return (t/t0)**(-eta)

def get_final_time(t0,min_speedIF, max_speedIF,eta):
    return t0*(min_speedIF/max_speedIF)**(-1/eta)
# def get_final_time(t0,speedIF, R0, n_H, L,eta,alpha):
#     base = speedIF*h*nu0 * alpha_dependence(alpha) * 4*pi*(R0)**2 * n_H*(1+chi)/L
#     return t0 * base**(-1/eta)
    # return t0/speedIF /h/nu0 * (alpha-1)/alpha * (4**-alpha-1)/(4**(1-alpha)-1) / 4*pi*(R0)**2 / n_H/(1*chi) * L
    
def get_skewer_length(speedIF, t0, t_end,eta):
    #t0,t_end must be in seconds
    if eta==1:
        t0*=(sec_per_yr*yr_per_Myr)
        t_end*=(sec_per_yr*yr_per_Myr)
        return speedIF*t0*np.log(t_end/t0)
    else:
        t0*=(sec_per_yr*yr_per_Myr)
        t_end*=(sec_per_yr*yr_per_Myr)
        return speedIF/(1-eta)*(t_end*(t_end/t0)**(-eta)-t0)
        # raise ValueError('This only works for eta=1 right now')
        
#t0=20
eta=2
#vIF_test_edges = [5.e+6, 8.e+7, 8.e+8, 5.e+9, 1.e+10]
#vIF_test_edges = [1.e+7, 2.e+9, 6.e+9, 1.e+10]
t0 = [40,20,20]
#vIF_test_edges = [1.e+7, 2e+8, 1.e+9, 6.e+9, 1.e+10]
#vIF_test_edges = [1.e+7, 8.e+8, 5.e+9, 1.e+10]
vIF_test_edges = [[7e+6,9e+8],[2e+8, 5.5e+9],[2.5e+9, 1.1e+10]]

stack = []
for alpha in alpha_bincenters:
    for i in range(len(vIF_test_edges)):  
    #for i in range(len(vIF_binedges_cms)-1):
        #max_speed = vIF_binedges_cms[i+1]
        #min_speed = vIF_binedges_cms[i]
        #max_speed = vIF_test_edges[i+1]*1.1
        #min_speed = vIF_test_edges[i]*0.5
        max_speed = vIF_test_edges[i][1]
        min_speed = vIF_test_edges[i][0]
        L = get_luminosity(max_speed,R0,n_H,alpha)
        # tmax = get_final_time(t0,min_speed, R0, n_H, L,eta,alpha)
        #if i==0:
        #    t_end = get_final_time(t0[i],min_speed, max_speed,eta)*2
        #    r_skew = get_skewer_length(max_speed,t0[i],t_end,eta)/cm_per_kpc * 1.5
        #else:
        #    t_end = get_final_time(t0[i],min_speed, max_speed,eta)*1.2
        t_end = get_final_time(t0[i],min_speed, max_speed,eta)*1.2
        r_skew = get_skewer_length(max_speed,t0[i],t_end,eta)/cm_per_kpc
        N_r = int(r_skew*2)
        stack.append([alpha,L,r_skew,N_r,R0/cm_per_kpc,t_end,t0[i]])
        # print(f"alpha={alpha:.1f}, vIF=({min_speed/1e5:.0e},{max_speed/1e5:.0e}) [km/s], L={L:.1e}[erg/s], "+
              # f"t_end={t_end:.1f} [Myr], r={r_skew:.1f}")
    # print("")
    
np.savetxt(fname="input_params/planeParallel_params.txt",X=stack,fmt=('%.1f','%.2e','%.1f','%d','%.1e','%.1f','%.1f'))
print("alpha, L[erg/s], r_skew[pkpc], N_r, R0[pkpc], t_end[Myr], t0[Myr]")
print(np.loadtxt(fname="input_params/planeParallel_params.txt",dtype=str))
