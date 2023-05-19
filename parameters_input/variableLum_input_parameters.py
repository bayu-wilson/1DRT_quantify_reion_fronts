import numpy as np

h=6.626068e-27 #erg s
h_eV=4.135667e-15 #eV s
nu0=13.6/h_eV
alpha=1.5
pi=np.pi

Omega_b = 0.048
z = 5.7
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rho_crit_z_b = rho_crit_z*Omega_b
m_H =  1.672622e-24 #g
Y = 0.24
n_H = (1-Y)*rho_crit_z_b/m_H
print(n_H)
n_H *=8 #BAYU MAY 16, 2023. Increasing mean density by 10. This should decrease IF width by 10!
print(n_H)
chi = 0.08


cm_per_kpc = 3.086e+21
cm_per_km  = 1e+5
logvIF_min,logvIF_max,nbins = 3.0,4.5,4

logvIF_bincenters = np.linspace(logvIF_min,logvIF_max,nbins)
delta_logvIF = (logvIF_max-logvIF_min)/(nbins-1)
logvIF_binedges = np.linspace(logvIF_min-delta_logvIF/2,logvIF_max+delta_logvIF/2,nbins+1)
vIF_binedges_cms = 10**logvIF_binedges * cm_per_km #km/s *cm/km = cm/s
vIF_bincenters_cms = 10**logvIF_bincenters * cm_per_km

R0 = 1e3*cm_per_kpc #pkpc to cm
R_buffer = 100*cm_per_kpc
sec_per_yr = 3.154e7
yr_per_Myr = 1e6

def get_ionizing_flux(L,r,R0): #photons/s/cm2
    return 1/h/nu0 * (alpha-1)/alpha * (4**-alpha-1)/(4**(1-alpha)-1) * L/4/pi/(R0+r)**2

def get_IF_speed(flux,n_H): #cms
    return flux/n_H/(1-chi)

def get_luminosity(speedIF,r,R0,n_H):
    return speedIF * h*nu0 * alpha/(alpha-1) * (4**(1-alpha)-1)/(4**-alpha-1) * 4*pi*(R0+r)**2 * n_H*(1+chi)

def get_skewer_length(speedIF,L,R0,n_H):
    radicand=1/h/nu0 * (alpha-1)/alpha * (4**-alpha-1)/(4**(1-alpha)-1) * L/4/pi/speedIF * 1/n_H/(1-chi)
    return np.sqrt(radicand)-R0

def get_time(L,R,R0,n_H):
    coeff = h*nu0 * alpha/(alpha-1) * (4**(1-alpha)-1)/(4**-alpha-1) * 4*pi/L * n_H*(1+chi)
    return coeff*1/3*((2*R0+R)**3-(2*R0)**3)/1.75 #1.75 is added via calibration!!!

def get_parameter_input(min_speed,max_speed,R0,n_H):
    lum_to_get_max_speed = get_luminosity(speedIF=max_speed,r=R_buffer,R0=R0,n_H=n_H) #erg/s
    len_to_get_min_speed = get_skewer_length(speedIF=min_speed,L=lum_to_get_max_speed,R0=R0,n_H=n_H)#/cm_per_kpc
    time_to_cross_len = get_time(L=lum_to_get_max_speed, R=len_to_get_min_speed, R0=R0,n_H=n_H)
    return lum_to_get_max_speed,len_to_get_min_speed,time_to_cross_len

def print_readable_units(min_speed,max_speed,lum_to_get_max_speed,len_to_get_min_speed,time_to_cross_len):
    print("Say we want to explore IF speeds from {:.2e} to {:.2e} km/s".format(min_speed/1e5,max_speed/1e5))
    print("In order for the maximum vIF to be reached at the beginning of the skewer,L={:.2e} erg/s".format(
        lum_to_get_max_speed))
    print("In order for flux to geometrically attentuate enough to reach the minimum speed, R_sim={:.2f} pkpc".format(
        len_to_get_min_speed/cm_per_kpc))
    print("The time it will take for the IF to travel accross the whole skewer length is t_sim={:.2f} Myr".format(
        time_to_cross_len/sec_per_yr/yr_per_Myr))


stack = []
for i in range(len(vIF_binedges_cms)-1):
    max_speed = vIF_binedges_cms[i+1]
    min_speed = vIF_binedges_cms[i]
    args = (min_speed,max_speed
           )+get_parameter_input(min_speed,max_speed,R0,n_H)
    L,R,t = get_parameter_input(min_speed,max_speed,R0,n_H)
    R+=R_buffer
    R/=cm_per_kpc
    t/=(sec_per_yr*yr_per_Myr)
    stack.append([L,R,R0/cm_per_kpc,t])
    print_readable_units(*args)
    print("")

print("Four columns: L, R, R0, t")
np.savetxt(fname="input_params/varyLum.txt",X=stack,fmt=('%.2e','%d','%.1e','%.1f'))
print(np.loadtxt(fname="input_params/varyLum.txt",dtype=str))

#stack = np.column_stack((skewer_temp,dotN_approx_arr,t_sim_arr))
#print(stack)
#np.savetxt(fname="input_params/flucRho.txt",X=stack,fmt=('%04d','%.5e','%.1f'))
#print(np.loadtxt(fname="input_params/flucRho.txt",dtype=str))

