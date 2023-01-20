import numpy as np

### ALL UNITS CGS UNLESS STATED OTHERWISE ###
omega_b = 0.048
m_H = 1.672622e-24 #g
Y = 0.24 #helium fraction
h = 6.626068e-27 #erg s
c = 2.998e10 #cm/s
k_B = 1.38065e-16
eV_cgs = 1.60218e-12
pi=np.pi
lambda_lya_cm = 1.21567e-5
kpc_to_cm = 3.086e21
chi = 0.08

#cross-section parameters from Verner et. al. 1996, Table 1
sigma0_H1 = 5.475e-14
Eth_H1 =1.360e1

z = 5.7
rho_crit_z = 8.688e-30*(1+z)**3 #g/cm3
rhoBaryons_0 = rho_crit_z*omega_b #g/cm3
