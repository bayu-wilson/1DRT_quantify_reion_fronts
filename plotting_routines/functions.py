import numpy as np
from constants import *

def get_inc_LyC_flux(Lum,spectral_index,R0_cm,t_Myr,t0_Myr,eta):
    if eta == 0: #no time dependence
        time_dependence = 1
    else:
        time_dependence = np.power((t_Myr+t0_Myr)/t0_Myr,-eta)
    if spectral_index == 0:
        alpha_dependence = 2*np.log(2)/3
    elif spectral_index == 1:
        alpha_dependence = 3/8/np.log(2)
    else:
        # pass
        alpha_dependence = (spectral_index-1)/spectral_index*(4**-spectral_index-1)/(4**(1-spectral_index)-1)
    return 1/h/nu0*alpha_dependence*Lum/4/np.pi/R0_cm**2 * time_dependence


def llambda(T, Eth):
    K2eV = k_B/eV_cgs
    eV2K = 1.0/K2eV
    xHI  = Eth*eV2K/T
    x    = 2*xHI
    return x

#Case A recombination cooling rate for HII from Hui & Gnedin 1996
def rcrA_H2(T):
    x = llambda(T, Eth_H1)
    return 1.778e-29*T*pow(x,1.965)/(pow(1. + pow(x/0.541,0.502),2.697))
#Case A recombination coefficient of HII from Hui & Gnedin 1996
def alphaA_H2(T):
    x = llambda(T, Eth_H1)
    return 1.269e-13*pow(x,1.503)/pow((1.0+pow((x/0.522),0.470)),1.923)


#Case B recombination coefficient of HII from Hui & Gnedin 1996
def alphaB_H2(T):
     x = llambda(T, Eth_H1)
     return 2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420)

def epsilonB_lya(T):
    T4 = T/1e+4
    return 0.686 - 0.106*np.log10(T4) - 0.009*T4**(-0.44)


#collisional excitation cooling rates from Cen 1992
def coll_ex_rate_H1(T):
        return 7.5e-19*pow((1 + pow(T/1e5, 0.5)), -1)*np.exp(-118348./T)

#// in units of [cm^3 / s]
# Uses polynomial fits from Giovanardi+1987, also used in Cantalupo+2008
def coll_ex_rate_H1_acc(T):
    hvHI    = 2.178717e-11  #in erg
    E_1s_2p = 0.75*hvHI
    E_1s_3s = (1.0-(1.0/9.0))*hvHI
    E_1s_3d = (1.0-(1.0/9.0))*hvHI
    E_1s_3p = (1.0-(1.0/9.0))*hvHI
    E_1s_4s = (1.0-(1.0/16.0))*hvHI
    E_1s_4p = (1.0-(1.0/16.0))*hvHI
    E_1s_4d = (1.0-(1.0/16.0))*hvHI
    E_1s_4f = (1.0-(1.0/16.0))*hvHI
    k_B = 1.38065e-16 # //Boltzmann constant (erg/K)

    T2 = pow(T,2.)
    T3 = pow(T,3.)
    Omega_1s_2p = 0.
    Omega_1s_3s = 0.
    Omega_1s_3p = 0. #//hmm?
    Omega_1s_3d = 0.
    Omega_1s_4p = 0. #//hmm?
    # Omega_1s_4s = 0 #//hmm?
    # Omega_1s_4d = 0 #//hmm?
    # Omega_1s_4f = 0 #//hmm?

    if (T > 5000.) & (T < 60000.): # //These coefficients good between 5000 and 72000K
        Omega_1s_2p = 0.3435 + 1.297e-5*T + 2.178e-12*T2 + 7.928e-17*T3
        Omega_1s_3s = 0.0625 - 1.299e-6*T + 2.666e-11*T2 - 1.596e-16*T3
        Omega_1s_3p = 0.0994 - 3.714e-7*T + 6.134e-11*T2 - 3.973e-16*T3
        Omega_1s_3d = 0.0503 + 7.514e-7*T - 2.826e-13*T2 - 1.098e-17*T3
        Omega_1s_4p = 1.527e-3+1.001e-6*T - 2.192e-12*T2 + 9.348e-18*T3
        #// Omega_1s_4s = 0 #//Can add these coefficients later for more precision
        #// Omega_1s_4d = 0
        #// Omega_1s_4f = 0
    elif (T >= 60000.) & (T < 500000.):
        Omega_1s_2p = 0.3162 + 1.472e-5*T - 8.275e-12*T2 - 8.794e-19*T3
        Omega_1s_3s = 0.03337 + 2.223e-7*T - 2.794e-13*T2 + 1.516e-19*T3
        Omega_1s_3p = 0.06985 + 2.538e-6*T - 8.729e-13*T2 - 1.291e-18*T3
        Omega_1s_3d = 0.05051 + 7.876e-7*T - 2.072e-12*T2 + 1.902e-18*T3
        Omega_1s_4p = 1.958e-3+9.525e-7*T - 9.668e-13*T2 + 4.807e-19*T3
        #// Omega_1s_4s = 0;
        #// Omega_1s_4d = 0;
        #// Omega_1s_4f = 0;
    else: #//At low temperatures the rates are small, and T should not get so hot
        Omega_1s_2p = 0.
        Omega_1s_3s = 0.
        Omega_1s_3p = 0.
        Omega_1s_3d = 0.
        Omega_1s_4p = 0.
        #// Omega_1s_4s = 0
        #// Omega_1s_4d = 0
        #// Omega_1s_4f = 0
    prefix = 8.629e-6/np.sqrt(T)
    q_1s_2p = prefix * Omega_1s_2p / 2.0 * np.exp(-E_1s_2p/(k_B*T))
    q_1s_3s = prefix * Omega_1s_3s / 2.0 * np.exp(-E_1s_3s/(k_B*T))
    q_1s_3p = prefix * Omega_1s_3p / 2.0 * np.exp(-E_1s_3p/(k_B*T))# //remove?
    q_1s_3d = prefix * Omega_1s_3d / 2.0 * np.exp(-E_1s_3d/(k_B*T))
    q_1s_4p = prefix * Omega_1s_4p / 2.0 * np.exp(-E_1s_4p/(k_B*T))  #//remove?
    #return q_1s_2p + q_1s_3s + q_1s_ 3d; //eq. 5 of Cantalupo+2008.  Can be improved.
    return q_1s_2p + q_1s_3s + q_1s_3d + q_1s_3p +q_1s_4p #//eq. 5 of Cantalupo+2008.  Can be improved.
    sigma0_H1 = 5.475e-14 #cm-2
    Eth_H1 = 1.360e+1 #erg
    Emax_H1  = 5e+4 #erg
    E0_H1 = 4.298e-1 #erg
    y0_H1 = 0  #not sure
    y1_H1 = 0 #not sure
    yw_H1 = 0 #not sure
    ya_H1 = 3.288e+1 #not sure
    P_H1 = 2.963 #not sure
    h_eV = 4.135667e-15 #eV s
    nu_0 = 13.6/h_eV
    def F_of_y(nu, Eth, Emax, E0, y0, y1, yw, ya, P):
        E = h_eV*nu
        x = E/E0 - y0
        y = pow(pow(x,2) + pow(y1,2),0.5)
        F = (pow(x-1,2) + pow(yw,2))*pow(y,0.5*P - 5.5)*pow(1 + pow(y/ya,0.5),-P)
        if (h_eV*nu >= Eth)&(h_eV*nu <= Emax):
            return F
        else:
            return 0
    def sigmapi_H1(nu):
        return sigma0_H1*F_of_y(nu, Eth_H1, Emax_H1, E0_H1, y0_H1, y1_H1, yw_H1, ya_H1, P_H1)
    # sigmapi_H1(nu_0)
