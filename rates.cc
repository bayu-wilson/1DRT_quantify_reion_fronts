#include <math.h>
#include "global_constants.h"
#include "rates.h"

using namespace g_constants;
using namespace rates;
double lambda(double T, double Eth)
{
    double K2eV = k_B/eV_cgs;
    double eV2K = 1.0/K2eV;
    double xHI  = Eth*eV2K/T;
    double x    = 2*xHI;
    return x;
}

double F_of_y(double nu, double Eth, double Emax, double E0, double y0, double y1, double yw, double ya, double P)
{
	double x, y;
	double E, F;

	E = h_eV*nu;
	x = E/E0 - y0;
	y = pow(pow(x,2) + pow(y1,2),0.5);
	F = (pow(x-1,2) + pow(yw,2))*pow(y,0.5*P - 5.5)*pow(1 + pow(y/ya,0.5),-P);

	if (h_eV*nu >= Eth && h_eV*nu <= Emax) {
		return F;
	}
	else  {
		return 0;
    }
}

//photoionization cross-sections

//photo-ionization cross-section of HI from Verner et. al. 1996
double sigmapi_H1(double nu) {
	return sigma0_H1*F_of_y(nu, Eth_H1, Emax_H1, E0_H1, y0_H1, y1_H1, yw_H1, ya_H1, P_H1);
}

//photo-ionization cross-section of HeI from Verner et. al. 1996
double sigmapi_He1(double nu) {
	return sigma0_He1*F_of_y(nu, Eth_He1, Emax_He1, E0_He1, y0_He1, y1_He1, yw_He1, ya_He1, P_He1);
}

//photo-ionization cross-section of HeII from Verner et. al. 1996
double sigmapi_He2(double nu) {
	return sigma0_He2*F_of_y(nu, Eth_He2, Emax_He2, E0_He2, y0_He2, y1_He2, yw_He2, ya_He2, P_He2);
}

//collisional ionization rates

//collisional ionization rate of HI from Hui & Gnedin 1996
double cic_H1(double T) {
	double x = lambda(T, Eth_H1);
	return 21.11*pow(T,-1.5)*exp(-x/2.)*pow(x,-1.089)/(pow(1. + pow(x/0.354,0.874),1.101));
}

//collisional ionization rate of HeI from Hui & Gnedin 1996
double cic_He1(double T) {
	double x = lambda(T, Eth_He1);
	return 32.38*pow(T,-1.5)*exp(-x/2.)*pow(x,-1.146)/(pow(1. + pow(x/0.416,0.987),1.056));
}

//collisional ionization rate of HeII from Hui & Gnedin 1996
double cic_He2(double T) {
	double x = lambda(T, Eth_He2);
	return 19.95*pow(T,-1.5)*exp(-x/2.)*pow(x,-1.089)/(pow(1. + pow(x/0.553,0.735),1.275));
}

//collisional ionization cooling rates

//collisional ionization cooling rate of HI from Hui & Gnedin 1996
double cicr_H1(double T) {
  return Eth_H1*ev_to_erg*cic_H1(T);
}

//collisional ionization cooling rate of HeI from Hui & Gnedin 1996
double cicr_He1(double T) {
	return Eth_He1*ev_to_erg*cic_He1(T);
}

//collisional ionization cooling rate of HeII from Hui & Gnedin 1996
double cicr_He2(double T) {
	return Eth_He2*ev_to_erg*cic_He2(T);
}

//recombination rates

//Case A recombination coefficient of HII from Hui & Gnedin 1996
double alphaA_H2(double T)
{
	double x = lambda(T, Eth_H1);
	return 1.269e-13*pow(x,1.503)/pow((1.0+pow((x/0.522),0.470)),1.923);
}

//Case B recombination coefficient of HII from Hui & Gnedin 1996
double alphaB_H2(double T)
{
	double x = lambda(T, Eth_H1);
	return 2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420);
}

//Case A recombination coefficient of HeII from Hui & Gnedin 1996
double alphaA_He2(double T)
{
	double x    = lambda(T, Eth_He1);
	return 3.0e-14*pow(x,0.654);
}

//Case B recombination coefficient of HeII from Hui & Gnedin 1996
double alphaB_He2(double T)
{
	double x    = lambda(T, Eth_He1);
	return 1.26e-14*pow(x,0.750);
}

//Case A recombination coefficient of HeIII from Hui & Gnedin 1996
double alphaA_He3(double T)
{
	double x    = lambda(T, Eth_He2);
	return 2.*1.269e-13*pow(x,1.503)/pow((1.0+pow((x/0.522),0.470)),1.923);
}

//Case B recombination coefficient of HeIII from Hui & Gnedin 1996
double alphaB_He3(double T)
{
	double x    = lambda(T, Eth_He2);
	return 2*2.753e-14*pow(x,1.500)/pow((1.0+pow((x/2.740),0.407)),2.2420);
}

//Dielectric recombination coefficient of He II from Hui & Gnedin 1996
double Dalpha_He2(double T)
{
	double x    = lambda(T, Eth_He2);
	return 1.90e-3 * pow(T, -1.5) * exp(-0.75*x/2)*(1 + 0.3*exp(-0.15*x/2));
}

//recombination cooling rates

//Case A recombination cooling rate for HII from Hui & Gnedin 1996
double rcrA_H2(double T)
{
	double x    = lambda(T, Eth_H1);
	return 1.778e-29*T*pow(x,1.965)/(pow(1. + pow(x/0.541,0.502),2.697));
}

//Case B recombination cooling rate for HII from Hui & Gnedin 1996
double rcrB_H2(double T)
{
	double x = lambda(T, Eth_H1);
	return 3.435e-30*T*pow(x,1.970)/(pow(1. + pow(x/2.250,0.376),3.720));
}

//Case A recombination cooling rate for HeII from Hui & Gnedin 1996
double rcrA_He2(double T)
{
	return k_B*T*alphaA_He2(T);
}

//Case B recombination cooling rate for HeII from Hui & Gnedin 1996
double rcrB_He2(double T)
{
	return k_B*T*alphaB_He2(T);
}

//Case A recombination cooling rate for HeIII from Hui & Gnedin 1996
double rcrA_He3(double T)
{
	double x    = lambda(T, Eth_He2);
	return 1.778e-29*T*pow(x,1.965)/(pow(1. + pow(x/0.541,0.502),2.697));
}

//Case B recombination cooling rate for HeIII from Hui & Gnedin 1996
double rcrB_He3(double T)
{
	double x = lambda(T, Eth_He2);
	return 3.435e-30*T*pow(x,1.970)/(pow(1. + pow(x/2.250,0.376),3.720));
}

//dielectric recombination cooling (He II).
double drcr_He2(double T)
{
	return 0.75*Eth_He2*ev_to_erg*Dalpha_He2(T);
}

//collisional excitation cooling rates from Cen 1992
double coll_ex_rate_H1(double T)
{
	return 7.5e-19*pow((1 + pow(T/1e5, 0.5)), -1)*exp(-118348./T);
}

// in units of [cm^3 / s]
//Uses polynomial fits from Giovanardi+1987, also used in Cantalupo+2008
double coll_ex_rate_H1_acc(double T)
{
  double T2 = pow(T,2.);
  double T3 = pow(T,3.);
  double Omega_1s_2p = 0.;
  double Omega_1s_3s = 0.;
  // double Omega_1s_3p = 0.; //hmm?
  double Omega_1s_3d = 0.;
  // double Omega_1s_4p = 0.;//hmm?
  // double Omega_1s_4s = 0;//hmm?
  // double Omega_1s_4d = 0;//hmm?
  // double Omega_1s_4f = 0;//hmm?

  if(T > 5000. && T < 60000.) //These coefficients good between 5000 and 72000K
    {
      Omega_1s_2p = 0.3435 + 1.297e-5*T + 2.178e-12*T2 + 7.928e-17*T3;
      Omega_1s_3s = 0.0625 - 1.299e-6*T + 2.666e-11*T2 - 1.596e-16*T3;
      // Omega_1s_3p = 0.0994 - 3.714e-7*T + 6.134e-11*T2 - 3.973e-16*T3;
      Omega_1s_3d = 0.0503 + 7.514e-7*T - 2.826e-13*T2 - 1.098e-17*T3;
      // Omega_1s_4p = 1.527e-3+1.001e-6*T - 2.192e-12*T2 + 9.348e-18*T3;
      // Omega_1s_4s = 0; //Can add these coefficients later for more precision
      // Omega_1s_4d = 0;
      // Omega_1s_4f = 0;
    }
  else if(T >= 60000. && T < 500000.)
    {
      Omega_1s_2p = 0.3162 + 1.472e-5*T - 8.275e-12*T2 - 8.794e-19*T3;
      Omega_1s_3s = 0.03337 + 2.223e-7*T - 2.794e-13*T2 + 1.516e-19*T3;
      // Omega_1s_3p = 0.06985 + 2.538e-6*T - 8.729e-13*T2 - 1.291e-18*T3;
      Omega_1s_3d = 0.05051 + 7.876e-7*T - 2.072e-12*T2 + 1.902e-18*T3;
      // Omega_1s_4p = 1.958e-3+9.525e-7*T - 9.668e-13*T2 + 4.807e-19*T3;
      // Omega_1s_4s = 0;
      // Omega_1s_4d = 0;
      // Omega_1s_4f = 0;
    }
  else //At low temperatures the rates are small, and T should not get so hot
    {
      Omega_1s_2p = 0.;
      Omega_1s_3s = 0.;
      // Omega_1s_3p = 0.;
      Omega_1s_3d = 0.;
      // Omega_1s_4p = 0.;
      // Omega_1s_4s = 0;
      // Omega_1s_4d = 0;
      // Omega_1s_4f = 0;
    }
  double prefix = 8.629e-6/sqrt(T);
  double q_1s_2p = prefix * Omega_1s_2p / 2.0 * exp(-E_1s_2p/(k_B*T));
  double q_1s_3s = prefix * Omega_1s_3s / 2.0 * exp(-E_1s_3s/(k_B*T));
  // double q_1s_3p = prefix * Omega_1s_3p / 2.0 * exp(-E_1s_3p/(k_B*T)); //remove?
  double q_1s_3d = prefix * Omega_1s_3d / 2.0 * exp(-E_1s_3d/(k_B*T));
  // double q_1s_4p = prefix * Omega_1s_4p / 2.0 * exp(-E_1s_4p/(k_B*T));  //remove?
  return q_1s_2p + q_1s_3s + q_1s_3d; //eq. 5 of Cantalupo+2008.  Can be improved.
  // return q_1s_2p + q_1s_3s + q_1s_3d + q_1s_3p +q_1s_4p; //eq. 5 of Cantalupo+2008.  Can be improved.

}

double coll_ex_rate_He1(double T) {
	return 9.10e-27*pow(T, -0.1687)*pow((1 + pow(T/1e5, 0.5)), -1)*exp(-13179./T);
}

double coll_ex_rate_He2(double T) {
	return 5.54e-17*pow(T, -0.397)*pow((1 + pow(T/1e5, 0.5)), -1)*exp(-473638./T);
}

//compton heating/cooling rate from Hui & Gnedin 1996
double compton_rate(double z) {
	return 5.6072e-36*pow(1 + z, 4.);
}
