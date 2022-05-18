#ifndef RATES_H
#define RATES_H

namespace rates
{
  //cross-section parameters from Verner et. al. 1996, Table 1
  inline constexpr double sigma0_H1  { 5.475e-14 };
  inline constexpr double sigma0_He1 { 9.492e-16 };
  inline constexpr double sigma0_He2 { 1.369e-14 };

  inline constexpr double Eth_H1     { 1.360e+1 };
  inline constexpr double Eth_He1    { 2.459e+1 };
  inline constexpr double Eth_He2    { 5.442e+1 };

  inline constexpr double Emax_H1    { 5e+4 };
  inline constexpr double Emax_He1   { 5e+4 };
  inline constexpr double Emax_He2   { 5e+4 };

  inline constexpr double E0_H1      { 4.298e-1 };
  inline constexpr double E0_He1     { 1.361e+1 };
  inline constexpr double E0_He2     { 1.720e0 };

  inline constexpr double y0_H1      { 0 };
  inline constexpr double y0_He1     { 4.434e-1 };
  inline constexpr double y0_He2     { 0 };

  inline constexpr double y1_H1      { 0 };
  inline constexpr double y1_He1     { 2.136 };
  inline constexpr double y1_He2     { 0 };

  inline constexpr double yw_H1      { 0 };
  inline constexpr double yw_He1     { 2.039 };
  inline constexpr double yw_He2     { 0 };

  inline constexpr double ya_H1      { 3.288e+1 };
  inline constexpr double ya_He1     { 1.469 };
  inline constexpr double ya_He2     { 3.288e+1 };

  inline constexpr double P_H1       { 2.963 };
  inline constexpr double P_He1      { 3.188 };
  inline constexpr double P_He2      { 2.963 };

  //for the new collisional excitation function: coll_ex_rate_H1_acc
  inline constexpr double hvHI     { 2.178717e-11 }; //in erg
  inline constexpr double E_1s_2p  { 0.75*hvHI };
  inline constexpr double E_1s_3s  { (1.0-(1.0/9.0))*hvHI };
  inline constexpr double E_1s_3d  { (1.0-(1.0/9.0))*hvHI };
  inline constexpr double E_1s_3p  { (1.0-(1.0/9.0))*hvHI };
  inline constexpr double E_1s_4s  { (1.0-(1.0/16.0))*hvHI };
  inline constexpr double E_1s_4p  { (1.0-(1.0/16.0))*hvHI };
  inline constexpr double E_1s_4d  { (1.0-(1.0/16.0))*hvHI };
  inline constexpr double E_1s_4f  { (1.0-(1.0/16.0))*hvHI };
}

double lambda(double T, double Eth);
double F_of_y(double nu, double Eth, double Emax, double E0, double y0, double y1, double yw, double ya, double P);

//photoionization cross-sections
//photo-ionization cross-section of HI from Verner et. al. 1996
double sigmapi_H1(double nu);
//photo-ionization cross-section of HeI from Verner et. al. 1996
double sigmapi_He1(double nu);
//photo-ionization cross-section of HeII from Verner et. al. 1996
double sigmapi_He2(double nu);

//collisional ionization rates
//collisional ionization rate of HI from Hui & Gnedin 1996
double cic_H1(double T);
//collisional ionization rate of HeI from Hui & Gnedin 1996
double cic_He1(double T);
//collisional ionization rate of HeII from Hui & Gnedin 1996
double cic_He2(double T);

//collisional ionization cooling rates
//collisional ionization cooling rate of HI from Hui & Gnedin 1996
double cicr_H1(double T);
//collisional ionization cooling rate of HeI from Hui & Gnedin 1996
double cicr_He1(double T);
//collisional ionization cooling rate of HeII from Hui & Gnedin 1996
double cicr_He2(double T);

//recombination rates
//Case A recombination coefficient of HII from Hui & Gnedin 1996
double alphaA_H2(double T);
//Case B recombination coefficient of HII from Hui & Gnedin 1996
double alphaB_H2(double T);
//Case A recombination coefficient of HeII from Hui & Gnedin 1996
double alphaA_He2(double T);
//Case B recombination coefficient of HeII from Hui & Gnedin 1996
double alphaB_He2(double T);
//Case A recombination coefficient of HeIII from Hui & Gnedin 1996
double alphaA_He3(double T);
//Case B recombination coefficient of HeIII from Hui & Gnedin 1996
double alphaB_He3(double T);
//Dielectric recombination coefficient of He II from Hui & Gnedin 1996
double Dalpha_He2(double T);

//recombination cooling rates
//Case A recombination cooling rate for HII from Hui & Gnedin 1996
double rcrA_H2(double T);
//Case B recombination cooling rate for HII from Hui & Gnedin 1996
double rcrB_H2(double T);
//Case A recombination cooling rate for HeII from Hui & Gnedin 1996
double rcrA_He2(double T);
//Case B recombination cooling rate for HeII from Hui & Gnedin 1996
double rcrB_He2(double T);
//Case A recombination cooling rate for HeIII from Hui & Gnedin 1996
double rcrA_He3(double T);
//Case B recombination cooling rate for HeIII from Hui & Gnedin 1996
double rcrB_He3(double T);
//dielectric recombination cooling (He II).
double drcr_He2(double T);
//collisional excitation cooling rates from Cen 1992
double coll_ex_rate_H1(double T);
////Uses polynomial fits from Giovanardi+1987, also used in Cantalupo+2008
double coll_ex_rate_H1_acc(double T);
double coll_ex_rate_He1(double T);
double coll_ex_rate_He2(double T);
//compton heating/cooling rate from Hui & Gnedin 1996
double compton_rate(double z);

#endif
