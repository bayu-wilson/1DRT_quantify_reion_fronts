#ifndef GLOBAL_CONSTANTS_H
#define GLOBAL_CONSTANTS_H

namespace g_constants
{
  //true/false
  inline constexpr bool TRUE { 1 };
  inline constexpr bool FALSE { 0 };

  // numbers
  inline constexpr double pi { 3.14159 };
  inline constexpr double e  { 2.71828 };

  //Fundamental constants
  inline constexpr double h { 6.626068e-27 }; //Planck's constant (erg s)
  inline constexpr double h_eV { 4.135667e-15 }; // Planck's constant (eV s)
  inline constexpr double year { 31557600 }; // Year (s)
  inline constexpr double k_B { 1.38065e-16 }; //Boltzmann constant (erg/K)
  inline constexpr double kpc { 3.08568e21 }; // Kiloparsec (cm)
  inline constexpr double c { 2.998e10 }; // Speed of light (cm/s)
  inline constexpr double sigma_T { 6.6524e-25 }; // Thompson scattering cross-section (cm^2)
  inline constexpr double m_H { 1.672622e-24 }; //Mass of hydrogen atom (g)
  inline constexpr double m_e { 9.10939e-28 }; //Mass of electron (g)
  inline constexpr double gam { 1.666666667 }; // Adiabatic coeff.
  inline constexpr double e_c { 4.803e-10 }; //electron charge (esu)
  inline constexpr double stef_boltz { 5.6704e-5 }; //stefan-boltzman constant in cgs
  inline constexpr double G_grav { 6.67259e-8 }; // Gravitational constant in cgs units
  inline constexpr double lambda_lya_cm { 1.21567e-5 }; //in cm //added 01/12/2022
  inline constexpr double lambda_HI_ion_cm { 0.912e-5 }; //in cm //added 06/07/2022
  inline constexpr double chi_He { 0.08 }; // 1+chi factor accounts for single-ionized helium

  //Cosmological parameters
  inline constexpr double Y { 0.24 }; // Mass ratio of helium
  inline constexpr double Omega_m { 0.3 }; // Density of matter in the universe
  inline constexpr double Omega_l { 0.7 }; // Dark energy density
  inline constexpr double Omega_b { 0.048 }; // Baryon density
  inline constexpr double h_const { 0.68 }; // H0 { h_const * 100 km/s/Mpc}

  //Spectral parameters
  inline constexpr double sigma0 { 6.304e-18 }; //cm^2

  //Units
  inline constexpr double ryd_to_ev { 13.6 };
  inline constexpr double pc_to_cm { 3.086e18 }; //86 not 68
  inline constexpr double kpc_to_cm { 3.086e21 };
  inline constexpr double mpc_to_km { 3.086e19 };
  inline constexpr double yr_to_s { 3.154e7 };
  inline constexpr double ev_to_erg { 1.6022e-12 };
  inline constexpr double eV_cgs { 1.60218e-12 };
}

#endif
