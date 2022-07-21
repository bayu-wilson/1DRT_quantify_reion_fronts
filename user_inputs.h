#ifndef USER_INPUTS_H
#define USER_INPUTS_H

#include "global_constants.h"
#include <string>

using g_constants:: TRUE;
using g_constants:: FALSE;
using g_constants::h_eV; // do it like this so namespace is less polluted. That can lead to neaming collisions.
using g_constants::k_B;
using g_constants::m_H;
using g_constants::Y;
using g_constants::Omega_b;
using g_constants::pi;
using g_constants::G_grav;
using g_constants::mpc_to_km;
using g_constants::h_const;

namespace user_inputs
{
  inline constexpr bool parallel     { TRUE }; //parallelization with open-mp

  inline constexpr char spectrum[] { "POWER" }; // "MONOCHROMATIC", "BLACKBODY", or "POWER".
  inline constexpr char rec_case[] { "A" }; //Use case "A" or "B" recombination coefficients

  inline constexpr bool FINITE_C         { FALSE }; //Finite speed of light?
  inline constexpr bool spherical        { TRUE }; //Spherical symmetry?  If FALSE use plane-parallel
  inline constexpr bool input_grid       { TRUE }; //User-defined input grid?
  inline constexpr bool input_source     { FALSE }; //User-defined source spectrum?
  inline constexpr bool save_initial_gas { TRUE }; //save initial grid

  inline constexpr bool recomb_cool       { TRUE }; //Recombination cooling fLag
  inline constexpr bool coll_ion          { TRUE }; //Collisional ionization flag
  inline constexpr bool coll_exc_cool	    { TRUE }; //Collisional excitation cooling flag
  inline constexpr bool sec_ion           { FALSE }; //secondary ionizations flag
  inline constexpr bool compton           { TRUE }; //compton heating/cooling flag
  inline constexpr bool temp_ev           { TRUE }; //temperature evolution flag
  inline constexpr bool hydro			        { FALSE }; //include hydrodynamics flag
  inline constexpr bool smooth_field      { TRUE }; //new BAYU July 13 2022
  inline constexpr bool correct_hardening { TRUE }; //new BAYU July 15 2022

  inline constexpr bool QSO_spec       { FALSE }; //CHRIS 05/16/22: output QSO spectrum?
  inline constexpr bool add_background { FALSE }; //CHRIS 05/17/22: add uniform background?
  inline constexpr bool LyC_opacity    { FALSE }; //Chris 06/16/22: make it an option to calculate LyC
  inline constexpr bool Lyn_series     { FALSE }; //Chris 06/14/22: include Lyn series in QSO spectrum?
  inline constexpr bool use_tau_GP     { FALSE }; //CHRIS 07/09/22: use Gunn-Petersen approximation to calculate lyman series
  inline constexpr int  NLyn	         {  39  }; //number of lyman series lines to use
  inline constexpr int  N_wav          { 1000 }; //number of wavelength bins
  inline constexpr double wav_start    { 1050.}; //start of the wavelength array
  inline constexpr double wav_end      { 800. }; //end of the wavelength array

  //input files (if user-defined)
  inline constexpr char grid_input[]      { "input_files/spec_xHeII1_003_ls_line0001.dat" };
  inline constexpr double R_start         { 5317 }; //kpc //choosing a starting point of skewer chunk
  inline constexpr char source_spectrum[] { "input_files/spectrum_hydro_10myr.txt" };

  inline constexpr double t_max { 200}; //Runtime (in Myr)
  //output files
  inline constexpr char gas_output[] { "output_files/gas_test_hydro_200myr.txt"};
  inline constexpr char initial_gas_output[]   { "output_files/gas_test_hydro_000myr.txt"};

  //grid sizes
  inline constexpr int N_r    { 408 }; //Number of spatial bins 10,000 //1360 //680... stepsize is 0.735 pkpc
  inline constexpr int N_nu   { 30 }; //Number of frequency bins 50

  //time stepping
  inline constexpr double tfrac	       { 0.1 }; //Maximum fraction by which any quantity can change in one time step
  inline constexpr double min_dt_early { 1e-8 }; //Minimum time step before 1 Myr
  inline constexpr double min_dt_late  { 1e-6 }; //Minimum time step after 1 Myr
  inline constexpr double min_frac     { 0. }; //Ignore fractional change constraint if f_i < min_frac

  inline constexpr double z        { 7.1 }; //redshift
  // inline constexpr double R0       { 1e4 }; //minimum radius of the RT, hydro, temperature and chemistry grids //10,000
  inline constexpr double R0       { 1e4 }; //minimum radius of the grids //keep large to minimize geometric attenuation //different than R_start
  inline constexpr double R        { R0+300. }; //Total radius for the spatial grid (in pkpc)
  inline constexpr double nu_min   { 13.6/h_eV }; //Minimum frequency (in Hz)
  inline constexpr double nu_max   { 13.6/h_eV*4 }; //Maximum frequency (in Hz)
  // inline constexpr double nu_max   { 9e3/h_eV }; //Maximum frequency
  inline constexpr double temp_0   { 5e2 }; //Initial gas temperature in K (if not input by user) 1e2 2.e4
  inline constexpr double T_source { 1.e5 }; //Source temperature (if blackbody)
  inline constexpr double alpha { 2.5 }; //spectral power law index if spectrum is a power law
  inline constexpr double Lum      { 4.e56*g_constants::h*nu_min*alpha/(alpha - 1) }; //for power law index L_MW = 1e43 erg/s
  inline constexpr double Gam0     { 1e-13 }; // uniform photoionization background
  inline constexpr double frontIF_fHI { 0.9 };
  inline constexpr double backIF_fHI { 1e-2 };
  inline constexpr double sigma_gauss { 3 }; //pkpc, sigma if field had gaussian smoothing applied
  inline constexpr double rho_crit_0 { 8.688e-30 }; // critical density at z=0 [g cm^-3] //3 Omega H0^2 / (8 pi // assume h=0.68
  inline constexpr double rho_0      { rho_crit_0*Omega_b*(1.+z)*(1.+z)*(1.+z) }; //initial baryonic density at redshift z
  inline constexpr double fH1_0      { 1. }; //Initial HI fraction
  inline constexpr double fHe1_0     { 1. }; //Initial HeI fraction
  inline constexpr double fH2_0      { 0. }; //Initial HII fraction
  inline constexpr double fHe2_0     { 0. }; //Initial HeII fraction
  inline constexpr double fHe3_0     { 0. }; //Initial HeIII fraction

  //leave these 0 for uniform initial abundances and density
  inline constexpr double rho_index  { 0 };
  inline constexpr double fH1_index  { 0 };
  inline constexpr double fHe1_index { 0 };
  inline constexpr double fH2_index  { 0 };
  inline constexpr double fHe2_index { 0 };
  inline constexpr double fHe3_index { 0 };

  //number of equally spaced time steps to output on-the-fly data.
  inline constexpr int N_output      { 250 };

  //folders
  inline constexpr char x_int_tables[] = "x_int_tables/";

  //variables for integration?? also referenced in secondary ionizations i think.
  inline constexpr int N_x { 14 };
  inline constexpr int N_e { 258 };
}

// inline constexpr char spec_output[] { "output_files/spectrum_hydro_100myr.txt"};
// inline constexpr char otf_output[] { "output_files/otf_test_hydro_100myr.txt"};
// inline constexpr char otf_output_bayu[] { "output_files/otf_bayu_100myr.txt"}; //added 21/12/10
// inline constexpr char otf_output_fred[] { "output_files/otf_fred_"}; //added 22/06/26
// inline constexpr double region_IF { 100 }; // region around the IF to check the smallest timescales. This is in units of kpc!!
// inline constexpr double Lum      = 4.e56*g_constants::h*nu_min*alpha/(alpha - 1); //for power law index L_MW = 1e43 erg/s
// inline constexpr double Lum      = 1e57*g_constants::h*nu_min*alpha/(alpha - 1); //for power law index
// inline constexpr double Lum      { 5e48*4.155*k_B*T_source}; //Luminosity in erg/s for blackbody
// inline constexpr double Lum      = 5e48*h*nu; //for monochromatic
// inline constexpr double rho_0      { 1e-3*m_H/(1 - Y) }; //uniform density (in cgs) //for Iliev 2009 test 1
// inline constexpr double rho_0      = Omega_b/Omega_m*rho_crit(1./(1. + z)); //For Fred Davies test

#endif
