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
using g_constants::h;
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
  inline constexpr bool input_grid       { FALSE }; //User-defined input grid? 
  inline constexpr bool input_source     { FALSE }; //User-defined source spectrum?
  inline constexpr bool save_initial_gas { TRUE }; //save initial grid

  inline constexpr bool recomb_cool       { TRUE }; //Recombination cooling fLag
  inline constexpr bool coll_ion          { TRUE }; //Collisional ionization flag
  inline constexpr bool coll_exc_cool	  { TRUE }; //Collisional excitation cooling flag
  inline constexpr bool sec_ion           { FALSE }; //secondary ionizations flag
  inline constexpr bool compton           { TRUE }; //compton heating/cooling flag
  inline constexpr bool temp_ev           { TRUE }; //temperature evolution flag
  inline constexpr bool hydro		  { FALSE }; //include hydrodynamics flag
  inline constexpr bool smooth_field      { FALSE }; //new BAYU July 13 2022
  inline constexpr bool correct_hardening { TRUE }; //new BAYU July 15 2022 

  //Input files if input_grid=True (user-defined) 
  inline constexpr char skewer[] {"0001"} ;
  inline constexpr char grid_input[] {"input_files/hydro_skewers/spec_xHeII1_007_ls_line0001.dat"};
  inline constexpr double R_start         { 50.}; //kpc //choosing a starting point of skewer chunk
  inline constexpr char source_spectrum[] { "input_files/spectrum_hydro_10myr.txt" };
  inline constexpr double sigma_gauss { 10. }; //pkpc, gaussian smoothing on inputed skewers if smooth_field=TRUE

  //Uniform density skewers when input_grid=False
  inline constexpr double scale_density_factor { 8. }; 

  //output files
  inline constexpr char gas_output[] {"output_files/gas_test_hydro_61.2myr.txt"};
  inline constexpr char initial_gas_output[]   { "output_files/gas_test_hydro_000myr.txt"};
  inline constexpr char otf_dir[] {"output_files/gasprops/ud_a=2.500_Lum=1.50e+46/"};

  //time stepping
  inline constexpr double tfrac	       { 0.1 }; //Maximum fraction by which any quantity can change in one time step
  inline constexpr double min_dt_early { 1.e-8 }; //Minimum time step before 1 Myr
  inline constexpr double min_dt_late  { 1.e-6 }; //Minimum time step after 1 Myr
  inline constexpr double min_frac     { 0. }; //Ignore fractional change constraint if f_i < min_frac

  inline constexpr int N_r { 2000 };
  inline constexpr int N_nu   { 50 }; //Number of frequency bins 50
  inline constexpr double t_max { 61.2 }; //Runtime (in Myr)
  inline constexpr double z        { 5.7 }; //redshift
  inline constexpr double R0       { 1.0e+03 }; //minimum radius of the grids, keep large to minimize geometric attenuation 
  inline constexpr double R        { R0+1219 };  //Total radius for the spatial grid (in pkpc) 
  inline constexpr double nu_min   { 13.6/h_eV }; //Minimum frequency (in Hz)
  inline constexpr double nu_max   { 13.6/h_eV*4 }; //Maximum frequency (in Hz)
  inline constexpr double temp_0   { 5.e1 }; //Initial gas temperature in K (if not input by user) 1e2 2.e4
  inline constexpr double T_source { 1.e5 }; //Source temperature (if blackbody)
  inline constexpr double alpha { 2.500 }; //spectral power law index if spectrum is a power law
  inline constexpr double Lum {  1.50e+46 }; 
  inline constexpr double frontIF_fHI { 0.99 };
  inline constexpr double backIF_fHI { 0.01 }; //where F_inc is defined. see io_funcs.cc:write_otf_spectrum and data_funcs.cc:calc_ifront_flux_method 
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
  inline constexpr int N_output      { 500 };

  //folders
  inline constexpr char x_int_tables[] = "x_int_tables/";

  //variables for integration?? also referenced in secondary ionizations i think.
  inline constexpr int N_x { 14 };
  inline constexpr int N_e { 258 };
}

#endif

