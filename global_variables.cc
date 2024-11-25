//this file is rather unique. when main.cc is run, it includes this file so
//I believe the preprocessor compiles this before compiling main.
#include <math.h>
#include "global_variables.h"

//for the hardening correction
int equil_index{0};
int equil_index_prev{0};
int ifront_index{0};
int ifront_index_prev{0};

//default differentials
double dr     			  = (R-R0)/N_r; //in kpc // got rid of N_r-1 because I don't think it's necessary.
double dlog_nu 				= 0;
double dt				      = 1e-4*t_max*yr_to_s*1e6/(N_output - 1);

//index of radial coordinate where IF is 10% neutral
// int index_rear_IF {};

void set_dlognu()  {
	if (N_nu > 1)  {
		dlog_nu = (log10(nu_max) - log10(nu_min))/(N_nu - 1); //in Hz
	}
}

//time
double t		= 0;
double t_step	= 0;
int step		= 0;

//grid variables
double r           [N_r]{}; //radial coordinate
double delta_r     [N_r - 1]{}; //radial differential
double nu          [N_nu]{}; //frequency
double delta_nu    [N_nu - 1]{}; //frequency differential

//hydrodynamics
double rho_total   [N_r]{}; //mass density
double rho_half    [N_r]{}; //half-time step
double rho_prev	   [N_r]{}; //previous time-step
double vel		   [N_r]{}; //gas velocity
double vel_prev    [N_r]{};
double p_total	   [N_r]{}; //pressure
double p_prev      [N_r]{};
double e_total     [N_r]{}; //specific energy density
double e_prev	   [N_r]{};
double phi		   [N_r]{}; //graviational potential
double phi_prev	   [N_r]{};
double cs		   [N_r]{}; //sound speed
//time derivatives
double drho_dt     [N_r]{};
double dvel_dt	   [N_r]{};
double de_dt       [N_r]{};
double dphi_dt     [N_r]{};

double rhov		   [N_r]{};
double rhov_prev   [N_r]{};
double Frhov	   [N_r]{};
double Fe		   [N_r]{};

//number densities
double nH          [N_r]{}; //number density of hydrogen
double nH_prev	   [N_r]{}; //previous time step
double nHe         [N_r]{}; //number density of helium
double nHe_prev    [N_r]{}; //previous time step
double nH1         [N_r]{}; //number density of neutral hydrogen
double nH1_prev    [N_r]{}; //previous time step
double nHe1        [N_r]{}; //number denstiy of neutral helium
double nHe1_prev   [N_r]{}; //previous time step
double nH2         [N_r]{}; //number density of ionized hydrogen
double nH2_prev    [N_r]{}; //previous time step
double nHe2        [N_r]{}; //number density of singly ionized helium
double nHe2_prev   [N_r]{}; //previous time step
double nHe3        [N_r]{}; //number density of doubly ionized helium
double nHe3_prev   [N_r]{}; //previous time step
double ne          [N_r]{}; //number density of free electrons
double ne_prev	   [N_r]{}; // previous time step
double n_tot       [N_r]{}; //total number density
double n_tot_prev  [N_r]{}; //previous time step

//fluxes
double FnH		   [N_r]{};
double FnHe        [N_r]{};
double FnH1        [N_r]{};
double FnH2        [N_r]{};
double FnHe1       [N_r]{};
double FnHe2       [N_r]{};
double FnHe3       [N_r]{};

//time derivatives
double dnH1_dt	   [N_r]{};
double dnH2_dt     [N_r]{};
double dnHe1_dt	   [N_r]{};
double dnHe2_dt    [N_r]{};
double dnHe3_dt    [N_r]{};
double dn_dt 	   [N_r]{};
double dne_dt	   [N_r]{};
double dT_dt	   [N_r]{};

//abundances
double f_H1        [N_r]{}; //HI fraction
double f_H1_prev   [N_r]{}; //TO CALCULATE IF SPEED BAYU: 04/13/23
double f_H1_step   [N_r]{};
double f_He1       [N_r]{}; //HeI fraction
double f_He1_step  [N_r]{};
double f_H2        [N_r]{}; //HII fraction
double f_H2_step   [N_r]{};
double f_He2       [N_r]{}; //HeII fraction
double f_He2_step  [N_r]{};
double f_He3       [N_r]{}; //HeIII fraction
double f_He3_step  [N_r]{};

//thermal parameters
double temp        [N_r]{}; //temperature
double temp_prev   [N_r]{}; //previous time step

//photoionization rates
double gamma_H1_tot  [N_r]{}; //photoionization rate of HI
double gamma_He1_tot [N_r]{}; //photoionization rate of HeI
double gamma_He2_tot [N_r]{}; //photionization rate of He2
double gamma_H1_prev [N_r]{}; //previous time step (HI)
double gamma_He1_prev[N_r]{}; //previous time step (He1)
double gamma_He2_prev[N_r]{}; //previous time step (HeII)

//secondary photoionization datas
double x_int_dom	[N_x][N_e + 1]{};
double f_ion		[N_x][N_e]{};
double f_heat		[N_x][N_e]{};
double f_exc		[N_x][N_e]{};
double n_lya		[N_x][N_e]{};
double n_ion_h1		[N_x][N_e]{};
double n_ion_he1	[N_x][N_e]{};
double n_ion_he2	[N_x][N_e]{};
double shull		[N_x][N_e]{};

double frac_ion     [N_r][N_nu]{};
double n_h1         [N_r][N_nu]{};
double n_he1        [N_r][N_nu]{};
double n_he2        [N_r][N_nu]{};
double frac_heat    [N_r][N_nu]{};

//recombination coefficients
double recomb_H2 	  [N_r]{};
double recomb_He2	  [N_r]{};
double recomb_He3	  [N_r]{};
double recomb_H2_prev [N_r]{};
double recomb_He2_prev[N_r]{};
double recomb_He3_prev[N_r]{};

//collisional ionization coefficients
double ci_H1		[N_r]{};
double ci_He1		[N_r]{};
double ci_He2		[N_r]{};

//heating/cooling rates
double heat_rate    [N_r]{}; //total heating rate
double cool_rate    [N_r]{}; //total cooling rate

//radiation energy
double u_nu        [N_r][N_nu]{};  //specific energy density of radiation

//rt equation
double gamma_nu_H1 [N_r][N_nu]{}; //spectral absorption coefficient of HI
double gamma_nu_He1[N_r][N_nu]{}; //spectral absorption coefficient of HeI
double gamma_nu_He2[N_r][N_nu]{}; //spectral absorption coefficient of HeII
double gamma_nu_tot[N_r][N_nu]{}; //total spectral absorption coefficient
double I_nu_prev   [N_r][N_nu]{}; //previous timestep
double I_nu        [N_r][N_nu]{}; //specific intensity
double I_nu_initial[N_nu]{};      //initial intensity of source cell

//Lya emissivity from collisional excitations
double j_lya [N_r]{};
double q_eff_lya[N_r]{};

//Measuring alpha for hardening
//double I_nu_avg_sk [N_r][N_nu]{};
//double alpha_meas [N_r]{};
