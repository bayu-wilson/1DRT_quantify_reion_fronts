#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H


#include "user_inputs.h"
#include "global_constants.h"


using user_inputs::N_e;
using user_inputs::N_x;
using user_inputs::N_r;
using user_inputs::N_nu;
using user_inputs::nu_min;
using user_inputs::nu_max;
using user_inputs::t_max;
using user_inputs::R;
using user_inputs::R0;
using user_inputs::N_output;
using g_constants::yr_to_s;


void set_dlognu();

extern int prev_index;

//default differentials
extern double dr        ; //in kpc
extern double dlog_nu   ;
extern double dt        ;

//index rear IF
extern int index_rear_IF ;

//time
extern double t	    	;
extern double t_step  ;
extern int step	      ;

//grid variables
extern double r           [N_r]; //radial coordinate
extern double delta_r     [N_r - 1]; //radial differential
extern double nu          [N_nu]; //frequency
extern double delta_nu    [N_nu - 1]; //frequency differential

//hydrodynamics
extern double rho_total   [N_r]; //mass density
extern double rho_half    [N_r]; //half-time step
extern double rho_prev	   [N_r]; //previous time-step
extern double vel		   [N_r]; //gas velocity
extern double vel_prev    [N_r];
extern double p_total	   [N_r]; //pressure
extern double p_prev      [N_r];
extern double e_total     [N_r]; //specific energy density
extern double e_prev	   [N_r];
extern double phi		   [N_r]; //graviational potential
extern double phi_prev	   [N_r];
extern double cs		   [N_r]; //sound speed
//time derivatives
extern double drho_dt     [N_r];
extern double dvel_dt	   [N_r];
extern double de_dt       [N_r];
extern double dphi_dt     [N_r];

extern double rhov		   [N_r];
extern double rhov_prev   [N_r];
extern double Frhov	   [N_r];
extern double Fe		   [N_r];

//number densities
extern double nH          [N_r]; //number density of hydrogen
extern double nH_prev	   [N_r]; //previous time step
extern double nHe         [N_r]; //number density of helium
extern double nHe_prev    [N_r]; //previous time step
extern double nH1         [N_r]; //number density of neutral hydrogen
extern double nH1_prev    [N_r]; //previous time step
extern double nHe1        [N_r]; //number denstiy of neutral helium
extern double nHe1_prev   [N_r]; //previous time step
extern double nH2         [N_r]; //number density of ionized hydrogen
extern double nH2_prev    [N_r]; //previous time step
extern double nHe2        [N_r]; //number density of singly ionized helium
extern double nHe2_prev   [N_r]; //previous time step
extern double nHe3        [N_r]; //number density of doubly ionized helium
extern double nHe3_prev   [N_r]; //previous time step
extern double ne          [N_r]; //number density of free electrons
extern double ne_prev	   [N_r]; // previous time step
extern double n_tot       [N_r]; //total number density
extern double n_tot_prev  [N_r]; //previous time step

//fluxes
extern double FnH		   [N_r];
extern double FnHe        [N_r];
extern double FnH1        [N_r];
extern double FnH2        [N_r];
extern double FnHe1       [N_r];
extern double FnHe2       [N_r];
extern double FnHe3       [N_r];

//time derivatives
extern double dnH1_dt	   [N_r];
extern double dnH2_dt     [N_r];
extern double dnHe1_dt	   [N_r];
extern double dnHe2_dt    [N_r];
extern double dnHe3_dt    [N_r];
extern double dn_dt 	   [N_r];
extern double dne_dt	   [N_r];
extern double dT_dt	   [N_r];

//abundances
extern double f_H1        [N_r]; //HI fraction
extern double f_H1_prev   [N_r]; //TO CALCULATE IF SPEED BAYU: 04/13/23
extern double f_H1_step   [N_r];
extern double f_He1       [N_r]; //HeI fraction
extern double f_He1_step  [N_r];
extern double f_H2        [N_r]; //HII fraction
extern double f_H2_step   [N_r];
extern double f_He2       [N_r]; //HeII fraction
extern double f_He2_step  [N_r];
extern double f_He3       [N_r]; //HeIII fraction
extern double f_He3_step  [N_r];

//thermal parameters
extern double temp        [N_r]; //temperature
extern double temp_prev   [N_r]; //previous time step

//photoionization rates
extern double gamma_H1_tot  [N_r]; //photoionization rate of HI
extern double gamma_He1_tot [N_r]; //photoionization rate of HeI
extern double gamma_He2_tot [N_r]; //photionization rate of He2
extern double gamma_H1_prev [N_r]; //previous time step (HI)
extern double gamma_He1_prev[N_r]; //previous time step (He1)
extern double gamma_He2_prev[N_r]; //previous time step (HeII)

//secondary photoionization datas
extern double x_int_dom	[N_x][N_e + 1];
extern double f_ion		[N_x][N_e];
extern double f_heat		[N_x][N_e];
extern double f_exc		[N_x][N_e];
extern double n_lya		[N_x][N_e];
extern double n_ion_h1		[N_x][N_e];
extern double n_ion_he1	[N_x][N_e];
extern double n_ion_he2	[N_x][N_e];
extern double shull		[N_x][N_e];

extern double frac_ion     [N_r][N_nu];
extern double n_h1         [N_r][N_nu];
extern double n_he1        [N_r][N_nu];
extern double n_he2        [N_r][N_nu];
extern double frac_heat    [N_r][N_nu];

//recombination coefficients
extern double recomb_H2 	  [N_r];
extern double recomb_He2	  [N_r];
extern double recomb_He3	  [N_r];
extern double recomb_H2_prev [N_r];
extern double recomb_He2_prev[N_r];
extern double recomb_He3_prev[N_r];

//collisional ionization coefficients
extern double ci_H1		[N_r];
extern double ci_He1		[N_r];
extern double ci_He2		[N_r];

//heating/cooling rates
extern double heat_rate    [N_r]; //total heating rate
extern double cool_rate    [N_r]; //total cooling rate

//radiation energy
extern double u_nu        [N_r][N_nu];  //specific energy density of radiation

//rt equation
extern double gamma_nu_H1 [N_r][N_nu]; //spectral absorption coefficient of HI
extern double gamma_nu_He1[N_r][N_nu]; //spectral absorption coefficient of HeI
extern double gamma_nu_He2[N_r][N_nu]; //spectral absorption coefficient of HeII
extern double gamma_nu_tot[N_r][N_nu]; //total spectral absorption coefficient
extern double I_nu_prev   [N_r][N_nu]; //previous timestep
extern double I_nu        [N_r][N_nu]; //specific intensity
extern double I_nu_initial[N_nu];      //initial intensity of source cell 
//Lya emissivity from collisional excitations
extern double j_lya [N_r];
extern double q_eff_lya[N_r];

#endif
