#ifndef DATA_FUNCS_H
#define DATA_FUNCS_H

//Volume-weighted average
double calc_vol_avg(double f[], int n);
//Mass-weighted average
double calc_mass_avg(double f[], double rho[], double m, int n);
//Calculate ionization front position and velocity
double* calc_ifront(double x[], double xprev[], double radius[], int n);
double* calc_ifront_avg(double x[], double xprev[], double radius[], int n);
double* calc_ifront_FlexRT();
double* calc_ifront_flux_method();
//Cakculate lya emissivity
void get_j_lya();
double sigmaHI_numeric(double nu_iter);
double avg_sigma_analytic(double alpha);
double f_alpha(double alpha, double numerical_soln);
double df_dalpha(double alpha);
double pred_alpha(double alpha_estimate, double numerical_soln);
double solve_spectral_index(int r_iter);

#endif
