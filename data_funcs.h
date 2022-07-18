#ifndef DATA_FUNCS_H
#define DATA_FUNCS_H

//Volume-weighted average
double calc_vol_avg(double f[], int n);
//Mass-weighted average
double calc_mass_avg(double f[], double rho[], double m, int n);
//Calculate ionization front position and velocity
double* calc_ifront(double x[], double xprev[], double radius[], int n);
double* calc_ifront_flux_method();
double* tau_bf(double wav_em[], double NHI[]);
void calc_mock_QSO_spec();
#endif
