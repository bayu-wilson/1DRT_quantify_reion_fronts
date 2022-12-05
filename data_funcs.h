#ifndef DATA_FUNCS_H
#define DATA_FUNCS_H

//using user_inputs::N_wav;

//Volume-weighted average
double calc_vol_avg(double f[], int n);
//Mass-weighted average
double calc_mass_avg(double f[], double rho[], double m, int n);
//Calculate ionization front position and velocity
double* calc_ifront(double x[], double xprev[], double radius[], int n);
double* calc_ifront_flux_method();
//Chris 7/18/22
//double sigma_HI_LyC(double lam); //Chris 7/18/22
//double* get_lyn_line_data(int series); //Chris 7/18/22
//double sigma_HI_Lyn(double lam, double T, int series, double NHI, double dv, double dlam, double lam_em, double scale); //Chris 7/18/22
//double tau_HI_Lyn_GP(int series, double a_init, double wav, double scale[N_r], double nHI[N_r]); //Chris 7/18/22
//double* tau_bf(double wav_em[N_wav], double v[N_r], double NHI[N_r], double scale[N_r], double nHI[N_r]); //Chris 7/18/22
//void calc_mock_QSO_spec(double is_length, double is_pos, int name_ind); //Chris 7/18/22
// double* tau_bf(double wav_em[], double NHI[]);
// void calc_mock_QSO_spec();
//Cakculate lya emissivity
void get_j_lya();
#endif
