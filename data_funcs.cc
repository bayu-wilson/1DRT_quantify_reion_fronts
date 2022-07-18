#include <math.h>
#include <fstream> //CHRIS: 05/16/22
#include <string> //CHRIS: 05/16/22
#include "global_variables.h"
#include "general_funcs.h"
#include "user_inputs.h"
#include "global_constants.h"
#include "data_funcs.h"
#include "io_funcs.h"
#include "cosmo_funcs.h"
#include "rates.h"
using namespace std; //CHRIS: 05/16/22

using g_constants::yr_to_s;
using user_inputs::N_output;
using user_inputs::t_max;
using user_inputs::z; //CHRIS: 05/16/22
using user_inputs::gas_output; //CHRIS: 05/16/22
using g_constants::pi;
using g_constants::c;
using g_constants::chi_He;
using g_constants::lambda_lya_cm;
using g_constants::h;
using g_constants::kpc; //CHRIS: 05/16/22
using g_constants::lambda_HI_ion_cm;

//Volume-weighted average
double calc_vol_avg(double f[], int n)
{
	double tot = 0;
	for (int i{ 1 }; i < n; i++)  {
		tot += f[i]*(pow(r[i],3) - pow(r[i-1],3))/pow(r[N_r-1],3);
	}
	return tot;
}
//Mass-weighted average
double calc_mass_avg(double f[], double rho[], double m, int n)
{
	double tot    = 0;
	double m_tot  = 0;
	for (int i{ 1 }; i < n; i++)
	{
		m_tot += m*rho[i]*4.*pi/3.*(pow(r[i],3) - pow(r[i-1],3));
		tot   += f[i]*m*rho[i]*4.*pi/3.*(pow(r[i],3) - pow(r[i-1],3));
  }
  return tot/m_tot;
}

//Calculate ionization front position and velocity
double* calc_ifront(double x[], double xprev[], double radius[], int n)
{
	double pos, pos_prev, vel_IF;
	double *ifront = (double*) malloc(sizeof(double) * 2);
	pos       = interpolate(x, radius, 0.5, n);
	pos_prev  = interpolate(xprev, radius, 0.5, n);
	vel_IF       = (pos - pos_prev)/(t_max*yr_to_s*1e6/(N_output - 1));
	*(ifront + 0) = pos;
	*(ifront + 1) = vel_IF;
	return ifront;
}

//Calculate ionization front incident photon flux, velocity, and neutral hydrogen number density
double* calc_ifront_flux_method()
{
	double I_div_nu_pair[2][N_nu]{};
	double inc_photon_flux_pair[2]{};
	double inc_photon_flux{};
	double f_H1_pair[2];
	double nH_boundary{};
	double vel_IF{};
	double *flux_vel = (double*) malloc(sizeof(double) * 3);
	// int index_rear_IF{};//DELETE

	double f_H1_IF_rear = 1e-1; //neutral hydrogen fraction at the rear of the IF
	for (int i=0;i<(N_r-1);i++) //find index of array at the location closest to 10% neutral
	{
		if ((f_H1[i]<f_H1_IF_rear) && (f_H1[i+1]>f_H1_IF_rear))
		{
			for (int j=0;j<N_nu;j++)
			{
				I_div_nu_pair[0][j] = I_nu[i][j]/nu[j]; //we track I divided by nu in order to do the integration in the next section
				I_div_nu_pair[1][j] = I_nu[i+1][j]/nu[j]; // int (I_nu / h / nu dnu ) = incident Flux ionizing photons
				f_H1_pair[0] = f_H1[i];
				f_H1_pair[1] = f_H1[i+1];
			}
			// index_rear_IF = i; //initialized in global_variables.cc. We track it in order to track the
												 // spectral evolution in io_funcs.cc: write_spectrum()
			break;
		}
	}
	inc_photon_flux_pair[0] = trapz_int(I_div_nu_pair[0],nu,N_nu)*4*pi/h;
	inc_photon_flux_pair[1] = trapz_int(I_div_nu_pair[1],nu,N_nu)*4*pi/h;
	inc_photon_flux = interpolate(f_H1_pair,inc_photon_flux_pair,f_H1_IF_rear,2);
	nH_boundary    = interpolate(f_H1, nH, 0.5, N_r); // number density neutral H at front boundary 50%
	vel_IF = (c*inc_photon_flux) / (inc_photon_flux + c*nH_boundary*(1+chi_He));
	*(flux_vel + 0) = inc_photon_flux;
	*(flux_vel + 1) = vel_IF;
	// *(flux_vel + 2) = nH[index_rear_IF];
	*(flux_vel + 2) = nH_boundary;
 	return flux_vel;
}

//CHRIS: 05/16/22
//Subroutine to calculate a mock QSO spectrum starting at r = R0
//Optional flag adds a neutral island along the sightline with a specified central position and width

double* tau_bf(double wav_em[], double NHI[])
{
	int sigma_ind = 3.;
	double lymanlimit = 911.267;
	double *tau_bf_out = (double*) malloc(sizeof(double) * N_r);

	*(tau_bf_out + 0) = 0.;
	for (int i=1; i<N_r; i++)
	{
		*(tau_bf_out + i) = 0.;
		double lam1 = wav_em[i];
		for (int j=0; j<i; j++)
		{
			double dlam = wav_em[i]*(wav_em[j]/wav_em[j+1] - 1.);
			*(tau_bf_out + i) += g_constants::sigma0*NHI[j]*pow(lam1 + dlam/2., sigma_ind)/pow(lymanlimit, sigma_ind);
			lam1 += dlam;
		}
	}
	return tau_bf_out;
}

void calc_mock_QSO_spec()
{
	double delta_v     [N_r - 1]{}; //velocity differential
	double wav         [N_r]{};
	double flux        [N_r]{};
	double NHI         [N_r-1]{};
	double a = 1./(1.+user_inputs::z);
	double lymanlimit = 911.76;

	wav[0] = lymanlimit;
	for (int i = 1; i<N_r; i++)
	{
		delta_v[i-1]   = (r[i] - r[i-1])*H(a);
		wav[i]         = wav[i-1]/(1.+delta_v[i-1]/g_constants::c);
		NHI[i-1]       = delta_r[i-1]*(nH1[i] + nH1[i-1])/2.;
	}

	double *tau_912 = tau_bf(wav, NHI);

	for (int i = 0; i<N_r; i++)
	{
		double tau = *(tau_912 + i);
		flux[i] = exp(-tau);
	}

	ofstream file;
	string s1  = "./output_files/QSO_spec_z=";
	string z_string = to_string(z);
	s1 = s1 + z_string;
	const char *spec_output = s1.c_str();
	file.open(spec_output, ios::out | ios::binary);

	for (int i=0; i<N_r; i++)
	{
		file.write((char*)&wav[i],  sizeof(double));
		file.write((char*)&flux[i], sizeof(double));
	}
	file.close();
}

void get_j_lya() {
	for (int i{ 0 }; i < N_r; i++) {
		j_lya[i]=coll_ex_rate_H1_acc(temp[i])*nH1[i]*ne[i]*h*c/lambda_HI_ion_cm/4/pi;
	}
}
