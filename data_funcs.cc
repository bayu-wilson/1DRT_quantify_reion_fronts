#include <math.h>
#include "global_variables.h"
#include "general_funcs.h"
#include "user_inputs.h"
#include "global_constants.h"
#include "data_funcs.h"
#include "io_funcs.h"

using g_constants::yr_to_s;
using user_inputs::N_output;
using user_inputs::t_max;
using g_constants::pi;
using g_constants::c;
using g_constants::chi_He;
using g_constants::lambda_lya_cm;
using g_constants::h;

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

	double f_H1_IF_rear = 0.1; //neutral hydrogen fraction at the rear of the IF
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
				index_rear_IF = i; //initialized in global_variables.cc. We track it in order to track the
													 // spectral evolution in io_funcs.cc: write_spectrum()
			}
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
	*(flux_vel + 2) = nH_boundary;
 	return flux_vel;
}
