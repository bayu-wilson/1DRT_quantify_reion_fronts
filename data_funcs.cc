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
double calc_vol_avg(double f[], int n)  {

	// int i;
	double tot = 0;

	for (int i{ 1 }; i < n; i++)  {
		tot += f[i]*(pow(r[i],3) - pow(r[i-1],3))/pow(r[N_r-1],3);
	}
	return tot;
}

//Mass-weighted average
double calc_mass_avg(double f[], double rho[], double m, int n)  {

	// int i;
	double tot    = 0;
	double m_tot  = 0;

	for (int i{ 1 }; i < n; i++)  {
		m_tot += m*rho[i]*4.*pi/3.*(pow(r[i],3) - pow(r[i-1],3));
		tot   += f[i]*m*rho[i]*4.*pi/3.*(pow(r[i],3) - pow(r[i-1],3));
    }
    return tot/m_tot;
}

//Calculate ionization front position and velocity
// double* calc_ifront(double x[], double xprev[], double radius[], double deltat, int n)  {
double* calc_ifront(double x[], double xprev[], double radius[], int n)
{
	// int i;
	double pos, pos_prev, vel_IF;
	double *ifront = (double*) malloc(sizeof(double) * 2);

	pos       = interpolate(x, radius, 0.5, n);
	pos_prev  = interpolate(xprev, radius, 0.5, n);
	vel_IF       = (pos - pos_prev)/(t_max*yr_to_s*1e6/(N_output - 1));

	*(ifront + 0) = pos;
	*(ifront + 1) = vel_IF;

	return ifront;
}

double* calc_ifront_flux_method()
{
	// double pos_IF_rear{};
	double I_div_nu_pair[2][N_nu]{};
	// double I_div_nu_below[N_nu]{};
	double inc_photon_flux_pair[2]{};
	double inc_photon_flux{};
	double f_H1_pair[2];
	double nH_boundary{};
	double vel_IF{};
	double *flux_vel = (double*) malloc(sizeof(double) * 3);


	//EDIT BAYU
	// index_rear_IF    		= interpolate(f_H1, r, 0.1, N_r);
	double f_H1_IF_rear = 0.1;
	for (int i=0;i<(N_r-1);i++) //get index of array at the location closest to 10% neutral
	{
		if ((f_H1[i]<f_H1_IF_rear) && (f_H1[i+1]>f_H1_IF_rear))
		{
			for (int j=0;j<N_nu;j++) //to make the following integration work
			{
				I_div_nu_pair[0][j] = I_nu[i][j]/nu[j]; // int (I_nu / h / nu dnu ) = incident Flux ionizing photons
				I_div_nu_pair[1][j] = I_nu[i+1][j]/nu[j];
				f_H1_pair[0] = f_H1[i];
				f_H1_pair[1] = f_H1[i+1];
				index_rear_IF = i;
			}
			break;
		}
	}
	inc_photon_flux_pair[0] = trapz_int(I_div_nu_pair[0],nu,N_nu)*4*pi/h;
	inc_photon_flux_pair[1] = trapz_int(I_div_nu_pair[1],nu,N_nu)*4*pi/h;
	inc_photon_flux = interpolate(f_H1_pair,inc_photon_flux_pair,f_H1_IF_rear,2);
	//EDIT BAYU

	//TESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTING
	// double I_div_nu[N_nu]{};
	// double pos_IF_rear{};
	// pos_IF_rear    		= interpolate(f_H1, r, 0.1, N_r);
	// for (int i=0;i<(N_r-1);i++) //get index of array at the location closest to 10% neutral
	// {
	// 	if ((r[i]<pos_IF_rear) && (r[i+1]>pos_IF_rear))
	// 	{
	// 		double a = abs(r[i]-pos_IF_rear);
	// 		double b = abs(r[i+1]-pos_IF_rear);
	// 		if (a<=b) {
	// 			index_rear_IF = i;
	// 			// write_spectrum(index);
	// 			break;
	// 		}
	// 		else {
	// 			index_rear_IF = i+1;
	// 			// write_spectrum(index);
	// 			break;
	// 		}
	// 	}
	// }
	// // printf("%f\n", I_nu[index][5]);
	// // printf("pos_IF_rear: %.1f index: %d N_r: %d \n", pos_IF_rear/g_constants::kpc_to_cm,index, N_r);
	// // index-=2;
	// for (int j=0;j<N_nu;j++) //to make the following integration work
	// {
	// 	I_div_nu[j] = I_nu[index_rear_IF][j]/nu[j]; // int (I_nu / h / nu dnu ) = incident Flux ionizing photons
	// }
	// inc_photon_flux = trapz_int(I_div_nu,nu,N_nu)*4*pi/h; //integrating over frequency.... I_nu = int j_nu dr
	//TESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTINGTESTING

	nH_boundary    = interpolate(f_H1, nH, 0.5, N_r); // number density neutral H at front boundary 50%
	vel_IF = (c*inc_photon_flux) / (inc_photon_flux + c*nH_boundary*(1+chi_He));
	// vel = inc_photon_flux/(nH_boundary*(1+chi_He));
	*(flux_vel + 0) = inc_photon_flux;
	*(flux_vel + 1) = vel_IF;
	*(flux_vel + 2) = nH_boundary;
 	return flux_vel;
}


// double calc_ifront_flux_method1(double neutral_HI_fraction[], double neutral_HI_number_density[], double radius[], int n)
// {
// 	double pos_IF_rear{};
// 	int index{};
// 	double I_div_nu[N_nu];
// 	double inc_photon_flux{};
// 	double nH1_boundary{};
// 	double vel{};
// 	pos_IF_rear    		= interpolate(neutral_HI_fraction, radius, 0.1, n);
// 	for (int i=0;i<n;i++) //get index of array at the location closest to 10% neutral
// 	{
// 		if ((radius[i]<pos_IF_rear) && (radius[i+1]>pos_IF_rear))
// 		{
// 			double a = abs(radius[i]-pos_IF_rear);
// 			double b = abs(radius[i+1]-pos_IF_rear);
// 			if (a<=b) {
// 				index = i;
// 			}
// 			else {
// 				index = i+1;
// 			}
// 		}
// 	}
//
// 	for (int j=0;j<N_nu;j++) //to make the following integration work
// 	{
// 		I_div_nu[j] = I_nu[index][j]/nu[j];
// 	}
// 	inc_photon_flux = trapz_int(I_div_nu,nu,N_nu)/pi/h; //integrating over frequency
// 	nH1_boundary    = interpolate(neutral_HI_fraction, neutral_HI_number_density, 0.5, n); // number density neutral H at front boundary 50%
// 	vel = (c*inc_photon_flux) / (inc_photon_flux + c*nH1_boundary*(1+chi_He));
// 	return vel;
// }
//
// double calc_ifront_flux_method2(double fH1_arr, double nH1_arr, double pos_arr, int n_arr)
// {
// 	double pos_IF_rear{};
// 	double pos_IF_half{};
// 	double regionIF{};
// 	pos_IF_rear    		= interpolate(fH1_arr, pos_arr, 0.1, n_arr); // position (pkpc) where fH1 is 10% (neutral neutral)
// 	pos_IF_half  			=	interpolate(fH1_arr, pos_arr, 0.5, n_arr);
// 	regionIF 					= 25*kpc_to_cm; //
// 	for (int i=0;i<n_arr;i++) //
// 		if ((pos_arr[i]>(locIF-regionIF))&&(pos_arr[i]<(locIF+regionIF)))
// 		{
// 			if ((pos_arr[i]<pos_IF_rear) && (pos_arr[i+1]>pos_IF_rear))
// 			{
// 				double rear_back = abs(pos_arr[i]-pos_IF_rear);
// 				double rear_front = abs(pos_arr[i+1]-pos_IF_rear);
// 				if (rear_back<=rear_front) {
// 					rear_index = i;
// 				}
// 				else {
// 					rear_index = i+1;
// 				}
// 			}
// 			if ((pos_arr[i]<pos_IF_half) && (pos_arr[i+1]>pos_IF_half))
// 			{
// 				double half_back  = abs(pos_arr[i]-pos_IF_half);
// 				double half_front =  abs(pos_arr[i+1]-pos_IF_half);
// 				if (half_back<=half_front) {
// 					half_index = i;
// 				}
// 				else {
// 					half_index = i+1;
// 				}
// 			}
// 		}
// 	}
// }

		// if ((radius[i]<pos_IF_rear) && (radius[i+1]>pos_IF_rear))
		// {
		// 	double a = abs(radius[i]-pos_IF_rear);
		// 	double b = abs(radius[i+1]-pos_IF_rear);
		// 	if (a<=b) {
		// 		index = i;
		// 	}
		// 	else {
		// 		index = i+1;
		// 	}
		// }


	// inc_photon_flux = trapz_int(spec_intensity/h/nu,r,N_r)/pi;


	// double flux_rear, nH1_boundary, vel;
	// // double ifront;//(double*) malloc(sizeof(double) * 2);
	// pos_IF_boundary   = interpolate(x, radius, 0.5, n); // location of the IF, front boundary
	// pos_IF_rear    		= interpolate(x, radius, 0.1, n); // location of IF, rear
	//
	// for (i=0;i<n;i++) //emissivity only behind IF rear
	// {
	// 	if (radius[i]>pos_IF_rear)
	// 	{
	// 		emissivity[i] = 0.;
	// 	}
	//
	// 	//
	// 	// {
	// 	// 	int index = i-1;
	// 	// 	return index;
	// 	// 	break;
	// 	}
	// }
	//
	//
	//
	// double flux = trapz_int(emissivity,r,N_r)*(lambda_lya_cm/h/c)/pi; //integrated lya flux
	//
	// flux_rear       = interpolate(x, emissivity, 0.1, n); // rear of IF, 10% neutral
	// nH1_boundary    = interpolate(x, emissivity, 0.5, n); // number density neutral H at front boundary 50%
	//
	// vel = (c * flux_rear) / (c*flux_rear + c*nH1_boundary*(1+chi_He));
