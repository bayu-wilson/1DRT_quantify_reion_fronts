#include "user_inputs.h"
#include "global_variables.h"
#include <stdlib.h>     /* abs */
#include "general_funcs.h" /* mind */
#include "rates.h"
#include "global_constants.h"
#include <math.h>

// #include <omp.h>
#include <stdio.h>

using user_inputs::parallel;
using user_inputs::tfrac;
using user_inputs::FINITE_C;
using user_inputs::hydro;
using user_inputs::temp_ev;
using user_inputs::min_frac;
using user_inputs::min_dt_late;
using user_inputs::min_dt_early;
using user_inputs::spherical; //Bayu Sept. 10 2023
using user_inputs::time_decay_index; //Bayu Sept. 11 2023
using user_inputs::time_naught_decay; //Bayu Sept. 11 2023
using g_constants::c;
using g_constants::pi;
using g_constants::kpc_to_cm;
using namespace rates;


void update_dt()
{
	double min_dt;
	// min_dt = min_dt_late; // Since we got rid of geometric attenuation, following commented section is no longer necessary
	min_dt = min_dt_early;//Bayu
	// //testing
	// if (t/yr_to_s/1e6 < 1.0)
  // {
	// 	min_dt = min_dt_early;
	// }
	// else
  // {
	// 	min_dt = min_dt_late;
	// }

	dt = t_max*yr_to_s*1e6/(N_output - 1);
	//dt = user_inputs::R/user_inputs::N_r*kpc_to_cm / c * 10; //s
	//printf("dt: %.3e \n",dt);
	//dt = 0.013045867681251322*yr_to_s*1e6; //Myr
	
	#pragma omp parallel for reduction(min:dt) if (parallel)
	for (int i=prev_index; i < N_r; i++) {
		// if ((r[i]>(pos_IF-region_IF*kpc_to_cm)) && (r[i]<(pos_IF+region_IF*kpc_to_cm)))
		// if ((r[i]>back_IF) && (r[i]<front_IF))
		{
			//chemical timescales
			if ( (nH1[i] > min_frac*nH[i]) && (abs(tfrac*nH1[i]/dnH1_dt[i]) >  min_dt*yr_to_s*1e6) ) {
				dt = mind(dt, abs(tfrac*nH1[i]/dnH1_dt[i]));
			}
			if ( (nH2[i] > min_frac*nH[i]) && (abs(tfrac*nH2[i]/dnH2_dt[i]) > min_dt*yr_to_s*1e6) ) {
				dt = mind(dt, abs(tfrac*nH2[i]/dnH2_dt[i]));
			}
			if ( (nHe1[i] > min_frac*nHe[i]) && (abs(tfrac*nHe1[i]/dnHe1_dt[i]) > min_dt*yr_to_s*1e6) ) {
				dt = mind(dt, abs(tfrac*nHe1[i]/dnHe1_dt[i]));
			}
			if ( (nHe2[i] > min_frac*nHe[i]) && (abs(tfrac*nHe2[i]/dnHe2_dt[i]) > min_dt*yr_to_s*1e6) ) {
				dt = mind(dt, abs(tfrac*nHe2[i]/dnHe2_dt[i]));
			}
			if ( (nHe3[i] > min_frac*nHe[i]) && (abs(tfrac*nHe3[i]/dnHe3_dt[i]) > min_dt*yr_to_s*1e6) ) {
				dt = mind(dt, abs(tfrac*nHe3[i]/dnHe3_dt[i]));
			}

			if ( (temp_ev == TRUE) && (abs(tfrac*temp[i]/dT_dt[i]) > min_dt*yr_to_s*1e6) )
	    {
				//thermal evolution timescale
				dt = mind(dt, abs(tfrac*temp[i]/dT_dt[i]));
			}

			if (hydro == TRUE) {
				dt = mind(dt, abs(tfrac*delta_r[0]/vel[i]));
				if ( abs(tfrac*rho_total[i]/drho_dt[i]) > min_dt*yr_to_s*1e6 ) {
					dt = mind(dt, abs(tfrac*rho_total[i]/drho_dt[i]));
				}
				if ( abs(tfrac*vel[i]/dvel_dt[i]) > min_dt*yr_to_s*1e6 ) {
					dt = mind(dt, abs(tfrac*vel[i]/dvel_dt[i]));
				}
				if ( abs(tfrac*e_total[i]/de_dt[i]) > min_dt*yr_to_s*1e6 ) {
					dt = mind(dt, abs(tfrac*e_total[i]/de_dt[i]));
				}
			}
		}
	}
}

void update_step()
{
	for (int i{ 1 }; i < N_r; i++)
	{
		f_H1_step[i]  = f_H1[i];
		f_H2_step[i]  = f_H2[i];
		f_He1_step[i] = f_He1[i];
		f_He2_step[i] = f_He2[i];
		f_He3_step[i] = f_He3[i];

		t_step        = t;
	}
}

void solve_planeparallel_rt()
{
	double time_naught_seconds = time_naught_decay*yr_to_s*1e6; //convert from Myr to seconds
	double time_dependence = ( pow(t+time_naught_seconds-dt,-time_decay_index)-pow(t+time_naught_seconds,-time_decay_index) \
                                 ) / pow(time_naught_seconds,-time_decay_index);
	//#pragma omp parallel for if (parallel)
        //for (int i=0; i < N_r; i++)
        //{
	#pragma omp parallel for if (parallel)
        for (int j=0; j < N_nu; j++)
        {
             	I_nu_prev[0][j] = I_nu[0][j]; //first remember previous
		I_nu[0][j] = I_nu_prev[0][j] - time_dependence*I_nu_initial[j]; //then change current
		//I_nu[i][j] = I_nu_prev[i][j]/1.1; // then change current
        }
	//}


	//printf("%.1e \t ", time_dependence);
	//printf("%.1e \n ", I_nu[1][0]);

	//double time_decay_factor{0};
	//double t_Myr = t/yr_to_s/1e6;
	//time_decay_factor = pow(1+t_Myr/time_naught_decay,-time_decay_index); 
	//printf("gamma: %.1e \t",time_decay_factor);
	//double time_naught_seconds = time_naught_decay*yr_to_s*1e6; //convert from Myr to seconds
	//double time_dependence =  time_decay_index * dt/time_naught_seconds * pow(1+t/time_naught_seconds,-time_decay_index-1);
	//double time_dependence = ( pow(t+time_naught_seconds-dt,-time_decay_index)-pow(t+time_naught_seconds,-time_decay_index) 
	//			 ) / pow(time_naught_seconds,-time_decay_index);
	//double time_dependence = pow(1+t/time_naught_seconds,-time_decay_index)
	//double time_dependence = 
	double optical_depth[N_nu]{};
	#pragma omp parallel for if (parallel)
        for (int j=0; j < N_nu; j++)
	{
		optical_depth[j]=0; 
		for (int i=1; i < N_r; i++)
		{
			//solve if we want artificial time-dependence
			//I_nu[i][j] = I_nu_prev[i][j]/1.1;//I_nu[0][j] * 9/1/yr_to_s/1e6 * dt;
			//I_nu[i][j] = I_nu_prev[i][j] - time_dependence*I_nu[0][j];
			//solve plane-parallel rt with trapezoidal integral of the exponent to get the optical depth. 
			optical_depth[j]+= 0.5 * (gamma_nu_tot[i][j]+gamma_nu_tot[i-1][j])*(r[i]-r[i-1]); //each step increases the optical depth 
			//optical_depth+=gamma_nu_tot[i][j]*(r[i]-r[i-1]); 
			I_nu[i][j] = I_nu[0][j] * exp(-optical_depth[j]); //* time_decay_factor;  		
			//I_nu_prev[i][j] = I_nu[i][j];
			//printf("%.1e \t",(gamma_nu_tot[i][j]+gamma_nu_tot[i-1][j]));
			//printf("%.1e \t",(r[i]-r[i-1]));
			//printf("%.1e \n",(gamma_nu_tot[i][j]+gamma_nu_tot[i-1][j])*(r[i]-r[i-1]));
			//printf("%.1e \n", 1/2 * (gamma_nu_tot[i][j]+gamma_nu_tot[i-1][j])*(r[i]-r[i-1]));
			//printf("i=%d, j=%d tau=%.1e\n",i,j,optical_depth);
			//fflush(stdout);
		}
		//printf("optical_depth: %.1e \n",optical_depth);
        	//fflush(stdout);
	}	
	//printf("%.1f \t ", time_decay_index);
        //printf("%.1e \t ", t/yr_to_s/1e6);
	//printf("%.1e \t ", dt/yr_to_s/1e6);
        //printf("%.1e \t ", time_naught_seconds);
	//printf("%.1e \t ", time_dependence);
	//printf("%.1e \t ",I_nu[0][0] * 9./1000./yr_to_s/1.e6*dt);
	//printf("%.1e \t ", I_nu_prev[1][0]/1.1);
	//printf("%.1e \n ", I_nu[1][0]);
        //printf("%.1e \n ", I_nu_prev[0][0]);
	//printf("I_nu: %.1e \n",I_nu[1][0]);
	//fflush(stdout);
}

void solve_spherical_rt()
{
	//BAYU: switched the j and i loops Jan 31, 2023 
	#pragma omp parallel for if (parallel)
	for (int j=0; j < N_nu; j++) //this is parallelizable
	{
		//for (int j=0; j < N_nu; j++) // cannot parallelize because there's an i and i-1
		for (int i=1; i < N_r; i++)
		{
			if (FINITE_C == FALSE) //you're probably using this one
			{
				double a = 1./delta_r[i-1] + gamma_nu_tot[i][j];
				double b = pow(r[i-1]/r[i],2.)/delta_r[i-1]*I_nu[i-1][j];
				I_nu[i][j] = b/a;
			}

			else if (FINITE_C == TRUE)
			{
				double a = 1./c/dt + 1./delta_r[i-1] + gamma_nu_tot[i][j];
				double b = 1./c/dt*I_nu[i][j] + pow(r[i-1]/r[i],2.)/delta_r[i-1]*I_nu[i-1][j];
				I_nu[i][j] = b/a;
			}
		}
	}

	if (FINITE_C == TRUE)
	{
		#pragma omp parallel for if (parallel)
		for (int i=0; i < N_r; i++)
		{
			for (int j=0; j < N_nu; j++)
			{
				I_nu_prev[i][j] = I_nu[i][j];
			}
		}
	}
}

void solve_rt()
{
        if (spherical){
                solve_spherical_rt();
        }
        else if (!spherical){ //using plane-parallel
                solve_planeparallel_rt();
        }
}


void update_gamma_nu()
{
	#pragma omp parallel for if (parallel)
	for (int i=0; i < N_r; i++)
	{
		for (int j=0; j < N_nu; j++)
		{
			// update absorption coefficient
			gamma_nu_H1[i][j]  = nH1[i]*sigmapi_H1(nu[j]);
			gamma_nu_He1[i][j] = nHe1[i]*sigmapi_He1(nu[j]);
			gamma_nu_He2[i][j] = nHe2[i]*sigmapi_He2(nu[j]);
			gamma_nu_tot[i][j] = gamma_nu_H1[i][j] + gamma_nu_He1[i][j] + gamma_nu_He2[i][j];
		}
	}
}

//update the energy density of radiation
void update_u_nu()
{
	#pragma omp parallel for if (parallel)
	for (int i=0; i < N_r; i++)
	{
		for (int j=0; j < N_nu; j++)
		{
			u_nu[i][j] = 4*pi/c*I_nu[i][j];
			//if (r[i] - (R0 + dr)*kpc_to_cm > c*t)  { //causal correction
			//	u_nu[i][j] = 0.;
			//}
		}

	}
}

//artificial attentuation in order to explore more IF speeds. time has units of seconds. eta is the power law of the attenuation function
//void artificial_attentuation(double t0, double eta_power)
//{
//	#pragma omp parallel for if (parallel)
//	for (int j=0; j<N_nu; j++)
//	{
//		I_nu[0][j] *= pow((t/t0),eta_power)*I_nu[0][j]; //
//	}
//}
