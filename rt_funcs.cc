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
using g_constants::c;
using g_constants::pi;
using g_constants::kpc_to_cm;
using namespace rates;

void update_dt()
{
	double min_dt;

	//testing
	if (t/yr_to_s/1e6 < 1.0)
  {
		min_dt = min_dt_early;
	}
	else
  {
		min_dt = min_dt_late;
	}
	// min_dt = min_dt_late; // Since we got rid of geometric attenuation, this is no longer necessary

	dt = t_max*yr_to_s*1e6/(N_output - 1);
// dt = min_dt*yr_to_s*1e6;

	//Adaptive time-stepping
	// double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);
	// double locIF = *(ifront_H1 + 0);
	double pos_IF    		= interpolate(f_H1, r, 0.5, N_r);
	double regionIF = (4*2.7+2)*kpc_to_cm;

	//This one seems to help save a substantial amount of time if (parallel).
	// should work in parallel if each of the cells don't depend on eachother.
	#pragma omp parallel for if (parallel)
	for (int i=0; i < N_r; i++)
  {
		if ((r[i]>(pos_IF-regionIF)) && (r[i]<(pos_IF+regionIF)))
		{
			//chemical timescales
			if ( (nH1[i] > min_frac*nH[i]) && (abs(tfrac*nH1[i]/dnH1_dt[i]) >  min_dt*yr_to_s*1e6) )
	    {
				dt = mind(dt, abs(tfrac*nH1[i]/dnH1_dt[i])); // TODO: THERE IS A PROBLEM WITH dnH1_dt!!!
				// printf("%.1e ",dnH1_dt[i]);
			}
			if ( (nH2[i] > min_frac*nH[i]) && (abs(tfrac*nH2[i]/dnH2_dt[i]) > min_dt*yr_to_s*1e6) )
	    {
				dt = mind(dt, abs(tfrac*nH2[i]/dnH2_dt[i]));
			}
			if ( (nHe1[i] > min_frac*nHe[i]) && (abs(tfrac*nHe1[i]/dnHe1_dt[i]) > min_dt*yr_to_s*1e6) )
	    {
				dt = mind(dt, abs(tfrac*nHe1[i]/dnHe1_dt[i]));
			}
			if ( (nHe2[i] > min_frac*nHe[i]) && (abs(tfrac*nHe2[i]/dnHe2_dt[i]) > min_dt*yr_to_s*1e6) )
	    {
				dt = mind(dt, abs(tfrac*nHe2[i]/dnHe2_dt[i]));
			}
			if ( (nHe3[i] > min_frac*nHe[i]) && (abs(tfrac*nHe3[i]/dnHe3_dt[i]) > min_dt*yr_to_s*1e6) )
	    {
				dt = mind(dt, abs(tfrac*nHe3[i]/dnHe3_dt[i]));
			}

			if ( (temp_ev == TRUE) && (abs(tfrac*temp[i]/dT_dt[i]) > min_dt*yr_to_s*1e6) )
	    {
				//thermal evolution timescale
				dt = mind(dt, abs(tfrac*temp[i]/dT_dt[i]));
			}

			if (hydro == TRUE)
			{
				dt = mind(dt, abs(tfrac*delta_r[0]/vel[i]));
				if ( abs(tfrac*rho_total[i]/drho_dt[i]) > min_dt*yr_to_s*1e6 )
				{
					dt = mind(dt, abs(tfrac*rho_total[i]/drho_dt[i]));
				}
				if ( abs(tfrac*vel[i]/dvel_dt[i]) > min_dt*yr_to_s*1e6 )
				{
					dt = mind(dt, abs(tfrac*vel[i]/dvel_dt[i]));
				}
				if ( abs(tfrac*e_total[i]/de_dt[i]) > min_dt*yr_to_s*1e6 )
				{
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

void solve_spherical_rt()
{
	//parallel is slower?? if (parallel)
	#pragma omp parallel for
	for (int i=1; i < N_r; i++)
	{
		for (int j=0; j < N_nu; j++)
		{
			if (FINITE_C == FALSE)
			{
				double a = 1./delta_r[i-1] + gamma_nu_tot[i][j];
				double b = pow(r[i-1]/r[i],2.)/delta_r[i-1]*I_nu[i-1][j];
				I_nu[i][j] = b/a;
				// printf("%.1e\n", gamma_nu_tot[i][j]);
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

void update_gamma_nu()
{
	// double dgamma; //NOTE: why was this initialized?
	//without parallelization it is 2.657e-03 seconds

	//collapse(2)
	//doesn't really do anything I think if (parallel)
	// #pragma omp parallel for if (parallel)
	#pragma omp parallel for
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
	//saves some time for if (parallel)
	#pragma omp parallel for
	for (int i=0; i < N_r; i++)
	{
		for (int j=0; j < N_nu; j++)
		{
			u_nu[i][j] = 4*pi/c*I_nu[i][j];
		}
	}
}
