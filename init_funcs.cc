#include <math.h>
#include <string.h> // for strcmp
#include <stdio.h>

#include "user_inputs.h"
#include "global_constants.h"
#include "global_variables.h"
#include "general_funcs.h" //like power_law
#include "rates.h" //for the rate equations and rates namespace
#include "io_funcs.h"

using g_constants::TRUE;
using g_constants::FALSE;
using g_constants::kpc_to_cm;
using g_constants::gam;
using g_constants::pi;
using namespace user_inputs;
using namespace rates;

void init_grid()
{
	//read in the spatial, frequency, and angular grids if available
	if (input_grid == TRUE)
  {
		read_grid();
    r[0] = r[0]*kpc_to_cm;
		for (int i{ 1 }; i < N_r; i++)
    {
			r[i] = r[i]*kpc_to_cm;
      delta_r[i-1] = r[i] - r[i-1];

			// rho_total[i] = rho_total[i]/10;s
		}
	}
	//if grids are not user-defined, generate a linear spatial grid and a
	//logarithmic frequency grid
	else
  {
		//spatial grid
		r[0] = R0*kpc_to_cm;
		for (int i{ 1 }; i < N_r; i++)
    {
			r[i]       = (R0 + i*dr)*kpc_to_cm;
			delta_r[i-1] = dr*kpc_to_cm;
		}

	}
}


void init_gas()
{
	if (input_grid == TRUE)
  {
		//do nothing if reading
	}
	else
  {
		for (int i{ 0 }; i < N_r; i++)
		{
			rho_total[i] = power_law(r[i], r[0], rho_0, rho_index); //rho_0 is in user_inputs.h, it is uniform density (in cgs)
			rho_prev[i]  = rho_total[i];

			f_H1[i] = power_law(r[i], r[0], fH1_0, fH1_index);

			if (f_H1[i] > 1.)
			{
				f_H1[i] = fH1_0;
			}
			f_H2[i] = 1. - f_H1[i];
			f_He1[i] = power_law(r[i], r[0], fHe1_0, fHe1_index);
			if (f_He1[i] > 1.)
			{
					f_He1[i] = fHe1_0;
			}
			f_He2[i] = power_law(r[i], r[0], fHe2_0, fHe2_index);
			if (f_He2[i] > 1. - f_He1[i])
			{
				f_He2[i] = 1. - f_He1[i];
			}
			f_He3[i]  = 1. - f_He1[i] - f_He2[i];
			temp[i]   = temp_0;
			nH[i]     = (1. - Y)*rho_total[i]/m_H;
			nHe[i]    = Y*rho_total[i]/(4.*m_H);
			nH2[i]    = f_H2[i]*nH[i];
			nHe2[i]   = f_He2[i]*nHe[i];
			nHe3[i]   = f_He3[i]*nHe[i];
			ne[i]     = nH2[i] + nHe2[i] + 2.*nHe3[i];
			n_tot[i]  = nH[i] + nHe[i] + ne[i];
			p_total[i] = k_B*temp[i]*n_tot[i]; //p is pressure
			//assume no initial radial velocity
			vel[i]       = 0;
		}
	}

	//compute other important quantities
	for (int i{ 0 }; i < N_r; i++)  {

		//number densities //NOTE: WHY IS THIS DONE TWICE?
		nH[i]         = (1. - Y)*rho_total[i]/m_H;
		nH_prev[i]    = nH[i];
		nHe[i]        = Y*rho_total[i]/(4.*m_H);
		nHe_prev[i]   = nHe[i];
		nH1[i]        = f_H1[i]*nH[i];
		nH1_prev[i]   = nH1[i];
		nHe1[i]       = f_He1[i]*nHe[i];
		nHe1_prev[i]  = f_He1[i]*nHe[i];
		nH2[i]        = f_H2[i]*nH[i];
		nH2_prev[i]   = nH2[i];
		nHe2[i]       = f_He2[i]*nHe[i];
		nHe2_prev[i]  = nHe2[i];
		nHe3[i]       = f_He3[i]*nHe[i];
		nHe3_prev[i]  = nHe3[i];
		ne[i]         = nH2[i] + nHe2[i] + 2.*nHe3[i];
		ne_prev[i]    = nH2[i] + nHe2[i] + 2.*nHe3[i];
		n_tot[i]      = nH[i] + nHe[i] + ne[i];
		n_tot_prev[i] = nH[i] + nHe[i] + ne[i];

		//recombination coefficients
		if (strcmp(rec_case, "A") == 0)  {
			recomb_H2[i]  = alphaA_H2(temp[i]);
			recomb_He2[i] = alphaA_He2(temp[i]);
			recomb_He3[i] = alphaA_He3(temp[i]);
		}
		else  {
			recomb_H2[i]  = alphaB_H2(temp[i]);
			recomb_He2[i] = alphaB_He2(temp[i]);
			recomb_He3[i] = alphaB_He3(temp[i]);
		}

		//compute first derivatives of the number densities
		dnH2_dt[i]  = gamma_H1_tot[i]*nH1[i] - recomb_H2[i]*ne[i]*nH2[i];
		dnH1_dt[i]  = -dnH2_dt[i];
		dnHe2_dt[i] = nHe1[i]*gamma_He1_tot[i] - recomb_He2[i]*ne[i]*nHe2[i]
						+ recomb_He3[i]*ne[i]*nHe3[i] - nHe2[i]*gamma_He2_tot[i];
		dnHe3_dt[i] = nHe2[i]*gamma_He2_tot[i] - recomb_He3[i]*ne[i]*nHe3[i];
		dnHe1_dt[i] = -dnHe2_dt[i] - dnHe3_dt[i];
		dne_dt[i]   = dnH2_dt[i] + dnHe2_dt[i] + 2.*dnHe3_dt[i];

		//hydrodynamics

		//previous
		rho_prev[i]  = rho_total[i];
		temp_prev[i] = temp[i];
		p_prev[i]    = p_total[i];
		vel_prev[i]  = vel[i];

		//thermal + kinetic energy
		e_total[i]   = p_total[i] / (gam - 1) + 0.5*rho_total[i]*pow(vel[i], 2.);
		e_prev[i]    = e_total[i];
		rhov[i]      = rho_total[i]*vel[i];
		rhov_prev[i] = rhov[i];
		cs[i]        = pow(gam*p_total[i]/rho_total[i], 0.5);
		//set the initial graviational potential function (which will be 0 everywhere in the
		//homogeneous density case).
		phi[i]       = 0.;

		//fluxes
		FnH[i]     = nH[i]*vel[i];
		FnHe[i]    = nHe[i]*vel[i];
		FnH1[i]    = nH1[i]*vel[i];
		FnH2[i]    = nH2[i]*vel[i];
		FnHe1[i]   = nHe1[i]*vel[i];
		FnHe2[i]   = nHe2[i]*vel[i];
		FnHe3[i]   = nHe3[i]*vel[i];

		Frhov[i]   = rhov[i]*vel[i];
		Fe[i]      = (e_total[i] + p_total[i])*vel[i];
	}
}

void init_intensity()
{
	double norm, phii, log_nu;
	double temp_nu_min = log10(1e11);
	double temp_nu_max = log10(1e18);
	double temp_nu[500];
	double phi_nu[500];

	//initialize intensity to be 0 everywhere except at r = 0.  Read in or generate source spectrum.
	if (input_source == TRUE)
	{
		read_source();
	}
	else
	{
		//frequency grid
		nu[0] = nu_min;
		for (int j{ 1 }; j < N_nu; j++)
		{
			log_nu          = log10(nu_min) + j*dlog_nu;
			nu[j]           = pow(10, log_nu);
			delta_nu[j - 1] = nu[j] - nu[j - 1];
		}

		if (strcmp(spectrum, "BLACKBODY") == 0) //PROBABLY USING THIS ONE
		{
			double temp_dlog_nu = (temp_nu_max - temp_nu_min)/(500. - 1.);

			for (int j{ 0 }; j < 500; j++)
			{
				temp_nu[j] = pow(10,temp_nu_min + j*temp_dlog_nu);
				phi_nu[j]  = b_nu(pow(10,temp_nu_min + j*temp_dlog_nu), T_source);
			}

			norm = trapz_int(phi_nu, temp_nu, 500);

			for (int j{ 0 }; j < N_nu; j++)
			{
				phii            = interpolate(temp_nu, phi_nu, nu[j], 500);
				I_nu_prev[0][j] = Lum/16./pow(pi,2)/pow(R0*kpc_to_cm,2)*phii/norm;
				I_nu[0][j]      = Lum/16./pow(pi,2)/pow(R0*kpc_to_cm,2)*phii/norm;
			}
		}

		else if (strcmp(spectrum, "POWER") == 0)
		{
			for (int j{ 0 }; j < N_nu; j++)
			{
				phii = (alpha - 1)/nu[0]*pow(nu[j]/nu[0],-alpha);
				I_nu_prev[0][j] = Lum/16./pow(pi,2)/pow(R0*kpc_to_cm,2)*phii;
				I_nu[0][j]      = Lum/16./pow(pi,2)/pow(R0*kpc_to_cm,2)*phii;
			}
		}

		else if (strcmp(spectrum, "MONOCHROMATIC") == 0)
		{
			I_nu_prev[0][0] = Lum/16./pow(pi,2)/pow(R0*kpc_to_cm,2);
			I_nu[0][0]      = Lum/16./pow(pi,2)/pow(R0*kpc_to_cm,2);
		}
	}

	//initialize to 0 outside the point source.
	for (int i{ 1 }; i < N_r; i++)
	{
		for (int j{ 0 }; j < N_nu; j++)
		{
			I_nu_prev[i][j] = 0;
			I_nu[i][j]      = 0;
		}
	}
}
