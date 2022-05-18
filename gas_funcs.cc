// #include "gas_funcs.h"
// #include "init_funcs.cc"
// #include "cosmo_funcs.h"
// #include "rates.h"
#include "global_variables.h"
#include "rates.h"
#include "global_constants.h"
#include "user_inputs.h"
#include <string.h>
#include "cosmo_funcs.h"
#include "general_funcs.h"
#include <math.h>
#include <stdio.h>


using namespace rates;
using g_constants::c;
using g_constants::h;
using g_constants::gam;
using user_inputs::compton;
using user_inputs::z;
using user_inputs::rec_case;
using user_inputs::hydro;
using user_inputs::temp_0;
using user_inputs::coll_exc_cool;
using user_inputs::recomb_cool;
using user_inputs::coll_ion;
using user_inputs::spectrum;
using user_inputs::parallel;

// int dt = dt_default;
//compute photoionization rates for HI, HeI, and HeII
void update_gamma()  {

	// int i,j;
	// double dgamma;
	// double f_h1;
	// double f_h1, f_he1, f_he2;
	// if (parallel)
	// reduction(+:sum,result)
	// gamma_H1_tot[i]  += cic_H1(temp[i])*ne[i];
	// gamma_He1_tot[i] += cic_He1(temp[i])*ne[i];
	// gamma_He2_tot[i] += cic_He2(temp[i])*ne[i];

	//if (parallel)
	// reduction(+:gamma_H1_tot[:N_r],gamma_He1_tot[:N_r],gamma_He2_tot[:N_r])
	//reduction(+:gamma_H1_tot[:N_r],gamma_He1_tot[:N_r],gamma_He2_tot[:N_r])
	// #pragma omp parallel for if (parallel)
	#pragma omp parallel for
	for (int i=0; i < N_r; i++)
	{
		gamma_H1_tot[i]  = 0;
		gamma_He1_tot[i] = 0;
		gamma_He2_tot[i] = 0;
		double dgamma;

		if (strcmp(spectrum, "MONOCHROMATIC") == 0)
		{
			gamma_H1_tot[i]  = sigmapi_H1(nu[0])*c*u_nu[i][0]/h/nu[0];
			gamma_He1_tot[i] = sigmapi_He1(nu[0])*c*u_nu[i][0]/h/nu[0];
			gamma_He2_tot[i] = sigmapi_He2(nu[0])*c*u_nu[i][0]/h/nu[0];
		}
		else
		{
			for (int j=0; j < N_nu - 1; j++)
			{
				//primary photoionization
				dgamma            = delta_nu[j]*sigmapi_H1(nu[j])*c*u_nu[i][j]/h/nu[j];
				dgamma           += delta_nu[j]*sigmapi_H1(nu[j + 1])*c*u_nu[i][j + 1]/h/nu[j + 1];
				gamma_H1_tot[i]  += dgamma/2.;

				dgamma            = delta_nu[j]*sigmapi_He1(nu[j])*c*u_nu[i][j]/h/nu[j];
				dgamma           += delta_nu[j]*sigmapi_He1(nu[j + 1])*c*u_nu[i][j + 1]/h/nu[j + 1];
				gamma_He1_tot[i] += dgamma/2.;

				dgamma            = delta_nu[j]*sigmapi_He2(nu[j])*c*u_nu[i][j]/h/nu[j];
				dgamma           += delta_nu[j]*sigmapi_He2(nu[j + 1])*c*u_nu[i][j + 1]/h/nu[j + 1];
				gamma_He2_tot[i] += dgamma/2.;
			}
		}
		//collisional ionization counts as an effective photoionization rate with coeff
		//cic x ne
		if ( (coll_ion == TRUE) && (temp[i] >= 1e4) && (temp[i] <= 1e9) )
		{
			gamma_H1_tot[i]  += cic_H1(temp[i])*ne[i];
			gamma_He1_tot[i] += cic_He1(temp[i])*ne[i];
			gamma_He2_tot[i] += cic_He2(temp[i])*ne[i];
		}
	}
}

void update_heat_cool()
{
	// double dheat;//, dcool; //dcool is not used
	// double f_h1_heat, f_he1_heat, f_he2_heat;
	// double f_h1_heat  = 1.;
	// double f_he1_heat = 1.;
	// double f_he2_heat = 1.;
	//
	//if (parallel)
	//reduction(+:heat_rate[:N_r],cool_rate[:N_r]) //for num_threads(2)
		// shared(N_r,heat_rate,cool_rate,dheat, frac_heat,nH1,nHe1,nHe2)
	// #pragma omp parallel for if (parallel)
		#pragma omp parallel for
		for (int i=0; i < N_r; i++)
		{
			// printf("%.2f \t", heat_rate[i]);
			heat_rate[i] = 0; // why is it set to zero each time???
			cool_rate[i] = 0;
			if (strcmp(spectrum, "MONOCHROMATIC") == 0)
			{
				//photoionization heating
				// double f_h1_heat  = 1.;
				// double f_he1_heat = 1.;
				// double f_he2_heat = 1.;
				heat_rate[i] += nH1[i]*frac_heat[i][0]*sigmapi_H1(nu[0])*c*u_nu[i][0]*(h_eV*nu[0] - Eth_H1)/h_eV/nu[0];
				heat_rate[i] += nHe1[i]*frac_heat[i][0]*sigmapi_He1(nu[0])*c*u_nu[i][0]*(h_eV*nu[0] - Eth_He1)/h_eV/nu[0];
				heat_rate[i] += nHe2[i]*frac_heat[i][0]*sigmapi_He2(nu[0])*c*u_nu[i][0]*(h_eV*nu[0] - Eth_He2)/h_eV/nu[0];
			}
			else
			{
				// #pragma omp parallel for if (parallel)
				// #pragma omp for  schedule(static)
				for (int j=0; j < N_nu - 1; j++)
				{
					//photoionization heating
					double f_h1_heat  = 1.;
					double f_he1_heat = 1.;
					double f_he2_heat = 1.;
					double dheat;
					// }
					//integrate over frequency to get the excess energy from ionizations that goes
					//into heating the gas.
					dheat         = nH1[i]*f_h1_heat*sigmapi_H1(nu[j])*u_nu[i][j]*(h_eV*nu[j] - Eth_H1)/nu[j];
					dheat        += nH1[i]*f_h1_heat*sigmapi_H1(nu[j+1])*u_nu[i][j+1]*(h_eV*nu[j+1] - Eth_H1)/nu[j+1];
					dheat        += nHe1[i]*f_he1_heat*sigmapi_He1(nu[j])*u_nu[i][j]*(h_eV*nu[j] - Eth_He1)/nu[j];
					dheat        += nHe1[i]*f_he1_heat*sigmapi_He1(nu[j+1])*u_nu[i][j+1]*(h_eV*nu[j+1] - Eth_He1)/nu[j+1];
					dheat        += nHe2[i]*f_he2_heat*sigmapi_He2(nu[j])*u_nu[i][j]*(h_eV*nu[j] - Eth_He2)/nu[j];
					dheat        += nHe2[i]*f_he2_heat*sigmapi_He2(nu[j+1])*u_nu[i][j+1]*(h_eV*nu[j+1] - Eth_He2)/nu[j+1];
					heat_rate[i] += delta_nu[j]*c*dheat/2./h_eV;
				}
			}
			if ( (coll_ion == TRUE) && (temp[i] >= 1e4) && (temp[i] <= 1e9) )
			{
				//collisional ionization cooling contributions from each ion
				cool_rate[i] += cicr_H1(temp[i])*ne[i]*nH1[i];
				cool_rate[i] += cicr_He1(temp[i])*ne[i]*nHe1[i];
				cool_rate[i] += cicr_He2(temp[i])*ne[i]*nHe2[i];
			}

			//case A or B recombination cooling
			if (recomb_cool == TRUE)  {
				if (strcmp(rec_case, "A") == 0)
				{
					cool_rate[i] += rcrA_H2(temp[i])*ne[i]*nH2[i];
					cool_rate[i] += rcrA_He2(temp[i])*ne[i]*nHe2[i];
					cool_rate[i] += rcrA_He3(temp[i])*ne[i]*nHe3[i];
				}
				else if (strcmp(rec_case, "B") == 0)
				{
					cool_rate[i] += rcrB_H2(temp[i])*ne[i]*nH2[i];
					cool_rate[i] += rcrB_He2(temp[i])*ne[i]*nHe2[i];
					cool_rate[i] += rcrB_He3(temp[i])*ne[i]*nHe3[i];
				}
				else
				{
					printf("Warning: Ignoring recombinations\n");
				}
			}
			//collisional excitation cooling contributions from HI, HeI, and HeII
			if (coll_exc_cool == TRUE)
			{
				cool_rate[i] += coll_ex_rate_H1(temp[i])*ne[i]*nH1[i];
				cool_rate[i] += coll_ex_rate_He1(temp[i])*pow(ne[i],2.)*nHe1[i];
				cool_rate[i] += coll_ex_rate_He2(temp[i])*ne[i]*nHe2[i];
			}
			//compton heating/cooling.  Counts as a negative heating rate if T > T_cmb.
			if (compton == TRUE)
			{
				//heat_rate[i] += compton_rate(z)*ne[i]/rho_total[i]*(2.726*(1 + z) - temp[i]);
				heat_rate[i] += compton_rate(z)*ne[i]*(2.726*(1 + z) - maxd(temp_0, temp[i]));
			}
		}
}

void solve_ion(double ne[])  {
	//solve the ionizing balance equation for H usin a 1st order implicit backwards-difference scheme

	#pragma omp parallel for if (parallel)
	for (int i=0; i < N_r; i++)  {
		if (nH1_prev[i] >= nH2_prev[i])  {
			//Solve for nH2 when there is less ionized gas
			nH2[i]	= (nH2_prev[i] + gamma_H1_tot[i]*nH[i]*dt)
					/(1. + gamma_H1_tot[i]*dt + recomb_H2[i]*ne[i]*dt);
			nH1[i] 	= nH[i] - nH2[i];
			// printf("%f %f %f %f %f\n", nH2[i],nH2_prev[i],gamma_H1_tot[i],recomb_H2[i],ne[i]);
		}
		else  {
			//solve for nH1 when there is less neutral gas
			nH1[i]	= (nH1_prev[i] + recomb_H2[i]*ne[i]*nH[i]*dt)
					/(1. + gamma_H1_tot[i]*dt + recomb_H2[i]*ne[i]*dt);
			nH2[i]  = nH[i] - nH1[i];
		}

		//solve the ionizing balance equations for He.  Solve the backwards difference equation for the
		//species that have the smaller abundances, and solve closing equation for the remaining species.
		if ( (nHe1_prev[i] <= nHe3_prev[i]) && (nHe2_prev[i] <= nHe3_prev[i]) )  {
			nHe1[i] = (nHe1_prev[i] + recomb_He2[i]*ne[i]*(nHe[i] - nHe3[i])*dt)
					/(1. + gamma_He1_tot[i]*dt + recomb_He2[i]*ne[i]*dt);
			nHe2[i] = (nHe2_prev[i] + gamma_He1_tot[i]*(nHe[i] - nHe3[i])*dt
				    + recomb_He3[i]*ne[i]*(nHe[i] - nHe1[i])*dt)
					/(1. + (gamma_He1_tot[i] + gamma_He2_tot[i])*dt + (recomb_He2[i]
					+ recomb_He3[i])*ne[i]*dt);
			nHe3[i] = nHe[i] - nHe1[i] - nHe2[i];
		}
		else if( (nHe1_prev[i] <= nHe2_prev[i]) && (nHe3_prev[i] <= nHe2_prev[i]) )  {
			nHe1[i] = (nHe1_prev[i] + recomb_He2[i]*ne[i]*(nHe[i] - nHe3[i])*dt)
					/(1. + gamma_He1_tot[i]*dt + recomb_He2[i]*ne[i]*dt);
			nHe3[i] = (nHe3_prev[i] + gamma_He2_tot[i]*(nHe[i] - nHe1[i])*dt)
					/(1. + gamma_He2_tot[i]*dt + recomb_He3[i]*ne[i]*dt);
			nHe2[i] = nHe[i] - nHe1[i] - nHe3[i];
		}
		else {
			nHe2[i] = (nHe2_prev[i] + gamma_He1_tot[i]*(nHe[i] - nHe3[i])*dt
					+ recomb_He3[i]*ne[i]*(nHe[i] - nHe1[i])*dt)
					/(1. + (gamma_He1_tot[i] + gamma_He2_tot[i])*dt + (recomb_He2[i]
					+ recomb_He3[i])*ne[i]*dt);
			nHe3[i] = (nHe3_prev[i] + gamma_He2_tot[i]*(nHe[i] - nHe1[i])*dt)
					/(1. + gamma_He2_tot[i]*dt + recomb_He3[i]*ne[i]*dt);
			nHe1[i] = nHe[i] - nHe2[i] - nHe3[i];
		}
		ne[i]    = nH2[i] + nHe2[i] + 2.*nHe3[i];
	}
}

void update_chem()  {
	// double a;//, b;
	//if (parallel)
	// #pragma omp parallel for if (parallel)
	#pragma omp parallel for
	for (int i=0; i < N_r; i++)  {

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

		//diaelectric recombination of He II.
		if ( (temp[i] >= 3e4) && (temp[i] <= 1e6) )  {
			recomb_He2[i] += Dalpha_He2(temp[i]);
		}
	}
	#pragma omp parallel for if (parallel)
	for (int j=0; j < 5; j++)  {
		solve_ion(ne);
	}
	// #pragma omp parallel for if (parallel)
	#pragma omp parallel for
	for (int i=0; i < N_r; i++)  {

		n_tot[i] = nH[i] + nHe[i] + ne[i];

		//compute first derivatives of the number densities
		//testing: if hydro is turned on, evaluate all the time derivatives numerically there rather
		//than getting them here analytically.  If hydro is off, evaluate them using the chemistry equations
		//as usual.

		if (hydro == FALSE)  {
			// printf("%.3e\n", dt);
			dnH1_dt[i]  = (nH1[i] - nH1_prev[i])/dt;
			dnH2_dt[i]  = (nH2[i] - nH2_prev[i])/dt;
			dnHe1_dt[i] = (nHe1[i] - nHe1_prev[i])/dt;
			dnHe2_dt[i] = (nHe2[i] - nHe2_prev[i])/dt;
			dnHe3_dt[i] = (nHe3[i] - nHe3_prev[i])/dt;
			dn_dt[i]    = (n_tot[i] - n_tot_prev[i])/dt;
			dne_dt[i]   = (ne[i] - ne_prev[i])/dt;

			f_H1[i]  = nH1[i]/nH[i];
			f_H2[i]  = nH2[i]/nH[i];
			f_He1[i] = nHe1[i]/nHe[i];
			f_He2[i] = nHe2[i]/nHe[i];
			f_He3[i] = nHe3[i]/nHe[i];//ionization states

			nH_prev[i]    = nH[i];
			nHe_prev[i]   = nHe[i];
			nH1_prev[i]   = nH1[i];
			nH2_prev[i]   = nH2[i];
			nHe1_prev[i]  = nHe1[i];
			nHe2_prev[i]  = nHe2[i];
			nHe3_prev[i]  = nHe3[i];
			ne_prev[i]    = ne[i];
			n_tot_prev[i] = n_tot[i];
		}
	}
}

void update_thermal()
{
	double a;
	a = 1/(1 + z);
	// #pragma omp parallel for if (parallel)
	#pragma omp parallel for
	for (int i=0; i < N_r; i++)
	{
		if (hydro == FALSE)
		{
			temp[i]	     = (temp[i] + 2./3./k_B/n_tot[i]*(heat_rate[i] - cool_rate[i])*dt)
						 /(1. + 2*H(a)*dt + 1./n_tot[i]*dne_dt[i]*dt);

			dT_dt[i]     = 2./3./k_B/n_tot[i]*(heat_rate[i] - cool_rate[i])  - 2*H(a)*temp[i]
					     - temp[i]/n_tot[i]*dne_dt[i];
			temp_prev[i] = temp[i];
		}
		else
		{
			//if hydro is on, extract the temperature directly from the energy density
			temp[i]      = p_total[i]/k_B/n_tot[i];
			dT_dt[i]     = (temp[i] - temp_prev[i])/dt;
			temp_prev[i] = temp[i];
		}
	}
}

double minmod(double a, double b)  {
	return 0.5*(sign(a) + sign(b))*mind(abs(a), abs(b));
}

void solve_TVD_1(double u[], double F[], double S[], double cf[])  {

	// int i;
	double uR[N_r], uL[N_r], w[N_r];
	double FR[N_r], FL[N_r];
	double FR_plus, FR_minus, FL_plus, FL_minus, F_plus, F_minus;

	for (int i{ 0 }; i < N_r; i++)  {
		w[i]  = F[i] / cf[i];
		uL[i] = (u[i] - w[i])/2.;
		uR[i] = (u[i] + w[i])/2.;
		FR[i] = cf[i]*uR[i];
		FL[i] = cf[i]*uL[i];
	}
	//Take 1/2 time step using first-order upwind scheme.
	for (int i{ 1 }; i < N_r - 1; i++)  {
		FR_plus  = FR[i];
		FR_plus  += minmod((FR[i] - FR[i-1])/2., (FR[i+1] - FR[i])/2.);
		FL_minus = FL[i];
		FL_minus += minmod(-(FL[i] - FL[i-1])/2., -(FL[i+1] - FL[i])/2.);
		FR_minus = FR[i-1];
		if (i >= 2)  {
			FR_minus += minmod((FR[i-1] - FR[i-2])/2., (FR[i] - FR[i-1])/2.);
		}
		FL_plus  = FL[i+1];
		if (i < N_r - 2)  {
			FL_plus  += minmod(-(FL[i+1] - FL[i])/2., -(FL[i+2] - FL[i+1])/2.);
		}

		F_plus   = FR_plus - FL_plus;
		F_minus  = FR_minus - FL_minus;

		u[i]     = u[i] - (F_plus - F_minus)/delta_r[i-1]*dt + S[i]*dt;
	}
}

void solve_conservation_eqn(double u[], double flux[], double S[])  {

	// int iter, i;

	double F[N_r], cf[N_r];

	//Relaxing TVD scheme

	//Multiply by r^2 and get the freezing speed
	for (int i{ 0 }; i < N_r; i++)  {
		u[i]  = pow(r[i],2.)*u[i]; //hydrodynamic variable
		F[i]  = pow(r[i],2.)*flux[i]; //flux
		cf[i] = abs(vel[i]) + cs[i]; //freezing speed
	}

	solve_TVD_1(u, F, S, cf);

	//divide by r^2 to recover the original variables
	for (int i{ 0 }; i < N_r; i++)  {
		u[i]  = u[i]/pow(r[i],2.);
	}
    //dumb boundary conditions
    u[0]     = u[1];
	u[N_r-1] = u[N_r-2];
}

void update_hydro()  {
	// int i,k;
	// double mb;
	double S[N_r];

	//solve the poisson equation
	//solve_poisson_eqn(phi);

	for (int i{ 0 }; i < N_r; i++)  {
		S[i] = 0.;

		//update the fluxes
		FnH1[i]    = nH1[i]*vel[i];
		FnH2[i]    = nH2[i]*vel[i];
		FnHe1[i]   = nHe1[i]*vel[i];
		FnHe2[i]   = nHe2[i]*vel[i];
		FnHe3[i]   = nHe3[i]*vel[i];
	}
	//FIRST
	//solve the continuity equations for the number densities
	solve_conservation_eqn(nH, FnH, S);
	solve_conservation_eqn(nHe, FnHe, S);
	solve_conservation_eqn(nH1, FnH1, S);
	solve_conservation_eqn(nH2, FnH2, S);
	solve_conservation_eqn(nHe1, FnHe1, S);
	solve_conservation_eqn(nHe2, FnHe2, S);
	solve_conservation_eqn(nHe3, FnHe3, S);

	S[0]       = -pow(r[0],2.)*(p_total[1] - p_total[0])/delta_r[0]
			   - rho_total[0]*pow(r[0],2.)*(phi[1] - phi[0])/delta_r[0];
	S[N_r - 1] = -pow(r[N_r-1],2.)*(p_total[N_r-1] - p_total[N_r-2])/delta_r[N_r-2]
			   - rho_total[N_r-1]*pow(r[N_r-1],2.)*(phi[N_r-1] - phi[N_r-2])/delta_r[N_r-2];
	for (int i{ 1 }; i < N_r - 1; i++)  {
		S[i]   = -pow(r[i],2.)*(p_total[i+1] - p_total[i-1])/(delta_r[i] + delta_r[i-1])
			   - rho_total[i]*pow(r[i],2.)*(phi[i+1] - phi[i-1])/(delta_r[i] + delta_r[i-1]);
	}

	solve_conservation_eqn(rhov, Frhov, S);

	for (int i{ 0 }; i < N_r; i++)  {
		S[i]  = pow(r[i],2.)*(heat_rate[i] - cool_rate[i]);
	}

	solve_conservation_eqn(e_total, Fe, S);

	for (int i{ 0 }; i < N_r; i++)  {

		rho_total[i] = m_H*nH[i] + 4.*m_H*nHe[i];
		vel[i]       = rhov[i]/rho_total[i];
		p_total[i]   = (gam - 1)*(e_total[i] - 0.5*rho_total[i]*pow(vel[i],2.));
		cs[i]        = pow(gam*p_total[i]/rho_total[i], 0.5);

		//update free electron and total number densities
		ne[i]        = nH2[i] + nHe2[i] + 2.*nHe3[i];
		n_tot[i]     = nH[i] + nHe[i] + ne[i];

		//apply the closing equation to the ionization state with the largest
		//number density
		if (nH1[i] >= nH2[i])  {
			nH1[i]       = nH[i] - nH2[i];
		}
		else  {
			nH2[i]       = nH[i] - nH1[i];
		}

		if ( (nHe3[i] > nHe1[i]) && (nHe3[i] > nHe2[i]) )  {
			nHe3[i]      = nHe[i] - nHe1[i] - nHe2[i];
		}
		else if ( (nHe2[i] > nHe1[i]) && (nHe2[i] > nHe3[i]) )  {
			nHe2[i]      = nHe[i] - nHe1[i] - nHe3[i];
		}
		else  {
			nHe1[i]      = nHe[i] - nHe2[i] - nHe3[i];
		}

		//update the fluxes
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

	// printf("%le\n",dt/yr_to_s/1e6);

	for (int i{ 0 }; i < N_r; i++)  {
		//update time derivatives
		// printf("%e\n", dt);
		drho_dt[i]  = (rho_total[i] - rho_prev[i])/dt;
		dvel_dt[i]  = (vel[i] - vel_prev[i])/dt;
		de_dt[i]    = (e_total[i] - e_prev[i])/dt;
		dphi_dt[i]  = (phi[i] - phi_prev[i])/dt;


		dnH1_dt[i]  = (nH1[i] - nH1_prev[i])/dt;
		dnH2_dt[i]  = (nH2[i] - nH2_prev[i])/dt;
		dnHe1_dt[i] = (nHe1[i] - nHe1_prev[i])/dt;
		dnHe2_dt[i] = (nHe2[i] - nHe2_prev[i])/dt;
		dnHe3_dt[i] = (nHe3[i] - nHe3_prev[i])/dt;
		dn_dt[i]    = (n_tot[i] - n_tot_prev[i])/dt;
		dne_dt[i]   = (ne[i] - ne_prev[i])/dt;


		//ionization states
		nH_prev[i]    = nH[i];
		nHe_prev[i]   = nHe[i];
		nH1_prev[i]   = nH1[i];
		nH2_prev[i]   = nH2[i];
		nHe1_prev[i]  = nHe1[i];
		nHe2_prev[i]  = nHe2[i];
		nHe3_prev[i]  = nHe3[i];
		ne_prev[i]    = ne[i];
		n_tot_prev[i] = n_tot[i];

		//Update previous time step
		rho_prev[i]  = rho_total[i];
		rhov_prev[i] = rhov[i];
		vel_prev[i]  = vel[i];
		p_prev[i]    = p_total[i];
		e_prev[i]    = e_total[i];
		phi_prev[i]  = phi[i];

		//update abundance fractions
		f_H1[i]  = nH1[i]/nH[i];
		f_H2[i]  = nH2[i]/nH[i];
		f_He1[i] = nHe1[i]/nHe[i];
		f_He2[i] = nHe2[i]/nHe[i];
		f_He3[i] = nHe3[i]/nHe[i];
	}
}
