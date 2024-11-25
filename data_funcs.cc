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
#include <limits> // std::numeric_limits //BAYU: 04/01/24

using g_constants::yr_to_s;
using user_inputs::N_output;
using user_inputs::t_max;
using user_inputs::frontIF_fHI; //BAYU 08/12/22
using user_inputs::backIF_fHI; //BAYU 08/12/22
using g_constants::pi;
using g_constants::c;
using g_constants::chi_He;
using g_constants::lambda_lya_cm;
using g_constants::h;
using g_constants::nu_HI;
using g_constants::lambda_HI_ion_cm;
using g_constants::kpc_to_cm; //BAYU 04/13/23
using g_constants::beta_sigmaHI; //Bayu 4/1/24

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

//Calculate average of ionization front speed //April 21, 2023 we found that xion average is equal to FlexRT method
double* calc_ifront_avg(double x[], double xprev[], double radius[], int n)
{
    	int N_int {60}; //number of ionization fractions in the discrete integral
    	double xion_list[N_int] = {0};
	double xion_iter = 0.001; //start at 5% ion fraction
       	double delta_xion = 0.998/(N_int-1); //step size going from 5% to 95% (90% range)
    	for (int i=0;i<N_int;i++){
        	xion_list[i]=xion_iter;
        	xion_iter+=delta_xion;
   	}	
	//{0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95};
    	//double xion_list[N_int] = {0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.75};
	double vIF_array[N_int] = {0};
	double em_array[N_int] = {0};
        double pos,pos_prev;
	double *ifront = (double*) malloc(sizeof(double) * 2);
	for (int i=0;i<N_int;i++) { // find vIF at every x_ion 
                pos       = interpolate(x, radius, xion_list[i], n);
                pos_prev  = interpolate(xprev, radius, xion_list[i], n);
                vIF_array[i] = (pos - pos_prev)/dt;
		em_array[i]  = interpolate(x, j_lya, xion_list[i], n); 
		//em_prev = interpolate(j_lya_prev, radius, xion_list[i], n);   
	}
	double numerator{0}; //for emissivity average
	double denominator{0}; //for emissivity average
	double vIF_ion{0};
	for (int k=1;k<N_int;k++) { // trapezoidal integral for emissivity average
                numerator += 0.5*(vIF_array[k-1]*em_array[k-1]+vIF_array[k]*em_array[k])*(xion_list[k]-xion_list[k-1]);
		denominator += 0.5*(em_array[k-1]+em_array[k])*(xion_list[k]-xion_list[k-1]);
       		vIF_ion += 0.5*(vIF_array[k-1]+vIF_array[k])*(xion_list[k]-xion_list[k-1]);
	}	
	//printf("Em: %.3e Em2: %.3e \t",em_array[5],em_array[6]);
	//fflush(stdout);
	*(ifront + 0) = vIF_ion/(xion_list[N_int-1]-xion_list[0]);
	*(ifront + 1) = numerator/denominator; //averaged IF speed
	return ifront;
}	

//Calculate ionization front velocity using the same method of FlexRT
double* calc_ifront_FlexRT()//double x[], double xprev[])
{
	double *ifront = (double*) malloc(sizeof(double) * 2);
	double vel_IF{0};
	double temp_IF{0};
	double denominator{0};
        for (int i=equil_index+1; i < N_r; i++) {
		vel_IF+=dr*(f_H1_prev[i]-f_H1[i])/dt; //pkpc/sec
		temp_IF += 0.5*(temp[i-1]+temp[i])*(f_H1[i]-f_H1[i-1]);
		denominator += f_H1[i]-f_H1[i-1];
	}
	temp_IF/=denominator;
	vel_IF*= kpc_to_cm/denominator; // cm/s
	*(ifront + 0) = vel_IF;
	*(ifront + 1) = temp_IF;
        return ifront;
}



//Calculate ionization front incident photon flux, velocity, and neutral hydrogen number density
double* calc_ifront_flux_method()
{
        double I_div_nu_pair[2][N_nu]{}; //will interpolate between pairs
        double F_inc_pair[2]{};
        double F_inc{0};
        //double fHI_pair[2];
        double vel_IF{0}; //vIF found via flux method
        double *otf_outputs = (double*) malloc(sizeof(double) * 8); //array of outputs
        //double three_avg{0}; //average of ne[i]*nH1[i]*q_eff_lya[i]
        double nH_avg{0}; //average over IF
        //double C_em{0}; //intensity/emissivity clumping factor
        int i{0}; //iterating over all 1D RT cells
        int numIF{0}; //number of cells within IF
        
	while ((f_H1[i]<frontIF_fHI)&(i+1<N_r)){ // iterating over the IF, if xHI less than front then...
                if (i>equil_index){ //equil_index tracks the location of photoionization equilibrium behind the IF
                        nH_avg+=nH[i];
                        numIF+=1;
                }
                i+=1;
        }

	for (int j=0;j<N_nu;j++) {
		I_div_nu_pair[0][j] = I_nu[equil_index+1][j]/nu[j];
		I_div_nu_pair[1][j] = I_nu[equil_index+2][j]/nu[j];
	}

	nH_avg/=numIF;
        double nH_center = interpolate(f_H1, nH, 0.5, N_r);
        double T_center =  interpolate(f_H1, temp, 0.5, N_r);
        double width_IF {interpolate(f_H1, r, 0.75, N_r)-interpolate(f_H1, r, 0.25, N_r)};
        double F_lya = 4*pi*lambda_lya_cm/h/c * trapz_int(j_lya,r,N_r);
        F_inc_pair[0] = trapz_int(I_div_nu_pair[0],nu,N_nu)*4*pi/h; //
        F_inc_pair[1] = trapz_int(I_div_nu_pair[1],nu,N_nu)*4*pi/h; //
	F_inc = (F_inc_pair[0]+F_inc_pair[1])/2;

	double gamma_loc =  interpolate(gamma_H1_tot, r, 1.047e-12, N_r);
        vel_IF = (c*F_inc) / (F_inc + c*nH_center*(1+chi_He));
        *(otf_outputs + 0) = F_lya; // photons/s/cm2
        *(otf_outputs + 1) = F_inc; // photons/s/cm2
        *(otf_outputs + 2) = vel_IF; //cm s-1 //flux method
        *(otf_outputs + 3) = T_center; //K
        *(otf_outputs + 4) = width_IF; //cm
        *(otf_outputs + 5) = nH_avg; //g cm-1
        *(otf_outputs + 6) = nH_center; //g cm-1
        *(otf_outputs + 7) = gamma_loc;
	//*(otf_outputs + 7) = C_em;
        return otf_outputs;
}

void get_j_lya() {
	for (int i{ 0 }; i < N_r; i++) {
                q_eff_lya[i] = coll_ex_rate_H1_acc(temp[i]); // cm3/s
                j_lya[i]=q_eff_lya[i]*nH1[i]*ne[i]*h*c/lambda_HI_ion_cm/4/pi; //erg/cm3/sr
	}
}

double sigmaHI_numeric(double nu_iter) {
	// Solves for HI photoionization cross section (cm^-2) as a function of frequency (Hz)
	// Note that sigmaHI0 is ignored
	double beta=beta_sigmaHI;
	return pow(nu_iter/nu_HI,-beta);
}

double avg_sigma_analytic(double alpha) {
	// Analytic solution to the photon-averaged sigmaHI0 
	// Note that sigmaHI0 is ignored
	double beta=beta_sigmaHI;
	double num = alpha * (1 - pow(4, -alpha - beta));
	double den = (alpha + beta) * (1 - pow(4, -alpha));
	return num / den;
}

double f_alpha(double alpha, double numerical_soln) {
	// f(alpha) = 0. this is the function we are finding the zeros of
	return  avg_sigma_analytic(alpha) - numerical_soln;
}

double df_dalpha(double alpha) {
	double beta=beta_sigmaHI;
	double A = (1 - pow(4, -beta - alpha)) / (1 - pow(4, -alpha)) / (beta + alpha);
	double B = alpha * (1 - pow(4, -beta - alpha)) / pow((1 - pow(4, -alpha)), 2) / (beta + alpha);
	double C = pow(4, -alpha) * alpha * log(4) * (1 - pow(4, -beta - alpha)) / pow((1 - pow(4, -alpha)), 2) / (beta + alpha);
	double D = alpha * log(4) * pow(4, -beta - alpha) / (1 - pow(4, -alpha)) / (beta + alpha);
	return A - B - C + D;
}

double pred_alpha(double alpha_estimate, double numerical_soln) {
	// Newton-Raphson method
	// alpha_iter is the iterator 
	// numerical_soln represents the solution for avg_sigma_numerical
	// f(alpha) is the function we are finding the zeros of. f(alpha) = avg_sigma_analytic(alpha) - avg_sigma_numeric
	int num_iter = 1000;
	for (int i = 0; i < num_iter; i++) {
		alpha_estimate = alpha_estimate - f_alpha(alpha_estimate, numerical_soln)/df_dalpha(alpha_estimate);
	}
	return alpha_estimate;
}

double solve_spectral_index(int r_iter) {
	// set analytic and numerical solution for the photon-averaged HI PI x-section equal to eachother and solve for alpha using Newton-Raphson method
	//double sigma_nu[N_nu] = {}; //note that this is divided by sigmaHI0
	// r_iter is the integer spatial iterator. 
	double I_nu_temp[N_nu] = {0};
	double I_nu_times_sigma_nu[N_nu] = {0};
	double avg_sigma_numerical {0};
	double alpha_guess {1.0};
	for (int j{0};j < N_nu;j++){
		I_nu_temp[j] = I_nu[r_iter][j] / nu[j]; //divide by h*nu (convert to photon intensity) but h's cancel out
		I_nu_times_sigma_nu[j] = I_nu_temp[j] * sigmaHI_numeric(nu[j]);
	}
	avg_sigma_numerical = trapz_int(I_nu_times_sigma_nu,nu,N_nu) / trapz_int(I_nu_temp,nu,N_nu); //this is a single number
	return pred_alpha(alpha_guess,avg_sigma_numerical);
}


