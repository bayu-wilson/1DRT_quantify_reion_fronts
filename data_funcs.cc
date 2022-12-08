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
//using namespace std; //CHRIS: 05/16/22

using g_constants::yr_to_s;
using user_inputs::N_output;
using user_inputs::t_max;
//using user_inputs::z; //CHRIS: 05/16/22
//using user_inputs::gas_output; //CHRIS: 05/16/22
//using user_inputs::LyC_opacity;
//using user_inputs::Lyn_series; //CHRIS 06/14/22
//using user_inputs::use_tau_GP; //CHRIS 07/09/22
using user_inputs::frontIF_fHI; //BAYU 08/12/22
using user_inputs::backIF_fHI; //BAYU 08/12/22
//using user_inputs::NLyn;
//using user_inputs::N_wav;
//using user_inputs::wav_start;
//using user_inputs::wav_end;
using g_constants::pi;
using g_constants::c;
using g_constants::chi_He;
using g_constants::lambda_lya_cm;
using g_constants::h;
using g_constants::kpc_to_cm; //CHRIS: 05/16/22
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
        double I_div_nu_pair[2][N_nu]{}; //will interpolate between pairs
        double F_inc_pair[2]{};
        double F_inc{0};
        double fHI_pair[2];
        double vel_IF{0}; //vIF found via flux method
        double *otf_outputs = (double*) malloc(sizeof(double) * 8); //array of outputs
        double three_avg{0}; //average of ne[i]*nH1[i]*q_eff_lya[i]
        double nH_avg{0}; //average over IF
        double C_em{0}; //intensity/emissivity clumping factor
        int i{0}; //iterating over all 1D RT cells
        int numIF{0}; //number of cells within IF
        
	while ((f_H1[i]<frontIF_fHI)&(i+1<N_r)){ // iterating over the IF, if xHI less than front then...
                if (i>prev_index){ //prev_index tracks the location of photoionization equilibrium behind the IF
                        three_avg+=ne[i]*nH1[i]*q_eff_lya[i];
                        nH_avg+=nH[i];
                        numIF+=1;
                        if ((f_H1[i]<backIF_fHI) && (f_H1[i+1]>backIF_fHI)) //find pair where xHI=0.1 (helps find Finc)
                                for (int j=0;j<N_nu;j++)
                                {
                                        I_div_nu_pair[0][j] = I_nu[i][j]/nu[j]; //we track I divided by nu in order to do the integration in the next section
                                        I_div_nu_pair[1][j] = I_nu[i+1][j]/nu[j]; // int (I_nu / h / nu dnu ) = incident Flux ionizing photons
                                        fHI_pair[0] = f_H1[i];
                                        fHI_pair[1] = f_H1[i+1];
                                }
                }
                i+=1;
        }
	nH_avg/=numIF;
        three_avg/=numIF;
        double nH_center = interpolate(f_H1, nH, 0.5, N_r);
        double ne_nH1_center = interpolate(f_H1, ne, 0.5, N_r) * interpolate(f_H1, nH1, 0.5, N_r);
        double T_center =  interpolate(f_H1, temp, 0.5, N_r);
        double q_T_center = coll_ex_rate_H1_acc(T_center);
        double width_IF {interpolate(f_H1, r, frontIF_fHI, N_r)-r[prev_index]};//interpolate(f_H1, r, backIF_fHI, N_r)};
        double F_lya = 4*pi*lambda_lya_cm/h/c * trapz_int(j_lya,r,N_r);
        C_em = three_avg/nH_avg/ne_nH1_center*nH_center/q_T_center;
        F_inc_pair[0] = trapz_int(I_div_nu_pair[0],nu,N_nu)*pi/h; //I'm not sure whether it's multiplied by pi or 4pi
        F_inc_pair[1] = trapz_int(I_div_nu_pair[1],nu,N_nu)*pi/h; //
        F_inc = interpolate(fHI_pair,F_inc_pair,backIF_fHI,2);
        vel_IF = (c*F_inc) / (F_inc + c*nH_center*(1+chi_He));
        *(otf_outputs + 0) = F_lya; // photons/s/cm2
        *(otf_outputs + 1) = F_inc; // photons/s/cm2
        *(otf_outputs + 2) = vel_IF; //cm s-1 //flux method
        *(otf_outputs + 3) = T_center; //K
        *(otf_outputs + 4) = width_IF; //cm
        *(otf_outputs + 5) = nH_avg; //g cm-1
        *(otf_outputs + 6) = nH_center; //g cm-1
        *(otf_outputs + 7) = C_em;
        return otf_outputs;
}

/*
//CHRIS: 05/16/22
//Subroutine to calculate a mock QSO spectrum starting at r = R0
//Optional flag adds a neutral island along the sightline with a specified central position and width
double sigma_HI_LyC(double lam)  {
	double sigma_ind = 2.75;
	double lymanlimit = 911.76;
	if (lam <= lymanlimit)  {
		return g_constants::sigma0*pow(lam, sigma_ind)/pow(lymanlimit, sigma_ind);
	}
	else  {
		return 0.;
	}
}

//This is the lyman series line data for the Lyn calculation.
double* get_lyn_line_data(int series)  {
	double lam0, Lambda_line, a_line, f_line;

	lam0 = 0.;
	Lambda_line = 0.;
	a_line = 0.;
	f_line = 0.;

	double *line_data = (double*) malloc(sizeof(double) * 4);
	switch(series)  {
                case 0: //Lyman-alpha
                lam0 = 1215.67e-8;
                Lambda_line = 6.2649e8;
                a_line = 4.6986e08;
                f_line = 0.41641;
                break;
                case 1: //Lyman-beta
                lam0 = 1025.72e-8;
                Lambda_line = 1.6725e8;
                a_line = 5.5751e7;
                f_line = 0.079142;
                break;
                case 2: //Lyman-gamma
                lam0 = 972.537e-8;
                Lambda_line = 6.8186e7;
                a_line = 1.2785e7;
                f_line = 0.029006;
                break;
                case 3: //Lyman-delta
                lam0 = 949.743e-8;
                Lambda_line = 3.4375e7;
                a_line = 4.1250e6;
                f_line = 0.013945;
                break;
                case 4: //Lyman-epsilon
                lam0 = 937.803e-8;
                Lambda_line = 1.9728e7;
                a_line = 1.6440e6;
                f_line = 7.8035e-3;
                break;
                case 5: //Lyman-zeta
                lam0 = 930.748e-8;
                Lambda_line = 1.1930e7;
                a_line = 7.5684e5;
                f_line = 4.8164e-3;
                break;
                case 6:
		lam0 = 926.2248e-8;
                Lambda_line = 0.;
                a_line = 3.870e+05;
                f_line = 3.185e-03;
                break;
                case 7:
                lam0 = 923.1494e-8;
                Lambda_line = 0.;
                a_line = 2.143e+05;
                f_line = 2.217e-03;
                break;
                case 8:
                lam0 = 920.9621e-8;
                Lambda_line = 0.;
                a_line = 1.263e+05;
                f_line = 1.606e-03;
                break;
                case 9:
                lam0 = 919.3505e-8;
                Lambda_line = 0.;
                a_line = 7.834e+04;
                f_line = 1.201e-03;
                break;
                case 10:
                lam0 = 918.128e-8;
                Lambda_line = 0.;
                a_line = 5.066e+04;
                f_line = 9.220e-04;
                break;
                case 11:
                lam0 = 917.179e-8;
                Lambda_line = 0.;
                a_line = 3.393e+04;
                f_line = 7.231e-04;
                break;
		case 12:
                lam0 = 916.428e-8;
                Lambda_line = 0.;
                a_line = 2.341e+04;
                f_line = 5.777e-04;
                break;
                case 13:
                lam0 = 915.822e-8;
                Lambda_line = 0.;
                a_line = 1.657e+04;
                f_line = 4.689e-04;
                break;
                case 14:
                lam0 = 915.328e-8;
                Lambda_line = 0.;
                a_line = 1.200e+04;
                f_line = 3.858e-04;
                break;
                case 15:
                lam0 = 914.918e-8;
                Lambda_line = 0.;
                a_line = 8.858e+03;
                f_line = 3.213e-04;
                break;
                case 16:
                lam0 = 914.575e-8;
                Lambda_line = 0.;
                a_line = 6.654e+03;
                f_line = 2.704e-04;
                break;
                case 17:
                lam0 = 914.285e-8;
                Lambda_line = 0.;
                a_line = 5.077e+03;
                f_line = 2.297e-04;
                break;
                case 18:
                lam0 = 914.037e-8;
                Lambda_line = 0.;
                a_line = 3.928e+03;
                f_line = 1.968e-04;
                break;
                case 19:
                lam0 = 913.824e-8;
                Lambda_line = 0.;
                a_line = 3.077e+03;
                f_line = 1.699e-04;
                break;
		case 20:
                lam0 = 913.640e-8;
                Lambda_line = 0.;
                a_line = 2.438e+03;
                f_line = 1.477e-04;
                break;
                case 21:
                lam0 = 913.479e-8;
                Lambda_line = 0.;
                a_line = 1.952e+03;
                f_line = 1.292e-04;
                break;
                case 22:
                lam0 = 913.338e-8;
                Lambda_line = 0.;
                a_line = 1.578e+03;
                f_line = 1.136e-04;
                break;
                case 23:
                lam0 = 913.213e-8;
                Lambda_line = 0.;
                a_line = 1.286e+03;
                f_line = 1.005e-04;
                break;
                case 24:
                lam0 = 913.103e-8;
                Lambda_line = 0.;
                a_line = 1.057e+03;
                f_line = 8.933e-05;
                break;
                case 25:
                lam0 = 913.004e-8;
                Lambda_line = 0.;
                a_line = 8.753e+02;
                f_line = 7.974e-05;
                break;
                case 26:
                lam0 = 912.917e-8;
                Lambda_line = 0.;
                a_line = 7.297e+02;
                f_line = 7.148e-05;
                break;
                case 27:
                lam0 = 912.8379e-8;
                Lambda_line = 0.;
                a_line = 6.122e+02;
                f_line = 6.432e-05;
                break;
		case 28:
                lam0 = 912.7667e-8;
                Lambda_line = 0.;
                a_line = 5.168e+02;
                f_line = 5.809e-05;
                break;
                case 29:
                lam0 = 912.7023e-8;
                Lambda_line = 0.;
                a_line = 4.386e+02;
                f_line = 5.264e-05;
                break;
                case 30:
                lam0 = 912.6438e-8;
                Lambda_line = 0.;
                a_line = 3.742e+02;
                f_line = 4.785e-05;
                break;
                case 31:
                lam0 = 912.5905e-8;
                Lambda_line = 0.;
                a_line = 3.208e+02;
                f_line = 4.362e-05;
                break;
                case 32:
                lam0 = 912.54197e-8;
                Lambda_line = 0.;
                a_line = 2.763e+02;
                f_line = 3.988e-05;
                break;
                case 33:
                lam0 = 912.497470e-8;
                Lambda_line = 0.;
                a_line = 2.390e+02;
                f_line = 3.655e-05;
                break;
                case 34:
                lam0 = 912.45663e-8;
                Lambda_line = 0.;
                a_line = 2.076e+02;
                f_line = 3.359e-05;
                break;
                case 35:
                lam0 = 912.41906e-8;
                Lambda_line = 0.;
                a_line = 1.810e+02;
                f_line = 3.093e-05;
                break;
                case 36:
                lam0 = 912.38442e-8;
                Lambda_line = 0.;
                a_line = 1.584e+02;
                f_line = 2.855e-05;
                break;
		case 37:
                lam0 = 912.35329e-8;
                Lambda_line = 0.;
                a_line = 1.391e+02;
                f_line = 2.641e-05;
                break;
                case 38:
                lam0 = 912.32366e-8;
                Lambda_line = 0.;
                a_line = 1.226e+02;
                f_line = 2.448e-05;
                break;
                default:
                printf("Invalid Ly-series line (valid: 0 through 5). Defaulting to Ly-alpha...\n");
                lam0 = 1215.67e-8;
                Lambda_line = 6.2649e8;
                f_line = 0.41641;
	}
	*(line_data + 0) = lam0;
	*(line_data + 1) = Lambda_line;
	*(line_data + 2) = a_line;
	*(line_data + 3) = f_line;

	return line_data;
}

//Calculate the lyman series cross-section using the voigt profile approximation from Tepper-Garcia+13.
double sigma_HI_Lyn(double lam, double T, int series, double NHI, double dv, double dlam, double lam_em, double scale)  {
	double b_doppler = sqrt(2.*g_constants::k_B*T/g_constants::m_H);
	double lam0 = 0., Lambda_line = 0., f_line = 0., a_line = 0.;

	double *line_data = get_lyn_line_data(series);
	lam0 	    = *(line_data + 0);
	Lambda_line = *(line_data + 1);
	a_line      = *(line_data + 2);
	f_line	    = *(line_data + 3);

	if (series >= 6)  {
		Lambda_line = a_line;
	}

	double delta_lam = b_doppler/g_constants::c*lam0;
	double a = lam0*lam0*Lambda_line/(4*pi*g_constants::c*delta_lam);
	double sigma_0 = a*4*pow(pi,1.5)*pow(g_constants::e_c,2.0)/(g_constants::m_e*g_constants::c)*f_line/Lambda_line;
	double x = (lam*1e-8 - lam0)/delta_lam; //this is calculated in velocity space in fred's code, wavelength space here (because easier)
	double Ho = exp(-(x*x));
    	double Q = 1.5/(x*x);
   	double x2 = x*x;
	double f_u = Ho - a/sqrt(pi)/x2*(Ho*Ho*(4*x2*x2+7*x2+4+Q)-Q-1);

	return f_u * sigma_0;
}

//Fluctuating Gunn-Petersen approximation for the Lyman series optical depth.
double tau_HI_Lyn_GP(int series, double a_init, double wav, double scale[N_r], double nHI[N_r])  {
	double lam0 = 0., f_line = 0.;

        double *line_data = get_lyn_line_data(series);
        lam0        = *(line_data + 0);
        f_line      = *(line_data + 3);

	if (wav*1e-8 <= lam0) {
		double a = a_init*lam0/(wav*1e-8);
		double HI_num = interpolate(scale, nHI, a, N_r);
		double tau_GP = pi*pow(g_constants::e_c,2.0)/(g_constants::m_e*g_constants::c)*f_line*lam0/H(a)*HI_num;
		return tau_GP;
	}
	else  {
		return 0.;
	}
}

//Calculate the opacity along the QSO sightline for a given set of wavelengths.
double* tau_bf(double wav_em[N_wav], double v[N_r], double NHI[N_r], double scale[N_r], double nHI[N_r])  {
	double *tau_bf_out = (double*) malloc(sizeof(double) * N_wav);

	#pragma omp parallel for
	for (int i=0; i<N_wav; i++)  {
		*(tau_bf_out + i) = 0.;
		double lam1 = wav_em[i];
		if ((LyC_opacity == TRUE) || ((Lyn_series == TRUE) && (use_tau_GP == FALSE)))  {
			for (int j=0; j<N_r-1; j++)  {
				double dlam = lam1*(v[j+1] - v[j])/g_constants::c;

				//LyC Opacity
				if (LyC_opacity == TRUE)  {
					*(tau_bf_out + i) += NHI[j]*sigma_HI_LyC(lam1 + dlam/2.);
				}

				//Ly-n series opacity if using the voigt profile fits
				if (Lyn_series == TRUE)  {
					if (use_tau_GP == FALSE)  {
						for (int n = 0; n < NLyn; n++)  {
							float column = NHI[j];
							*(tau_bf_out + i) += column*sigma_HI_Lyn(lam1 + dlam/2., temp[j], n, NHI[j], v[j+1]-v[j], dlam, wav_em[i], scale[j]);
						}
					}
				}
				lam1 += dlam;
			}
		}

		//Ly-n series opacity if using the Gunn-Petersen approximation.
		if ((use_tau_GP == TRUE) && (Lyn_series == TRUE))  {
			for (int n = 0; n < NLyn; n++)  {
                                *(tau_bf_out + i) += tau_HI_Lyn_GP(n, scale[0], wav_em[i], scale, nHI);
                        }
		}
	}
	return tau_bf_out;

}

void calc_mock_QSO_spec(double is_length, double is_pos, int name_ind)  {
	//int N_wav = 1000;
	double delta_v     [N_r - 1]{}; //velocity differential
	double wav         [N_wav]{};
	double vel   	   [N_r]{};
	double scale       [N_r]{};
	double flux        [N_wav]{};
	double nH1_spec    [N_r]{};
	double NHI         [N_r-1]{};
	double a = 1./(1.+user_inputs::z);

	//nHI for column density calculation
	#pragma omp parallel for
	for (int i = 0; i<N_r; i++)  {
		nH1_spec[i] = nH1[i];
		if ((r[i] >= kpc_to_cm*(is_pos  - is_length/2.)) && (r[i] <= kpc_to_cm*(is_pos + is_length/2.)))  {
			nH1_spec[i] = nH[i];
		}
	}

	//emitted wavelength spectrum
	wav[0] = wav_start;
	double dwav = (wav_start - wav_end)/N_wav;
	for (int i = 1; i<N_wav; i++)  {
		wav[i] = wav[i-1] - dwav;
	}

	//Sorry this section is super messy.  The evolution of nHI along the line of sight is hard-coded here - which needs to be fixed.
	//The 5.5 scaling should reproduce tau ~ (1+z)^4 assumed in George's paper.
	vel[0] = 0.;
	scale[0] = a;
	for (int i = 1; i<N_r; i++)  {
		delta_v[i-1]   = (r[i] - r[i-1])*H(scale[i-1]); //get velocity bins for hubble flow
		vel[i] 	       = vel[i-1] + delta_v[i-1];
		scale[i]       = scale[i-1]*(1.+delta_v[i-1]/g_constants::c); //get scale factor along the line of sight
		NHI[i-1]       = delta_r[i-1]*(nH1_spec[i] + nH1_spec[i-1])/2.*pow(scale[0]/scale[i],5.5); //account for redshift evolution of nHI
		nH1[i] *= pow(scale[0]/scale[i],5.5);
	}

	double *tau_spec = tau_bf(wav, vel, NHI, scale, nH1); //get the opacity over the spectrum

	//calculate flux
	for (int i = 0; i<N_wav; i++)  {
		double tau = *(tau_spec + i);
		flux[i] = exp(-tau);
	}

	//write file
	ofstream file;
	string s1  = "./output_files/QSO_spec_z=";
	string z_string = to_string(z);
	string ind_string = to_string(name_ind);
	s1 = s1 + z_string + "_" + ind_string;

	if (LyC_opacity == TRUE)  {
		s1 = s1 + "_LyC";
	}

	if (Lyn_series == TRUE)  {
		s1 = s1 + "_Lyn"+to_string(NLyn);
	}

	if (use_tau_GP == TRUE)  {
		s1 = s1 + "_GP";
	}

	const char *spec_output = s1.c_str();
	file.open(spec_output, ios::out | ios::binary);

	for (int i=0; i<N_wav; i++)  {
		file.write((char*)&wav[i],  sizeof(double));
		file.write((char*)&flux[i], sizeof(double));
		printf("wav: %le\n", wav[i]);
		printf("flux: %le\n", flux[i]);
	}
	file.close();
}
*/

void get_j_lya() {
	for (int i{ 0 }; i < N_r; i++) {
                q_eff_lya[i] = coll_ex_rate_H1_acc(temp[i]); // cm3/s
                j_lya[i]=q_eff_lya[i]*nH1[i]*ne[i]*h*c/lambda_HI_ion_cm/4/pi; //erg/cm3/sr
	}
}

