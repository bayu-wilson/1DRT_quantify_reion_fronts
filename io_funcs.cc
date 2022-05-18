#include <stdio.h>
#include "io_funcs.h"
#include "user_inputs.h"
#include "global_constants.h"
#include "global_variables.h"
#include "data_funcs.h"
#include "rates.h"
#include "general_funcs.h"
#include <iostream>
#include <vector>

using user_inputs::grid;
using user_inputs::source_spectrum;
using user_inputs::otf_output;
using user_inputs::otf_output_bayu;
using user_inputs::gas_output;
using user_inputs::initial_gas_output;
using user_inputs::spec_output;
using g_constants::yr_to_s;
using g_constants::kpc_to_cm;
using g_constants::h; //added 01/12/2022
using g_constants::lambda_lya_cm; //added 01/12/2022
using g_constants::c; // cm per sec //added 01/12/2022
using g_constants::pi;

void get_j_lya()
{
	for (int i{ 0 }; i < N_r; i++)
	{
		j_lya[i]=coll_ex_rate_H1_acc(temp[i])*nH1[i]*ne[i]*h*c/lambda_lya_cm/4/pi;
		// printf("%le\n", j_lya[i]);
	}
	// printf("%le\n", trapz_int(j_lya,r,N_r));
}
void read_grid()
{
	// int i, m;
	FILE *file = NULL;
	char hd[999];
	float trash{};

	file = fopen(grid, "r");
	if (file == NULL)
  {
		printf("Input grid not found\n");
		// return -1;
	}
	else
  {
		fscanf(file, "%s %s\n", hd, hd); //skip header
		fscanf(file, "%s", hd);
		fscanf(file, "%le\n", &t); //get time from input file
		t *= yr_to_s*1e6;
		for (int m{ 0 }; m < 12; m++)
    { //skip header
			fscanf(file, "%s", hd);
		}
		// fscanf(file, "%s\n", hd);
		for (int i{ 0 }; i < N_r; i++)
    {
			fscanf(file, "%le", &r[i]);
			fscanf(file, "%le", &(rho_total[i]));
			fscanf(file, "%le", &p_total[i]);
			fscanf(file, "%le", &temp[i]);
			fscanf(file, "%le", &vel[i]);
			fscanf(file, "%le", &nH1[i]);
			fscanf(file, "%le", &f_H1[i]);
			fscanf(file, "%le", &f_He1[i]);
			fscanf(file, "%le", &f_H2[i]);
			fscanf(file, "%le", &f_He2[i]);
			fscanf(file, "%le", &f_He3[i]);
			fscanf(file, "%f\n", &trash);
			// printf("%e\t", r[i]);
		}
	}
}

void read_source()
{

	// int i,j;
	FILE *file = NULL;
	char hd[999];

	file = fopen(source_spectrum, "r");
	if (file == NULL)  {
		printf("Source spectrum not found\n");
	}
	else  {
		fscanf(file, "%s %s\n", hd, hd); //skip header
		for (int j{ 0 }; j < N_nu; j++)  {
			fscanf(file, "%le", &nu[j]);
			fscanf(file, "%le\n", &I_nu[0][j]);
			I_nu_prev[0][j] = I_nu[0][j];
		}
		for (int j{ 0 }; j < N_nu - 1; j++)  {
			delta_nu[j] = nu[j+1] - nu[j];
		}
	}
}

void make_output()
{
	FILE *file = NULL;
	file = fopen(spec_output, "w"); //filenames defined in user_inputs.h
	fprintf(file, "Source spectrum\n");
	fclose(file);

	file = fopen(gas_output, "w");
	fprintf(file, "Gas data\n");
	fclose(file);

	file = fopen(otf_output, "w");
	fprintf(file, "On-the-fly output\n");
	fprintf(file, "\n");
	fclose(file);

	file = fopen(otf_output_bayu, "w");
	fprintf(file, "step_number\t");
	fprintf(file, "time\t");
	fprintf(file, "timestep\t");
	fprintf(file, "I_lya\t");
	fprintf(file, "inc_photon_flux\t");
	fprintf(file, "nH_boundary\t");
	fprintf(file, "H1_IF_x\t");
	fprintf(file, "H1_IF_v\t");
	fprintf(file, "HI_IF_v_flxMthd\t");
	fprintf(file, "He1_IF_x\t");
	fprintf(file, "He1_IF_v\t");
	fprintf(file, "He3_IF_x\t");
	fprintf(file, "He3_IF_v\t");
	fprintf(file, "width_IF\n");
	// fprintf(file, "On-the-fly output\n");
	// fprintf(file, "\n");
	fclose(file);

	// if (bool_initial_gas)
	// {
	// 	file = fopen(initial_gas_output, "w");
	// 	// fprintf(file, "On-the-fly output\n");
	// 	// fprintf(file, "\n");
	// 	fclose(file);
	// }

}

void write_spectrum()
{
	FILE *file = NULL;

	file = fopen(spec_output, "a");
	for (int j{ 0 }; j < N_nu; j++)
	{
		if (index_rear_IF<=N_r)
		{
			fprintf(file, "%le \t", nu[j]);
			fprintf(file, "%le \n", I_nu[index_rear_IF][j]);
		}
		else
		{
			printf("Something is wrong\n");
		}
	}
}

void write_gas(bool bool_initial_gas)
{
	FILE *file = NULL;
	if (bool_initial_gas)	{
		file = fopen(initial_gas_output, "w");
		fprintf(file, "Gas data\n");
		fprintf(file, "Time: %le\n", t/yr_to_s/1e6);
		// file = fopen(initial_gas_output, "a");
	}
	else {
		file = fopen(gas_output, "a");
		fprintf(file, "Time: %le\n", t/yr_to_s/1e6);
	}
	fprintf(file, "radius \t");
	fprintf(file, "rho \t");
	fprintf(file, "P \t");
	fprintf(file, "T \t");
	fprintf(file, "vel \t");
	fprintf(file, "nH1 \t");
	fprintf(file, "fH1 \t");
	fprintf(file, "fHeI \t");
	fprintf(file, "fHII \t");
	fprintf(file, "fHeII \t");
	fprintf(file, "fHeIII \t");
	fprintf(file, "q_eff_lya \n");
	for (int i{ 0 }; i < N_r; i++)  {
		fprintf(file, "%le \t", r[i]/kpc_to_cm);
		fprintf(file, "%le \t", rho_total[i]);
		fprintf(file, "%le \t", p_total[i]);
		fprintf(file, "%le \t", temp[i]);
		fprintf(file, "%le \t", vel[i]);
		fprintf(file, "%le \t", nH1[i]);
		fprintf(file, "%le \t", f_H1[i]);
		fprintf(file, "%le \t", f_He1[i]);
		fprintf(file, "%le \t", f_H2[i]);
		fprintf(file, "%le \t", f_He2[i]);
		fprintf(file, "%le \t", f_He3[i]);
		fprintf(file, "%le \n", coll_ex_rate_H1_acc(temp[i]));
	}
	fclose(file);
}
void write_otf()
{
	FILE *file = NULL;

	file = fopen(otf_output, "a");

	fprintf(file, "Step Number    : %d\n", step);
	fprintf(file, "Time [Myr]     : %le\n", t/yr_to_s/1e6);
	fprintf(file, "Time Step [Myr]: %le\n", dt/yr_to_s/1e6);
	fprintf(file, "\n");

	double fH1_avg_vol   = calc_vol_avg(f_H1, N_r);
	double fH2_avg_vol   = calc_vol_avg(f_H2, N_r);
	double fHe1_avg_vol  = calc_vol_avg(f_He1, N_r);
	double fHe2_avg_vol  = calc_vol_avg(f_He2, N_r);
	double fHe3_avg_vol  = calc_vol_avg(f_He3, N_r);

	double fH1_avg_mass  = calc_mass_avg(f_H1, nH, m_H,  N_r);
	double fH2_avg_mass  = calc_mass_avg(f_H2, nH, m_H,  N_r);
	double fHe1_avg_mass = calc_mass_avg(f_He1, nHe, 4.*m_H, N_r);
	double fHe2_avg_mass = calc_mass_avg(f_He2, nHe, 4.*m_H, N_r);
	double fHe3_avg_mass = calc_mass_avg(f_He3, nHe, 4.*m_H, N_r);

	fprintf(file, "Global Averages\n");
	fprintf(file, "fH1  : %le    %le\n", fH1_avg_vol, fH1_avg_mass);
	fprintf(file, "fH2  : %le    %le\n", fH2_avg_vol, fH2_avg_mass);
	fprintf(file, "fHe1 : %le    %le\n", fHe1_avg_vol, fHe1_avg_mass);
	fprintf(file, "fHe2 : %le    %le\n", fHe2_avg_vol, fHe2_avg_mass);
	fprintf(file, "fHe3 : %le    %le\n", fHe3_avg_vol, fHe3_avg_mass);
	fprintf(file, "\n");

	// printf("computing ifronts\n");
	// double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, t - t_step, N_r);
	// double *ifront_He1 = calc_ifront(f_He1, f_He1_step, r, t - t_step, N_r);
	// double *ifront_He3 = calc_ifront(f_He3, f_He3_step, r, t - t_step, N_r);
	double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);
	double *ifront_He1 = calc_ifront(f_He1, f_He1_step, r, N_r);
	double *ifront_He3 = calc_ifront(f_He3, f_He3_step, r, N_r);

	// printf("writing I-fronts\n");
	// fprintf(file, "I-fronts\n");
	fprintf(file, "H1   : %le    %le\n", *(ifront_H1 + 0), *(ifront_H1 + 1));
	fprintf(file, "He1  : %le    %le\n", *(ifront_He1 + 0), *(ifront_He1 + 1));
	fprintf(file, "He3  : %le    %le\n", *(ifront_He3 + 0), *(ifront_He3 + 1));
	fprintf(file, "\n");

	fclose(file);
}

void write_otf_bayu()
{
	FILE *file = NULL;
	file = fopen(otf_output_bayu, "a");

	double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);
	double *ifront_He1 = calc_ifront(f_He1, f_He1_step, r, N_r);
	double *ifront_He3 = calc_ifront(f_He3, f_He3_step, r, N_r);

	// int IF_index {};
	double locIF = *(ifront_H1 + 0);  //location of ÎF (50% neutral)
	double regionIF = (3*2.7+2)*kpc_to_cm; // half-width of region containing all the emissivity for the IF, 5*sigma+2 pkpc
	std::vector<double> r_IF{}; // vector containing the IF locations
	std::vector<double> j_lya_IF{}; // vector containing the IF emissivities
	for (int i=0;i<N_r;i++) //populating the vectors using only the defined IF region
	{
		if ((r[i]>(locIF-regionIF))&&(r[i]<(locIF+regionIF))) //only use region inside IF
		{
			r_IF.push_back(r[i]);
			j_lya_IF.push_back(j_lya[i]);
			// IF_index = i;
			// printf("%d\n", i);
		}
	}
	int numIF {static_cast<int>(r_IF.size())}; //size function of vector doesn't result in an integer directly, "lu" I think
	// double I_lya = trapz_int(j_lya,r,N_r);
	double I_lya = trapz_int(&j_lya_IF[0],&r_IF[0],numIF); //integrated lya intensity (not flux yet)
	double *iflux_vel  = calc_ifront_flux_method(); //uses flux method from Treion paper (D'Aloisio et al. 2019) to find the incI & velocity
	double inc_photon_flux = *(iflux_vel + 0);
	double vIF_H1 = *(iflux_vel + 1);
	double nH_boundary = *(iflux_vel + 2);

	double width_IF_measured = interpolate(f_H1, r, 0.9, N_r)-interpolate(f_H1, r, 0.1, N_r);

	// remove filtering by setting fHI=0 behind the IF
	// for (int idx{0};idx<IF_index;idx++)
	// {
	// 	// f_H1[idx] = 0.;//1.e-10;
	//
	// 	// nH1[idx] = 0.;
	// 	// nH1_prev[idx] = 0.;
	// 	// nH2[idx] = nH[idx];
	// 	// nH2_prev[idx] = nH[idx];
	// 	// // f_He1[idx] = 1.e-10;
	// 	// nHe1[idx] = 0.;
	// 	// nHe1_prev[idx] = 0.;
	// 	// // f_H2[idx] = 1.-1.e-10;
	// 	gamma_H1_tot[idx] = 1.e-9;
	// 	gamma_H1_prev[idx] = 1.e-9;
	// 	// printf("%.3e\n", gamma_H1_tot[idx]);
	// }
	// if (step==0) //this is the first line
	// {
	// 	fprintf(file, "step_number\t");
	// 	fprintf(file, "time\t");
	// 	fprintf(file, "timestep\t");
	// 	fprintf(file, "I_lya\t");
	// 	fprintf(file, "inc_photon_flux\t");
	// 	fprintf(file, "nH_boundary\t");
	// 	fprintf(file, "H1_IF_x\t");
	// 	fprintf(file, "H1_IF_v\t");
	// 	fprintf(file, "HI_IF_v_flxMthd\t");
	// 	fprintf(file, "He1_IF_x\t");
	// 	fprintf(file, "He1_IF_v\t");
	// 	fprintf(file, "He3_IF_x\t");
	// 	fprintf(file, "He3_IF_v\n");
	// }
	fprintf(file, "%d \t", step);
	fprintf(file, "%le \t", t/yr_to_s/1e6); //Myr
	fprintf(file, "%le \t", dt/yr_to_s/1e6); //Myr
	fprintf(file, "%le \t", I_lya);
	fprintf(file, "%le \t", inc_photon_flux);
	fprintf(file, "%le \t", nH_boundary);
	fprintf(file, "%le \t", *(ifront_H1 + 0));
	fprintf(file, "%le \t", *(ifront_H1 + 1));
	fprintf(file, "%le \t", vIF_H1);
	fprintf(file, "%le \t", *(ifront_He1 + 0));
	fprintf(file, "%le \t", *(ifront_He1 + 1));
	fprintf(file, "%le \t", *(ifront_He3 + 0));
	fprintf(file, "%le \t", *(ifront_He3 + 1));
	fprintf(file, "%le \t", width_IF_measured);
	fprintf(file, "\n");
	fclose(file);
}

void loading_bar(double number, double max_number)
{
      int barWidth = 70;

      float progress = number/max_number;
      std::cout << "[";
      int pos = barWidth * progress;
      for (int i{ 0 }; i < barWidth; ++i) {
          if (i < pos) std::cout << "=";
          else if (i == pos) std::cout << ">";
          else std::cout << " ";
      }
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
}
