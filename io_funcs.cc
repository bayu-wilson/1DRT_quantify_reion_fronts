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

using namespace user_inputs;
using g_constants::yr_to_s;
using g_constants::kpc_to_cm;
using g_constants::h; //added 01/12/2022
using g_constants::lambda_lya_cm; //added 01/12/2022
using g_constants::lambda_HI_ion_cm; //added 06/07/2022
using g_constants::c; // cm per sec //added 01/12/2022
using g_constants::pi;

//Write gas properties at each location in the skewer. Used only once to save initial properties
void write_gas(char output_path[]){
	FILE *file = NULL;
	file = fopen(output_path, "w");
	for (int i{ 0 }; i < N_r; i++){
		fprintf(file, "%le \t", r[i]/kpc_to_cm);
		fprintf(file, "%le \t", rho_total[i]);
		fprintf(file, "%le \t", ne[i]);
		fprintf(file, "%le \t", f_H1_prev[i]);
		fprintf(file, "%le \t", n_tot[i]);
		fprintf(file, "%le \t", temp[i]);
		fprintf(file, "%le \t", f_H1[i]);
		fprintf(file, "%le \t", f_He1[i]);
		fprintf(file, "%le \t", f_He2[i]);
		fprintf(file, "%le \t", gamma_H1_tot[i]);
		fprintf(file, "%le \n", coll_ex_rate_H1_acc(temp[i]));
	}
	fclose(file);
}

//Read inhomogeneous density fields from Matt Mcquinn's simulations.
void read_grid_mmc(){
	FILE *file = NULL;
	char hd[999]; //for un-needed header information
	double trash{};
	double Delta_cmv{};
	double R_cmv{};
	double T_cmv{};
	std::vector<double> R_mmc{};
	std::vector<double> rho_temp{};
	std::vector<double> rho_mmc{};
	std::vector<double> T_temp{};
	std::vector<double> T_mmc{};
	file = fopen(grid_input, "r");
  fgets(hd,1000,file);
  fgets(hd,1000,file);
  fgets(hd,1000,file);
  fgets(hd,1000,file);
	int i{0};
  while (fscanf(file,"%le %le %le %le %le %le %le %le %le %le\n",&trash,&trash,&trash,&trash,&trash,&trash,&Delta_cmv,&T_cmv,&R_cmv,&trash)!=EOF){
		double R_temp = R_cmv/(1+z)/h_const; //pkpc //converting from R_comoving/h to R_proper (no little h) 
		double buffer {5}; //Skip length 'R_start' pkpc from the front of the skewer (sometimes near halos). Add small 5 pkpc buffer to each side
		if ((R_temp>R_start-buffer)&&(R_temp<R_start+(R-R0+buffer))){ //extract the desired chunk of the skewer. length is R-R0
			R_mmc.push_back(R_temp-R_start);
			rho_temp.push_back(Delta_cmv*rho_0);
			T_temp.push_back(T_cmv);
		}
    i++;
  }

	int size_mmc {static_cast<int>(R_mmc.size())}; //size of extracted skewer
	if (smooth_field){ //this is probably turned off
		for (int i=0;i<size_mmc;i++){
			rho_mmc.push_back(smooth_gaussian(&R_mmc[0],&rho_temp[0],R_mmc[i],sigma_gauss,size_mmc));
			T_mmc.push_back(smooth_gaussian(&R_mmc[0],&T_temp[0],R_mmc[i],sigma_gauss,size_mmc));
		}
	}
	else{
		rho_mmc=rho_temp;
		T_mmc=T_temp;
	}

	//interpolate Matt's skewers onto our grid and then initialize all gas properties
	for (int i=0;i<N_r;i++){
		rho_total[i] = interpolate(&R_mmc[0],&rho_mmc[0],(r[i]/kpc_to_cm-R0),size_mmc);
		rho_prev[i]  = rho_total[i];
		temp[i] = interpolate(&R_mmc[0],&T_mmc[0],(r[i]/kpc_to_cm-R0),size_mmc);
		f_H1[i] = power_law(r[i], r[0], fH1_0, fH1_index);
		if (f_H1[i] > 1.){
			f_H1[i] = fH1_0;
		}
		f_H2[i] = 1. - f_H1[i];
		f_He1[i] = power_law(r[i], r[0], fHe1_0, fHe1_index);
		if (f_He1[i] > 1.){
			f_He1[i] = fHe1_0;
		}
		f_He2[i] = power_law(r[i], r[0], fHe2_0, fHe2_index);
		if (f_He2[i] > 1. - f_He1[i]){
			f_He2[i] = 1. - f_He1[i];
		}
		f_He3[i]  = 1. - f_He1[i] - f_He2[i];
		nH[i]     = (1. - Y)*rho_total[i]/m_H;
		nHe[i]    = Y*rho_total[i]/(4.*m_H);
		nH2[i]    = f_H2[i]*nH[i];
		nHe2[i]   = f_He2[i]*nHe[i];
		nHe3[i]   = f_He3[i]*nHe[i];
		ne[i]     = nH2[i] + nHe2[i] + 2.*nHe3[i];
		n_tot[i]  = nH[i] + nHe[i] + ne[i];
		p_total[i] = k_B*temp[i]*n_tot[i]; //p is pressure
		// //assume no initial radial velocity
		vel[i]       = 0;
	}
	fclose(file);
}

void read_source(){
	FILE *file = NULL;
	char hd[999];

	file = fopen(source_spectrum, "r");
	if (file == NULL){
		printf("Source spectrum not found\n");
	}
	else{
		fscanf(file, "%s %s\n", hd, hd); //skip header
		for (int j{ 0 }; j < N_nu; j++){
			fscanf(file, "%le", &nu[j]);
			fscanf(file, "%le\n", &I_nu[0][j]);
			I_nu_prev[0][j] = I_nu[0][j];
		}
		for (int j{ 0 }; j < N_nu - 1; j++){
			delta_nu[j] = nu[j+1] - nu[j];
		}
	}
}


void write_Inu(char output_string[]) { //this should be done for each skewer.
        FILE *file = NULL;
        file = fopen(output_string, "w");
	for (int i=0; i< N_r; i++) {
		for (int j=0; j<N_nu; j++) {
			fprintf(file, "%le \t", I_nu[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void write_otf_spectrum(char output_string[]){
	FILE *file = NULL;
	file = fopen(output_string, "w");
	for (int j{ 0 }; j < N_nu; j++){
		int index_rear_IF = find_index(f_H1,backIF_fHI,N_r); // BAYU CHANGED FEB. 6 ///backIF_fHI,N_r);
		//int index_rear_IF = find_index(f_H1,0.5,N_r); // boost test July 24, 2023
		if (index_rear_IF<=N_r){
			fprintf(file, "%le \t", nu[j]);
			fprintf(file, "%le \n", I_nu[index_rear_IF][j]);
		}
		else{
			printf("Something is wrong\n");
		}
	}
	fclose(file);
}

void write_otf_fred(char output_string[]){
	FILE *file = NULL;
	file = fopen(output_string, "w");
	double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);
	// double *ifront_He1 = calc_ifront(f_He1, f_He1_step, r, N_r);
	// double *ifront_He3 = calc_ifront(f_He3, f_He3_step, r, N_r);	
	double *FlexRT_avgs = calc_ifront_FlexRT();
	double *vIF_avgs = calc_ifront_avg(f_H1, f_H1_prev, r, N_r); 
	double vIF_avg_xion = *(vIF_avgs+0);
	double vIF_FlexRT = *(FlexRT_avgs+0);
	double temp_avg_IF = *(FlexRT_avgs+1);
	double *otf_outputs  = calc_ifront_flux_method(); //uses flux method from Treion paper (D'Aloisio et al. 2019) to find the incI & velocity
	double F_lya = *(otf_outputs + 0); // photons/s/cm2
        double F_inc = *(otf_outputs + 1); // photons/s/cm2
        //double vIF_H1 = *(otf_outputs + 2); //cm s-1 //flux method
        double T_center = *(otf_outputs + 3); //K
        double width_IF = *(otf_outputs + 4); //cm
        //double nH_avg = *(otf_outputs + 5); //g cm-1
        double nH_center = *(otf_outputs + 6); //g cm-1
        //double gamma_loc = *(otf_outputs + 7);

	fprintf(file, "%d \t", step);
	fprintf(file, "%le \t", t/yr_to_s/1e6); //Myr
	fprintf(file, "%le \t", dt/yr_to_s/1e6); //Myr
	fprintf(file, "%le \t", F_lya); //photons/s/cm2
	fprintf(file, "%le \t", F_inc); //incident photon flux // #/s/cm2
	fprintf(file, "%le \t", *(ifront_H1 + 0)); // IF loc
	fprintf(file, "%le \t", *(ifront_H1 + 1)); //IF speed (finite differencing)
	fprintf(file, "%le \t", vIF_FlexRT); //IF speed (flexRT method)
	fprintf(file, "%le \t", temp_avg_IF); //average temperature in IF
	fprintf(file, "%le \t", T_center);
	fprintf(file, "%le \t", width_IF);
        fprintf(file, "%le \t", nH_center);	
	fprintf(file, "%le \n", vIF_avg_xion); //IF speed over xion // BAYU, 10/30/2024. replacing gamma_lov with avg vIF

	for (int i{ 0 }; i < N_r; i++){
		fprintf(file, "%le \t", r[i]/kpc_to_cm);
		fprintf(file, "%le \t", rho_total[i]);
		fprintf(file, "%le \t", ne[i]);
		fprintf(file, "%le \t", f_H1_prev[i]);
		fprintf(file, "%le \t", n_tot[i]);
		fprintf(file, "%le \t", temp[i]);
		fprintf(file, "%le \t", f_H1[i]);
		fprintf(file, "%le \t", f_He1[i]);
		fprintf(file, "%le \t", f_He2[i]);
		fprintf(file, "%le \t", gamma_H1_tot[i]);
		fprintf(file, "%le \n", coll_ex_rate_H1_acc(temp[i]));
	}
	fclose(file);
}

void loading_bar(double number, double max_number){
	int barWidth = 70;
	float progress = number/max_number;
	std::cout << "[";
	int pos = barWidth * progress;
	for (int i{ 0 }; i < barWidth; ++i){
	  if (i < pos) std::cout << "=";
	  else if (i == pos) std::cout << ">";
	  else std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();
}

