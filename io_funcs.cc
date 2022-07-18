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


void write_gas(char output_path[]){
	FILE *file = NULL;
	// if (bool_initial_gas){
	// 	file = fopen(initial_gas_output, "w");
	// 	fprintf(file, "Gas data\n");
	// 	fprintf(file, "Time: %le\n", t/yr_to_s/1e6);
	// }
	// else{
	// 	file = fopen(gas_output, "a");
	// 	fprintf(file, "Time: %le\n", t/yr_to_s/1e6);
	// }
	file = fopen(output_path, "w");
	fprintf(file, "Time: %le\n", t/yr_to_s/1e6);
	fprintf(file, "radius \t");
	fprintf(file, "rho \t");
	fprintf(file, "ne \t");
	fprintf(file, "ne_prev \t");
	fprintf(file, "dne_dt \t");
	fprintf(file, "n_tot \t");
	fprintf(file, "dt \t");
	fprintf(file, "T \t");
	fprintf(file, "vel \t");
	fprintf(file, "nH1 \t");
	fprintf(file, "fH1 \t");
	fprintf(file, "fHeI \t");
	fprintf(file, "fHII \t");
	fprintf(file, "fHeII \t");
	fprintf(file, "fHeIII \t");
	fprintf(file, "j_lya \n");
	for (int i{ 0 }; i < N_r; i++){
		fprintf(file, "%le \t", r[i]/kpc_to_cm);
		fprintf(file, "%le \t", rho_total[i]);
		fprintf(file, "%le \t", ne[i]);
		fprintf(file, "%le \t", ne_prev[i]);
		fprintf(file, "%le \t", dne_dt[i]);
		fprintf(file, "%le \t", n_tot[i]);
		fprintf(file, "%le \t", dt/yr_to_s/1e6);
		fprintf(file, "%le \t", temp[i]);
		fprintf(file, "%le \t", vel[i]);
		fprintf(file, "%le \t", nH1[i]);
		fprintf(file, "%le \t", f_H1[i]);
		fprintf(file, "%le \t", f_He1[i]);
		fprintf(file, "%le \t", f_H2[i]);
		fprintf(file, "%le \t", f_He2[i]);
		fprintf(file, "%le \t", f_He3[i]);
		fprintf(file, "%le \n", j_lya[i]);//coll_ex_rate_H1_acc(temp[i]));
	}
	fclose(file);
}

void read_grid_mmc(){
	FILE *file = NULL;
	char hd[999]; //for un-needed header information
	double trash{};
	double z_mmc{};
	double Delta_cmv{};
	double R_cmv{};
	double T_cmv{};
	// double R_start{};
	std::vector<double> R_mmc{};
	std::vector<double> rho_temp{};
	std::vector<double> rho_mmc{};
	std::vector<double> T_temp{};
	std::vector<double> T_mmc{};
	file = fopen(grid_input, "r");
  fscanf(file, "%s %s %le\n", hd, hd, &z_mmc); //get the redshift of simulation
  fgets(hd,1000,file);
  fgets(hd,1000,file);
  fgets(hd,1000,file);
  fgets(hd,1000,file);
	int i{0};
  while (fscanf(file,"%le %le %le %le %le %le %le %le %le %le\n",&trash,&trash,&trash,&trash,&trash,&trash,&Delta_cmv,&T_cmv,&R_cmv,&trash)!=EOF){
		double R_temp = R_cmv/(1+z_mmc)/h_const; //pkpc
		double buffer {5};
		if ((R_temp>R_start-buffer)&&(R_temp<R_start+(R-R0+buffer))){
			R_mmc.push_back(R_temp-R_start);
			rho_temp.push_back(Delta_cmv*rho_0);
			T_temp.push_back(T_cmv);
		}
    i++;
  }

	int size_mmc {static_cast<int>(R_mmc.size())};
	if (smooth_field){
		for (int i=0;i<N_r;i++){
			rho_mmc.push_back(smooth_gaussian(&R_mmc[0],&rho_temp[0],R_mmc[i],sigma_gauss,size_mmc));
			T_mmc.push_back(smooth_gaussian(&R_mmc[0],&T_temp[0],R_mmc[i],sigma_gauss,size_mmc));
		}
	}
	else{
		rho_mmc=rho_temp;
		T_mmc=T_temp;
	}

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

void make_output(){ //not being used as of 7/15/22
	FILE *file = NULL;
	// file = fopen(spec_output, "w"); //filenames defined in user_inputs.h
	// fprintf(file, "Source spectrum\n");
	// fclose(file);
	// file = fopen(gas_output, "w");
	// fprintf(file, "Gas data\n");
	// fclose(file);
	// file = fopen(otf_output, "w");
	// fprintf(file, "On-the-fly output\n");
	// fprintf(file, "\n");
	// fclose(file);
	// file = fopen(otf_output_bayu, "w");
	// fprintf(file, "step_number\t");
	// fprintf(file, "time\t");
	// fprintf(file, "timestep\t");
	// fprintf(file, "I_lya\t");
	// fprintf(file, "inc_photon_flux\t");
	// fprintf(file, "nH_boundary\t");
	// fprintf(file, "H1_IF_x\t");
	// fprintf(file, "H1_IF_v\t");
	// fprintf(file, "HI_IF_v_flxMthd\t");
	// fprintf(file, "He1_IF_x\t");
	// fprintf(file, "He1_IF_v\t");
	// fprintf(file, "He3_IF_x\t");
	// fprintf(file, "He3_IF_v\t");
	// fprintf(file, "width_IF\n");
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

void write_otf_spectrum(char output_string[]){
	FILE *file = NULL;
	file = fopen(output_string, "w");
	for (int j{ 0 }; j < N_nu; j++){
		int index_rear_IF = find_index(f_H1,backIF_fHI,N_r);
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

	// double locIF = *(ifront_H1 + 0);  //location of ÎF (50% neutral)
	std::vector<double> r_IF{}; // vector containing the IF locations
	std::vector<double> j_lya_IF{}; // vector containing the IF emissivities
	std::vector<double> nH_IF{}; // vector containing the Hydrogen number densities
	// double frontIF_loc{interpolate(f_H1, r, frontIF_fHI, N_r)};
	// double backIF_loc{interpolate(f_H1, r, backIF_fHI, N_r)};
	// double width_IF_measured{frontIF_loc-backIF_loc};
	for (int i=0;i<N_r;i++){ //populating the vectors using only the defined IF region
		// if ((r[i]>backIF_loc)&&(r[i]<frontIF_loc)) //only use region inside IF //BAYU TEST
		// {
			r_IF.push_back(r[i]);
			j_lya_IF.push_back(j_lya[i]);
			nH_IF.push_back(nH[i]);
		// }
	}
	int numIF {static_cast<int>(r_IF.size())}; //size function of vector doesn't result in an integer directly, "lu" I think
	double I_lya = trapz_int(&j_lya_IF[0],&r_IF[0],numIF); //integrated lya intensity (not flux yet)
	double *iflux_vel  = calc_ifront_flux_method(); //uses flux method from Treion paper (D'Aloisio et al. 2019) to find the incI & velocity
	double inc_photon_flux = *(iflux_vel + 0);
	double vIF_H1 = *(iflux_vel + 1);
	double nH_avg = *(iflux_vel + 2);
	// double nH_avg = average(&nH_IF[0], static_cast<int>(nH_IF.size())); //BAYU CHANGE THIS BACK

	// fprintf(file, "%d \t", step);
	fprintf(file, "%d \t", step);
	fprintf(file, "%le \t", t/yr_to_s/1e6); //Myr
	fprintf(file, "%le \t", dt/yr_to_s/1e6); //Myr
	fprintf(file, "%le \t", I_lya); //lya intensity
	fprintf(file, "%le \t", inc_photon_flux); //incident photon flux
	fprintf(file, "%le \t", nH_avg); // average nH in the IF
	fprintf(file, "%le \t", *(ifront_H1 + 0)); // IF loc
	fprintf(file, "%le \t", *(ifront_H1 + 1)); //IF speed (finite differencing)
	fprintf(file, "%le \t", vIF_H1); //IF speed (flux method)
	// fprintf(file, "%le \n", width_IF_measured); // width of the IF (inner 99%)
	fprintf(file, "%le \n", interpolate(f_H1, temp, 0.5, N_r));
	for (int i{ 0 }; i < N_r; i++){
		fprintf(file, "%le \t", r[i]/kpc_to_cm);
		fprintf(file, "%le \t", rho_total[i]);
		fprintf(file, "%le \t", ne[i]);
		fprintf(file, "%le \t", ne_prev[i]);
		fprintf(file, "%le \t", n_tot[i]);
		fprintf(file, "%le \t", temp[i]);
		fprintf(file, "%le \t", f_H1[i]);
		fprintf(file, "%le \t", f_He1[i]);
		fprintf(file, "%le \t", f_He2[i]);
		fprintf(file, "%le \t", f_He3[i]);
		fprintf(file, "%le \n", coll_ex_rate_H1_acc(temp[i]));
		// fprintf(file, "%le \n", j_lya[i]);//coll_ex_rate_H1_acc(temp[i]));
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

// void read_grid() {
// 	FILE *file = NULL;
// 	char hd[999];
// 	float trash{};
//
// 	file = fopen(grid_input, "r");
// 	if (file == NULL) {
// 		printf("Input grid not found\n");
// 	}
// 	else {
// 		fscanf(file, "%s %s\n", hd, hd); //skip header
// 		fscanf(file, "%s", hd);
// 		fscanf(file, "%le\n", &t); //get time from input file
// 		t *= yr_to_s*1e6;
// 		for (int m{ 0 }; m < 16; m++) { //skip header
// 			fscanf(file, "%s", hd);
// 			// printf("%s\n", hd);
// 		}
// 		for (int i{ 0 }; i < N_r; i++)
//     {
// 			fscanf(file, "%le", &r[i]);
// 			fscanf(file, "%le", &(rho_total[i]));
// 			// fscanf(file, "%le", &p_total[i]);
// 			fscanf(file, "%f", &trash);
// 			fscanf(file, "%f", &trash);
// 			fscanf(file, "%f", &trash);
// 			fscanf(file, "%f", &trash);
// 			fscanf(file, "%f", &trash);
// 			fscanf(file, "%le", &temp[i]);
// 			fscanf(file, "%le", &vel[i]);
// 			fscanf(file, "%le", &nH1[i]);
// 			fscanf(file, "%le", &f_H1[i]);
// 			fscanf(file, "%le", &f_He1[i]);
// 			fscanf(file, "%le", &f_H2[i]);
// 			fscanf(file, "%le", &f_He2[i]);
// 			fscanf(file, "%le", &f_He3[i]);
// 			fscanf(file, "%f\n", &trash);
// 			// printf("%.3f\n", nH1[i]);
// 		}
// 	}
// }

//
// void write_spectrum()
// {
// 	FILE *file = NULL;
// 	file = fopen(spec_output, "a");
// 	for (int j{ 0 }; j < N_nu; j++)
// 	{
// 		int index_rear_IF = find_index(f_H1,backIF_fHI,N_r);
// 		if (index_rear_IF<=N_r)
// 		{
// 			fprintf(file, "%le \t", nu[j]);
// 			fprintf(file, "%le \n", I_nu[index_rear_IF][j]);
// 		}
// 		else {
// 			printf("Something is wrong\n");
// 		}
// 	}
// 	fclose(file);
// }
//
//
// void write_otf()
// {
// 	FILE *file = NULL;
// 	file = fopen(otf_output, "a");
//
// 	fprintf(file, "Step Number    : %d\n", step);
// 	fprintf(file, "Time [Myr]     : %le\n", t/yr_to_s/1e6);
// 	fprintf(file, "Time Step [Myr]: %le\n", dt/yr_to_s/1e6);
// 	fprintf(file, "\n");
//
// 	double fH1_avg_vol   = calc_vol_avg(f_H1, N_r);
// 	double fH2_avg_vol   = calc_vol_avg(f_H2, N_r);
// 	double fHe1_avg_vol  = calc_vol_avg(f_He1, N_r);
// 	double fHe2_avg_vol  = calc_vol_avg(f_He2, N_r);
// 	double fHe3_avg_vol  = calc_vol_avg(f_He3, N_r);
//
// 	double fH1_avg_mass  = calc_mass_avg(f_H1, nH, m_H,  N_r);
// 	double fH2_avg_mass  = calc_mass_avg(f_H2, nH, m_H,  N_r);
// 	double fHe1_avg_mass = calc_mass_avg(f_He1, nHe, 4.*m_H, N_r);
// 	double fHe2_avg_mass = calc_mass_avg(f_He2, nHe, 4.*m_H, N_r);
// 	double fHe3_avg_mass = calc_mass_avg(f_He3, nHe, 4.*m_H, N_r);
//
// 	fprintf(file, "Global Averages\n");
// 	fprintf(file, "fH1  : %le    %le\n", fH1_avg_vol, fH1_avg_mass);
// 	fprintf(file, "fH2  : %le    %le\n", fH2_avg_vol, fH2_avg_mass);
// 	fprintf(file, "fHe1 : %le    %le\n", fHe1_avg_vol, fHe1_avg_mass);
// 	fprintf(file, "fHe2 : %le    %le\n", fHe2_avg_vol, fHe2_avg_mass);
// 	fprintf(file, "fHe3 : %le    %le\n", fHe3_avg_vol, fHe3_avg_mass);
// 	fprintf(file, "\n");
//
// 	double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);
// 	double *ifront_He1 = calc_ifront(f_He1, f_He1_step, r, N_r);
// 	double *ifront_He3 = calc_ifront(f_He3, f_He3_step, r, N_r);
//
// 	fprintf(file, "H1   : %le    %le\n", *(ifront_H1 + 0), *(ifront_H1 + 1));
// 	fprintf(file, "He1  : %le    %le\n", *(ifront_He1 + 0), *(ifront_He1 + 1));
// 	fprintf(file, "He3  : %le    %le\n", *(ifront_He3 + 0), *(ifront_He3 + 1));
// 	fprintf(file, "\n");
// 	fclose(file);
// }
//
// void write_otf_bayu()
// {
// 	FILE *file = NULL;
// 	file = fopen(otf_output_bayu, "a");
//
// 	double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);
// 	double *ifront_He1 = calc_ifront(f_He1, f_He1_step, r, N_r);
// 	double *ifront_He3 = calc_ifront(f_He3, f_He3_step, r, N_r);
//
// 	double locIF = *(ifront_H1 + 0);  //location of ÎF (50% neutral)
// 	// double regionIF = 10 *kpc_to_cm;//(3*2.7+2)*kpc_to_cm; // half-width of region containing all the emissivity for the IF, 5*sigma+2 pkpc
// 	std::vector<double> r_IF{}; // vector containing the IF locations
// 	std::vector<double> j_lya_IF{}; // vector containing the IF emissivities
// 	std::vector<double> nH_IF{}; // vector containing the IF emissivities
// 	for (int i=0;i<N_r;i++) //populating the vectors using only the defined IF region
// 	{
// 		if ((r[i]>(locIF-region_IF*kpc_to_cm))&&(r[i]<(locIF+region_IF*kpc_to_cm))) //only use region inside IF
// 		{
// 			r_IF.push_back(r[i]);
// 			j_lya_IF.push_back(j_lya[i]);
// 			// nH_IF.push_back(nH[i]);
// 		}
// 		if ((r[i]>(locIF-2*kpc_to_cm))&&(r[i]<(locIF+2*kpc_to_cm))) //only use region inside IF
// 		{
// 			nH_IF.push_back(nH[i]);
// 		}
//
// 	}
// 	int numIF {static_cast<int>(r_IF.size())}; //size function of vector doesn't result in an integer directly, "lu" I think
// 	double I_lya = trapz_int(&j_lya_IF[0],&r_IF[0],numIF); //integrated lya intensity (not flux yet)
// 	double *iflux_vel  = calc_ifront_flux_method(); //uses flux method from Treion paper (D'Aloisio et al. 2019) to find the incI & velocity
// 	double inc_photon_flux = *(iflux_vel + 0);
// 	double vIF_H1 = *(iflux_vel + 1);
// 	double nH_avg_temp = average(&nH_IF[0], static_cast<int>(nH_IF.size()));
// 	// double nH_avg_temp = *(iflux_vel + 2); //BAYU 22/05/25
//
// 	double width_IF_measured = interpolate(f_H1, r, 0.9, N_r)-interpolate(f_H1, r, 0.1, N_r); //crude estimation
//
// 	fprintf(file, "%d \t", step);
// 	fprintf(file, "%le \t", t/yr_to_s/1e6); //Myr
// 	fprintf(file, "%le \t", dt/yr_to_s/1e6); //Myr
// 	fprintf(file, "%le \t", I_lya);
// 	fprintf(file, "%le \t", inc_photon_flux);
// 	fprintf(file, "%le \t", nH_avg_temp);
// 	fprintf(file, "%le \t", *(ifront_H1 + 0));
// 	fprintf(file, "%le \t", *(ifront_H1 + 1));
// 	fprintf(file, "%le \t", vIF_H1);
// 	fprintf(file, "%le \t", *(ifront_He1 + 0));
// 	fprintf(file, "%le \t", *(ifront_He1 + 1));
// 	fprintf(file, "%le \t", *(ifront_He3 + 0));
// 	fprintf(file, "%le \t", *(ifront_He3 + 1));
// 	fprintf(file, "%le \t", width_IF_measured);
// 	fprintf(file, "\n");
// 	fclose(file);
// }

// grid_mmc = i-1;
// double rho_spline[]{};
// double T_spline[]{};
// spline(R_mmc,rho_mmc,grid_mmc,1e30,1e30,rho_spline);
// spline(R_mmc,T_mmc,grid_mmc,1e30,1e30,T_spline);

// double skSect{0.515}; //skewer section
// int len = sizeof(R_mmc)/sizeof(R_mmc[0]);
// double R_max{static_cast<double>(max_array(R_mmc,len))}; //max los position
// double R_start{R_max*skSect}; //starting point along the los
// std::vector<double> rho_subset{};
// std::vector<double> T_subset{};
// // double rho_subset[]{}; // density in the desired skewer section
// // double T_subset[]{}; // density in the desired skewer section
//
// // printf("%f \n", R_mmc[-1]);
// for (int j=0;((R_mmc[j]>R_start)&(R_mmc[j]<R_start+(R-R0)));j++){
// // while ((R_mmc[j]>R_start)&(R_mmc[j]<R_start+(R-R0))) {
// 	rho_subset.push_back(rho_mmc[j]);
// 	T_subset.push_back(T_mmc[j]);
// }
// for (int i=0;i<N_r;i++){
// 	rho_total[i] = interpolate(&rho_subset[0],r,R_mmc[i],N_r);
// 	temp[i] = interpolate(&T_subset[0],r,R_mmc[i],N_r);
// }



// std::vector<double> r_IF{}; // vector containing the IF locations &j_lya_IF[0]
// std::vector<double> j_lya_IF{}; // vector containing the IF emissivities
// std::vector<double> nH_IF{}; // vector containing the Hydrogen number densities
// // double frontIF_loc{interpolate(f_H1, r, frontIF_fHI, N_r)};
// // double backIF_loc{interpolate(f_H1, r, backIF_fHI, N_r)};
// // double width_IF_measured{frontIF_loc-backIF_loc};
// for (int i=0;i<N_r;i++) //populating the vectors using only the defined IF region
// {
// 	// if ((r[i]>backIF_loc)&&(r[i]<frontIF_loc)) //only use region inside IF //BAYU TEST
// 	{
// 		r_IF.push_back(r[i]);
// 		j_lya_IF.push_back(j_lya[i]);
// 		nH_IF.push_back(nH[i]);
// 	}
