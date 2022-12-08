// headers from this project
#include "gas_funcs.h"
#include "general_funcs.h"
#include "global_variables.h"
#include "io_funcs.h"
#include "init_funcs.h"
#include <omp.h>
#include "rt_funcs.h"
#include "user_inputs.h"
#include "data_funcs.h" //CHRIS 05/16/22

// Standard library headers
#include <stdio.h>
#include <time.h>

using user_inputs::hydro;
using user_inputs::temp_ev;
using user_inputs::input_grid;
//using user_inputs::QSO_spec; //CHRIS 05/16/22
using user_inputs::save_initial_gas;
using user_inputs::initial_gas_output;
using user_inputs::gas_output;
using user_inputs::correct_hardening;
using user_inputs::alpha;
using user_inputs::skewer;
using user_inputs::dotN;

using g_constants::kpc_to_cm; //debug december 2022

int main() {
  double start_all, end_all; //to measure performance
  int nProcessors=omp_get_num_procs();
  char outputstring[300];
  omp_set_num_threads(nProcessors);
  printf("omp_get_num_procs(): %d\n", nProcessors); //depends on the machine this runs on
  start_all = omp_get_wtime();

  set_dlognu(); // from global_variables.cc, frequency step-size in logspace
  init_grid();  //from init_funcs.cc, generate a linear spatial grid
  init_gas(); //from init_funcs.cc, initialize grid with gas distribution and gas properties like ionization fraction
  if (save_initial_gas){
    sprintf(outputstring, initial_gas_output);
    write_gas(outputstring); //from io_funcs.cc, saves the initial grid
  }
  init_intensity(); //from init_funcs.cc, initialize the source (probably BB) in the first grid-cell
  update_gamma_nu();  //from rt_funcs.cc, initializes opacity (absorption coefficient)
  // make_output(); //from io_funcs.cc, makes new data files
  int otf_step = 1; //counting on the fly files starting at 1
  double t_tot = t_max*yr_to_s*1e6; //t_max (Myr) defined in user_inputs.h. t_tot is in seconds
  
  while (t < t_tot) //starting at t=0 to t_tot
  {
	solve_spherical_rt(); //from rt_funcs.cc, results in the the intensity at each location and frequencies (spectrum) per location
	update_u_nu(); //radiation energy density at each location and frequency //from rt_funcs.cc
	update_gamma(); //photoionization rate. calc at each location //from gas_funcs.cc
	update_chem(); //chemistry equations. calc at each location //from gas_funcs.cc

    //update heating and cooling functions
 	if ( (hydro) || (temp_ev) ) { update_heat_cool(); 
		}//from gas_funcs.cc
    //update hydrodynamics
	if (hydro) { update_hydro(); 
		}	//from gas_funcs.cc

    //update heating/cooling rates and then temperature
	if (temp_ev) { update_thermal(); 
		} //from gas_funcs.c

    //update attenuation coefficient
    update_gamma_nu(); //from rt_funcs.cc
    
    // on-the-fly outputs
    if ((absd(remainder_chris(t, t_tot/(N_output - 1))) < dt)&&(input_grid)) {
//double *ifront_H1  = calc_ifront(f_H1, f_H1_step, r, N_r);//*(ifront_H1 + 0)
	   
//printf("otf step: %d    t_frac: %.3f IF_loc: %.3f\n",otf_step,t/t_tot,(*(ifront_H1 + 0))/kpc_to_cm-R0);
//fflush(stdout);
      //Absolute value of the remainder between time and other chunk of time
      get_j_lya(); //from data_funcs.cc

      //printf("output_files/gasprops/sk%s_a=%.3f/n%d_gasprops.txt\n",skewer,alpha,otf_step);
      sprintf(outputstring, "output_files/gasprops/sk%s_a=%.3f/n%d_gasprops.txt", skewer, alpha, otf_step);
      write_otf_fred(outputstring); //otf
      //sprintf(outputstring, "output_files/incident_spectra/sk%s_a=%.3f/n%d_spectrum.txt", skewer, alpha, otf_step);
      //write_otf_spectrum(outputstring);
      otf_step+=1;
      update_step(); //from rt_funcs.cc
    }
    else if ((absd(remainder_chris(t, t_tot/(N_output - 1))) < dt)&&(!input_grid)){
      get_j_lya(); //from data_funcs.cc
      sprintf(outputstring, "output_files/gasprops/dotN%.5e_a=%.3f/n%d_gasprops.txt", dotN, alpha, otf_step);
      write_otf_fred(outputstring); //otf
      //sprintf(outputstring, "output_files/incident_spectra/rampup_a=%.3f/n%d_spectrum.txt", alpha, otf_step);
      //write_otf_spectrum(outputstring);
      otf_step+=1;
      update_step(); //from rt_funcs.cc
    }

    //loading_bar(t, t_tot); //from io_funcs.cc
    t += dt;
    //update time steps using Chris' conservative method
    update_dt(); //from rt_funcs.cc
    step += 1; //enumerating the time-steps
  }


  //loading_bar(t_tot, t_tot); //done!

  // write_spectrum(); //from io_funcs.cc, writing the data to files
  // bool_initial_gas = FALSE;
  sprintf(outputstring, gas_output);
  write_gas(outputstring); //from io_funcs.cc, writing the data to files

  end_all   = omp_get_wtime();
  printf("\nProgram complete, real-time: \t%.3e seconds\n",end_all-start_all);
  //printf("`cat Logfile.log` for more information\n");
  return 0;
}

