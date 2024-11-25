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
using user_inputs::save_initial_gas;
using user_inputs::initial_gas_output;
using user_inputs::gas_output;
using user_inputs::correct_hardening;
using user_inputs::alpha;
using user_inputs::skewer;
using user_inputs::Lum;
using user_inputs::otf_dir;

int main() {
  double start_all, end_all; //to measure performance
  char outputstring[300];
  printf("Begin 1dRT. On-the-fly outputs can be found here: /expanse/lustre/projects/uot171/bwils033/master_1d_rt/%s\n", otf_dir);
  start_all = omp_get_wtime();

  set_dlognu(); // from global_variables.cc, frequency step-size in logspace
  init_grid();  //from init_funcs.cc, generate a linear spatial grid
  init_gas(); //from init_funcs.cc, initialize grid with gas distribution and gas properties like ionization fraction
  if (save_initial_gas){
    sprintf(outputstring, initial_gas_output);
    write_gas(outputstring); //from io_funcs.cc, saves the initial grid
  }
  init_intensity(); //from init_funcs.cc, initialize the source in the first grid-cell
  update_gamma_nu();  //from rt_funcs.cc, initializes opacity (absorption coefficient)
  int otf_step = 1; //counting on the fly files starting at 1
  double t_tot = t_max*yr_to_s*1e6; //t_max (Myr) defined in user_inputs.h. t_tot is in seconds

  while (t < t_tot) { //starting at t=0 to t_tot //seconds!
	solve_rt(); //from rt_funcs.cc, results in the the intensity at each location and frequency
	update_u_nu(); //radiation energy density at each location and frequency //from rt_funcs.cc
	update_gamma(); //photoionization rate. calc at each location //from gas_funcs.cc
	update_chem(); //chemistry equations. calc at each location //from gas_funcs.cc

    	//update heating and cooling functions
 	if ( (hydro) || (temp_ev) ) { 
		update_heat_cool(); } //from gas_funcs.cc
    	//update hydrodynamics
	if (hydro) { 
		update_hydro(); } //from gas_funcs.cc

    	//update heating/cooling rates and then temperature
	if (temp_ev) { 
		update_thermal(); } //from gas_funcs.c

   	//update attenuation coefficient
   	update_gamma_nu(); //from rt_funcs.cc
    
    	// on-the-fly outputs
    	if ((absd(remainder_chris(t, t_tot/(N_output - 1))) < dt)) { 	   
		//Absolute value of the remainder between time and other chunk of time
     		get_j_lya(); //from data_funcs.cc

      		sprintf(outputstring, "%s/n%d_gasprops.txt", otf_dir, otf_step); //saving otf outputs in the otf directory
      		write_otf_fred(outputstring);
      		otf_step+=1;
      		update_step(); //from rt_funcs.cc
    	}
	
        ifront_index_prev = ifront_index; //index of I-front location. it is used at the end of program to check how far the I-front progressed	
	ifront_index = max(find_index(f_H1,0.5,N_r),ifront_index_prev);
   	t += dt;
	//update time steps using Chris' conservative method
   	update_dt(); //from rt_funcs.cc
	step += 1; //enumerating the time-steps
  }

  //save final I_nu matrix. The incident spectrum is saved in every grid cell as the I-front passes through
  sprintf(outputstring, "%s/Inu_output.txt", otf_dir);
  write_Inu(outputstring);
 
  sprintf(outputstring, gas_output); // saving the final gas file. this is good for debugging
  write_gas(outputstring); //from io_funcs.cc, writing the data to files

  end_all   = omp_get_wtime();
  
  printf("skewer completion: %.2f \n",static_cast<float>(ifront_index)/N_r);
  printf("equilibrium check: %.2f \n",static_cast<float>(equil_index)/N_r);
  printf("Program complete, real-time: \t%.3e seconds\n\n",end_all-start_all);
  
  return 0;
}

