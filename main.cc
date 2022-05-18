// headers from this project
#include "gas_funcs.h"
#include "general_funcs.h"
#include "global_variables.h"
#include "io_funcs.h"
#include "init_funcs.h"
#include <omp.h>
#include "rt_funcs.h"
#include "user_inputs.h"
// #include "global_constants.h"

// Debugging via plog (3rd party library headers)
#include "plog/Logger.h"
#include "plog/Log.h"
#include "plog/Initializers/RollingFileInitializer.h"

// Standard library headers
#include <stdio.h>
#include <time.h>

using user_inputs::hydro;
using user_inputs::temp_ev;
// using g_constants::c;
// using g_constants::kpc_to_cm;
// using user_inputs::temp_gas_output;

int main()
{
plog::init(plog::debug, "Logfile.log"); //initilizing the logfile
PLOGI.printf("main() called");
  clock_t start, end; //to measure performance
  double start_all, end_all;
	double cpu_time_used;
  int nProcessors=omp_get_num_procs();
  omp_set_num_threads(nProcessors);
  printf("omp_get_num_procs(): %d\n", nProcessors);
  start_all = omp_get_wtime();

  set_dlognu(); // from global_variables.cc, frequency step-size in logspace
  init_grid();  //from init_funcs.cc, generate a linear spatial grid
  init_gas(); //from init_funcs.cc, initialize grid with gas distribution and gas properties like ionization fraction
  bool bool_initial_gas { FALSE }; //only turn on if you need to make a new one
  write_gas(bool_initial_gas); //from io_funcs.cc, saves the initial grid
  init_intensity(); //from init_funcs.cc, initialize the source (probably BB) in the first grid-cell
  update_gamma_nu();  //from rt_funcs.cc, initializes opacity (absorption coefficient)
  //NOTE: why is this function initialized twice???

  make_output(); //from io_funcs.cc, makes new data files

  double t_tot = t_max*yr_to_s*1e6; //t_max (Myr) defined in user_inputs.h. t_tot is in seconds

  // temp[0]=2e4;
PLOGI.printf("------------ BEGINNING FOR LOOP ------------\n");
  while (t < t_tot) //starting at t=0 to t_tot
  {
    // for (int i=0;i<N_r;i++)
    // {
    //   // if (temp[i] != temp[i])
    //   {
    //     printf("%.2e ",temp[i]);
    //     // printf("step: %d dt: %.1e heat_rate:%.3f \t",step,dt,temp[i]);
    //   }
    //   // break;
    // }
PLOGI.printf("step:%03d t:%.1e dt:%.1e",step,t,dt);
    if (FALSE)//(step==10) //choosing one in case something weird happens at the 0'th step
    //only doing this once to see performance of each function in the loop. this is a global variable
    {
start = clock();
      solve_spherical_rt(); //see comments for each function in the else statement
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("solve_spherical_rt(), runtime: %.3e seconds",cpu_time_used);

start = clock();
      update_u_nu();
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_u_nu(), runtime: \t%.3e seconds",cpu_time_used);

start = clock();
      update_gamma();
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_gamma(), runtime: \t%.3e seconds",cpu_time_used);

start = clock();
      update_chem();
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_chem(), runtime: \t%.3e seconds",cpu_time_used);

start = clock();
      if ( (hydro == TRUE) || (temp_ev == TRUE) )
      {
        update_heat_cool();
      }
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_heat_cool(), runtime: %.3e seconds",cpu_time_used);

start = clock();
      if (temp_ev == TRUE)
      {
        update_thermal();
      }
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_thermal(), runtime: \t%.3e seconds",cpu_time_used);

start = clock();
      update_gamma_nu();
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_gamma_nu(), runtime: \t%.3e seconds",cpu_time_used);


      if (absd(remainder_chris(t, t_tot/(N_output - 1))) < dt)
      {
        get_j_lya(); //from io_funcs.cc
        write_otf(); //from io_funcs.cc
        write_otf_bayu();
        update_step(); //from rt_funcs.cc
      }
      loading_bar(t, t_tot);

      t += dt;
start = clock();
      update_dt();
end = clock();
cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
PLOGI.printf("update_dt(), runtime: \t%.3e seconds",cpu_time_used);
      step += 1;
    }
    else
    {
      solve_spherical_rt(); //from rt_funcs.cc, results in the the intensity at each location and frequencies (spectrum) per location
      //from gas_funcs.cc, update the radiation energy density and photoionization rates at each location and frequency
  		update_u_nu(); //energy //from rt_funcs.cc
  		update_gamma(); //photoionization rates //from gas_funcs.cc
  		update_chem(); //chemistry equations //from gas_funcs.cc

      //update heating and cooling functions //hydro is FALSE for Bayu's project, but temp_ev is TRUE
    		if ( (hydro == TRUE) || (temp_ev == TRUE) )
        {
    			update_heat_cool(); //from gas_funcs.cc
    		}

  		if (hydro == TRUE)  //hydrodynamics = False for Bayu's project
      {
  			update_hydro(); //from gas_funcs.cc
  		}

      //update heating/cooling rates and then temperature
      if (temp_ev == TRUE)
      {
        update_thermal(); //from gas_funcs.cc
      }

      //update attenuation coefficient
      update_gamma_nu(); //from rt_funcs.cc
      //NOTE: why is this function called twice?

      //I want to add a function where I find the IF location and then set then region behind the IF to be 0.
      if (FALSE)
    {
      double IF_loc{};
      double testme{};
      IF_loc = interpolate(f_H1, r, 0.5, N_r);
      // printf("%.2e\n", IF_loc);
      for (int j=0;r[j]<IF_loc;j++)
      {
        //dne_dt[i]*n_tot[i]*(R-R0)*g_constants::kpc_to_cm/N_r/g_constants::c;
        testme= absd(dne_dt[j]*n_tot[j]*(R-R0)*g_constants::kpc_to_cm/N_r/g_constants::c);
        if (testme<1e-11) // if less than threshold, then photo-io eq. set to zero
        {
          // printf("%.3e\n", testme);
          nH1[j] = 1.e-10*0;
          f_H1[j] = 1.e-10*0;
          f_H2[j] = 1.;
          nH1_prev[j] = 1.e-10*0;
          nH2[j] = nH[j];
          nH2_prev[j] = nH[j];
          nHe1[j] = 1.e-10*0;
          nHe1_prev[j] = 1.e-10*0;
          f_He1[j]=1e-10*0;
        }
      }
    }
      // FILE *file = NULL;
      // file = fopen("testme.txt", "w");
      // for (int i{0};i<N_r;i++)
      // {
      //   fprintf(file, "%le \t", testme[i]);
      // }
      // fclose(file);
      // printf("%.3e \t",max(testme));
      // printf("%.3e \t",min(testme));
      // printf("%.3e \n",mean(testme));

      if (FALSE)//((step>0)&(TRUE))
      {
        int prev_index{0};
        double IF_rear_loc{};
        IF_rear_loc = interpolate(f_H1, r, 0.5, N_r);
        for (int j=prev_index;r[j]<IF_rear_loc-10*g_constants::kpc_to_cm;j++)
        {
          nH1[j] = 1.e-10*0;
          f_H1[j] = 1.e-10*0;
          f_H2[j] = 1.;
          nH1_prev[j] = 1.e-10*0;
          nH2[j] = nH[j];
          nH2_prev[j] = nH[j];
          nHe1[j] = 1.e-10*0;
          nHe1_prev[j] = 1.e-10*0;
          f_He1[j]=1e-10*0;
          prev_index = j;
          // 	// // f_H2[idx] = 1.-1.e-10;
        	// gamma_H1_tot[j] = 1.e+3;
      	  // gamma_H1_prev[j] = 1.e+3;
        }
      }
      if (FALSE)
      {
        int prev_index{0};
        // for (int j=prev_index;j<N_r-1;j++)
        double IF_rear_loc{};
        IF_rear_loc = interpolate(f_H1, r, 0.5, N_r);
        for (int j=prev_index;r[j]<IF_rear_loc+1.5*g_constants::kpc_to_cm;j++) // only use if temp ev is OFF
        {
          temp[j] = 2.e4;
          temp_prev[j]=2.e4;
          prev_index = j;
        }
      }

      // for (int j{0};j<N_r;j++)
      // {
      //   if (r[j]>IF_rear_loc+10*g_constants::kpc_to_cm)
      //   {
      //     temp[j] = 100;
      //   }
      //   else
      //   {
      //     temp[j] = 2e4;
      //   }
      // }
      // for (int j{0};r[j]<IF_rear_loc+15*g_constants::kpc_to_cm;j++)
      // {
      //   temp[j] = 2e4;
      // }
    // }

      if (absd(remainder_chris(t, t_tot/(N_output - 1))) < dt) //Absolute value of the remainder between time and other chunk of time
      {
        // printf("Writing on-the-fly output\n");
        get_j_lya(); //from io_funcs.cc
        write_otf(); //from io_funcs.cc
        write_otf_bayu();
        update_step(); //from rt_funcs.cc
        // printf("%d\n", step);
      }

      loading_bar(t, t_tot); //from io_funcs.cc
      t += dt;

      //update time steps using Chris' conservative method
      update_dt(); //from rt_funcs.cc
      // printf("step=%d,nHI[0]=%f\n",step,::nH1[0]);
      step += 1; //enumerating the time-steps
    }
  }
  PLOGI.printf("------------ ENDING FOR LOOP ------------\n");

  loading_bar(t_tot, t_tot); //done!

  write_spectrum(); //from io_funcs.cc, writing the data to files
  // printf("%s\n", gas_output);
  bool_initial_gas = FALSE;
  write_gas(bool_initial_gas);

  end_all   = omp_get_wtime();
  // cpu_time_used = ((double) (end_all - start_all)) / CLOCKS_PER_SEC;
  PLOGI.printf("Program complete, real-time: \t%.3e seconds",end_all-start_all);
  printf("\nProgram complete, real-time: \t%.3e seconds\n",end_all-start_all);
  printf("`cat Logfile.log` for more information\n");
  return 0;
}

// make clean
// make
// ./1d_radtransfer.x
