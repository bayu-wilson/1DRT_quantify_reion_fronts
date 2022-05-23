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
using user_inputs::QSO_spec; //CHRIS 05/16/22

int main()
{
  double start_all, end_all;//to measure performance
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

  make_output(); //from io_funcs.cc, makes new data files

  double t_tot = t_max*yr_to_s*1e6; //t_max (Myr) defined in user_inputs.h. t_tot is in seconds
  while (t < t_tot) //starting at t=0 to t_tot
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

    //I want to add a function where I find the IF location and then set then region behind the IF to be 0.
    if (FALSE) //sophisticated filtering correction from Bolton & Haehnelt 2007
    {          // based on location of photo-ionization equilibrium
      double IF_loc{};
      double testme{};
      IF_loc = interpolate(f_H1, r, 0.5, N_r);
      for (int j=0;r[j]<IF_loc;j++)
      {
        //dne_dt[i]*n_tot[i]*(R-R0)*g_constants::kpc_to_cm/N_r/g_constants::c;
        testme= absd(dne_dt[j]*n_tot[j]*(R-R0)*g_constants::kpc_to_cm/N_r/g_constants::c);
        if (testme<1e-11) // if less than threshold, then photo-io eq. set to zero
        {
          nH1[j] = 0;
          f_H1[j] = 0;
          f_H2[j] = 1.;
          nH1_prev[j] = 0;
          nH2[j] = nH[j];
          nH2_prev[j] = nH[j];
          nHe1[j] = 0;
          nHe1_prev[j] = 0;
          f_He1[j]=0;
        }
      }
    }
    if (FALSE) // crude filtering correction
    {
      int prev_index{0};
      double IF_rear_loc{};
      IF_rear_loc = interpolate(f_H1, r, 0.5, N_r); //cm
      for (int j=prev_index;r[j]<IF_rear_loc-10*g_constants::kpc_to_cm;j++)
      {
        nH1[j] = 0;
        f_H1[j] = 0;
        f_H2[j] = 1.;
        nH1_prev[j] = 0;
        nH2[j] = nH[j];
        nH2_prev[j] = nH[j];
        nHe1[j] = 0;
        nHe1_prev[j] = 0;
        f_He1[j]=0;
        prev_index = j;
      }
    }
    if (FALSE) //only use if temp_ev is FALSE. Creates step function where temp is hot behind IF and cold in front of it
    {
      int prev_index{0};
      double IF_rear_loc{};
      IF_rear_loc = interpolate(f_H1, r, 0.5, N_r);
      for (int j=prev_index;r[j]<IF_rear_loc+1.5*g_constants::kpc_to_cm;j++) // only use if temp ev is OFF
      {
        temp[j] = 2.e4;
        temp_prev[j]=2.e4;
        prev_index = j;
      }
    }

    if (absd(remainder_chris(t, t_tot/(N_output - 1))) < dt) //Absolute value of the remainder between time and other chunk of time
    {
      get_j_lya(); //from io_funcs.cc
      write_otf_bayu();
      update_step(); //from rt_funcs.cc
    }

    loading_bar(t, t_tot); //from io_funcs.cc
    t += dt;

    //update time steps using Chris' conservative method
    update_dt(); //from rt_funcs.cc
    step += 1; //enumerating the time-steps
  }

  loading_bar(t_tot, t_tot); //done!

  write_spectrum(); //from io_funcs.cc, writing the data to files
  bool_initial_gas = FALSE;
  write_gas(bool_initial_gas);

  if (QSO_spec == TRUE)  { //CHRIS 05/16/22: output QSO spectrum to output_files
  calc_mock_QSO_spec();
  }

  end_all   = omp_get_wtime();
  printf("\nProgram complete, real-time: \t%.3e seconds\n",end_all-start_all);
  printf("`cat Logfile.log` for more information\n");
  return 0;
}
