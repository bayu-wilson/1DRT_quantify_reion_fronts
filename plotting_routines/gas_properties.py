import numpy as np
import matplotlib.pyplot as plt
import rw_functions as rw
import options as opt
import pandas as pd
import matplotlib.ticker as ticker


# home_dir = "/Volumes/Extreme SSD/anson/benchmark_optimize_1d_rt_code/"#os.getcwd()
# sim_tot_times = [10,20,100,200,500,700,900]
sim_tot_times = [100]#[10,20,100,200,500,700]
colors = ["red","orange","gold","green","blue","pink","purple"]
# sim_tot_time=200
rows = 2
cols = 3
fig,ax=plt.subplots(rows,cols,sharex=True)
# fig.set_size_inches(6,6)
lambda_lya_cm = 1.21567e-5
c = 2.998e10
h = 6.626068e-27
pi = np.pi

for i,sim_tot_time in enumerate(sim_tot_times):
    data_path=opt.output_dir+"gas_test_hydro_{}myr.txt".format(str(sim_tot_time))
    radius, rho, p, temp, vel, nh1, fh1, fhe1, fh2, fhe2, fhe3, q_eff_lya = rw.read_gas(data_path)
    # initial_grid = pd.read_csv("../input_files/rho_fluc_line0000.txt",skiprows=2,delim_whitespace=True)
    # initial_nh1 = initial_grid["nH1"]

    #(nh1/fh1) #np.ones_like(temp)*1e-4
    gas_prop_list= [temp/1e4, (nh1/fh1)*1e4, np.log10(fh1), np.log10(fhe1), q_eff_lya, np.log10(q_eff_lya*nh1**2*h*c/lambda_lya_cm/4/pi)]
    gas_prop_list_str= [r"T$_4$ (K)",r"$10^4 n_H$", r"log($f_{HI}$)", r"log($f_{HeI}$)", r"$q^{eff}_{Ly\alpha}$", r"log($j_{Ly\alpha}$)"]
    counter=0

    # ax[0][1].plot(radius,initial_nh1,alpha=0.5, color='gray')


    for j in range(rows):
        for k in range(cols):
            # if gas_prop_list_str[counter] != "q_eff_lya":
            r_relative = radius-radius[0]
            # mask = (r_relative>250)&(r_relative<270)
            mask = (r_relative>0)&(r_relative<300)
            ax[j][k].plot(r_relative[mask],gas_prop_list[counter][mask],color=colors[i],label=str(sim_tot_times[i])+"Myr")
            # ax[j][k].set_xlim(230,242)
            ax[j][k].tick_params(bottom=True, top=True, left=True, right=False, direction="in")
            ax[j][k].xaxis.set_major_locator(ticker.MultipleLocator(100))

            ax[j][k].text(0.3,0.5,gas_prop_list_str[counter],transform=ax[j][k].transAxes)
            counter+=1
            # else:
            #     ax[j][k].plot(radius/box_length,gas_prop_list[counter]*nh1**2*h*c/lambda_lya_cm/4/pi
            #     ,color=colors[i],label=str(sim_tot_times[i])+"Myr")
            #     ax[j][k].text(0.6,0.7,"j_lya",transform=ax[j][k].transAxes)
            #     counter+=1

    # ax[0][1].plot(radius/box_length,gas_prop_list[2],color=colors[2]
    # ax[0][0].set_yscale('log')
    # ax[0][1].set_yscale('log')
    ax[0][1].legend(loc='upper center',ncol=len(sim_tot_times),bbox_to_anchor=(0.5, 1.2))
    # ax[1][0].set_yscale('log')
    # ax[0][2].set_yscale('log')
    # ax[1][2].set_yscale('log')

    ax[1][1].set_xlabel(r"L$_{box}$ kpc")
# fig.savefig(opt.figure_dir+"gas_properties.png",dpi=250)
# fig.savefig("/Users/bayuwilson/Desktop/flucRho_50myr_N53.pdf",dpi=250)
# fig.savefig("/Users/bayuwilson/Desktop/test_adaptiveTime_fiducial.pdf",dpi=250)
plt.subplots_adjust(hspace=0)
# fig.tight_layout()
fig.savefig("/Users/bayuwilson/Desktop/figures_IF/sample_skewers_constT/six_panels.pdf")
# plt.show()
plt.close()
