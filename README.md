# 1D Radiative Transfer code

Written originally by Christopher Cain and then modified, improved, and maintained by Bayu Wilson.


## Using the 1D-RT code

### 1) Clone the repository 
```
git clone git@github.com:bayu-wilson/master_1d_rt.git
```
### 2) Gather input data (if necessary)
If you would like density skewers from another simulation, the input data should be placed in the ```input_files/``` directory. The code is set up to read inhomogeneous density skewers used in [Davies et al. 2016](https://arxiv.org/abs/1409.0855). Contact [bwils033@ucr.edu](mailto:bwils033@ucr.edu) to request these skewers.

### 3) Set up user-input parameters
User-input parameters can be edited in ```user_inputs.h```

### 4) Parameter space range to be be explored
The rate of incident ionizing photon productiion and simuluation time are scaled such that the ionization fronts will reach a certain speed over a sufficient length of skewer. This is calculated for each skewer in ```parameters_input/variableSkewer_input_parameters.py```. 

```
cd parameters_input/
module load anaconda3/2020.11
python variableSkewer_input_parameters.py
```
A table of skewer number, ionizing photon rate, and simulation time are located ```input_params/flucRho.txt```.

Additionally, the spectral indices data file is ```input_params/spectral_indicies.txt``` which can be adjusted manually.

### 5) Submit batch job to run the code and explore parameter space
Edit ```go_skewers_1dRT.pbs``` to explore the desired parameter space as well as check user-input values. This code will automatically produce output files located here: ```output_files/gasprops/sk[X]_a=[X]/``` where `[X]` is placeholder for skewer number and spectral index alpha, respectively.
```
cd pbs_scripts/
sbatch go_skewers_1dRT.pbs
```

### 6) On-the-fly outputs
Here we read in the output files to produce a table of on-the-fly outputs (like IF speed and incident flux at the IF location). This is saved in ```results/[otf directory]/otf.csv``` where the otf directory is defined within ```results/array_results.py```. 
```
cd pbs_scripts/
sbatch go_resultArrays.pbs
```

### 7) Plotting routines and interpolation table
We plot the flux ratio parameter space as well as an interpolation table that can be used given alpha and vIF. The interpolation tables are located in ```results/interp_tables/```
```
cd plotting_routines/
module load anaconda3/2020.11
python Frat_parameter_space.py
```

### 8) Other plotting routines
```
python flux_ratio.py
python plot_coll_exc_coef.py
python profileIF.py
python emissivity_dependencies.py
```





