Purpose of each pbs file

go_varyLum_uniform_density.pbs
Parameter space study using uniform density skewers where a time dependent, plane parallel source luminosity is used to explore different I-front speeds. The input parameters are from parameters_input/input_params/varyLum.txt.  

go_pp_fd.pbs
Parameter space study using inhomogeneous density skewers where a time independent, plane parallel source luminosity is used. Different I-front speeds are explored using a few different source luminosities as well as the density fluctuations (I-fronts slow down through over-densities). The input parameters are from parameters_input/input_params/fd_April_params.txt.

go_resultArrays.pbs
For data post-processing and organization. It is used to organize the on-the-fly results into a single data table.

fd_cell_size_convergence.pbs
For the inhomogeneous density field, cell size convergence test in the appendix.

ud_cell_size_convergence.pbs
For the uniform density field, cell size convergence test in the appendix.

go_quickTest.pbs  
A diagnostic test to see if the code is functioning correctly. Input parameters are fixed. Not exploring parameter space at all. 


