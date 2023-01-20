import numpy as np

#IF speed
min_logvIF = 2.75
max_logvIF = 4.25
nbins_logvIF = 6
delta_bin_logvIF = (max_logvIF-min_logvIF)/(nbins_logvIF-1)
bincenters_logvIF = np.linspace(min_logvIF,max_logvIF,nbins_logvIF)
binedges_logvIF = np.linspace(min_logvIF-delta_bin_logvIF/2,max_logvIF+delta_bin_logvIF/2,nbins_logvIF+1)

#Spectral Index
bincenters_alpha = np.loadtxt("../parameters_input/input_params/spectral_indices.txt")
min_alpha = np.min(bincenters_alpha)
max_alpha = np.max(bincenters_alpha)
nbins_alpha = len(bincenters_alpha)
delta_bin_alpha = (max_alpha-min_alpha)/(nbins_alpha-1)
binedges_alpha = np.linspace(min_alpha-delta_bin_alpha/2,max_alpha+delta_bin_alpha/2,nbins_alpha+1)

#min_alpha = -1.0
#max_alpha = 2.5
#nbins_alpha = 8
#delta_bin_alpha = (max_alpha-min_alpha)/(nbins_alpha-1)
#bincenters_alpha = np.linspace(min_alpha,max_alpha,nbins_alpha)
#bincenters_alpha = np.array(["{:.3f}".format(i) for i in bincenters_alpha])
#binedges_alpha = np.linspace(min_alpha-delta_bin_alpha/2,max_alpha+delta_bin_alpha/2,nbins_alpha+1)

#mask = bincenters_alpha ==0
#bincenters_alpha[mask] = 0.001
# binedges_alpha = np.linspace(min_alpha-delta_bin_alpha/2,max_alpha+delta_bin_alpha/2,nbins_alpha+1)
