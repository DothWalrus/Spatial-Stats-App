import numpy as np
import csv
import scipy.io as sio
import scipy.spatial
import scipy.stats
import matplotlib.pyplot as plt
import PIL.Image
import PIL.ImageFilter
from math_research.stats_pipeline_package.spatial_randomness_functions import spatial_randomness_functions as srf

class Private(object):
    pass
p = Private()
p.srf = srf()

nuclei_data = sio.loadmat('#34 a2-1 20xthresh0_2.mat')
# print(nuclei_data['Selection'].shape, np.min(nuclei_data['Selection']))

roi_data = nuclei_data['Selection']
roi_data_img = PIL.Image.fromarray(roi_data*255,'L')

mask, border = p.srf.compute_roi_mask_and_border(roi_data_img, 150, 1000)

# roi_array = p.srf.cartesian_to_array(roi_data,(mask.shape))

# np.mean(np.min(cdist(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(border)), axis=1))
# compute_mean_nuclei_center_dist_to_border(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(border))

cell_ep_dist_mean, cell_ep_dist_data = p.srf.compute_mean_nuclei_center_dist_to_border(p.srf.array_to_cartesian(nuclei_data['Nuclei_Dendritic']), p.srf.array_to_cartesian(border))
cell_nuc_dist_mean, cell_nuc_dist_data = p.srf.compute_mean_min_nuclei_center_dist(p.srf.array_to_cartesian(nuclei_data['Nuclei_Dendritic']))

# np.set_printoptions(threshold=np.inf)
# print(data['Nuclei_Centers'])
nuclei_centers = nuclei_data['Nuclei_Dendritic']

# nuclei_cent_img = PIL.Image.fromarray(nuclei_centers*255,'L')
# nuclei_cent_img.show()
#
# mask_img = PIL.Image.fromarray(mask*128,'L')
# mask_img.show()
#
mask_and_nuc = mask*127+nuclei_centers*128
mask_and_nuc_img = PIL.Image.fromarray(mask_and_nuc)
mask_and_nuc_img.show(title='Real Data')

# print("data", cell_dist_data, len(cell_dist_data), cell_dist_mean, max(cell_dist_data))

#for loop is for X^2 test
# x = []
# for i in range(1,100):

sim_ep_dist = p.srf.perform_border_bootstrap_test('#34 a2-1 20xthresh0_2.mat', roi_data_img, 150, 1, 1000, plot_num=0)
sim_nuc_dist = p.srf.perform_nuclei_bootstrap_test('#34 a2-1 20xthresh0_2.mat', roi_data_img, 150, 1, 1000, plot_num=0)

#Real Data
cell_ep_dist_data = np.sort(cell_ep_dist_data, axis=None) #distance to epithelium
cell_nuc_dist_data = np.sort(cell_nuc_dist_data, axis=None) #distance to nearest neighbor

#simulated Data
sim_ep_dist = np.array(sim_ep_dist)
sim_ep_dist = np.sort(sim_ep_dist, axis=None) #distance to epithelium
sim_nuc_dist = np.array(sim_nuc_dist)
sim_nuc_dist = np.sort(sim_nuc_dist, axis=None) #distance to nearest neighbor
# print(sim_ep_dist)

    # #Some X^2 test information
    ep_dist_chi_square_arr =  np.array([cell_ep_dist_data+1, sim_ep_dist+1])
    ep_dist_chi_2, ep_dist_p_val, ep_dist_d_o_f, ep_expected_frequencies = scipy.stats.chi2_contingency(ep_dist_chi_square_arr)
    # print(ep_dist_chi_2, ep_dist_p_val, ep_dist_d_o_f)

    nuc_dist_chi_square_arr =  np.array([cell_nuc_dist_data+1, sim_nuc_dist+1])
    nuc_dist_chi_2, nuc_dist_p_val, nuc_dist_d_o_f, nuc_expected_frequencies = scipy.stats.chi2_contingency(nuc_dist_chi_square_arr)
    # print(nuc_dist_chi_2, nuc_dist_p_val, nuc_dist_d_o_f)
    x.append(nuc_dist_p_val)

    plt.hist(x, label='p_value from X^2 Test of Nuclei to Nearest Neighbor')
    plt.title('Distribution of p-values from X^2 Test of Nuclei to Nearest Neighbor')
    plt.show()

bins = np.linspace(0,125)

plt.hist(cell_ep_dist_data, density=True, bins=bins,label='Real Data', alpha=0.5)
plt.hist(sim_ep_dist, density=True, bins=bins, label='Simulated Data', alpha=0.5)
plt.title('Distance of Each DC Nuclei to Epithelium')
plt.ylabel('Relative Counts')
plt.xlabel('Distance to Epithelium (pxl)')
plt.xlim(0,125)
plt.legend()
plt.show()

plt.hist(cell_nuc_dist_data, density=True, bins=bins,label='Real Data', alpha=0.5)
plt.hist(sim_nuc_dist, density=True, bins=bins, label='Simulated Data', alpha=0.5)
plt.title('Distance of Each DC Nuclei to Nearest Neighbor')
plt.ylabel('Relative Counts')
plt.xlabel('Distance to Nearest Neigbor (pxl)')
plt.xlim(0,125)
plt.legend()
plt.show()

print(np.mean(cell_nuc_dist_data))