import numpy as np
import csv
import scipy.io as sio
import scipy.spatial
# from scipy.misc import toimage
# from scipy.sparse import csr_matrix, csc_matrix
import matplotlib.pyplot as plt
# import matplotlib.image as mpimg
# import seaborn
# import PIL
# from PIL import Image, ImageFilter
import PIL.Image
import PIL.ImageFilter
# from numpy.linalg import norm
from scipy.spatial.distance import cdist
import skimage.measure
# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

def compute_roi_mask_and_border(image_path, intensity_threshold, blob_count_threshold):
    img = PIL.Image.open(image_path)
    img = np.array(np.array(img, dtype='uint8') > intensity_threshold, dtype='uint8') * 255
    blobs_labels = skimage.measure.label(img, connectivity=1.5)
    blob_counts = np.bincount(blobs_labels.reshape((blobs_labels.shape[0] * blobs_labels.shape[1],)))
    # Need to make sure the [1:] is working as an index to filter the cell border
    roi_mask = np.array(np.isin(blobs_labels, np.where(blob_counts > blob_count_threshold)[0][1:]) == False, dtype='uint8')
    border_mask = np.array(PIL.Image.fromarray(np.array(roi_mask, dtype='uint8')).filter(PIL.ImageFilter.FIND_EDGES))
    # Consider returning edge in cartesian format
    return roi_mask, border_mask

def array_to_cartesian(matrix):
    points = []
    return np.array(np.where(matrix)).transpose()

def cartesian_to_array(points, shape):
    array = np.zeros(shape)
    for point in points:
        array[point[0], point[1]] = 1
    return np.array(array, dtype='uint8')

def compute_mean_nuclei_center_dist_to_border(nuclei_centers, border_points):
    distances = scipy.spatial.distance.cdist(nuclei_centers, border_points)
    min_dist = np.min(distances, axis=1) #array of minimum distances to epithelium layer from each nuclei center
    mean_dist = np.mean(min_dist) #mean distance to epithelium from all cells
    return mean_dist, min_dist

# def compute_mean_min_nuclei_center_dist(nuclei_centers):
#     distances = cdist(nuclei_centers, nuclei_centers)
#     distances[distances == 0] = np.inf
#     return np.mean(np.min(distances, axis=1)), np.min(distances, axis=1)

def border_bootstrap_simulation(nuclei_centers_array, roi_mask, border_mask, num_sim):
    # Count number of nuclei_centers
    num_centers = np.sum(nuclei_centers_array > 0)
    # Convert ROI mask to points
    roi_mask_points = array_to_cartesian(roi_mask)
    # Define Border Points From Mask
    border_points = array_to_cartesian(border_mask)
    # Begin Boostrap Procedure
    spatial_statistics = np.zeros(num_sim)
    dist_list = []

    for i in range(num_sim):
        # Randomly sample points from roi_mask_points without replacement
        roi_mask_point_indices = np.random.choice(len(roi_mask_points), size=num_centers, replace=False)
        bootstrap_centers = roi_mask_points[roi_mask_point_indices, :]

        # bootstrap_centers_array = cartesian_to_array(bootstrap_centers, nuclei_centers_array.shape)
        # sim_img = PIL.Image.fromarray(bootstrap_centers*156,'L')
        # print(bootstrap_centers, bootstrap_centers_array.shape)
        # sim_img.show()
        #
        # roi_mask_img = PIL.Image.fromarray(roi_mask*156,'L')
        # roi_mask_img.show()
        #
        # roi_border_img = PIL.Image.fromarray(border_mask*156,'L')
        # roi_border_img.show()
        #
        # sim_and_border = PIL.Image.fromarray(bootstrap_centers_array*256+border_mask*158)
        # sim_and_border.show()
        #
        # sim_and_mask = PIL.Image.fromarray(roi_mask*128+bootstrap_centers_array*127)
        # sim_and_mask.show()

        mean_dist, dist = compute_mean_nuclei_center_dist_to_border(bootstrap_centers, border_points)
        spatial_statistics[i] = mean_dist
        dist_list.append(dist)
    return spatial_statistics, dist_list

# def nuclei_distance_bootstrap_simulation(nuclei_centers_array, roi_mask, num_sim):
#     # Count number of nuclei_centers
#     num_centers = np.sum(nuclei_centers_array > 0)
#     # Convert ROI mask to points
#     roi_mask_points = array_to_cartesian(roi_mask)
#     # Begin Boostrap Procedure
#     spatial_statistics = np.zeros(num_sim)
#     bootstrap_centers_array = []
#     for sim_index in range(num_sim):
#         # Randomly sample points from roi_mask_points without replacement
#         roi_mask_point_indices = np.random.choice(len(roi_mask_points), size=num_centers, replace=False)
#         bootstrap_centers = roi_mask_points[roi_mask_point_indices, :]
#         bootstrap_centers_array[sim_index] = bootstrap_centers
#         spatial_statistics[sim_index] = compute_mean_min_nuclei_center_dist(bootstrap_centers)
#     return spatial_statistics, bootstrap_centers_array

def perform_border_bootstrap_test(data_path, reference_image_path, intensity_threshold, blob_count_threshold, num_sim):
    data = sio.loadmat(data_path)
    nuclei_centers = data['Nuclei_Centers']
    # roi_mask, border_mask = compute_roi_mask_and_border(sample_image_path, intensity_threshold, blob_count_threshold)
    roi_mask, border_mask = compute_roi_mask_and_border(reference_image_path, intensity_threshold, blob_count_threshold)
    bootstrap_results, sim_ep_dist = border_bootstrap_simulation(nuclei_centers, roi_mask, border_mask, num_sim, plot_num = 1)
    # observed_spatial_statistic = compute_mean_nuclei_center_dist_to_border(array_to_cartesian(nuclei_centers), array_to_cartesian(border_mask))

    # plt.hist(sim_ep_dist)
    # plt.show()

    # p_val = np.sum(bootstrap_results > observed_spatial_statistic) / len(bootstrap_results)
    # if p_val == 1.0:
    #     print('Observed is farther from border than observed')
    #     p_val = np.sum(bootstrap_results < observed_spatial_statistic) / len(bootstrap_results)
    # ax = seaborn.distplot(bootstrap_results)
    #plt.axvline(observed_spatial_statistic)
    return sim_ep_dist #p_val, observed_spatial_statistic

# def perform_nuclei_bootstrap_test(data_path, reference_image_path, intensity_threshold, blob_count_threshold, num_sim):
#     data = sio.loadmat(data_path)
#     nuclei_centers = data['Nuclei_Centers']
#     # roi_mask, border_mask = compute_roi_mask_and_border(sample_image_path, intensity_threshold, blob_count_threshold)
#     roi_mask, border_mask = compute_roi_mask_and_border(reference_image_path, intensity_threshold, blob_count_threshold)
#     bootstrap_results = nuclei_distance_bootstrap_simulation(nuclei_centers, roi_mask, border_mask, num_sim)
#     observed_spatial_statistic = compute_mean_min_nuclei_center_dist(array_to_cartesian(nuclei_centers))
#     p_val = np.sum(bootstrap_results > observed_spatial_statistic) / len(bootstrap_results)
#     if p_val == 1.0:
#         print('Observed has higher mean nuclei distance than observed')
#         p_val = np.sum(bootstrap_results < observed_spatial_statistic) / len(bootstrap_results)
#     # ax = seaborn.distplot(bootstrap_results)
#     plt.axvline(observed_spatial_statistic)
#     return p_val, observed_spatial_statistic

# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

data = sio.loadmat('34a2_1_med_0std.mat')

roi_path = 'Selection.csv'
with open(roi_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    # get header from first row
    headers = next(reader)
    # get all the rows as a list
    roi_data = list(reader)
    # transform data into numpy array
    roi_data = np.array(roi_data).astype(int)
    # print(max(roi_data[:,2]))
    roi_data = roi_data[:,:2]

# print(roi_data,roi_data.shape)
# roi_data_img = PIL.Image.fromarray(roi_data,'L')
# roi_data_img.show()




mask, border = compute_roi_mask_and_border('34_a2_1_20x_1_nuclei_theresholded.tif', 150, 1)

roi_array = cartesian_to_array(roi_data,(mask.shape))

# np.mean(np.min(cdist(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(border)), axis=1))
# compute_mean_nuclei_center_dist_to_border(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(border))

cell_dist_mean, cell_dist_data = compute_mean_nuclei_center_dist_to_border(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(border))

# np.set_printoptions(threshold=np.inf)
# # print(data['Nuclei_Centers'])
nuclei_centers = data['Nuclei_Centers']

# nuclei_cent_img = PIL.Image.fromarray(nuclei_centers*255,'L')
# nuclei_cent_img.show()
#
# mask_img = PIL.Image.fromarray(mask*128,'L')
# mask_img.show()
#
# mask_and_nuc = mask*128+nuclei_centers*127
# mask_and_nuc_img = PIL.Image.fromarray(mask_and_nuc)
# mask_and_nuc_img.show()




# print("data", cell_dist_data, len(cell_dist_data), cell_dist_mean, max(cell_dist_data))
#
sim_ep_dist = perform_border_bootstrap_test('34a2_1_med_0std.mat', '34_a2_1_20x_1_nuclei_theresholded.tif', 150, 1000, 1)

cell_dist_data = np.sort(cell_dist_data, axis=None)

sim_ep_dist = np.array(sim_ep_dist)
sim_ep_dist = np.sort(sim_ep_dist, axis=None)

# print(sim_ep_dist)

bins = np.linspace(0,50)

plt.hist(cell_dist_data, density=True, bins=bins,label='Real Data', alpha=0.5)
plt.hist(sim_ep_dist, density=True, bins=bins, label='Simulated Data', alpha=0.5)
plt.title('Distance of DC Nuclei to Epithelium')
plt.ylabel('Relative Counts')
plt.xlabel('Distance to Epithelium (pxl)')
plt.xlim(0,30)
plt.legend()
# plt.show()

##################################################################################################################################################################

# mat_file = askopenfilename(filetypes = ('*.mat'))
# root = Tk()
# root.mat_file =  filedialog.askopenfilename(initialdir = "/home/dothwalrus/Downloads",title = "Select file",filetypes = (("mat files","*.mat"),("all files","*.*")))
# print (root.mat_file)
# data = sio.loadmat(root.mat_file)
#
# root.tif_file = filedialog.askopenfilename(initialdir = "/home/dothwalrus/Downloads",title = "Select file",filetypes = (("tif files","*.tif"),("all files","*.*")))
# test1, test2 = compute_roi_mask_and_border(root.tif_file, 150, 1000)
#
# np.mean(np.min(cdist(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(test2)), axis=1))
# compute_mean_nuclei_center_dist_to_border(array_to_cartesian(data['Nuclei_Centers']), array_to_cartesian(test2))

# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

# paths_dict = {'#34 a2-1': ['34a2_1_med_0std.mat', '#34 a2-1 20x.tif'],
#               '#34 a2-2': ['34a2_2_med_0std.mat', '#34 a2-2 20x.tif'],
#               '#34 a3-2': ['34a3_2_med_0std.mat', '#34 a3-2 20x.tif'],
#               '#34 a5-2': ['34a5_2_med_0std.mat', '#34 a5-2 20x.tif'],
#               '#34 a6-1': ['34a6_1_med_0std.mat', '#34 a6-1 20x.tif'],
#               '#34 a7-1': ['34a7_1_med_0std.mat', '#34 a7-1 20x.tif'],
#               '#34 a9-2': ['34a9_2_med_0std.mat', '#34 a9-2 20x.tif'],
#               '#34 b1-1': ['34b1_1_med_0std.mat', '#34 b1-1 20x.tif'],
#               '#34 b2-1': ['34b2_1_med_0std.mat', '#34 b2-1 20x.tif'],
#               '#34 b3-1': ['34b3_1_med_0std.mat', '#34 b3-1 20x.tif'],
#               '#34 b4-2': ['34b4_2_med_0std.mat', '#34 b4-2 20x.tif'],
#               '#34 b5-2': ['34b5_2_med_0std.mat', '#34 b5-2 20x.tif'],
#               '#34 b7-3': ['34b7_3_med_0std.mat', '#34 b7-3 20x.tif'],
#               '#34 b8-1': ['34b8_1_med_0std.mat', '#34 b8-1 20x-1.tif'],
#               '#34 b9-1': ['34b9_1_med_0std.mat', '#34 b9-1 20x.tif'],
#               '#34 c4-2': ['34c4_2_med_0std.mat', '#34 c4-2 20x.tif']
#               }

# for key in paths_dict.keys():
#     try:
#         data = sio.loadmat('New images March 2019/' + key + '/' + paths_dict[key][0])
#         img = Image.open('New images March 2019/' + key + '/' + paths_dict[key][1])
#     except:
#         print('Had an error opening ' + key)

# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

# sample_data_path = 'New images March 2019/#34 a2-1/34a2_1_med_0std.mat'
# sample_image_path = 'New images March 2019/#34 a2-1/#34 a2-1 20x.tif'
# sample_intensity_threshold = 150
# sample_blob_count_threshold = 1000



# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

# for key in paths_dict.keys():
#     print(key)
#     data_path = 'New images March 2019/' + key + '/' + paths_dict[key][0]
#     reference_image_path = 'New images March 2019/' + key + '/' + paths_dict[key][1]
#     print(perform_border_bootstrap_test(data_path, reference_image_path, 150, 1000, 100))
#
# not_infected_border_dists = [25.260582691751985, 21.708321290575178, 18.90475177182989, 18.3503460194175,
#                              19.195422872995987, 16.64837500141929, 19.539147247238763, 18.985181205617867]
# infected_border_dists = [21.23510223634364, 19.712141898846927, 19.918259141696986, 19.426236210638088,
#                          19.100607167426983, 19.442076821424074, 20.156568303602125, 20.01597163099225]
#
# not_infected_border_dists = [25.260582691751985, 21.708321290575178, 18.90475177182989, 18.3503460194175,
#                              19.195422872995987, 16.64837500141929, 19.539147247238763, 18.985181205617867]
# infected_border_dists = [21.23510223634364, 19.712141898846927, 19.918259141696986, 19.426236210638088,
#                          19.100607167426983, 19.442076821424074, 20.156568303602125, 20.01597163099225]
#
# fig, ax = plt.subplots(dpi=1000)
#
# sns.distplot(not_infected_border_dists, ax=ax, label='Non-Infected')
# sns.distplot(infected_border_dists, ax=ax, label='Infected')
# ax.set_title('Density Plot of Mean Nuclei \n Distance to Cell Edge', fontsize=18)
# ax.set_xlabel('Pixels', fontsize=18)
# ax.legend(fontsize=16)
# ax.set_yticks([])

# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

# for key in paths_dict.keys():
#     print(key)
#     data_path = 'New images March 2019/' + key + '/' + paths_dict[key][0]
#     reference_image_path = 'New images March 2019/' + key + '/' + paths_dict[key][1]
#     print(perform_nuclei_bootstrap_test(data_path, reference_image_path, 150, 1000, 100))

# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

# non_infected_mean_min = [18.11330373, 17.71030505, 18.76769729, 19.00663107, 18.70095769, 17.85960056, 19.09293695,
#                          18.2627691]
# infected_mean_min = [18.61139906, 18.13983507, 18.37084686, 17.53309838, 18.46291703, 17.88457097, 18.36523828,
#                      17.99851943]
#
# fig, ax = plt.subplots(dpi=1000)
#
# sns.distplot(non_infected_mean_min, ax=ax, label='Non-Infected')
# sns.distplot(infected_mean_min, ax=ax, label='Infected')
# ax.set_title('Density Plot of Mean-Min Inter-Nuclei Distance', fontsize=18)
# ax.set_xlabel('Pixels', fontsize=18)
# ax.legend(fontsize=14)
# ax.set_yticks([])
#
# # % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

# np.mean(non_infected_mean_min), np.mean(infected_mean_min)
#
# np.std(non_infected_mean_min), np.std(infected_mean_min)