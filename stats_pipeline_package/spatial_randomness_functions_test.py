import numpy as np
import scipy.io as sio
import scipy.spatial
import matplotlib.pyplot as plt
import PIL.Image
import PIL.ImageFilter
import PIL.ImageOps
from scipy.spatial.distance import cdist
import skimage.measure

# % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

class spatial_randomness_functions(object):

    def compute_roi_mask_and_border(self, img, intensity_threshold, blob_count_threshold):
        # img = PIL.Image.open(image_path)
        img = PIL.ImageOps.invert(img) # images need to be negative for some reason
        img = np.array(np.array(img, dtype='uint8') > intensity_threshold, dtype='uint8') * 255
        blobs_labels = skimage.measure.label(img, connectivity=1.5)
        blob_counts = np.bincount(blobs_labels.reshape((blobs_labels.shape[0] * blobs_labels.shape[1],)))
        # Need to make sure the [1:] is working as an index to filter the cell border
        roi_mask = np.array(np.isin(blobs_labels, np.where(blob_counts > blob_count_threshold)[0][1:]) == False, dtype='uint8')
        border_mask = np.array(PIL.Image.fromarray(np.array(roi_mask, dtype='uint8')).filter(PIL.ImageFilter.FIND_EDGES))
        # Consider returning edge in cartesian format
        return roi_mask, border_mask

    def array_to_cartesian(self, matrix):
        points = []
        return np.array(np.where(matrix)).transpose()

    def cartesian_to_array(self, points, shape):
        array = np.zeros(shape)
        for point in points:
            array[point[0], point[1]] = 1
        return np.array(array, dtype='uint8')

    def compute_mean_nuclei_center_dist_to_border(self, nuclei_centers, border_points):
        distances = scipy.spatial.distance.cdist(nuclei_centers, border_points,'euclidean')
        # print('Total cells',np.sum(nuclei_centers>0))
        # print('border Points', np.size(border_points),border_points.shape)
        # print('nuclei centers',np.size(nuclei_centers),nuclei_centers.shape)
        # distances = np.sort(distances,axis=None)
        # min_dist = np.min(distances)
        min_dist = np.min(distances, axis=1) #array of minimum distances to epithelium layer from each nuclei center
        mean_dist = np.mean(min_dist) #mean distance to epithelium from all cells
        # print('size of array containg all pairwise distances',np.size(distances))
        # print('in compute_meam_nuc_ep_dist: min dist to ep',np.size(min_dist))
        return mean_dist, min_dist

    def compute_mean_min_nuclei_center_dist(self, nuclei_centers):
        distances = cdist(nuclei_centers, nuclei_centers,'euclidean')
        # print('in compute_mean_min_nuclei_center_dist:')
        # print('shape of distances:',distances.shape)
        distances[distances == 0] = np.inf
        min_dist = np.min(distances, axis=1) #distance to closest nuclei
        mean_dist = np.mean(np.min(distances, axis=1)) #average distance to nearest nuclei
        return mean_dist, min_dist

    def border_bootstrap_simulation(self, nuclei_centers_array, roi_mask, border_mask, num_sim, plot_num = 0):
        # Count number of nuclei_centers
        num_centers = np.sum(nuclei_centers_array > 0)
        # Convert ROI mask to points
        roi_mask_points = self.array_to_cartesian(roi_mask)
        # Define Border Points From Mask
        border_points = self.array_to_cartesian(border_mask)
        # Begin Boostrap Procedure
        spatial_statistics = np.zeros(num_sim)
        # cell_count = num_centers // 2 # the division by two is because num_centers is a list of xy indices and thus twice as large
        # dist_list = np.zeros(cell_count*num_sim)
        dist_list = []
        plot_num_array = np.linspace(0,plot_num, num=plot_num)

        for i in range(num_sim):
            # j = i*(i+cell_count)
            # Randomly sample points from roi_mask_points without replacement
            roi_mask_point_indices = np.random.choice(len(roi_mask_points), size=num_centers, replace=False)
            bootstrap_centers = roi_mask_points[roi_mask_point_indices,:]

            if i in plot_num_array:
                bootstrap_centers_array = self.cartesian_to_array(bootstrap_centers, nuclei_centers_array.shape)
                sim_and_mask = PIL.Image.fromarray(roi_mask*128+bootstrap_centers_array*127)
                title = 'Simulated Data Trial {}'.format(i)
                sim_and_mask.show(title=title)

            # print('size bootstrap centers before',np.size(bootstrap_centers))

            mean_dist, dist = self.compute_mean_nuclei_center_dist_to_border(bootstrap_centers, border_points)

            # print('size bootstrap centers after', np.size(bootstrap_centers))
            spatial_statistics[i] = mean_dist
            # print('size of dist', np.size(dist), dist.shape)
            # dist_list[j:j+cell_count+1] = dist
            dist_list.append(dist)
            # print('i = ',i)

        dist_list = np.sort(dist_list,axis=None)
        print('shape of epithelium dist_list', dist_list.shape)
        return spatial_statistics, dist_list

    def nuclei_distance_bootstrap_simulation(self, nuclei_centers_array, roi_mask, num_sim, plot_num=0):
        # Count number of nuclei_centers
        num_centers = np.sum(nuclei_centers_array > 0)
        # Convert ROI mask to points
        roi_mask_points = self.array_to_cartesian(roi_mask)
        # Begin Boostrap Procedure
        spatial_statistics = np.zeros(num_sim)
        # bootstrap_centers_array = np.zeros()
        dist_list = []
        plot_num_array = np.linspace(0, plot_num, num=plot_num)
        for i in range(num_sim):
            # Randomly sample points from roi_mask_points without replacement
            roi_mask_point_indices = np.random.choice(len(roi_mask_points), size=num_centers, replace=False)
            bootstrap_centers = roi_mask_points[roi_mask_point_indices]

            if i in plot_num_array:
                bootstrap_centers_array = self.cartesian_to_array(bootstrap_centers, nuclei_centers_array.shape)
                # sim_and_mask = PIL.Image.fromarray(bootstrap_centers_array * 255) # used for figure generation
                sim_and_mask = PIL.Image.fromarray(roi_mask*128+bootstrap_centers_array*127)
                title = 'Simulated Data Trial {}'.format(i)
                sim_and_mask.show(title=title)

            # bootstrap_centers_array[i] = bootstrap_centers
            mean_dist, dist = self.compute_mean_min_nuclei_center_dist(bootstrap_centers)
            spatial_statistics[i] = mean_dist
            dist_list.append(dist)
        dist_list = np.sort(dist_list,axis=None)
        print('shape of nuclei dist_list', dist_list.shape)
        return spatial_statistics, dist_list #bootstrap_centers_array

    def perform_border_bootstrap_test(self, nuclei_centers_array, reference_image, intensity_threshold, blob_count_threshold, num_sim, plot_num=0):
        # data = sio.loadmat(data_path)
        nuclei_centers = nuclei_centers_array
        roi_mask, border_mask = self.compute_roi_mask_and_border(reference_image, intensity_threshold, blob_count_threshold)
        bootstrap_results, sim_ep_dist = self.border_bootstrap_simulation(nuclei_centers, roi_mask, border_mask, num_sim, plot_num)

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

    def perform_nuclei_bootstrap_test(self, nuclei_centers_array, reference_image,intensity_threshold, blob_count_threshold, num_sim, plot_num=0):
        # data = sio.loadmat(data_path)
        nuclei_centers = nuclei_centers_array
        roi_mask, border_mask = self.compute_roi_mask_and_border(reference_image, intensity_threshold, blob_count_threshold)
        bootstrap_results, sim_nuc_dist = self.nuclei_distance_bootstrap_simulation(nuclei_centers, roi_mask, num_sim, plot_num)

        # observed_spatial_statistic = self.compute_mean_min_nuclei_center_dist(self.array_to_cartesian(nuclei_centers))
        # p_val = np.sum(bootstrap_results > observed_spatial_statistic) / len(bootstrap_results)
        # if p_val == 1.0:
        #     print('Observed has higher mean nuclei distance than observed')
        #     p_val = np.sum(bootstrap_results < observed_spatial_statistic) / len(bootstrap_results)
        # # ax = seaborn.distplot(bootstrap_results)
        # plt.axvline(observed_spatial_statistic)
        return sim_nuc_dist #p_val, observed_spatial_statistic



