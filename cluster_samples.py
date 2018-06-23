import numpy as np
from sklearn.cluster import SpectralClustering


# Calculate the final correlation matrix
def calculate_corr_mat(temp_distance_mat):
    temp_mat = temp_distance_mat - np.diag(np.diag(temp_distance_mat))
    row_sum = temp_mat.sum(axis=1).reshape(-1, 1)
    row_mean = row_sum / ((temp_mat.shape[1]) - 1)
    col_sum = temp_mat.sum(axis=0)
    col_mean = col_sum / ((temp_mat.shape[0]) - 1)
    w_mat = np.exp(-(temp_distance_mat ** 2) / (0.3 * (row_mean + col_mean + temp_distance_mat) / 3))
    return w_mat


# Clustering the data with spectral clustering
def do_clustering(cluster_num, mrna_corr_mat, mirna_corr_mat, mrna_corr_weight, sample_id_list):
    mrna_distance_mat = 1 - mrna_corr_mat
    mrna_normal_mat = calculate_corr_mat(mrna_distance_mat)
    mirna_distance_mat = 1 - mirna_corr_mat
    mirna_normal_mat = calculate_corr_mat(mirna_distance_mat)

    a = mrna_corr_weight
    normal_mat = a * mrna_normal_mat + (1 - a) * mirna_normal_mat

    cluster = SpectralClustering(n_clusters=cluster_num, affinity='precomputed', n_init=100)
    cluster.fit(normal_mat)
    predict_label = cluster.labels_

    sample_id_col = np.array(["SampleID"])
    sample_id_col = np.hstack((sample_id_col, sample_id_list))
    clustering_result = sample_id_col.reshape(-1, 1)
    label_col = np.array(["Label"])
    predict_label.astype(str)
    label_col = np.hstack((label_col, predict_label))
    label_col = label_col.reshape(-1, 1)
    clustering_result = np.hstack((clustering_result, label_col))
    return normal_mat, clustering_result
