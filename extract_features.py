import numpy as np
from sklearn.decomposition import PCA


# Calculate the average of the list
def calculate_list_avg(lst):
    if len(lst) == 0:
        avg_list = 0.0
    else:
        avg_list = sum(lst) / len(lst)
    return avg_list


# Extract the information for each sample
def extract_msg(mrna_exp_mat, tf_exp_mat, mirna_exp_mat, mrna_id_list, tf_id_list, mirna_id_list,
                mrna_to_mrna_dict, tf_to_mrna_dict, mirna_to_mrna_dict_for_mrna,
                mirna_to_mrna_dict_for_mirna, mirna_to_tf_dict, tf_to_mirna_dict):

    mrna_num = len(mrna_id_list)
    tf_num = len(tf_id_list)
    mirna_num = len(mirna_id_list)
    sample_num = mrna_exp_mat.shape[1]

    mrna_feature_mat = np.zeros((mrna_num, sample_num))
    mrna_to_mrna_feature_mat = np.zeros((mrna_num, sample_num))
    tf_to_mrna_feature_mat = np.zeros((mrna_num, sample_num))
    mirna_to_mrna_feature_mat_for_mrna = np.zeros((mrna_num, sample_num))

    mirna_feature_mat = np.zeros((mirna_num, sample_num))
    mirna_to_mrna_feature_mat_for_mirna = np.zeros((mirna_num, sample_num))
    mirna_to_tf_feature_mat = np.zeros((mirna_num, sample_num))
    tf_to_mirna_feature_mat = np.zeros((mirna_num, sample_num))

    # extract the useful information for each sample
    for sample_index in range(sample_num):
        mrna_index = 0
        mirna_index = 0
        # mRNA/TF/miRNA expression data
        # Format：{ID：exp}
        mrna_id_exp_dict = {}
        tf_id_exp_dict = {}
        mirna_id_exp_dict = {}

        # Read the mRNA expression data save in the dictionary
        for i in range(mrna_num):
            mrna_id = mrna_id_list[i]
            mrna_exp = float(mrna_exp_mat[i][sample_index])
            mrna_id_exp_dict[mrna_id] = mrna_exp
        for i in range(tf_num):
            tf_id = tf_id_list[i]
            tf_exp = float(tf_exp_mat[i][sample_index])
            tf_id_exp_dict[tf_id] = tf_exp
        for i in range(mirna_num):
            mirna_id = mirna_id_list[i]
            mirna_exp = float(mirna_exp_mat[i][sample_index])
            mirna_id_exp_dict[mirna_id] = mirna_exp

        # mRNA feature matrix
        for mrna in mrna_id_list:
            mrna_exp = mrna_id_exp_dict[mrna]
            mrna_to_mrna_exp_list = []
            tf_to_mrna_exp_list = []
            mirna_to_mrna_exp_list_for_mrna = []
            for i in mrna_to_mrna_dict[mrna]:
                mrna_to_mrna_exp_list.append(mrna_id_exp_dict[i])
            for i in tf_to_mrna_dict[mrna]:
                tf_to_mrna_exp_list.append(tf_id_exp_dict[i])
            for i in mirna_to_mrna_dict_for_mrna[mrna]:
                mirna_to_mrna_exp_list_for_mrna.append(mirna_id_exp_dict[i])

            # calculate the average of the list
            avg_mrna_to_mrna_exp = calculate_list_avg(mrna_to_mrna_exp_list)
            avg_tf_to_mrna_exp = calculate_list_avg(tf_to_mrna_exp_list)
            avg_mirna_to_mrna_exp_for_mrna = calculate_list_avg(mirna_to_mrna_exp_list_for_mrna)

            mrna_feature_mat[mrna_index, sample_index] = mrna_exp
            mrna_to_mrna_feature_mat[mrna_index, sample_index] = avg_mrna_to_mrna_exp
            tf_to_mrna_feature_mat[mrna_index, sample_index] = avg_tf_to_mrna_exp
            mirna_to_mrna_feature_mat_for_mrna[mrna_index, sample_index] = avg_mirna_to_mrna_exp_for_mrna

            mrna_index += 1

        # mRNA feature matrix
        for mirna in mirna_id_list:
            mirna_to_mrna_exp_list_for_mirna = []
            mirna_to_tf_exp_list = []
            tf_to_mirna_exp_list = []
            mirna_exp = mirna_id_exp_dict[mirna]
            for i in mirna_to_mrna_dict_for_mirna[mirna]:
                mirna_to_mrna_exp_list_for_mirna.append(mrna_id_exp_dict[i])
            for i in mirna_to_tf_dict[mirna]:
                mirna_to_tf_exp_list.append(tf_id_exp_dict[i])
            for i in tf_to_mirna_dict[mirna]:
                tf_to_mirna_exp_list.append(tf_id_exp_dict[i])

            # calculate the average of the list
            avg_mirna_to_mrna_exp_for_mirna = calculate_list_avg(mirna_to_mrna_exp_list_for_mirna)
            avg_mirna_to_tf_exp = calculate_list_avg(mirna_to_tf_exp_list)
            avg_tf_to_mirna_exp = calculate_list_avg(tf_to_mirna_exp_list)

            mirna_feature_mat[mirna_index, sample_index] = mirna_exp
            mirna_to_mrna_feature_mat_for_mirna[mirna_index, sample_index] = avg_mirna_to_mrna_exp_for_mirna
            mirna_to_tf_feature_mat[mirna_index, sample_index] = avg_mirna_to_tf_exp
            tf_to_mirna_feature_mat[mirna_index, sample_index] = avg_tf_to_mirna_exp

            mirna_index += 1

    return mrna_feature_mat, mrna_to_mrna_feature_mat, tf_to_mrna_feature_mat, mirna_to_mrna_feature_mat_for_mrna, \
           mirna_feature_mat, mirna_to_mrna_feature_mat_for_mirna, mirna_to_tf_feature_mat, tf_to_mirna_feature_mat


# Use PCA to reduce dimension
def get_dim(total_ratio, temp_mat):
    pca = PCA(n_components=total_ratio, svd_solver='full')
    pca.fit_transform(temp_mat)
    main_dim = pca.n_components_
    return main_dim


# Use PCA to reduce dimension
def reduce_dim(dim, temp_mat):
    pca = PCA(n_components=dim)
    reduce_dim_mat = pca.fit_transform(temp_mat)
    return reduce_dim_mat.T
