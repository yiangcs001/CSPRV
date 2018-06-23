import math
import numpy as np


# Calculate the correlation between matrix use RV2
def rv2_corr(temp_x, temp_y):
    X = np.dot(temp_x, temp_x.T)
    _X = np.diag(np.diag(X))
    X_X = X - _X
    Y = np.dot(temp_y, temp_y.T)
    _Y = np.diag(np.diag(Y))
    Y_Y = Y - _Y
    RV2_corr = (np.trace(np.dot(X_X, Y_Y))) / math.sqrt(
        np.dot(np.trace(np.dot(X_X, X_X)), np.trace(np.dot(Y_Y, Y_Y))))
    return RV2_corr


def calculate_corr_mat(temp_mat_1, temp_mat_2, temp_mat_3, temp_mat_4):
    sample_num = temp_mat_1.shape[1]

    sample_feature_dict = {}
    corr_mat = np.zeros((sample_num, sample_num))

    # extract information for each sample
    for sample_index in range(sample_num):
        col_1 = temp_mat_1[:, sample_index].reshape((-1, 1))
        msg_for_each_sample = col_1
        col_2 = temp_mat_2[:, sample_index].reshape((-1, 1))
        msg_for_each_sample = np.hstack((msg_for_each_sample, col_2))
        col_3 = temp_mat_3[:, sample_index].reshape((-1, 1))
        msg_for_each_sample = np.hstack((msg_for_each_sample, col_3))
        col_4 = temp_mat_4[:, sample_index].reshape(-1, 1)
        msg_for_each_sample = np.hstack((msg_for_each_sample, col_4))

        sample_feature_dict[sample_index] = msg_for_each_sample

    for i in range(sample_num):
        for j in range(sample_num):
            x = sample_feature_dict[i]
            y = sample_feature_dict[j]

            corr_mat[i, j] = rv2_corr(x, y)
            corr_mat[j, i] = rv2_corr(x, y)

    return corr_mat
