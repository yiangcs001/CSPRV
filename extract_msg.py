import numpy as np
from sklearn import preprocessing


# Extract mRNA/TF/miRNA expression message
def get_exp_msg(exp_file_path):
    temp_mat = np.loadtxt(fname=exp_file_path, dtype=bytes, delimiter='\t').astype(str)
    # id list
    id_list = list(temp_mat[1:, 0])
    # z_score
    exp_mat = temp_mat[1:, 1:].astype(float)
    normal_mat = preprocessing.scale(exp_mat, axis=1)
    return id_list, normal_mat


def get_sample_id_list(file_path):
    temp_mat = np.loadtxt(fname=file_path, dtype=bytes, delimiter='\t').astype(str)
    # sample id list
    id_list = list(temp_mat[0, 1:])
    return id_list


# Extract interaction message
# Format：{ID:[ID1, ID2...]}
def get_interaction_dict(mrna_id_list, tf_id_list, mirna_id_list, network_file_path):
    index = 0
    is_first_line = True
    mrna_to_mrna_dict = {}
    tf_to_mrna_dict = {}
    mirna_to_mrna_dict_for_mrna = {}
    mirna_to_mrna_dict_for_mirna = {}
    mirna_to_tf_dict = {}
    tf_to_mirna_dict = {}

    mrna_num = len(mrna_id_list)
    tf_num = len(tf_id_list)
    mirna_num = len(mirna_id_list)
    id_list = mrna_id_list + tf_id_list + mirna_id_list

    # initialize the dict
    for i in range(mrna_num):
        mrna_to_mrna_dict[id_list[i]] = []
        tf_to_mrna_dict[id_list[i]] = []
        mirna_to_mrna_dict_for_mrna[id_list[i]] = []
    for i in range(mrna_num + tf_num, mrna_num + tf_num + mirna_num):
        mirna_to_mrna_dict_for_mirna[id_list[i]] = []
        mirna_to_tf_dict[id_list[i]] = []
        tf_to_mirna_dict[id_list[i]] = []

    with open(network_file_path) as f:
        for line in f:
            if is_first_line:
                is_first_line = False
            else:
                temp_list = line.strip().split('\t')
                if index < mrna_num:
                    # mRNA : mRNA -> mRNA
                    for i in range(mrna_num):
                        if temp_list[i + 1] == '1':
                            mrna_to_mrna_dict[id_list[i]].append(id_list[index])
                elif index < mrna_num + tf_num:
                    # mRNA : TF -> mRNA
                    for i in range(mrna_num):
                        if temp_list[i + 1] == '1':
                            tf_to_mrna_dict[id_list[i]].append(id_list[index])
                    # miRNA : TF -> miRNA
                    for i in range(mrna_num + tf_num, mrna_num + tf_num + mirna_num):
                        if temp_list[i + 1] == '1':
                            tf_to_mirna_dict[id_list[i]].append(id_list[index])
                else:
                    # mRNA : miRNA -> mRNA
                    for i in range(mrna_num):
                        if temp_list[i + 1] == '1':
                            mirna_to_mrna_dict_for_mrna[id_list[i]].append(id_list[index])
                    # miRNA : miRNA -> mRNA
                    for i in range(mrna_num):
                        if temp_list[i + 1] == '1':
                            mirna_to_mrna_dict_for_mirna[id_list[index]].append(id_list[i])
                    # miRNA : miRNA -> TF
                    for i in range(mrna_num, mrna_num + tf_num):
                        if temp_list[i + 1] == '1':
                            mirna_to_tf_dict[id_list[index]].append(id_list[i])
                index += 1
    f.close()
    return mrna_to_mrna_dict, tf_to_mrna_dict, mirna_to_mrna_dict_for_mrna, \
           mirna_to_mrna_dict_for_mirna, mirna_to_tf_dict, tf_to_mirna_dict
