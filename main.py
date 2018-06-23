import os
import sys
import getopt
import numpy as np

import extract_msg
import extract_features
import calculate_corr_matrix
import cluster_samples


# Get input data from console
def get_args():
    mrnaFile = ""   # -m : mRNA file path
    tfFile = ""     # -t : TF file path
    mirnaFile = ""  # -r ： miRNA file path
    netFile = ""    # -n : network file path
    saveFile = ""   # -s : save directory

    clusterNum = 0  # -c : cluster number
    pcaRatio = 0    # -p : PCA ratio
    weight = 0      # -w : the weight of mRNA correlation weight

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:t:r:n:s:c:p:w:",
                                   ["help", "mrnafile=", "tffile=", "mirnafile=", "netfile=",
                                    "savefile=", "clusternum=", "pcaratio=", "weight="])
    except getopt.GetoptError:
        print()
        sys.exit()

    for opt, arg in opts:
        if (opt not in ["-h", "--help", "-m", "--mrnafile", "-t", "--tffile", "-r", "--mirnafile", "-n", "--netfile",
                        "-s", "--savefile", "-c", "--clusternum", "-p", "--pcaratio", "-w", "--weight"]):
            print("Error:Input data does not match.")
            sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Try: main.py -m <mRNA File> -t <TF File> -r <miRNA File> -n <Network File> "
                  "-s <Save Directory> -c <Number of Clusters> -p <PCA Ratio> -w <Weight of mRNA>")
            print()
            print("or: main.py --mrnafile <mRNA File> --tffile <TF File> --mirnafile <miRNA File> "
                  "--netfile <Network File> --clusternum <Number of Clusters> "
                  "--savefile <Save File> --pcaratio <PCA Ratio> --weight <Weight of mRNA>")
            sys.exit()
        elif opt in ("-m", "--mranfile"):
            mrnaFile = arg
            continue
        elif opt in ("-t", "--tffile"):
            tfFile = arg
            continue
        elif opt in ("-r", "--mirnafile"):
            mirnaFile = arg
            continue
        elif opt in ("-n", "--netfile"):
            netFile = arg
            continue
        elif opt in ("-s", "--savefile"):
            saveFile = arg
            continue
        elif opt in ("-c", "--clusternum"):
            clusterNum = int(arg)
            continue
        elif opt in ("-p", "--pcaratio"):
            pcaRatio = float(arg)
            continue
        elif opt in ("-w", "--weight"):
            weight = float(arg)
            continue

    return mrnaFile, tfFile, mirnaFile, netFile, saveFile, clusterNum, pcaRatio, weight


def main():
    # get the input data
    mrna_file_path, tf_file_path, mirna_file_path, net_file_path, save_file_path,\
        cluster_num, pca_ratio, weight = get_args()

    print("Calculating...")

    mrna_id_list, mrna_exp_mat = extract_msg.get_exp_msg(mrna_file_path)
    tf_id_list, tf_exp_mat = extract_msg.get_exp_msg(tf_file_path)
    mirna_id_list, mirna_exp_mat = extract_msg.get_exp_msg(mirna_file_path)

    sample_id_list = extract_msg.get_sample_id_list(mrna_file_path)

    mrna_to_mrna_dict, tf_to_mrna_dict, mirna_to_mrna_dict_for_mrna, \
        mirna_to_mrna_dict_for_mirna, mirna_to_tf_dict, tf_to_mirna_dict = \
        extract_msg.get_interaction_dict(mrna_id_list, tf_id_list, mirna_id_list, net_file_path)

    mrna_feature_mat, mrna_to_mrna_feature_mat, tf_to_mrna_feature_mat, mirna_to_mrna_feature_mat_for_mrna, \
        mirna_feature_mat, mirna_to_mrna_feature_mat_for_mirna, mirna_to_tf_feature_mat, tf_to_mirna_feature_mat = \
        extract_features.extract_msg(mrna_exp_mat, tf_exp_mat, mirna_exp_mat,
                                          mrna_id_list, tf_id_list, mirna_id_list,
                                          mrna_to_mrna_dict, tf_to_mrna_dict, mirna_to_mrna_dict_for_mrna,
                                          mirna_to_mrna_dict_for_mirna, mirna_to_tf_dict, tf_to_mirna_dict)

    # PCA for mRNA
    mrna_dim = extract_features.get_dim(pca_ratio, mrna_feature_mat.T)
    mrna_to_mrna_dim = extract_features.get_dim(pca_ratio, mrna_to_mrna_feature_mat.T)
    tf_to_mrna_dim = extract_features.get_dim(pca_ratio, tf_to_mrna_feature_mat.T)
    mirna_to_mrna_dim_for_mrna = extract_features.get_dim(pca_ratio, mirna_to_mrna_feature_mat_for_mrna.T)
    mrna_dim = max([mrna_dim, mrna_to_mrna_dim, tf_to_mrna_dim, mirna_to_mrna_dim_for_mrna])

    pca_mrna_mat = extract_features.reduce_dim(mrna_dim, mrna_feature_mat.T)
    pca_mrna_to_mrna_mat = extract_features.reduce_dim(mrna_dim, mrna_to_mrna_feature_mat.T)
    pca_tf_to_mrna_mat = extract_features.reduce_dim(mrna_dim, tf_to_mrna_feature_mat.T)
    pca_mirna_to_mrna_mat_for_mrna = extract_features.reduce_dim(mrna_dim, mirna_to_mrna_feature_mat_for_mrna.T)

    # PCA for miRNA
    mirna_dim = extract_features.get_dim(pca_ratio, mirna_feature_mat.T)
    mirna_to_mrna_dim_for_mirna = extract_features.get_dim(pca_ratio, mirna_to_mrna_feature_mat_for_mirna.T)
    mirna_to_tf_dim = extract_features.get_dim(pca_ratio, mirna_to_tf_feature_mat.T)
    tf_to_mirna_dim = extract_features.get_dim(pca_ratio, tf_to_mirna_feature_mat.T)

    mirna_dim = max([mirna_dim, mirna_to_mrna_dim_for_mirna, mirna_to_tf_dim, tf_to_mirna_dim])
    pca_mirna_mat = extract_features.reduce_dim(mirna_dim, mirna_feature_mat.T)
    pca_mirna_to_mrna_mat_for_mirna = extract_features.reduce_dim(mirna_dim, mirna_to_mrna_feature_mat_for_mirna.T)
    pca_mirna_to_tf_mat = extract_features.reduce_dim(mirna_dim, mirna_to_tf_feature_mat.T)
    pca_tf_to_mirna_mat = extract_features.reduce_dim(mirna_dim, tf_to_mirna_feature_mat.T)

    mrna_corr_mat = calculate_corr_matrix.calculate_corr_mat(pca_mrna_mat, pca_mrna_to_mrna_mat,
                                                             pca_tf_to_mrna_mat, pca_mirna_to_mrna_mat_for_mrna)
    mirna_corr_mat = calculate_corr_matrix.calculate_corr_mat(pca_mirna_mat, pca_mirna_to_mrna_mat_for_mirna,
                                                              pca_mirna_to_tf_mat, pca_tf_to_mirna_mat)

    corr_mat, clustering_result = cluster_samples.do_clustering(cluster_num, mrna_corr_mat, mirna_corr_mat,
                                                                weight, sample_id_list)

    # save result
    if os.path.exists(save_file_path):
        pass
    else:
        os.makedirs(save_file_path)

    corr_path = os.path.join(save_file_path, "CorrMatrix.txt")
    clustering_result_path = os.path.join(save_file_path, "Result.txt")

    np.savetxt(corr_path, corr_mat, fmt='%0.3f', delimiter='\t')
    np.savetxt(clustering_result_path, clustering_result, fmt='%s', delimiter='\t')

    print("Done!")

    return 0


if __name__ == "__main__":
    main()
