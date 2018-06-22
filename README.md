
# Improvement of cancer subtype prediction by incorporating transcriptome expression data and heterogeneous biological network #

Requirements
----------
1. python 3.6.1

2. numpy (1.12.1)

3. scikit-learn (0.18.1)

Installation
----------
The source code can be directly called from python.

Usage
----------
python main.py

	-h|--help: get help

	-m|--mrnafile: the mRNA expression data file path

	-t|--tffile: the TF expression data file path

	-r|--mirnafile: the mRNA expression data file path

	-n|--netfile: the interaction network file path

	-s|--savefile: the directory to save the results

	-c|--clusternum: the number of clusters

	-p|--pcaratio: the PCA component ratio

	-w|--weight: the weight of mRNA feature weight

Example
----------
python main.py -m mRNA_Expression_Data.txt -t TF_Expression_Data.txt -r miRNA_Expression_Data.txt -n Interaction_Network_Data.txt -s Result_Directory -c 3 -p 0.95 -w 0.3

**Input**

1. mRNA_Expression_Data.txt

2. TF_Expression_Data.txt

3. miRNA_Expression_Data.txt

4. Interaction_Network_Data.txt


**Output**  

'Patient Subtypes Labels': the first column represents patients' ID, and the second column represents patients' label.  
'Correlation Matrix': the patients' correlation matrix.

----------
Copyright and License Information
----------
Copyright (C) 2018 Northwestern Polytechnical University, Xi’an, China. Yang Guo(gyang@mail.nwpu.edu.cn) and Yang Qi(yang.qi@mail.nwpu.edu.cn)
