This program outputs two tabulated text files and several folders. Each folder corresponds to a cluster; it contains: the fasta file 
containing the sequences, the kmer counts of the fasta, the eigen values file and the pca file obtained from the count file. 

The first tabulated text file is named cluster_param.txt, where param is the chain of carachters of the main arguments (kmer length, 
epsilon, delta epsilon, minpoints, dimpca). Its header and content is the following:

cluster_name	epsilon	father_size	son1_size	son2_size

Epsilon  is either a positive float (the real epsilon) or a negative integer that indicates an error code: -1 signifies no sub-cluster 
was found in the cluster_name, and -2 indicates that one cluster was detected in cluster_name; this subcluster is named and put in a 
folder that will not be visited by the algorithm. 

The second tabulated text file is sequence_parameters, where parameters is written as described above. It contains the names of the 
sequences and the respective name of the last cluster it belonged to. 



Next and problems, 28/05/2021 : 

.2 clusters are not visible in sequence file, source of bug ?
Must put cluster of epsilon_code = -1 in a folder/fasta but must not be visited, and disclose size.

louis-mael.gueguen@etu.univ-lyon1.fr
