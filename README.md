This program outputs 3 tabulated files and several folders. Each folder corresponds to a cluster; it contains: the fasta file 
containing the sequences, the kmer counts of the fasta, the eigen values file and the pca file obtained from the count file. 

The first tabulated text file is named cluster_param.txt, where param is the chain of carachters of the main arguments (kmer length, 
epsilon, delta epsilon, minpoints, dimpca). Its header and content is the following:

cluster_name	epsilon	father_size	son1_size	son2_size

Epsilon  is either a positive float (the real epsilon) or a negative integer that indicates an error code: -1 signifies no sub-cluster 
was found in the cluster_name, and -2 indicates that one cluster was detected in cluster_name; this subcluster is named and put in a 
folder that will not be visited by the algorithm. 

The second tabulated file is sequence_parameters, where parameters is written as described above. It contains the names of the 
sequences and the respective name of the last cluster it belonged to. Exemple given:

sequenceXXXXXXXX cluster
sequenceYYYYYYYY cluster.1.2.1

The third file is a summary of sequence_parameters, conveniently named sequence_summary.txt . It contains the number of sequences 
remaining in a cluster. In the case of a leaf, the number corresponds to the number of sequence in the leaf. In case of a branch, the 
number corresponds to the sequences fromn the custer that were not assigned to a sub-cluster; these are orphans.  


###louis-mael.gueguen@etu.univ-lyon1.fr
