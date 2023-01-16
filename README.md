# Summary

The Recursive Binary Divisor clusters sequences using a density-based algorithm. For this, the script uses the mnhn-tree-tools (Haschka T. 2021) in a specific way. It first operates the count of kmers, then applies a pca on these counts. The projection is used as a support for dbscan: distance between each point is calculated; if enough points are within distance epsilon of each-others, the density threshold is reached and a cluster is found. The density parameter epsilon is reduced until exactly 2 clusters are found. The sequences forming each cluster are separated, the others are left out. The steps count/pca/projection/density-based-clustering are repeated on each precendently found clusters. The final result is a tree. Its leefs are the clusters that are visualised.

# Installation

* Dependencies:

Make sure you have a GCC compiler version >= 4.9.2, pip and the R language installed.

For MNHN-Tree-Tools:

```
sudo apt-get install git build-essential libpng-dev libsdl2-dev liblapack-dev libopenmpi-dev libpocl-dev ocl-icd-opencl-dev pocl-opencl-icd

git clone https://github.com/Louis-MG/MNHN-Tree-Tools.git

cd MNHN-Tree-Tools
mkdir bin
make all
cd bin

# make MNHN-Tree-Tools available from any folder temporarily
export PATH=$PATH:$PWD
```

For RBD:

```
#if you dont have pip installed already (UNIX):

pip install regex shutils numpy subprocess argparse textwrap
```

For rbd_pca.r, open R in the console and paste:

```
install.packages('ggplot2', 'seqinr')
```

# Usage

```
python3 rbd.py -h
python3 seq_highlight.py -h
```

e.g. :
``` 
python3 rbd.py -f test.fasta -e 0.5 -d 0.01 -o output_rbd -g
Rscript rbd_pca.r output_rbd
python3 seq_highlight.py -s seq_of_interest.txt -f test.fst -o seq_highlight.txt -r sequence_parameters.txt
```

# Output files format

The rbd program outputs 3 tabulated files and several folders. Each folder corresponds to a cluster; it contains: the fasta file 
containing the sequences, the kmer counts of the fasta, the eigen values file and the pca file obtained from the count file. 

* The first tabulated text file is named cluster_parameters.txt, where parameters is the chain of carachters of the main arguments (kmer length, 
epsilon, delta epsilon, minpoints, dimpca). Its header and content is the following:

|cluster_name|epsilon|father_size|son1_size|son2_size|
|:----------:|:-----:|:---------:|:-------:|:-------:|
|cluster|1.2|1000|400|500|
|cluster.1|1.5|400|200|100|	
|cluster.1.1|-1|200|NONE|NONE|
|cluster.1.2|-5|100|NONE|NONE|

Epsilon  is either a positive float (the real epsilon) or a negative integer that indicates an error code: -1 signifies no or just one sub-cluster was found.
 
* The second tabulated file is sequence_parameters, where parameters is written as described above. It contains the names of the 
sequences and the respective name of the last cluster it belonged to. Exemple given:

|sequences_names|cluster_name|
|:-------------:|:-----:|
|sequenceXXXXXXX|cluster|
|sequenceYYYYYYY|cluster.1.2.1|
|sequenceZZZZZZZ|cluster.2.1|
|sequenceAAAAAAA|cluster.1.0|

* The third file is a summary of sequence_parameters, conveniently named sequence_summary.txt . It contains the number of sequences 
remaining in a cluster. In the case of a leaf, the number corresponds to the number of sequence in the leaf. In case of a branch, the 
number corresponds to the sequences fromn the custer that were not assigned to a sub-cluster; these are orphans and are denominated with a .0 added ot the cluster name.  

|Sequences left|cluster_name|
|:-:|:--------:|
|240|cluster.0|
|45|cluster.1.0|
|30|cluster.1.2|
|12|cluster.2|

* The R script yields 3 pca plots of the sequences for each parent cluster (clusters that are divided by rbd). Dimensions 1 and 2, 2 and 3, 1 and 3 are used. 3 additionnal plots are produced, showing final clusters in colors only. They are placed at the output folder's root. 
* The seq_highlight.py script yields a tabulated file that contains the ids of the sequences of interest and their last cluster belonging assigned by rbd.

|sequence_id|cluster_name|
|:---------:|:----------:|
|sequenceXXXXXXX|cluster|
|sequenceUUUUUUU|cluster.1.1.2|
|sequenceOOOOOOO|cluster.1.1.1.1|


###louis-mael.gueguen@crchudequebec.ulaval.ca
