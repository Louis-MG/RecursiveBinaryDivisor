# Summary

The Recursive Binary Divisor recursiverly splits clusters in two until it cannot annymore. For this, the  script uses the mnhn-tree-tools (Haschka T. 2021) in a specific way. It first operates the count of kmers, the applies a pca on these counts. The projection is used as a support for dbscan. The epsilon paramter of dbscan is increased until a binary division occures. Then it adds the resulting clusters to a list. The tool visits and expands the list at the same time.  
# Installation

* Dependencies:

Make sure you have a GCC compiler version >= 4.9.2

For MNHN-Tree-Tools:

```
sudo apt-get install git build-essential libpng-dev libsdl2-dev liblapack-dev libopenmpi-dev libpocl-dev ocl-icd-opencl-dev pocl-opencl-icd

git clone https://github.com/Louis-MG/MNHN-Tree-Tools.git

cd MNHN-Tree-Tools
mkdir bin
make all
cd bin

# make MNHN-Tree-Tools available from any folder ( this is temporary, you
# may modify your .bashrc and similar files to make this permanent
export PATH=$PATH:$PWD
```

For RBD:

```
#if you dont have pip installed already (UNIX):

curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py 

#in case of installation troubles, refer to https://pip.pypa.io/en/stable/installing/
#python packages required:

pip install regex
pip install shutils
pip install numpy
pip install subprocess
pip install argparse
pip install textwrap

```

For rbd_pca.r:

```
#if your dont have r installed already, refer to https://linuxize.com/post/how-to-install-r-on-ubuntu-20-04/
#open an R session in the console and type:
install.packages('ggplot2', 'seqinr')

#or install them through the CLI following the tutorial https://linuxize.com/post/how-to-install-r-on-ubuntu-18-04/

```

# Usage

e.g. :
``` 
python3 rbd.py -f test.fasta -e 0.5 -d 0.01 -o output_rbd -g
Rscript rbd_pca.r output_rbd
python3 seq_highlight.py -s seq_of_interest.txt -f test.fst -o seq_highlight.txt -r sequence_parameters.txt

#reminder: sequence_parameters.txt is one of the 3 outputs from rbd (hence the long option name --rbd). See output section below for more details. 
#for more options, you can use:

python3 rbd.py -h
python3 seq_highlight.py -h 
```

# Results

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

* The R script yields 3 pca plots of the sequences for each parent cluster (clusters that are divided by rbd). Dimensions 1 and 2, 2 and 3 are used. 3 additionnal plots are produced, showing final clusters in colors only. They are placed at the output folder's root. 
* The seq_highlight.py script yields a tabulated file that contains the ids of the sequences of interest and their last cluster belonging assigned by rbd.

|sequence_id|cluster_name|
|:---------:|:----------:|
|sequenceXXXXXXX|cluster|
|sequenceUUUUUUU|cluster.1.1.2|
|sequenceOOOOOOO|cluster.1.1.1.1|


###louis-mael.gueguen@etu.univ-lyon1.fr
