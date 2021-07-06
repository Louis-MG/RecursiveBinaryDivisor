#!/usr/bin/env python3
#pyhton script to bi-cluster the sequences in fasta file using mnhn-tree-tools

import argparse
import os
import shutil
import subprocess
import numpy as np
import sys
from textwrap import dedent
import re

## arguments
#recap argument: fasta, kmer, epsilon, delta epsilon, dimpca, threads, verbose (2 levels), minpoints,  output, help

parser = argparse.ArgumentParser(prog = 'script.py', formatter_class = argparse.RawDescriptionHelpFormatter, description = dedent('''\
			Process fasta file, kmer length, epsilon arguments. PCA dimensions and numer of threads optionnal.'

			This program outputs three tabulated text files and several folders. Each folder corresponds to a cluster; it contains: 
			 - the fasta file containing the sequences
			 - the kmer counts of the fasta
			 - the eigen values file and the pca file obtained from the count file. 
			The first tabulated text file is named cluster_param.txt, where param is the chain of characters of the main arguments (kmer length, epsilon, delta epsilon, minpoints, dimpca). 
			Its header and content is the following:
			
			cluster_name	epsilon	father_size	son1_size	son2_size	
			
			Epsilon  is either a positive float (the real epsilon) or a negative integer that indicates an error code: 
			
			 - -1 signifies no sub-cluster was found in the cluster_name
			 - -2 indicates that one cluster was detected in cluster_name; this subcluster is named and put in a folder that will not be visited by the algorithm.
			
			The second tabulated text file is sequence_parameters.txt, where parameters is written as described above. It contains the names of the sequences and the respective name of the last cluster it belonged to.

			>sequence_name cluster_name

			The third tabulated text file is sequence_summary.txt. It contains the number of sequences in each cluster after the complete distribution. If the cluster is a leaf, the number is the number of sequences in the cluster. If the cluster is a branch the number is the number of sequences that are orphans of the next division. 
			For 1000 sequences :			

			200 cluster
			500 cluster.1
			300 cluster.2
			'''))

parser.add_argument('--fasta_file', '-f', dest = 'fastafile', action = 'store', type = str, required = True, 
				help = 'Name of the fasta file containing the sequences to cluster during the different steps.')

parser.add_argument('--epsilon', '-e', dest = 'epsilon', action = 'store', type = float, required = True, 
						help = 'Minium value of epsilon used for the different calculation steps. Float.') 

parser.add_argument('--delta_epsilon', '-d', dest = 'delta', action = 'store', type = float, default = 0.001,
						help = 'Minium difference between each epsilon of the dichotomy: when difference is inferior to delta, loop stops (default: 0.001).')

parser.add_argument('--growth', '-g', dest = 'growth', action = 'store', default = False, type = bool, 
						help = 'Gives directionnality of epsilon: increase or decrease')

parser.add_argument('--min_points', '-m', dest = 'minpoints', action = 'store', default = 3, type = int,
			help = 'Minimun number of sequences to find to form a cluster during the clustering steps (default: 3).')

parser.add_argument('--dim_pca', '-p', dest = 'dimpca', action = 'store', default = 7, type = int, 
				help = 'Number of dimensions to use from the pca during the calculation steps (default: 7).')

parser.add_argument('--kmer_len', '-k', dest = 'kmer', action = 'store', type = int, default = 5, 
		help = 'Length of the kmers processed by mnhn-tree-tools from the fastafile for the different calculation steps. Default: 5')

parser.add_argument('--threads', '-n', dest = 'threads', action = 'store', default = 4, type = int, 
			help = 'Number of threads used by mnhn-tree-tools for the different calculation steps (default: 4).') 

parser.add_argument('--verbose', '-v', dest = 'verbose', action = 'count', default = 0, 
						help = 'Shows progression messages in standard output.')

parser.add_argument('--output', '-o', dest = 'output', action = 'store', type = str, required = True, 
							help = 'Specifies output file name if needed.')

parser.add_argument('--minsize', '-s', dest = 'minsize', action = 'store', type = int, default = 50, 
							help = 'Specifies minimum size of a cluster to try to subcluster it.')

args = parser.parse_args()

## subfunctions

def create_dir(path, verbose) :
	"""
	Checks if folder already exists or not to create it. If yes, asks for permission to overwrite: otherwise it creates it.
	path: path-like object
	Command usage: create_folder(./temp)
	"""
	path = str(path)
	if os.path.isdir(path) == True :
		answer  = input("The directory {} already exists, overwrite ? [O/n]".format(path))
		good_answer = False
		while good_answer == False : #loops until user responds either O (capitalised for to be conservative) or n.
			if answer == "O" :
				shutil.rmtree(path) 
				os.mkdir(path)
				good_answer = True
			elif answer == "n" :
				good_answer = False
				print("Directory not overwitten, choose another name.")
				sys.exit()
			else :
				answer = input("Sorry, but '{}'is not a valid answer.\nThe directory {} already exists, overwrite ? [O/n]\n".format(answer, answer))
	else :
		os.mkdir(path)
		if args.verbose >= 1:
        		print("The directory has been created.")

def count_files(path) :
	"""
	Count files in a directory. Returns the number of files otherwise.
	path: path-like object
	Command usage: count_files(./cluster1)
	"""
	path = str(path)
	files = [name for name in os.listdir(path) if os.path.isfile(os.path.join(path,name))]
	return len(files)

def iter_epsilon(epsilon, delta, dimpca, growth, minpoints, verbose) :
	"""
	Creates two new clusters from one by a dichotomic search of the right epsilon value. Creates two new fodlers for the new clusters. each name takes the prant's name and adds 1 and 2 respectively.
	epsilon: value of epsilon desired for the split.
	delta: value of the incremental value.
	growth: boolean for growth or degrowth/ increase or decrease of the epsilon parameter.
	minpoints: minimum number of points (which are pca projections of sequences) in the epsilon radius required to create a new cluster.
	verbose: (optionnal) more information about the process in the standard output, aka the console.
	Command usage: iter_epsilon(0.8, 0.01, 7, True, 3, 1)
	"""
	curr_dir = os.getcwd().split("/")[-1]
	fasta = [i for i in os.listdir("./") if i.startswith("fasta")][0] #selects the fasta file in the current directory
	pca = [i for i in os.listdir("./") if i.endswith(".pca")][0] #selects counts.pca
	with open(fasta, 'r') as f :
		lines = f.readlines()
		sequences = [i.rstrip('\n') for i in lines if i.startswith(">")]
	command = " ".join(["grep", "-c", "'>'", fasta]) #this line and the following count the number of sequences in  children files to store in in output file
	parent = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines=True).rstrip()
	out_loop = False
	while out_loop != True : #allow to loop until a return
		if verbose >= 1 :
			print(epsilon)
		if os.path.exists("cluster") :
			shutil.rmtree('cluster')#del the directory itself too, dont want that
		os.mkdir("cluster")
		command = " ".join(['cluster_dbscan_pca', fasta, str(pca), str(dimpca), str(epsilon), str(minpoints), 'cluster/fastaCL', 'cluster/ev'])
		if verbose >= 2 :
			print(command)
		p = subprocess.Popen(command, shell = True, stderr = subprocess.PIPE) #finds clusters
		p.wait()
		test = count_files('cluster') #counts clusters
		if test > 4 :
			if args.growth :
				epsilon += args.delta #increases epsilon if more than 2 clusters
				shutil.rmtree("cluster") #cleans things up
			else : 
				shutil.rmtree("cluster") #if epsilon decreases, this is a dead-end and the program ends
				if verbose >= 1 :
					print("Clustering of {} yielded more than 2 clusters.".format(curr_dir))
				with open("../cluster_"+parameters+".txt", "a") as f :
					f.writelines([str(curr_dir), "\t", "-1", "\t", parent, "\t", "NONE", "\t", "NONE", "\n"]) #error code
				sys.exit()
		if test < 4 :
			if test <= 0 :
				if verbose >= 1 :
					print('Clustering of {} yielded 0 cluster'.format(curr_dir))
				with open("../cluster_"+parameters+".txt", "a") as f:
					f.writelines([str(curr_dir), "\t", "-1", "\t", parent, "\t","NONE","\t","NONE", "\n"]) #error code $
				shutil.rmtree("cluster") #cleans things up
				return([])
			if args.growth :
				if verbose >= 1:
					print("Clustering of {} yielded less than 2 clusters.".format(curr_dir))
				with open("../cluster_"+parameters+".txt", "a") as f:
                        	        f.writelines([str(curr_dir), "\t", "-1", "\t", parent, "\t","NONE","\t","NONE", "\n"]) #error code and no sons so no sizes
				shutil.rmtree("cluster") #cleans things up
				return [] #returns list of length = 0
			else :
				epsilon -= args.delta #decreases epsilon if less than 2 clusters
				shutil.rmtree("cluster")
		elif test == 4 : #if 2 clusters are found, saved as .1 and .2 and added to the list for further split
			if verbose == True :
				print("Clustering of {} done, with epsilon = {}.".format(curr_dir, epsilon))
			command = "grep -c '>' cluster/fastaCL-000" #this line and the following count the number of sequences in  children files to store in in output file
			son1 = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines=True).rstrip()
			command = "grep -c '>' cluster/fastaCL-001"
			son2 = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines=True).rstrip()
			if verbose >= 2 :
				print("Parent with size {} has been clustered in {}  of size {} and {} of size {}".format(parent, curr_dir+".1", son1, curr_dir+".2", son2))
			with open("../cluster_"+parameters+".txt", "a") as f:
				f.writelines([str(curr_dir), "\t", str(epsilon), "\t", parent, "\t", son1, "\t", son2, "\n"])
			dir_1 = curr_dir.join(["../", ".1"]) #names of the next directories
			dir_2 = curr_dir.join(["../", ".2"])
			os.mkdir(dir_1) #creates the output directories of the new clusters
			os.mkdir(dir_2)
			if verbose >= 2 :
				print("Moving files of the new cluster to their directory ...")
			os.rename("cluster/fastaCL-000", "cluster/fastaCL1.fst")
			shutil.move('cluster/fastaCL1.fst', dir_1) #moves the files to the next directories
			shutil.move('cluster/ev-000', dir_1)
			os.rename("cluster/fastaCL-001","cluster/fastaCL2.fst")
			shutil.move('cluster/fastaCL2.fst', dir_2)
			shutil.move('cluster/ev-001', dir_2)
			shutil.rmtree("cluster") #cleans things up
			return [dir_1.lstrip("../"), dir_2.lstrip("../")] #returns the two directories created, which contains respectively the newly created cluster 1 and 2 from the previous cluster

def extract_names(source, verbose) :
	"""
	Extracts the names of the sequences from the parent fasta file located in source. Index is given by the output of extract_kmer().
	source: parent folder where to find the parent fasta file. Path.
	verbose = boolean
	Command usage: extract_name(/path/to/parent/folder, 1)
	"""
	if verbose >= 1 :
		print("Extracting names of the sequences in the fasta file ...")
	old_fasta = [i for i in os.listdir("../"+source) if i.endswith(".fst")][0] #selects the parent fasta file
	curr_fasta = [i for i in os.listdir("./") if i.endswith(".fst")][0] #selects the child fasta file
	old_fasta = "../" + source+ "/" + old_fasta 
	with open(old_fasta, "r") as f :
		lines = f.readlines()
		names  = [i.split(' ')[0].rstrip()+"\n" for i in lines if i.startswith(">")] #list of sequences names from parent file, looks like ">ESHIH49767"
	f.close()
	with open(curr_fasta, "r") as f :
		lines = f.readlines()
		index = [int(i.lstrip(">sequence_")) for i in lines if i.startswith(">")] #index of sequences from child fasta file we want in parent fasta file
		n = 0
		if verbose >= 1:
			print("Writing names of the sequences ...")
		for j in range(len(lines)) : #goes through lines of the new fasta file where we rename the seuquences
			if lines[j].startswith(">") : #if line is a sequence id
				#print("ligne j",lines[j].rstrip())
				#print("index[n]", index[n])			#these three lines are for debugging
				#print('nom de la sequence', names[index[n]])
				lines[j] = names[index[n]] #changes sequence name at line j by the name at position [index[n]] which is the n_th position of the index, itself the  number of the sequence extracted from ">sequence_X". At each iteration of sequence name the n is incremented (+1)
				n += 1
	with open(curr_fasta, "w") as f :
		f.writelines(lines) #replaces the lines by new ones

def extract_kmer(source, verbose) :
	"""
	Extracts the kmers corresponding to the sequences in the cluster's fasta file from the primary kmer file. Avoids repeating a demanding computation task.
	The extraction is done prior to clustering and pca.
	source: parent folder where to find the parent kmer file.
	Command usage: extract_kmer(/path/to/parent/folder, 1)
	"""
	if verbose >= 1 :
		print("Extracting counts of the sequences in the fasta file ...")
	file = [i for i in os.listdir("./") if i.endswith(".fst")][0] #gets the name of the fasta file
	with open(file, "r") as f :
		lines = f.readlines()
		index = [int(i.lstrip(">sequence_")) for i in lines if i.startswith(">")]# ">sequence_22" --> [22]
	f.close()
	counts = np.genfromtxt("../"+str(source)+"/counts.kmer", dtype = str) #generates a table of counts from the kmer file in source folder using numpy
	sub_counts = counts[index, :] #selects the rows corresponding to our sequences
	np.savetxt("./counts.kmer", sub_counts, fmt = '%s', delimiter = "\t") #saves in a new counts.kmer file
	if verbose >= 2 :
		print("Sub-counts saved in counts.kmer")
	
## main

create_dir(args.output, args.verbose) #create output dir
p = subprocess.Popen(['cp', args.fastafile, args.output], stderr = subprocess.PIPE)
p.wait()
if args.verbose >= 1 :
	print("Going to {} directory".format(args.output))
os.chdir(args.output)
fasta = [i for i in os.listdir("./") if i.endswith(".fst")][0]
parameters = "_".join([str(args.kmer), str(args.epsilon), str(args.delta),str(args.minpoints), str(args.dimpca)]) #name for file containing ckuster name and its corresponding epsilon value ("-" if cluster is a leaf)
with open("cluster_"+parameters+'.txt', "w") as f :
	f.writelines(["cluster_name\tepsilon\tfather_size\tchild1_size\tchild2_size\n"])
os.mkdir("cluster.1") #creates the directories for the two first clusters
os.mkdir("cluster.2")
folders = ["cluster.1", "cluster.2"] #initialises the list of cluster directories to visit
other_folders = [] #initiates the list of clusters yielded by a single cluster in main loop
if args.verbose >= 1 :
	print("Counting kmers and saving in a primary counts file ...")
command = " ".join(['fasta2kmer', fasta, str(args.kmer), str(args.threads), '0', '>', './counts.kmer'])
p = subprocess.Popen(command, shell = True, stderr = subprocess.PIPE)
p.wait()
command = " ".join(['kmer2pca', 'counts.kmer', 'counts.pca', 'counts.ev', str(args.dimpca), str(args.threads)]) 
p = subprocess.Popen(command, shell = True, stderr=subprocess.PIPE)
p.wait()
epsilon = args.epsilon
command = " ".join(["grep", "-c", "'>'", fasta]) #this line and the following count the number of sequences in  children files to store in in output file
parent = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines = True).rstrip()
out_loop = 0
curr_dir = os.getcwd().split("/")[-1]
while out_loop != 4 :
	if args.verbose >= 1 :
		print(epsilon)
	if os.path.exists('cluster') : 
		shutil.rmtree('cluster')
	os.mkdir('cluster')
	command = " ".join(['cluster_dbscan_pca', fasta, "counts.pca", str(args.dimpca), str(epsilon), str(args.minpoints), 'cluster/fastaCL', 'cluster/ev'])
	p = subprocess.Popen(command, shell = True, stderr = subprocess.PIPE)
	p.wait()
	out_loop = count_files('cluster')
	if out_loop > 4 :
		if args.growth : 
			epsilon += args.delta
		else :
			print("Clustering of {} yielded more than two clusters.".format(args.fastafile))
			shutil.rmtree("cluster")
			sys.exit() 
	elif out_loop < 4 : #if less than 2 clusters thus 4 files found:
		if args.growth :
			print("Clustering of {} yielded less than 2 clusters.".format(args.fasta))
			with open("./cluster_"+parameters+".txt", "a") as f:
				f.writelines(['cluster', "\t", str(epsilon), "\t", parent, "\t", "NONE", "\t", "NONE", "\n"])
			sys.exit()
		else :
			epsilon -= args.delta
	elif out_loop == 4 :
		if args.verbose >= 1 :
			print("Clustering of {} done, with epsilon = {}.".format(curr_dir, epsilon))
		command = "grep -c '>' cluster/fastaCL-000" #this line and the following count the number of sequences in  children files to store in in output file
		son1 = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines=True).rstrip()
		command = "grep -c '>' cluster/fastaCL-001"
		son2 = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines=True).rstrip()
		if args.verbose >= 2 :
			print("Parent with size {} has been clustered in {}  of size {} and {} of size {}".format(parent, curr_dir+".1", son1, curr_dir+".2", son2))
		with open("./cluster_"+parameters+".txt", "a") as f:
			f.writelines(['cluster', "\t", str(epsilon), "\t", parent, "\t", son1, "\t", son2,"\n"])

os.rename("cluster/fastaCL-000", "cluster/fastaCL1.fst")
shutil.move("cluster/fastaCL1.fst", "cluster.1") #move files to iniate first iterations
shutil.move("cluster/ev-000", "cluster.1")
os.rename("cluster/fastaCL-001", "cluster/fastaCL2.fst")
shutil.move("cluster/fastaCL2.fst", "cluster.2")
shutil.move("cluster/ev-001", "cluster.2")

files = [fasta, "counts.kmer", "counts.pca", "counts.ev"] #moves the files to the folder of the first 
for i in files :
	shutil.move(i, "cluster")

with open("cluster/"+fasta, "r") as f : #saves the names of the sequences
	lines = f.readlines()
	reference = [i.rstrip('\n').split(' ')[0]+"\t"+"cluster\n" for i in lines if i.startswith(">")]
	if args.verbose >= 2 :
		print("Reading sequences names ...")

with open("sequence_"+parameters+".txt", "w") as f : # writes names of the sequences in the table referencing sequences and their respective last branch level
	if args.verbose >= 2 :
		print("Writing reference sequences names in the sequence_{}.txt".format(parameters))
	f.writelines(reference)

leaf = [] #list of leavesn clusters to update sequences_parameters.txt

# loop that goes through all the directories of different clusters to find new clusters
for i in folders :
	if args.verbose >= 1 :
		print("Going to cluster directory {}.".format(i))
	os.chdir(i)
	if i.endswith(".1") :
		source = re.sub(r"\.1$", "", i)
	elif i.endswith(".2"): 
		source = re.sub(r"\.2$", "", i)
	extract_kmer(source, args.verbose) #takes parent counts.kmer, extracts kmers into a new counts.kmer 
	extract_names(source, args.verbose) #extracts real sequences's names and inject them 
	with open([i for i in os.listdir("./") if i.startswith("fasta")][0], "r") as f :
		lines = f.readlines()
		sequences = [i for i in lines if i.startswith(">")] #list of the sequences names in the fasta of the currently visited folder
	if len(sequences) > args.minsize :
		command = " ".join(['kmer2pca', 'counts.kmer', 'counts.pca', 'counts.ev', str(args.dimpca), str(args.threads)])
		if args.verbose >= 2 :
			print(command)
		p = subprocess.Popen(command, shell = True, stderr=subprocess.PIPE)#produces the pca file
		p.wait()
		a = iter_epsilon(epsilon = args.epsilon, delta = args.delta, growth = True, dimpca = args.dimpca, minpoints = args.minpoints, verbose = args.verbose)
		if len(a) == 2 :
			folders.extend(a) #adds the new folders to the list so they are visited too
		else :
			leaf.append(i)
	else :
		leaf.append(i)
	if args.verbose >= 2 :
		print("Updating the last cluster belonging of sequences ...")
	table = np.genfromtxt('../sequence_'+parameters+".txt", dtype = str) #numpy table of elements, n lines and 2 columns 
	for j in range(len(table[:,0])) :
		if table[j,0]+"\n" in sequences : #if the sequences name of the j line is in the sequences list:
			table[j,1] = i
	np.savetxt("../sequence_"+parameters+".txt", table, fmt = '%s', delimiter = "\t")
	os.chdir('../')

table = np.genfromtxt('./sequence_'+parameters+".txt", dtype = str) #numpy table of elements, n lines and 2 columns 
for j in range(len(table[:,1])) :
	if table[j,1] not in leaf : #if the sequences name of the j line is in the sequences list:
		table[j,1] = table[j,1]+'.0'
np.savetxt("./sequence_"+parameters+".txt", table, fmt = '%s', delimiter = "\t") #file is updated with .0 added to oprhan sequences

if args.verbose >= 2 :
	print("Saving summary of sequence_{}.txt".format(parameters))

command = "awk '{x=2; print $x}' "+"{}".format("sequence_"+parameters+".txt")+" | sort | uniq -c"
output = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines = True)

with open("sequence_summary.txt", "w") as f:
	f.writelines(output)

#TODO: change output for orphans
##### louis-mael.gueguen@etu.univ-lyon1.fr alpha10.05.2021
