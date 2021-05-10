#!/usr/bin/env python3
#pyhton script to cluster by succesive dychotomy the sequences in fasta file using mnhn-tree-tools

##FAIRE DU CODE BEAU ET FACILE A RELIRE

import argparse
import os
import shutil
import subprocess
import math
import numpy as np

## arguments
#recap argument: fasta, kmer, epsilon min/max (range), dimpca, threads, verbose, minpoints, help, output

parser = argparse.ArgumentParser(description = 'Process fasta file, kmer length, epsilon arguments. PCA dimensions and numer of threads optionnal.')

parser.add_argument('--fasta_file', '-f', dest = 'fastafile', action = 'store', type = str, required = True,
					help = 'Name of the fasta file containing the sequences to cluster during the different steps.')

parser.add_argument('--epsilon_min' , '-e1', dest = 'epsilonmin', action = 'store', type = int, required = True, 
					help = 'minium value of epsilon used for the different calculation steps.') 

parser.add_argument('--epsilon_max' , '-e2', dest = 'epsilonmax', action = 'store', type = int, required = True,
                                        help = 'maximum value of epsilon used for the different calculation steps.') 

parser.add_argument('--min_points', '-m', dest = 'minpoints', action = 'store', default = 3, type = int, 
					help = 'Minimun number of sequences to find to form a cluster during the clustering steps (default: 3).')

parser.add_argument('--dim_pca', '-d', dest = 'dimpca', action = 'store', default = 5, type = int, 
					help = 'Number of dimensions to use from the pca during the calculation steps (default: 5).')

parser.add_argument('--kmer_len', '-k', dest = 'kmer', action = 'store', type = int, required = True, 
					help = 'Length of the kmers processed by mnhn-tree-tools from the fastafile for the different calculation steps.')

parser.add_argument('--threads', '-n', dest = 'threads', action = 'store', default = 4, type = int, 
					help = 'Number of threads used by mnhn-tree-tools for the different calculation steps (default: 4).') 

parser.add_argument('--verbose', '-v', dest = 'verbose', action = 'store_true', default = False, type = 'string',
					help = 'Shows progression messages in standard output.')

parser.add_argument('--output', '-o', dest = 'output', action = 'store', type = 'string', required = True,
					help = 'Specify output file name if needed.')

args = parser.parse_args()

## subfunctions

def create_dir(path, verbose) :
	"""
	Checks if folder already exists or not to create it. If yes, aksks for permission to overwrite: otherwise it creates it.
	path: path-like object
	Command usage: create_folder(./temp)
	"""
	path = str(path)
	if os.path.isdir(path) == True :
		answer  = input("The directory {path} already exists, overwrite ? [O/n]")
		good_answer = False
		while good_answer == False :
			if answer == "O" :
				shutil.rmtree(path) 
				os.mkdir(path)
				good_answer = True
			elif answer == "n" :
 				good_answer = True
			else :
	                	answer = input("Sorry, but '{answer}'is not a valid answer.\nThe directory {path} already exists, overwrite ? [O/n]")
	else :
		os.mkdir(path)
	if args.verbose :
        	print("The directory  has been created.")

def count_files(path) :
	"""
	Count files in a directory. Returns True if there are two files, the number of files otherwise.
	path: path-like object
	Command usage: count_files(./cluster1)
	"""
	path = str(path)
	files = [name in name for os.listdir(path) if os.path.isfile(name) == True]
	if len(files) == 2 :
		return True
	else :
		return len(files)

def iter_epsilon(args.epsilonmin, args.epsilonmax, args.minpoints, args.verbose) :
	"""
	Creates two new clusters from one by a dichotomic search of the right epsilon value. Creates two new fodlers for the new clusters. each name takes the prant's name and adds 1 and 2 respectively.
	args.epsilonmin: minimum value of epsilon desired for the dichotomic search.
	args.epsilonmax: maximum value of epsilon desired for the dichotomic search.
	args.minpoints: minimum number of points (which are pca projections of sequences) in the epsilon radius required to create a new cluster.
	args.verbose: (optionnal) more information about the process in the standard output, aka the console.
	"""
	eps1 = args.epsilonmin
	eps2 = args.epsilonmax
	test = False
	epsilon = round((eps1+eps2)/2) #starts dichotomy in the middle of the given range of epsilon values
	curr_dir = os.getcwd().split("/")[-1]
	fasta = [i for i in os.listdir("./") if i.startswith("fastaCL")][0]
	pca = [i for i in os.listdir("./") if i.endswith(".pca")][0]
	while test != True :
	        shutil.rmtree('cluster')#del the directory itself too, dont want that
        	subprocess.Popen(['cluster_dbscan_PCA', fasta, pca, args.epsilon, args.minpoints, 'cluster/fastaCL']), 'cluster/ev'])], shell = True,stdout = subprocess.PIPE)
        	test = count_files('cluster')
        	if epsilon == args.epsilonmin or epsilon == args.epsilonmax :
			shutil.rmtree('cluster')
			os.chdir("../")
			if args.verbose :
				print("The cluster is a leaf, leaving directory ...")
		elif test < 2 :
                	eps1 = args.epsilonmin
                	eps2 = epsilon
			epsilon = math.floor((eps1+eps2)/2) #necessary to avoid case eps1 = 1 and eps2 = 2, round((1+2)/2) = 2, would result in infinite loop
        	elif test > 2 :
                	eps1 = epsilon
                	eps2 = args.epsilonmax
                        epsilon = math.ceil((eps1+eps2)/2) #same case with eps1 = 9 and eps2 = 10
		else :
			if args.verbose == True :
				print("Clustering of {curr_dir} done, with epsilon = {epsilon}.")
			with open(parameters.join(["../",".txt"]), "a") as f: 
				f.writelines([str(curr_dir), "\t", str(epsilon), "\n"])
			dir_1 = curr_dir.join(["../", ".1"]) #names of the next directories
			dir_2 = curr_dir.join(["../", ".2"])
			os.mkdir(dir_1) #creates the output directories of the new clusters
			os.mkdir(dir_2)
			if args.verbose :
				print("Moving files of the new cluster to their directory ...")
			shutil.move('cluster/fastaCL1', dir_1) #moves the files to the next directories
			shutil.move('cluster/ev1', dir_1)
			shutil.move('cluster/fastaCL2', dir_2)
                        shutil.move('cluster/ev2', dir_2)
			os.remove('cluster') #cleans things up 
	os.chdir("../")
	return dir_1.lstrip("../"), dir_2.lstrip("../") #returns the two directories created, which contains respectively the newly created cluster 1 and 2 from the previous cluster


def extract_names(source, args.verbose) :
	"""
	Extracts the names of the sequences from the parent fasta file located in source. Index is given by the outptu of extract_kmer().
	source: parent folder where to find the parent fasta file. Path.
	args.verbose = boolean
	Command usage: extract_name(/path/to/parent/folder, [1, 2, 3, 20], True)
	"""
	if args.verbose  :
                print("Extracting names of the sequences in the fasta file ...")
	old_fasta = [i for i in os.listdir(source) if i.endswith(".fst")][0] #selects the parent fasta file
	curr_fasta = [i for i in os.listdir("./") if i.endswith(".fst")][0] #selects the child fasta file
	with open(old_fasta, "r") as f :
                lines = f.readlines()
                names  = [i for i in lines if i.startswith(">")] #list of sequences nqmes from parent file, looks like ">ESHIH49767"
        f.close()
	if args.verbose :
		print("Writing names of the sequences ...")
	with open(curr_fasta, "r+") as f :
		lines = f.readlines()
		index = [int(i.lstrip(">sequence_")) for i in lines if i.startswith(">")] #index of sequences from child fasta file we want in parent fasta file
		n = 0
		for j in range(len(lines)) : #goes through lines by index value
			if j.startswith(">") : #if line is a sequence id
				n += 1 #counts indexes of sequences ids
				if n in index : #if the index of the sequence is in the list of the ones to get
					lines[j] = names[n] #changes sequence name corresponding to index
	f.close()

def extract_kmer(source, args.verbose) :
	"""
	Extracts the kmers corresponding to the sequences in the cluster's fasta file from the primary kmer file. Avoids repeating a demanding computation task.
	The extraction is done prior to clustering and pca.
	source: parent folder where to find the parent kmer file.
	Command usage: extract_kmer(/path/to/parent/folder, True)
	"""
	if args.verbose  :
		print("Extracting counts of the sequences in the fasta file ...")
	file = [i in os.listdir("./") if i.endswith(".fst")][0] #gets the name of the fasta file 
	with open(file, "r") as f :
		lines = f.readlines()
		index = [int(i.lstrip(">sequence_")) for i in lines if i.startswith(">")]# ">sequence_22" --> [22]
	f.close()
	counts = np.genfromtxt(str(source).join(["counts.kmer"], dtype = str) #generates a table of counts using numpy
	sub_counts = counts[index, :] #selects the rows corresponding to our sequences
	np.savetxt("./counts.kmer", sub_counts, fmt = '%s', delimiter = "\t") #saves in a new counts.kmer file
	if args.verbose :
		print("Sub-counts saved in counts.kmer")

## main

create_dir(args.ouput, args.verbose) #create output dir
if args.verbosity :
	print("Going to {args.output} directory")
os.chdir(args.output)
parameters = ".".join(str(args.kmer), str(args.epsilonmin), str(args,epsilonmax), str(args.minpoints), str(args.dimpca), "txt") #name for file containing ckuster name and its corresponding epsilon value ("-" if cluster is a leaf)
open(parameters, "w") #creates parameters text file
os.mkdir("cluster1") #creates the directories for the two first clusters
os.mkdir("cluster2")
folders = ["cluster1", "cluster2"] #initialises the list of cluster directories to visit
if os.path.isfile("./counts.kmer") :
	#checks if the counts are already done, to avoid repeating a demanding calculation step
	if args.verbose :
		print("File counts.kmer already exists, skipping counting step.")	
else :
	if args.verbose :
                print("Counting kmers and saving in a primary counts file ...")
	subprocess.Popen(['fasta2kmer', args.fastafile, args.kmer, args.threads, '0', '>', 'counts.kmer'], shell = True, stdout = subprocess.PIPE)

subprocess.Popen(['kmer2pca', 'counts.kmer', 'counts.pca', 'counts.ev', args.kmer, args.threads], shell = True, stdout = subprocess.PIPE)
#subprocess to excecute the command to get the two first clusters
source = os.getcwd().split("/")[-1] #gets the name of the parent folder
for i in folders :
	os.chdir(i)
	extract_kmer(source, args.verbose) #takes parent counts.kmer, extracts kmers into a new counts.kmer 
	extract_names(source, args.verbose) #extracts real sequences's names and inject them in new fasta file
	subprocess(['kmer2pca', 'counts.kmer', 'counts.pca', 'counts.ev', args.kmer, args.threads], stdout = subprocess.PIPE)#produces the pca file
	a, b = iter_epsilon(epsilonmin, epsilonmax, minpoints, output, verbose)
	folders.append(a, b) #adds the new folders to the list so they are visited too
	source = a.rstrip(".1")


##### louis-mael.gueguen@etu.univ-lyon1.fr alpha10.05.2021
