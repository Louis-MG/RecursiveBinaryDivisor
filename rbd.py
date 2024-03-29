#!/usr/bin/env python3
#pyhton script to bi-cluster the sequences in fasta file using mnhn-tree-tools

import argparse
import os
import numpy as np
import shutil
import subprocess
import sys
from textwrap import dedent
import re
import errno

## arguments
#recap argument: fasta, kmer, epsilon, delta epsilon, growth,  dimpca, threads, verbose (2 levels), minpoints,  output, help

parser = argparse.ArgumentParser(prog = 'script.py', formatter_class = argparse.RawDescriptionHelpFormatter, description = dedent('''\
			Process fasta file, kmer length, epsilon arguments. PCA dimensions and number of threads optionnal.'
			'''))

parser.add_argument('--fasta_file', '-f', dest = 'fastafile', action = 'store', type = str, required = True, 
				help = 'Name of the fasta file containing the sequences to cluster during the different steps.')

parser.add_argument('--epsilon', '-e', dest = 'epsilon', action = 'store', type = float, required = True, 
						help = 'Minium value of epsilon used for the different calculation steps. Float.') 

parser.add_argument('--delta_epsilon', '-d', dest = 'delta', action = 'store', type = float, default = 0.001,
						help = 'Minium difference between each epsilon of the dichotomy: when difference is inferior to delta, loop stops (default: 0.001).')

parser.add_argument('--growth', '-g', dest = 'growth', action = 'count', default = 0, 
						help = 'Gives directionnality of epsilon: increase or decrease. Default: 0.')

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
							help = 'Specifies minimum size of a cluster to try to subcluster it. Default: 50.')

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
	Counts and returns the number of files in a directory.
	path: path-like object
	Command usage: count_files(./cluster1)
	"""
	path = str(path)
	files = [name for name in os.listdir(path) if os.path.isfile(os.path.join(path,name))]
	return len(files)

def iter_epsilon(epsilon, delta, dimpca, growth, minpoints, verbose) :
	"""
	Creates two new clusters from one by search of the right epsilon value. Creates two new folders for the new clusters. Each name takes the parent's name and adds .1 and .2 respectively.
	epsilon: base value used for the split.
	delta: value of the increment or decrement
	growth: 1 or 0 for increase or decrease of the epsilon parameter.
	minpoints: minimum number of points (which are pca projections of sequences) in the epsilon radius required to create a new cluster.
	verbose: (optional) more information about the process in the standard output. 3 levels.
	Command usage: iter_epsilon(0.8, 0.01, 7, 1, 3, 1)
	"""
	curr_dir = os.getcwd().split("/")[-1]
	fasta = [i for i in os.listdir("./") if i.startswith("fasta")][0] #selects the fasta file in the current directory
	pca = [i for i in os.listdir("./") if i.endswith(".pca")][0] #selects counts.pca
	with open(fasta, 'r') as f :
		lines = f.readlines() #stores lines of the fasta file
		sequences = [i.rstrip('\n') for i in lines if i.startswith(">")] #stores sequences from the lines
	command = " ".join(["grep", "-c", "'>'", fasta]) #this line and the following count the number of sequences in  children files to store in in output file
	parent = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines=True).rstrip()
	test = 0
	while test != 4 : #will loop until 4 files are found (success, 2 clusters found) or 0 clusters are found (failed)
		if verbose >= 1 :
			print(epsilon)
		if os.path.exists("cluster") :
			shutil.rmtree('cluster')#del the directory from the previous try if it exists
		os.mkdir("cluster")
		command = " ".join(['cluster_dbscan_pca', fasta, str(pca), str(dimpca), str(epsilon), str(minpoints), 'cluster/fastaCL', 'cluster/ev'])
		if verbose >= 2 :
			print(command)
		p = subprocess.Popen(command, shell = True, stderr = subprocess.PIPE) #finds clusters
		p.wait()
		test = count_files('cluster') #counts clusters
		if test == 0: #if 0 cluster found
			if verbose >= 1 :
				print("Clustering of {} yielded 0 cluster.".format(curr_dir))
			with open("../cluster_"+parameters+".txt", "a") as f :
				f.writelines([str(curr_dir), "\t", "0", "\t", parent, "\t", "NONE", "\t", "NONE", "\n"])
			return([])
		elif test > 4 : #if more than 3 clusters are found (2 files per cluster)
			if args.growth == 1 :
				epsilon += args.delta #increases epsilon if more than 2 clusters
				shutil.rmtree("cluster") #cleans things up
			else : 
				shutil.rmtree("cluster") #if epsilon decreases, this is a dead-end and the program stops
				if verbose >= 1 :
					print("Clustering of {} yielded more than 2 clusters.".format(curr_dir))
				with open("../cluster_"+parameters+".txt", "a") as f :
					f.writelines([str(curr_dir), "\t", str(-test/2), "\t", parent, "\t", "NONE", "\t", "NONE", "\n"])
				return([])
		elif test < 4 : #if less than 2 clusters found
			if args.growth == 1 :
				if verbose >= 1:
					print("Clustering of {} yielded less than 2 clusters.".format(curr_dir)) #if growth, dead-end so program stops
				with open("../cluster_"+parameters+".txt", "a") as f:
					f.writelines([str(curr_dir), "\t", str(-test/2), "\t", parent, "\t","NONE","\t","NONE", "\n"]) 
				shutil.rmtree("cluster") #cleans things up
				return([]) #returns list of length = 0
			else :
				epsilon -= args.delta #decreases epsilon if less than 2 clusters
				shutil.rmtree("cluster")
				if epsilon <= 0 :
					with open("../cluster_"+parameters+".txt", "a") as f:
						f.writelines([str(curr_dir), "\t", "0", "\t", parent, "\t","NONE","\t","NONE", "\n"])
					return([])
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
			##################################################################################################### mettre ca dans une fonciton clean() ou arrange()
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
	Extracts the names of the sequences from the parent fasta file located in source.
	source: parent folder where to find the parent fasta file. Is a path.
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

##########################################################################################recommenter la section, ca veut rien dire
## main

if os.path.isfile(args.fastafile) == False :
	raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.fastafile)

create_dir(args.output, args.verbose) #create output dir
p = subprocess.Popen(['cp', args.fastafile, args.output], stderr = subprocess.PIPE) #copies the fasta file in the output dir
p.wait()
os.rename(args.output+'/'+args.fastafile.split('/')[-1], args.output+'/'+'cluster.fst')
if args.verbose >= 1 :
	print("Going to {} directory".format(args.output))
os.chdir(args.output) #goes to output dir
fasta = [i for i in os.listdir("./") if i.endswith(".fst")][0]
if args.growth :
	parameters = "_".join([str(args.kmer), str(args.epsilon), str(args.delta),str(args.minpoints), str(args.dimpca), str(args.minsize), 'g']) #name for file containing ckuster name and its corresponding epsilon value ("-" if cluster is a leaf)
else :
	parameters = "_".join([str(args.kmer), str(args.epsilon), str(args.delta),str(args.minpoints), str(args.dimpca), str(args.minsize)])
with open("cluster_"+parameters+'.txt', "w") as f : #prepares output file
	f.writelines(["cluster_name\tepsilon\tfather_size\tchild1_size\tchild2_size\n"])
os.mkdir("cluster.1") #creates the directories for the two first clusters
os.mkdir("cluster.2")
folders = ["cluster.1", "cluster.2"] #initialises the list of cluster directories to visit
if args.verbose >= 1 :
	print("Counting kmers and saving in a primary counts file ...")
command = " ".join(['fasta2kmer', fasta, str(args.kmer), str(args.threads), '0', '>', './counts.kmer']) #counts kmers in sequences
p = subprocess.Popen(command, shell = True, stderr = subprocess.PIPE)
p.wait()
command = " ".join(['kmer2pca', 'counts.kmer', 'counts.pca', 'counts.ev', str(args.dimpca), str(args.threads)]) #PCA from kmer matrix
p = subprocess.Popen(command, shell = True, stderr=subprocess.PIPE)
p.wait()
epsilon = args.epsilon
command = " ".join(["grep", "-c", "'>'", fasta]) #counts the number of sequences in  children files to store counts in output file
parent = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines = True).rstrip()
test = 0
curr_dir = os.getcwd().split("/")[-1]
while test != 4 :
	if args.verbose >= 1 :
		print(epsilon)
	if os.path.exists('cluster') : 
		shutil.rmtree('cluster')
	os.mkdir('cluster')
	command = " ".join(['cluster_dbscan_pca', fasta, "counts.pca", str(args.dimpca), str(epsilon), str(args.minpoints), 'cluster/fastaCL', 'cluster/ev'])
	p = subprocess.Popen(command, shell = True, stderr = subprocess.PIPE)
	p.wait()
	test = count_files('cluster')
	if test > 4 : #if more than 2 clusters found
		if args.growth == 1 : #if growth of epsilon
			epsilon += args.delta #increase epsilon
		else :
			print("Clustering of {} yielded more than two clusters.".format(args.fastafile))
			with open("./cluster_"+parameters+".txt", "a") as f:
				f.writelines([str(curr_dir), "\t", str(-test/2), "\t", parent, "\t","NONE","\t","NONE", "\n"])
			shutil.rmtree("cluster") #if not growth, then dead-end reached, stops program
			sys.exit() 
	elif test < 4 : #if less than 2 clusters found:
		if args.growth == 1 : #if growth
			print("Clustering of {} yielded less than 2 clusters.".format(args.fastafile)) #dead-end reached, stops program
			with open("./cluster_"+parameters+".txt", "a") as f:
				f.writelines(['cluster', "\t", str(-test/2), "\t", parent, "\t", "NONE", "\t", "NONE", "\n"])
			sys.exit()
		else :
			epsilon -= args.delta
			if epsilon <= 0 :
				with open("./cluster_"+parameters+".txt", "a") as f:
					f.writelines(['cluster', "\t", "0", "\t", parent, "\t", "NONE", "\t", "NONE", "\n"])
	elif test == 4 :
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
	os.chdir(i) #goes to folder
	if i.endswith(".1") :
		source = re.sub(r"\.1$", "", i) #finds parent folder name
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
		a = iter_epsilon(epsilon = args.epsilon, delta = args.delta, growth = args.growth, dimpca = args.dimpca, minpoints = args.minpoints, verbose = args.verbose)
		if len(a) == 2 :
			folders.extend(a) #adds the new folders to the list so they are visited too
		else :
			leaf.append(i)
	else :
		leaf.append(i)
	if args.verbose >= 2 :
		print("Updating the last cluster belonging of sequences ...")
	with open('../sequence_'+parameters+".txt", 'r') as f:
		table = f.readlines()
	for j in range(len(table)) : #for each line of the file
		line = table[j].split('\t') #store line splitted on tab to separate the columns
		if line[0]+'\n' in sequences : #if the sequence id of the line is in the list of seq_id from the cluster we explore
			line[1] = i+'\n' #update the cluster belonging
		table[j] = line[0]+'\t'+line[1] #restore the original line format
	with open('../sequence_'+parameters+".txt", 'w') as f:
		f.writelines(table) #save
	os.chdir('../')

#following block of code adds a .0 to the cluster name where orphans sequences were at last
with open('./sequence_'+parameters+".txt", 'r') as f:
	table = f.readlines() 
for j in range(len(table)) :
	line = table[j].split() #store line splitted on tab to separate the columns
	if line[1] not in leaf :
		line[1] = line[1].rstrip()+'.0' #if the cluster of the sequence id is not in the leaves list, .0 added
	table[j] = line[0]+'\t'+line[1]+'\n' #restores line format
with open('./sequence_'+parameters+".txt", 'w') as f:
	f.writelines(table)

if args.verbose >= 2 :
	print("Saving summary of sequence_{}.txt".format(parameters))

command = "awk '{x=2; print $x}' "+"{}".format("sequence_"+parameters+".txt")+" | sort | uniq -c" #creates third output file
output = subprocess.check_output(command, shell = True, stderr = subprocess.PIPE, universal_newlines = True)

#the following block of code is used to give a bit of statistical output to the user, % of orphans
total = 0
orphans = 0
for i in output.strip().split('\n') :
	total += int(i.strip().split(' ')[0])
	if i.endswith('.0') :
		orphans += int(i.strip().split(' ')[0])

with open("sequence_summary.txt", "w") as f:
	f.writelines(output)

if args.verbose >= 1 :
	print('There are {} final clusters in the classification tree (they are not divided).'.format(len(leaf)))
	#idea: add % of the clusters that are taged 0 or -test/2
	print('The proportion of orphans is {}.'.format(orphans/total))

###louis-mael.gueguen@crchudequebec.ulaval.ca v1.2 29.07.2021
