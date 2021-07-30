#!/user/bin/env python3
import argparse

parser = argparse.ArgumentParser(prog = 'script.py', formatter_class = argparse.RawDescriptionHelpFormatter, description = '''\
Takes a file of one-line sequences as input, sequences_parameters.txt from rbd, the original fasta file, returns a file with ids and last cluster belonging (assigned by rbd) from a reference fasta file as output for each match of each sequences of interest. The order follows the one of the file sequence_of_interest.txt
			''')

parser.add_argument('--seq', '-s', dest = 'seq', action = 'store', type = str, required = True, help = 'Input file with one sequence per 									line')

parser.add_argument('--fasta', '-f', dest = 'fasta', action = 'store', type = str, required = True, help = 'Input reference fasta file')

parser.add_argument('--output', '-o', dest = 'output', action = 'store', type = str, required = True, help = 'Name of the output file')

parser.add_argument('--rbd', '-r', dest = 'rbd', action = 'store', type = str, required = True, help = 'Input file which is the sequence_parameters.txt file from rbd output')

args = parser.parse_args()

with open(args.seq, 'r') as f:
	lines = f.readlines() #reads the sequences from input file.

seq_bank = [] #list with all the seq in one line

with open(args.fasta, 'r') as f:
	temp = f.readlines()
	for i in range(len(temp)) :
		if temp[i].startswith('>') : #creates a list with ids followed by the corresponding sequences: 
			seq_bank.append(temp[i]) #adds the id of the seq
			seq_bank.append(temp[i+1].rstrip()+temp[i+2].rstrip()+temp[i+3]) #adds the sequences in one string

seq = []
for i in range(len(lines)):
	for j in range(len(seq_bank)):
		if lines[i] == seq_bank[j] :
			seq.append(seq_bank[j-1]) #seq contains the sequence references of the sequences
print('Number of ids corresponding to the sequences given : {}'.format(len(seq)))

with open(args.rbd, 'r') as f :
	cluster_ref = f.readlines() #reads the sequences cluster assignement from rbd

out = []
for i in range(len(seq)):
	for j in range(len(cluster_ref)):
		if seq[i].split(' ')[0] == cluster_ref[j].split('\t')[0]:
			out.append(seq[i].rstrip().split(' ')[0]+'\t'+cluster_ref[j].split('\t')[1])
#out[] contains the ids of the sequences that corresponds to the ones passed in the input file, and their last cluster belonging

with open(args.output, 'w') as f :
	f.writelines(out)
