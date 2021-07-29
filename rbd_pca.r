#!/usr/bin/env Rscript

library('seqinr')
library('ggplot2')

defaultW <- getOption("warn") 
options(warn = -1)

args = commandArgs()

setwd(args[6])
directories = dir("./", pattern = '.1$')
directories = directories[!sapply(directories, is.null)]
#the global startegy is to list the directories that ends with .1. For each dir, take the root where there are the coordinates, the .2, adn get the figure 

getplot <- function(dir1) {
	dir2 = gsub('.1$', '.2', dir1) #gets the .2 folder
	parent_dir = gsub('.1$', '', dir1) #gets parent directory
	names_cluster1 = getName(read.fasta(paste(dir1, '/', 'fastaCL1.fst', sep = ''))) #gets the sequence names in fastaCL1.fst 
	names_cluster2 = getName(read.fasta(paste(dir2, '/', 'fastaCL2.fst', sep = ''))) #gets the sequence names in fastaCL2.fst
	names_parent = getName(read.fasta(paste(parent_dir, '/', list.files(parent_dir, pattern = '.fst$')[1], sep = ''))) #gets the sequence names of the parent cluster
        groups = c(rep('cluster.0', length(names_parent))) #initialises the group belonging
	df = read.table(paste(parent_dir, '/', list.files(parent_dir, pattern = '.pca$')[1], sep = ''), sep = '\t') #gets coordinates
	dfone = data.frame(names_parent, df[,1:3], stringsAsFactors = FALSE) #adds the sequence ids and selects 3 first dimensions for coo
	df = data.frame(dfone, groups, stringsAsFactors = FALSE) #adds the groups
        colnames(df) = c('seq_names', 'coo1', 'coo2', 'coo3', 'groups') #names the columns
	df$groups[which(df$seq_names %in% names_cluster1)] <- 'cluster.1' #changes values in column groups for sequences of cluster.1
        df$groups[which(df$seq_names %in% names_cluster2)] <- 'cluster.2' #same for cluster.1	
	#figures :
	g = ggplot(df, aes(coo1, coo2, colour = groups))+geom_jitter(alpha = 0.3)+xlab('First principal component')+ylab('Second principal component')+labs(color='CLuster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size=20))
	ggsave(paste(parent_dir, '/', parent_dir, '_plotPCA12', '.png', sep = ''), width = 16, height = 9, dpi = 300)
	g = ggplot(df, aes(coo3, coo2, colour = groups))+geom_jitter(alpha = 0.3)+xlab('Third principal component')+ylab('Second principal component')+labs(colour = 'Cluster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size=20))
	ggsave(paste(parent_dir, '/', parent_dir, '_plotPCA32', '.png', sep = ''), width = 16, height = 9, dpi = 300)
	g = ggplot(df, aes(coo1, coo3, colour = groups))+geom_jitter(alpha = 0.3)+xlab('First principal component')+ylab('Third principal component')+labs(colour = 'Cluster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size=20))
	ggsave(paste(parent_dir, '/', parent_dir, '_plotPCA13', '.png', sep = ''), width = 16, height = 9, dpi = 300)
}

invisible(lapply(directories, getplot)) #applies function to list of folders and silences the output

directories = dir()[!file_test('-f', dir())] #gives the list of folders in the current dir
directories = directories[!sapply(directories, is.null)]
leaves = directories[lapply(lapply(directories, list.files), length) < 7] #gets leaves of the classification tree
leaves = leaves[!sapply(leaves, is.null)]
df = read.table("cluster/counts.pca", sep = '\t') #coordinates of the sequences in the pca projection
sequences = read.table(dir()[file_test('-f', dir())][2], sep = '\t') #gets sequences names and belonging from output file of rbd
df = data.frame(sequences, df[,1:3], stringsAsFactors = FALSE)
colnames(df) = c('seqnames', 'groups', 'coo1', 'coo2', 'coo3')
df$groups <- as.character(df$groups)
df$groups[which(!(df$groups %in% leaves))] <- "orphans"
g = ggplot(df, aes(coo1, coo2, colour = groups))+geom_jitter(alpha = 0.3)+xlab('First principal component')+ylab('Second principal component')+labs(colour = 'Cluster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size = 20))
ggsave('leaves_plotPCA12.png', width = 16, height = 9, dpi = 300)
g = ggplot(df, aes(coo3, coo2, colour = groups))+geom_jitter(alpha = 0.3)+xlab('Third principal component')+ylab('Second principal component')+labs(colour = 'Cluster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size = 20))
ggsave('leaves_plotPCA32.png', width = 16, height = 9, dpi = 300)
g = ggplot(df, aes(coo1, coo3, colour = groups))+geom_jitter(alpha = 0.3)+xlab('First principal component')+ylab('Third principal component')+labs(colour = 'Cluster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size = 20))
ggsave('leaves_plotPCA13.png', width = 16, height = 9, dpi = 300)

options(warn = defaultW) #suppresses warnings
