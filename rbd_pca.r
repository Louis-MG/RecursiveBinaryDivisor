#!/usr/bin/env Rscript

library('seqinr')
library('ggplot2')


args = commandArgs()

setwd(args[6])
directories = dir("./", pattern = '.1$') 

#the global startegy is to list the directories that ends with .1. For each dir, take the root, the .2, do the trick and get the figure 

getplot <- function(dir1) {
	dir2 = gsub('.1$', '.2', dir1)
	parent_dir = gsub('.1$', '', dir1)
	files = list.files(dir1) #gets all the file names in the dir .1
	names_cluster1 = getName(read.fasta(paste(dir1, '/', tail(files, n=1), sep = ''))) #gets the sequence names in fastaCL1.fst 
	files = list.files(dir2) #gets all the file names in the dir .2
	names_cluster2 = getName(read.fasta(paste(dir2, '/', tail(files, n=1), sep = ''))) #gets the sequence names in fastaCL2.fst
	files = list.files(parent_dir)
	names_parent = getName(read.fasta(paste(parent_dir, '/', tail(files, n=1), sep = ''))) #gets the sequence names of the parent cluster
        groups = c(rep('0', length(names_parent)))
	df = read.table(paste(parent_dir, '/', files[3], sep = ''), sep = '\t')
	df = data.frame(names_parent, df[,1:2])
	df = data.frame(df, groups, stringAsFactor = FALSE)
        colnames(df) = c('seq_names', 'coo1', 'coo2', 'coo3', 'groups')
	df$groups[which(df$seq_names %in% names_cluster1)] = '1' #changes values in column groups for sequences of cluster.1
        df$groups[which(df$seq_names %in% names_cluster2)] = '2' #same for cluster.1
	g = ggplot(df[,2:3], aes(df$coo1, df$coo2, colour = df$groups))+geom_jitter(alpha = 0.3)+xlab('First principal component')+ylab('Second principal component')+labs(color='CLuster belonging of the sequences')+theme(axis.text = element_text(size=15), axis.title = element_text(size=20))
	ggsave(paste(parent_dir, '/', 'plotPCA', parent_dir,'.png', sep = ''), width = 16, height = 9, dpi = 300)
}

lapply(directories, getplot)

