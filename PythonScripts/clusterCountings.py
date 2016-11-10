#! /usr/bin/python
#-*- coding: utf-8 -*-

import argparse
import sys
import os

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required = True, metavar = "input counts.tsv Seqcluster", help = "CSV file, counts.tsv of Seqcluster")
	args = parser.parse_args()
	
	input_counts = args.input
	
	if os.path.exists(args.input):
		inPut = args.input
	else:
		print "No such file : " + args.input
		sys.exit(1)
	
	i = open(input_counts, "r")

	header = i.readline()
	header = header.split("\t")
	samples_names = header[3:]
	
	for n in samples_names:
		exec("file_" + n.strip() + " = open('" + n.strip() + "_Clustercount.tsv','w')")
	
	for line in i:
		l = line.split("\t")
		rawCluster = l[2]
		rawCounts = l[3:]
		noCluster = str(l[0]).strip()
		
		if rawCluster.strip() != "|":
			ClusterAnnot = rawCluster.split("|")[1]
			ClusterAnnot = ClusterAnnot.split(";")
			ann = []
			linkAnnot = "+"
			
			if len(ClusterAnnot) > 1:
				for annot in ClusterAnnot:
					a = annot.split("::")[1]
					a = a.split(",")
					a = sorted(set(a))
					ann = ann + a
				ann = sorted(set(ann))
				ann = [a.strip('"') for a in ann]
				finalAnnot = linkAnnot.join(ann)
			
			else:
				annot = rawCluster.strip().split("|")[1]
				annot = annot.split("::")[1]
				annot = annot.split(",")
				annot = sorted(set(annot))
				ann = ann + annot
				ann = [a.strip('"') for a in ann]
				finalAnnot = linkAnnot.join(ann)
				
			for n,c in zip(samples_names,rawCounts):
				exec("file_" + n.strip() + ".write(noCluster + ':' + finalAnnot + '\t' + str(" + c.strip() + "))")
				exec("file_" + n.strip() + ".write('\\n')")
			
		else:
			pass
	
	i.close()
	for n in samples_names:
		exec("file_" + n.strip() + ".close()")		
