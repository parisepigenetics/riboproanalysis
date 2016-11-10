#! /usr/bin/python
#-*- coding: utf-8 -*-

import argparse
import sys
import os

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required = True, metavar = "input annotations", help = "Gtf file to filter")
	parser.add_argument("-o", "--output", required = True, metavar = "cds range")
	args = parser.parse_args()
	
	input_annotations = args.input
	output_cds_range = args.output
	
	if os.path.exists(args.input):
		inPut = args.input
	else:
		print "No such file : " + args.input
		sys.exit(1)
	
	i = open(input_annotations, "r")
	o = open(output_cds_range, "w")
	
	dict_transcript = {}
	
	for line in i:
		if not line.startswith("#"):
			line = line.split('\t')
			
			if line[2] == "transcript":
				
				attributes  = line[8].split(';')
				l_attributes = [info.strip() for info in attributes]

				for info in l_attributes:
					if info.find("transcript_id") != -1:
						info_transcript_id = info
				
				transcript_id = info_transcript_id[15:-1]
				
				if transcript_id not in dict_transcript:
					#dict_transcript[transcript_id] = [{},{}]
					dict_transcript[transcript_id] = [{},{},{"CDS_starts":[]},{"CDS_stops":[]}]
					
			elif line[2] == "exon":

				attributes  = line[8].split(';')
				l_attributes = [info.strip() for info in attributes]
				
				strand = str(line[6])
				
				if strand == "+":
					exon_start = int(line[3])
					exon_stop = int(line[4])
				elif strand == "-":
					exon_start = int(line[4])
					exon_stop = int(line[3])
				
				exon_lgt = abs(exon_stop - exon_start) + 1
				
				for info in l_attributes:
					if info.find("transcript_id") != -1:
						info_transcript_id = info
					elif info.find("exon_number") != -1:	
						info_exon_number = info
							
				transcript_id = info_transcript_id[15:-1]
				exon_number = int(info_exon_number[13:-1])

				if transcript_id in dict_transcript:
					dict_transcript[transcript_id][0][exon_number] = (exon_lgt, exon_start)
				
			elif line[2] == "CDS":

				strand = str(line[6])
				
				if strand == "+":
					cds_start = int(line[3])
					cds_stop = int(line[4])
				elif strand == "-":
					cds_start = int(line[4])
					cds_stop = int(line[3])
				
				cds_lgt = abs(cds_stop - cds_start) + 1
				
				attributes  = line[8].split(';')
				l_attributes = [info.strip() for info in attributes]
				
				for info in l_attributes:
					if info.find("transcript_id") != -1:
						info_transcript_id = info
					elif info.find("exon_number") != -1:
						info_exon_number = info
				
				transcript_id = info_transcript_id[15:-1]
				exon_number = int(info_exon_number[13:-1])
				
				cds_id = "CDS_" + str(cds_start) + "_" + str(cds_stop) + "_" + transcript_id
				
				if transcript_id in dict_transcript and cds_id not in dict_transcript[transcript_id][1]:
					dict_transcript[transcript_id][1][cds_id] = (cds_lgt, exon_number, cds_start)
					
	for k in dict_transcript.keys():
		first_exon_transcript = min(dict_transcript[k][0].keys())
		
		if len(dict_transcript[k][1].keys()) == 0:
			del dict_transcript[k]
			continue
			
		for k_cds in dict_transcript[k][1].keys():
			exon_start_cds = dict_transcript[k][1][k_cds][1]
			
			if exon_start_cds == first_exon_transcript:
				if dict_transcript[k][1][k_cds][2] == dict_transcript[k][0][first_exon_transcript][1]:
					start = 0
				else:
					start = abs(dict_transcript[k][1][k_cds][2] - dict_transcript[k][0][first_exon_transcript][1])
			else:
				if dict_transcript[k][0][exon_start_cds][1] == dict_transcript[k][1][k_cds][2]:
					start = 0
					for x in range(1, exon_start_cds):
						#start = start + (dict_transcript[k][0][x][0] - 1)
						start = start + dict_transcript[k][0][x][0]
				else:
					start = abs(dict_transcript[k][1][k_cds][2] - dict_transcript[k][0][exon_start_cds][1])
					for x in range(1, exon_start_cds):
						#start = start + (dict_transcript[k][0][x][0] - 1)
						start = start + dict_transcript[k][0][x][0]
						
			stop = start + dict_transcript[k][1][k_cds][0] + 1
			
			dict_transcript[k][2]["CDS_starts"].append(start)
			dict_transcript[k][3]["CDS_stops"].append(stop)
			
		o.write(k + "\t" + str(min(dict_transcript[k][2]["CDS_starts"])) + "\t" + str(max(dict_transcript[k][3]["CDS_stops"])) + "\n")
				
	i.close()
	o.close()
