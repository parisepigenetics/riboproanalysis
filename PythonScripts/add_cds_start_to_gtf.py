#! /usr/bin/python
#-*- coding: utf-8 -*-

import argparse
import sys
import os

# FIXME If a python script does not have functions which can be imported as a module in another program, there is no need for this __name == "__main__" trick.
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required = True, metavar = "input annotations", help = "Gtf file to filter")
	parser.add_argument("-o", "--output", required = True, metavar = "Gtf file with cds range")
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

	features_of_interest = ["exon","CDS","start_codon","stop_codon"]
	n = 0
	for line in i:
		if not line.startswith("#"):
			line = line.strip().split('\t')

			if line[2] == "transcript":
				attributes  = line[8].split(';')
				l_attributes = [info.strip() for info in attributes]

				for info in l_attributes:
					if info.find("transcript_id") != -1:
						info_transcript_id = info

				transcript_id = info_transcript_id[15:-1]

				if transcript_id not in dict_transcript:
					dict_transcript[transcript_id] = [{},{},{},{},{},{}] # Dict transcript : key = transcript_id, value = list with 0 = dict exons lgt, 1 = dict CDSs lgt, 2 = dict cdsStart & cdsEnd, 3 = strand +/-, 4 = startExon1 & stopExon1, 5 = startCDS1 & stopCDS1
					dict_transcript[transcript_id][3]["strand"] = line[6].strip()

			elif line[2] == "exon":

				attributes  = line[8].split(';')
				l_attributes = [info.strip() for info in attributes]

				exon_start = int(line[3])
				exon_stop = int(line[4])

				exon_lgt = abs(exon_stop - exon_start) + 1

				for info in l_attributes:
					if info.find("transcript_id") != -1:
						info_transcript_id = info
					elif info.find("exon_number") != -1:
						info_exon_number = info

				transcript_id = info_transcript_id[15:-1]
				exon_number = int(info_exon_number[13:-1])

				if transcript_id in dict_transcript:
					dict_transcript[transcript_id][0][exon_number] = exon_lgt
					if exon_number == 1:
						dict_transcript[transcript_id][4]["startGenoExon1"] = int(line[3])
						dict_transcript[transcript_id][4]["stopGenoExon1"] = int(line[4])

			elif line[2] == "CDS":

				cds_start = int(line[3])
				cds_stop = int(line[4])

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

				if transcript_id in dict_transcript:
					dict_transcript[transcript_id][1][exon_number] = cds_lgt
					if exon_number == 1:
						dict_transcript[transcript_id][5]["startGenoCDS1"] = int(line[3])
						dict_transcript[transcript_id][5]["stopGenoCDS1"] = int(line[4])
		n = n + 1
		#print n

	for k in dict_transcript.keys():
		transcript_lgt = sum(dict_transcript[k][0].values())
		cds_length = sum(dict_transcript[k][1].values())

		if (cds_length == 0) or (cds_length % 3 != 0):
			del dict_transcript[k]
			continue

		#assert cds_length % 3 == 0, "CDS length of " + k + " not divisible by 3 !"

		exon_number_stop = max(dict_transcript[k][1].keys())

		#stop = start + cds_length

		if len(dict_transcript[k][0]) == len(dict_transcript[k][1]):
			if dict_transcript[k][0][1] == dict_transcript[k][1][1]:
					start = 0
			else:
				if dict_transcript[k][4]["startGenoExon1"] == dict_transcript[k][5]["startGenoCDS1"]:
					start = 0
				else:
					start = dict_transcript[k][0][1] - dict_transcript[k][1][1]
		else:
			exon_number_start = min(dict_transcript[k][1].keys())

			start = 0
			for x in range(1,exon_number_start):
				start += dict_transcript[k][0][x]
			#print k
			#print dict_transcript

			d_cdsStart_exon = dict_transcript[k][0][exon_number_start] - dict_transcript[k][1][exon_number_start]
			start += d_cdsStart_exon

		stop = start + cds_length

		dict_transcript[k][2]["cdsStart"] = start
		dict_transcript[k][2]["cdsEnd"] = stop

		#print str(dict_transcript[k][2]["cdsStart"]) + "###"
		#print str(dict_transcript[k][2]["cdsEnd"]) + "###"
		#if stop == transcript_lgt:
		#	del dict_transcript[k]
		#	continue

		#o.write(k + "\t" + str(start) + "\t" + str(stop) + "\n")
		#print(dict_transcript)
	#print(dict_transcript)
	i.seek(0)

	n = 0
	for line in i:
		if not line.startswith("#"):
			line = line.strip()
			l = line.strip().split("\t")

			if l[2] in features_of_interest:
				attributes  = l[8].split(';')
				l_attributes = [info.strip() for info in attributes]

				for info in l_attributes:
					if info.find("transcript_id") != -1:
						info_transcript_id = info

				transcript_id = info_transcript_id[15:-1]

				if transcript_id in dict_transcript.keys():
					line = line + " cds_start " + '"' + str(dict_transcript[transcript_id][2]["cdsStart"]) + '"' + "; cds_end " + '"' + str(dict_transcript[transcript_id][2]["cdsEnd"]) + '"' + ";"
					o.write(line + "\n")
		#n = n + 1
		#print(str(n))
	i.close()
	o.close()
