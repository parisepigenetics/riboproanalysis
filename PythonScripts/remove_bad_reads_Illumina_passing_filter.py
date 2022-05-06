#!/usr/bin/env python

import sys
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument("infileh", nargs='?', type = argparse.FileType('r'), metavar = "input_fastq", help = "A fastq file with all the reads.")
parser.add_argument("outfileh", nargs='?', type = argparse.FileType('w'), metavar = "output_fastq", default = sys.stdout, help = "A fastq file with reads passed the quality control")
args = parser.parse_args()

infileh = args.infileh
outfileh = args.outfileh

nbKept = 0
nbReads = 0

sys.stderr.write("START TIME: " + str(time.strftime("%Y-%m-%d %H:%M:%S")) + "\n")
while True:
	line = infileh.readline().strip()
	if line == "":
		break
	res = line.split(":")[-3]
	if res == "N":
		outfileh.write(line + "\n" + infileh.readline() + infileh.readline() + infileh.readline())
		nbKept += 1
	else:	
	  infileh.readline()
	  infileh.readline()
	  infileh.readline()
	nbReads += 1

infileh.close()
outfileh.close()

## output some stats to stdout
perct = float(nbKept)/float(nbReads) * 100
perct = "{0:.2f}".format(perct)
sys.stderr.write("TOTAL COUNT: " + str(nbKept) + " reads kept from " + str(nbReads) + " total (" + perct + "%)." + "\n")
sys.stderr.write("INFO: reads have been written to " + outfileh.name + " in FASTQ format." + "\n")
sys.stderr.write("DATE OF THE END OF THE STEP: " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
