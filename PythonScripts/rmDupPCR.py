#!/usr/bin/envpython

## !!!! IMPORTANT this script needs as stdin the stdout of the following command...  awk '{ i=(NR-1) % 4; tab[i]=$0 ; if (i==3) { print tab[1]"\t"tab[0]"\t"tab[3]"\t"tab[2]} }' | sort -T /tmp

import sys
import argparse
import time

def computeSequenceQuality(qualityLine):
	qualityScore = 0
	for char in qualityLine:
		qualityScore += ord(char)
	return qualityScore

parser = argparse.ArgumentParser()
parser.add_argument("output", metavar = "output.fastq", help = "Fastq where PCR duplicates have been removed", default = sys.stdout)
args = parser.parse_args()

outPut = args.output
o = open(outPut,"w")
# The input MUST always be the output of the awk command above! Nothing else is required, i.e. NOT a fastq file.
entree = sys.stdin

# Initialisations
seq = ""
quality = 0

nbReads = 0
nbKept  = 0

print "DATE OF THE STARTING : " + str(time.strftime("%Y-%m-%d %H:%M:%S"))
for line in entree:
	line = line.strip()
	s = line.split("\t")
	if len(s) != 4:
		sys.stderr.write("WARNING : fastq file does not have a multiple of 4 lines. Please check your file !")
		break
	seqCurr = s[0]
	readIdCurr = s[1]
	thirdLineCurr = s[3]
	qualityStrCurr = s[2]
	qualityCurr = computeSequenceQuality(qualityStrCurr)
	nbReads += 1
	if nbReads%1000000 == 0 :
		print "COUNT : " + str(nbReads) + " processed reads (" + str(nbKept)+" kept)."
	if seq == seqCurr:
		if qualityCurr > quality:
			quality = qualityCurr
		readId = readIdCurr
		thirdLine = thirdLineCurr
		qualityStr = qualityStrCurr
	else:
		if seq:
			o.write(readId + "\n")
			o.write(seq + "\n")
			o.write(thirdLine + "\n")
			o.write(qualityStr + "\n")
		seq = seqCurr
		quality = qualityCurr
		readId = readIdCurr
		thirdLine = thirdLineCurr
		qualityStr = qualityStrCurr
		nbKept += 1
if seq:
	o.write(readId + "\n")
	o.write(seq + "\n")
	o.write(thirdLine + "\n")
	o.write(qualityStr + "\n")
o.close()

## output some stats to stdout
try:
	perct = float(nbKept)/float(nbReads) * 100
	perct = "{0:.2f}".format(perct)
except ZeroDivisionError:
	sys.stderr.write("Cannot read sequences !")
	sys.exit(1)

print "TOTAL COUNT : "+str(nbKept)+" reads kept from "+str(nbReads)+" total ("+str(perct)+"%)."
print "DATE OF THE END OF THE STEP : " + str(time.strftime("%Y-%m-%d %H:%M:%S"))
