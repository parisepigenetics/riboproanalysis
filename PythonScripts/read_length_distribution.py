#!/usr/bin/env python

'''Python script to generate a histogram plot from a single column integer number file (or STDIN).'''

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs = "?", help = "Input single column number file, or empty for STDIN", default = sys.stdin)
parser.add_argument("output", help = "Read length distribution figure.")
parser.add_argument("-e", "--experiment", dest = "experiment", help = "A sigle word describing the experiment.", required = True)
args = parser.parse_args()

# Directly read the data.
readLenghts = np.loadtxt(args.input, dtype = int)

nBins = max(readLenghts) - min(readLenghts)

plt.hist(readLenghts, bins = nBins)

plt.title('Reads length distribution from ' + args.experiment)
plt.xlabel("Read lengths")
plt.ylabel('Number of reads')
plt.grid(True)
plt.xlim(min(readLenghts) - 1, max(readLenghts) + 1)
plt.savefig(args.output)
